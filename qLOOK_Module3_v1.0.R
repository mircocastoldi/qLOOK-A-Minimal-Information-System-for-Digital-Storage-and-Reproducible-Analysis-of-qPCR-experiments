# ============ License ======================================================
# qLOOK (qPCR-LOg-boOK)
# Copyright (c) 2025 [Mirco Castoldi]
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
# Description: Modular R scripts for extraction, normalization, statistical
# analysis, and visualization of qPCR data from Thermo Fisher/ABI *.EDS* files.
# ==============================================================================

# ---------------------------------------------------------------
# Author: MIRCO CASTOLDI (for RStudio)
# GraphPad-ready export script + Boxplots + Heatmap + Enhanced T-tests & ANOVA
# Script name: qLOOK_Module3_v1.0.R
# Script version: 1.0
#  - Logs p-value cutoff and all output files (including plots)
#  - Logs script name and version for traceability
# October 2025
# Version: v1.0 (14-10-2025)
# ---------------------------------------------------------------

# ---- Script metadata ----
script_name <- "qLOOK_Module3_v1.0.R"
script_version <- "v1.0 (14-10-2025)"

# ---- Automatically install & load required packages silently ----
required_pkgs <- c("openxlsx", "tidyr", "dplyr", "ggplot2", "rstatix", "scales")
for(pkg in required_pkgs){
  if(!requireNamespace(pkg, quietly = TRUE)){
    install.packages(pkg, quiet = TRUE, dependencies = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# ---- Choose Excel file (remembers last folder) ----
choose_file <- function(){
  last_dir_file <- file.path(tempdir(), "last_used_dir.txt")
  init_dir <- if (file.exists(last_dir_file)) readLines(last_dir_file, warn = FALSE)[1] else getwd()
  
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()){
    f <- rstudioapi::selectFile(caption = "Select Excel file", filter = "xlsx", path = init_dir)
  } else if (.Platform$OS.type == "windows"){
    f <- utils::choose.files(default = file.path(init_dir, "*.xlsx"), caption = "Select Excel file",
                             filters = matrix(c("Excel", "*.xlsx"), 1, 2, byrow = TRUE))
  } else {
    f <- readline("Enter path to Excel file: ")
  }
  
  if(length(f) == 0 || f == "") stop("No file selected")
  
  writeLines(dirname(f), last_dir_file)
  normalizePath(f, mustWork = TRUE)
}

# ---- Main function ----
make_graphpad_file <- function(file_path){
  # Read fold change sheet
  df <- openxlsx::read.xlsx(file_path, sheet = "FoldChange_2^-DDCq")
  colnames(df)[1] <- "Sample"
  df_long <- tidyr::pivot_longer(df, cols = -Sample, names_to = "Gene", values_to = "Value")
  
  gene_list <- list()
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  folder <- dirname(file_path)
  wb_anova <- openxlsx::createWorkbook()
  wb_ttest <- openxlsx::createWorkbook()
  
  # Create plot folder
  plot_folder <- file.path(folder, paste0("qLOOK_Plots_", ts))
  if(!dir.exists(plot_folder)) dir.create(plot_folder)
  
  # ---- Heatmap ----
  heatmap_file <- file.path(plot_folder, paste0("Heatmap_", ts, ".png"))
  df_heatmap <- df_long %>%
    group_by(Gene, Sample) %>%
    summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")
  
  p_heatmap <- ggplot(df_heatmap, aes(x = Sample, y = Gene, fill = mean_value)) +
    geom_tile(color = "gray70") +
    scale_fill_gradient2(
      low = "orange2", mid = "gray90", high = "red4",
      midpoint = median(df_heatmap$mean_value, na.rm = TRUE),
      name = "Expression"
    ) +
    labs(title = "Gene Expression Heatmap", x = "Sample", y = "Gene") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))
  
  ggsave(heatmap_file, plot = p_heatmap, width = 9, height = 7, dpi = 300)
  
  # ---- Ask for p-value cutoff ----
  pval_input <- utils::winDialogString("Enter p-value cutoff (e.g., 0.05):", default = "0.05")
  p_cutoff <- as.numeric(trimws(pval_input))
  if(is.na(p_cutoff) || p_cutoff <= 0 || p_cutoff > 1) stop("Invalid p-value cutoff entered.")
  message("Using p-value cutoff: ", p_cutoff)
  
  # ---- Process each gene ----
  for(g in unique(df_long$Gene)){
    tmp <- df_long %>% filter(Gene == g)
    
    # GraphPad-style table
    tmp_wide <- tmp %>% 
      group_by(Sample) %>%
      mutate(rep_id = row_number()) %>%
      pivot_wider(names_from = Sample, values_from = Value) %>%
      select(-rep_id)
    gene_list[[g]] <- tmp_wide
    
    # Boxplot
    boxplot_file <- file.path(plot_folder, paste0(g, "_boxplot.png"))
    p <- ggplot(tmp, aes(x = Sample, y = Value, fill = Sample)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
      labs(title = paste("Fold Change (2^-ΔΔCq):", g),
           y = "Fold Change (2^-ΔΔCq)", x = "Sample") +
      theme_bw(base_size = 14) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(boxplot_file, plot = p, width = 7, height = 5, dpi = 300)
    
    # ANOVA
    anova_res <- aov(Value ~ Sample, data = tmp)
    anova_df <- as.data.frame(summary(anova_res)[[1]])
    colnames(anova_df) <- c("Df", "Sum_Sq", "Mean_Sq", "F_value", "Pr(>F)")
    openxlsx::addWorksheet(wb_anova, g)
    openxlsx::writeData(wb_anova, g, anova_df)
    
    # t-tests
    stats_summary <- tmp %>% group_by(Sample) %>%
      summarise(mean = mean(Value, na.rm = TRUE),
                sd = sd(Value, na.rm = TRUE), .groups = "drop")
    
    ttest_df <- tmp %>%
      rstatix::pairwise_t_test(Value ~ Sample, p.adjust.method = "bonferroni") %>%
      as.data.frame() %>%
      left_join(stats_summary, by = c("group1" = "Sample")) %>%
      rename(mean_group1 = mean, sd_group1 = sd) %>%
      left_join(stats_summary, by = c("group2" = "Sample")) %>%
      rename(mean_group2 = mean, sd_group2 = sd) %>%
      mutate(
        fold_change = mean_group2 / mean_group1,
        direction = case_when(
          fold_change > 1.1 ~ "UP",
          fold_change < 0.9 ~ "DOWN",
          TRUE ~ "UNCH"
        )
      ) %>%
      select(group1, group2, p, p.adj, p.signif,
             mean_group1, sd_group1, mean_group2, sd_group2,
             fold_change, direction)
    
    filtered_df <- ttest_df %>% filter(p <= p_cutoff)
    
    openxlsx::addWorksheet(wb_ttest, g)
    openxlsx::writeData(wb_ttest, g, ttest_df)
    
    filtered_name <- paste0(g, "_p<", gsub("\\.", "", as.character(p_cutoff)))
    openxlsx::addWorksheet(wb_ttest, substr(filtered_name, 1, 31))
    openxlsx::writeData(wb_ttest, substr(filtered_name, 1, 31), filtered_df)
  }
  
  # ---- Save outputs ----
  out_xlsx <- file.path(folder, paste0("qLOOK_GraphPad_", ts, ".xlsx"))
  anova_file <- file.path(folder, paste0("qLOOK_ANOVA_", ts, ".xlsx"))
  ttest_file <- file.path(folder, paste0("qLOOK_TTEST_", ts, ".xlsx"))
  
  openxlsx::write.xlsx(gene_list, out_xlsx)
  openxlsx::saveWorkbook(wb_anova, anova_file, overwrite = TRUE)
  openxlsx::saveWorkbook(wb_ttest, ttest_file, overwrite = TRUE)
  
  # ---- Collect all .png plots ----
  png_files <- list.files(plot_folder, pattern = "\\.png$", full.names = FALSE)
  if (length(png_files) == 0) png_files <- "(none)"
  
  # ---- Log info ----
  log_file <- file.path(folder, "qLOOK_Summary.txt")
  log_entry <- paste0(
    "===== Statistical Analysis Log =====", "\n",
    "\n--- ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ---\n",
    "Script: ", script_name, " (", script_version, ")\n",
    "Input Excel: ", basename(file_path), "\n",
    "p-value cutoff: ", p_cutoff, "\n",
    "Created files:\n",
    "   - ", basename(out_xlsx), "\n",
    "   - ", basename(anova_file), "\n",
    "   - ", basename(ttest_file), "\n",
    "Plots folder: ", basename(plot_folder), "\n",
    "Plots:\n   - ", paste(png_files, collapse = "\n   - "), "\n", "\n"
  )
  cat(log_entry, file = log_file, append = TRUE)
  
  # ---- Messages ----
  message("✅ GraphPad-ready file saved: ", out_xlsx)
  message("✅ Boxplots + Heatmap saved in: ", plot_folder)
  message("✅ ANOVA results saved: ", anova_file)
  message("✅ T-test results saved: ", ttest_file)
  message("✅ Logged results + p-value cutoff + script info to: ", log_file)
  
  return(list(graphpad = out_xlsx, anova = anova_file,
              ttest = ttest_file, heatmap = heatmap_file))
}

# ---- Run ----
file <- choose_file()
make_graphpad_file(file)
