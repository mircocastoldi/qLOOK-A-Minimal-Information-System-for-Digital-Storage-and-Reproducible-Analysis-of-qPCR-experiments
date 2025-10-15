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
# qPCR normalization + ΔΔCt + 2^-ΔΔCt + per-gene boxplots
# Script name: qLOOK_Module2_v1.0.R
# Script version: 1.0
# Adds logging to qLOOK_Summary.txt and remembers last folder
# October 2025
# Version: v1.0 (14-10-2025)
# ---------------------------------------------------------------


# ---- Script metadata ----
script_name <- "qLOOK_Module2_v1.0.R"
script_version <- "v1.0 (14-10-2025)"


# ---- 1) Load/install required packages ----
required_pkgs <- c("openxlsx", "dplyr", "tidyr", "ggplot2")
for(pkg in required_pkgs){
  if(!requireNamespace(pkg, quietly = TRUE)){
    install.packages(pkg, quiet = TRUE, dependencies = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# ---- 2) Remember last used folder ----
last_folder_file <- file.path(tempdir(), "qLOOK_last_folder.txt")
if (file.exists(last_folder_file)) {
  last_folder <- readLines(last_folder_file, warn = FALSE)[1]
  if (!dir.exists(last_folder)) last_folder <- getwd()
} else {
  last_folder <- getwd()
}

# ---- 3) Choose Excel file interactively ----
choose_file <- function(start_folder = getwd()){
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()){
    f <- rstudioapi::selectFile(caption = "Select raw qPCR Excel file", path = start_folder, filter = "xlsx")
  } else if (.Platform$OS.type == "windows"){
    f <- utils::choose.files(default = file.path(start_folder, "*.xlsx"),
                             caption = "Select raw qPCR Excel file",
                             filters = matrix(c("Excel", "*.xlsx"), 1, 2, byrow=TRUE))
  } else {
    f <- readline("Enter path to raw qPCR Excel file: ")
  }
  if(length(f)==0 || f=="") stop("No file selected")
  normalizePath(f, mustWork = TRUE)
}

# ---- 4) Normalize qPCR data (ΔCt) ----
normalize_qpcr <- function(df, ref_genes){
  df[, -1] <- lapply(df[, -1], as.numeric)
  missing_genes <- setdiff(ref_genes, names(df))
  if(length(missing_genes) > 0){
    stop(paste("Reference genes missing:", paste(missing_genes, collapse=", ")))
  }
  ref_values <- df[, ref_genes, drop = FALSE]
  geo_mean <- apply(ref_values, 1, function(x) exp(mean(log(x[x > 0]))))
  norm_df <- df
  for(g in names(df)[-1]){
    norm_df[[g]] <- df[[g]] - geo_mean
  }
  return(norm_df)
}

# ---- 5) Calculate ΔΔCt using calibrator ----
calculate_ddcq <- function(norm_df, calibrator_sample){
  df <- norm_df
  df[,1] <- trimws(as.character(df[,1]))
  df[, -1] <- lapply(df[, -1], as.numeric)
  calib_rows <- which(df[,1] %in% calibrator_sample)
  if(length(calib_rows)==0) stop("No calibrator sample found.")
  calib_mean <- colMeans(df[calib_rows, -1, drop=FALSE], na.rm=TRUE)
  ddq_df <- df
  for(g in names(df)[-1]){
    ddq_df[[g]] <- df[[g]] - calib_mean[g]
  }
  fold_df <- ddq_df
  for(g in names(ddq_df)[-1]){
    fold_df[[g]] <- 2^(-ddq_df[[g]])
  }
  list(ddc_df = df, ddq_df = ddq_df, fold_df = fold_df)
}

# ---- 6) Generate one boxplot per gene ----
plot_boxplots <- function(fold_df, folder){
  long_df <- tidyr::pivot_longer(fold_df, cols=-1, names_to="Gene", values_to="FoldChange")
  colnames(long_df)[1] <- "Sample"
  plot_folder <- file.path(folder, paste0("qLOOK_Norm_Plots_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  if(!dir.exists(plot_folder)) dir.create(plot_folder)
  genes <- unique(long_df$Gene)
  for(g in genes){
    plot_data <- long_df %>% dplyr::filter(Gene == g)
    p <- ggplot(plot_data, aes(x=Sample, y=FoldChange, fill=Sample)) +
      geom_boxplot(outlier.shape=NA, alpha=0.7) +
      geom_jitter(width=0.2, size=1.5, color="black") +
      labs(title=paste0("Relative Expression 2^-ΔΔCq - ", g),
           y="2^-ΔΔCq", x="Sample") +
      theme_minimal(base_size=12) +
      theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none")
    ggsave(file.path(plot_folder, paste0(g, "_boxplot.png")), plot=p, width=6, height=4)
  }
  message("✅ Boxplots saved in folder: ", plot_folder)
  return(plot_folder)
}

# ---- 7) Main workflow ----
file <- choose_file(last_folder)
raw_df <- openxlsx::read.xlsx(file)
folder <- dirname(file)
writeLines(folder, last_folder_file)  # remember this folder
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ---- 8) Ask for reference genes ----
ref_genes_input <- utils::winDialogString(
  "Enter reference gene(s) for normalization, comma-separated", 
  default = "RefGene"
)
ref_genes <- strsplit(ref_genes_input, ",")[[1]]
ref_genes <- trimws(ref_genes)
cat("Reference genes selected:", paste(ref_genes, collapse=", "), "\n")

normalized_df <- normalize_qpcr(raw_df, ref_genes)
normalized_file <- file.path(folder, paste0("qLOOK_Norm_", ts, ".xlsx"))
openxlsx::write.xlsx(normalized_df, normalized_file)
message("✅ Normalized ΔCt saved: ", normalized_file)

# ---- 9) Ask for calibrator samples ----
calibrator_input <- utils::winDialogString(
  "Enter calibrator sample(s) for relative quantities, comma-separated",
  default = raw_df[[1]][1]
)
calibrator_sample <- strsplit(calibrator_input, ",")[[1]]
calibrator_sample <- trimws(calibrator_sample)
cat("Calibrator samples selected:", paste(calibrator_sample, collapse=", "), "\n")

results <- calculate_ddcq(normalized_df, calibrator_sample)
express_file <- file.path(folder, paste0("qLOOK_Express_", ts, ".xlsx"))
openxlsx::write.xlsx(
  list("DDCq"=results$ddc_df,
       "DeltaDeltaCq"=results$ddq_df,
       "FoldChange_2^-DDCq"=results$fold_df),
  express_file
)
message("✅ ΔΔCt and 2^-ΔΔCt saved: ", express_file)

# ---- 10) Generate boxplots per gene ----
plot_folder <- plot_boxplots(results$fold_df, folder)

# ---- 11) Logging step ----
log_file <- file.path(folder, "qLOOK_Summary.txt")
log_entry <- paste0(
  "===== qPCR Normalization and Relative Expression Log =====", "\n",
  "\n--- ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ---\n",
  "Script: ", script_name, " (", script_version, ")\n",
  "Input File: ", basename(file), "\n",
  "Reference Gene(s): ", paste(ref_genes, collapse=", "), "\n",
  "Calibrator Sample(s): ", paste(calibrator_sample, collapse=", "), "\n",
  "Created Files:\n  - ", basename(normalized_file), "\n  - ",
  basename(express_file), "\n  - Boxplots Folder: ", basename(plot_folder), "\n",
  "\n"
)

if (file.exists(log_file)) {
  existing <- readLines(log_file, warn = FALSE)
  updated <- c(existing, "", log_entry)
  writeLines(updated, log_file)
} else {
  writeLines(log_entry, log_file)
}
message("🧾 Log updated in: ", log_file)

message("\n🎉 Analysis completed successfully!\n")
