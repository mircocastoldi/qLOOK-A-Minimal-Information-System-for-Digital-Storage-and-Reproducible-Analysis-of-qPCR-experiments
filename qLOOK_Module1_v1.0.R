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


# ============================================================
# Author: MIRCO CASTOLDI (for RStudio)
# qLOOK Integrated Script - qLOOK_Module1_v1.0.R
# This script integrates the TXT and JSON processing modules to handle both types of .EDS files.
# It detects the type per .EDS file, processes accordingly, stacks the results, generates qLOOK_Data,
# and then runs the reference gene analysis. Now logs whether each .EDS was TXT or JSON.
# Reference gene identified using the script qLOOK_RefGene_v1.R
# Version: v1.0 (13-10-2025)
# ============================================================

# ---- Load all required packages from the original scripts ----
# (Unchanged from original)
required_packages <- c("readr", "dplyr", "openxlsx", "tcltk", "stringr", "fs", "tidyr",
                       "archive", "jsonlite", "writexl", "rstudioapi", "readxl")
for(pkg in required_packages){
  if(!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# ---- Folder selection popup (consistent with originals) ----
# (Unchanged from original)
choose_folder <- function(){
  if (.Platform$OS.type == "windows") {
    fld <- utils::choose.dir(default = getwd(), caption = "Select folder containing .EDS files")
  } else if (requireNamespace("tcltk", quietly = TRUE)) {
    fld <- tcltk::tk_choose.dir(default = getwd(), caption = "Select folder containing .EDS files")
  } else {
    fld <- readline("Enter path to folder containing .EDS files: ")
  }
  if(is.na(fld) || fld == "") stop("No folder selected. Exiting.")
  normalizePath(fld, winslash = "/", mustWork = TRUE)
}
folder_path <- choose_folder()

# ---- Global timestamp ----
# (Unchanged from original)
timestamp_global <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ---- Vectors to track files and types ----
created_files <- character(0)
all_dfs <- list()  # To collect processed data frames from all .EDS
eds_file_types <- list()  # To track whether each .EDS was TXT or JSON

# ---- Helper from TXT script: find column name ----
# (Unchanged from original)
find_col <- function(cols, patterns){
  cols_norm <- tolower(gsub("\\s+", "", cols))
  for(p in patterns){
    idx <- which(cols_norm == tolower(gsub("\\s+","",p)))
    if(length(idx)) return(cols[idx[1]])
  }
  # fallback: try partial match
  for(p in patterns){
    idx <- grep(tolower(p), cols_norm, fixed = TRUE)
    if(length(idx)) return(cols[idx[1]])
  }
  return(NA_character_)
}

# ---- Function to process a single TXT-type .EDS ----
# (Unchanged from original)
process_txt_eds <- function(eds_file, temp_dir) {
  # Recursively locate analysis_result.txt (case-insensitive)
  all_extracted_files <- fs::dir_ls(temp_dir, recurse = TRUE, type = "file")
  analysis_candidates <- all_extracted_files[grepl("(?i)analysis_result\\.txt$", fs::path_file(all_extracted_files), perl = TRUE)]
  if(length(analysis_candidates) == 0){
    warning("analysis_result.txt not found inside ", eds_file)
    return(NULL)
  }
  analysis_file <- analysis_candidates[1]
  
  # Read the file and remove the first row (suppress parsing warnings)
  data <- tryCatch(
    suppressWarnings(readr::read_tsv(analysis_file, skip = 1, col_types = readr::cols(.default = "c"))),
    error = function(e){
      warning("Failed reading analysis_result.txt in ", eds_file, " — skipping. (", e$message, ")")
      return(NULL)
    }
  )
  if(is.null(data)) return(NULL)
  
  # Select relevant columns and rename
  cols <- names(data)
  well_col     <- find_col(cols, c("Well"))
  sample_col   <- find_col(cols, c("SampleName","Sample","SampleName"))
  detector_col <- find_col(cols, c("Detector","Target","Gene"))
  ct_col       <- find_col(cols, c("Ct","Cq","CtMean"))
  
  if(any(is.na(c(well_col,sample_col,detector_col,ct_col)))){
    warning("Required columns not found in ", fs::path_file(eds_file), ". Found columns: ", paste(cols, collapse = ", "), " — skipping.")
    return(NULL)
  }
  
  # Keep only rows where 'Well' is numeric (suppress coercion warnings)
  data_sel <- data %>%
    dplyr::select(all_of(c(well_col, sample_col, detector_col, ct_col)))
  names(data_sel) <- c("Well","Sample","Gene","Ct")
  
  # keep rows with numeric Well
  data_sel <- suppressWarnings(data_sel %>% dplyr::filter(!is.na(as.numeric(Well))))
  
  # Add 'Type' column between 'Well' and 'Sample'
  data_sel <- data_sel %>%
    dplyr::mutate(Type = "UNKN") %>%
    dplyr::relocate(Type, .after = Well) %>%
    dplyr::mutate(Quantity = NA_real_) %>%
    dplyr::relocate(Quantity, .after = Ct)
  
  # Clean Sample, Gene, Ct; remove "NTC"
  data_sel <- data_sel %>%
    dplyr::mutate(
      # Replace empty or NA Sample/Gene
      Sample = dplyr::case_when(is.na(Sample) | trimws(Sample) == "" ~ "MISS", TRUE ~ Sample),
      Gene   = dplyr::case_when(is.na(Gene)   | trimws(Gene)   == "" ~ "NoValue", TRUE ~ Gene),
      # Numeric conversion and cleaning of Ct
      Ct_tmp = suppressWarnings(as.numeric(Ct)),
      Ct = dplyr::case_when(
        is.na(Ct_tmp) ~ 36,
        Ct_tmp > 36  ~ 36,
        TRUE         ~ Ct_tmp
      )
    ) %>%
    dplyr::select(-Ct_tmp) %>%
    
    # --- Normalize Sample for robust NTC filtering ---
    dplyr::mutate(Sample_clean = toupper(trimws(gsub("_", "", Sample)))) %>%
    dplyr::filter(!grepl("^NTC", Sample_clean)) %>%   # remove NTC, ntc, NTC_1, etc.
    dplyr::select(-Sample_clean) %>%                  # drop helper column
    
    # --- Final cleanup ---
    dplyr::mutate(
      Sample = stringr::str_replace_all(Sample, "\\s+", "_"),
      Gene   = stringr::str_replace_all(Gene, "\\s+", "_")
    )
  
  return(data_sel)
}

# ---- Function to process a single JSON-type .EDS ----
# (Unchanged from original)
process_json_eds <- function(eds_file, temp_dir) {
  # Recursively locate relative_quantification_result.json (case-insensitive)
  all_extracted_files <- fs::dir_ls(temp_dir, recurse = TRUE, type = "file")
  json_candidates <- all_extracted_files[grepl("(?i)relative_quantification_result\\.json$", fs::path_file(all_extracted_files), perl = TRUE)]
  if(length(json_candidates) == 0){
    warning("relative_quantification_result.json not found inside ", eds_file)
    return(NULL)
  }
  json_file_path <- json_candidates[1]
  
  # Read and clean JSON
  raw_json <- readLines(json_file_path, warn = FALSE)
  raw_json <- paste(raw_json, collapse = "\n")
  raw_json_clean <- gsub("\\bNaN\\b", "null", raw_json)
  raw_json_clean <- gsub("\\bInfinity\\b", "null", raw_json_clean)
  raw_json_clean <- gsub("\\-Infinity\\b", "null", raw_json_clean)
  
  # Parse and flatten JSON
  parsed_json <- fromJSON(raw_json_clean, flatten = TRUE)
  
  # Extract data frame
  if (is.data.frame(parsed_json)) {
    data <- parsed_json
  } else if (is.list(parsed_json)) {
    list_lengths <- sapply(parsed_json, function(x) is.list(x) && length(x) > 1)
    if (any(list_lengths)) {
      data <- parsed_json[[which.max(list_lengths)]]
    } else {
      warning(paste("Unable to find tabular data in JSON inside", basename(eds_file), "- skipping."))
      return(NULL)
    }
  } else {
    warning(paste("Unsupported JSON structure inside", basename(eds_file), "- skipping."))
    return(NULL)
  }
  
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  
  # Remove reactions column if present
  if ("reactions" %in% names(data)) {
    data$reactions <- NULL
  }
  
  # === Format output dataframe as per user specs ===
  required_cols <- c("sampleName", "targetName", "cqMean")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    warning(paste("Missing columns in JSON data:", paste(missing_cols, collapse = ", "), "- skipping file."))
    return(NULL)
  }
  
  n <- nrow(data)
  
  # Clean and prepare Ct column
  Ct <- suppressWarnings(as.numeric(data$cqMean))
  Ct[is.na(Ct) | Ct > 36] <- 36
  
  # Build new dataframe
  formatted_df <- data.frame(
    well = integer(n),
    Type = character(n),
    Sample = as.character(data$sampleName),
    Gene = as.character(data$targetName),
    Ct = Ct,
    Quantity = character(n),
    stringsAsFactors = FALSE
  )
  
  # Fill 'well' and 'Type' columns:
  current_well <- 1
  formatted_df$well[1] <- current_well
  formatted_df$Type[1] <- "UNKN"
  
  for (i in 2:n) {
    if (formatted_df$Sample[i] == formatted_df$Sample[i-1]) {
      current_well <- current_well + 1
    } else {
      current_well <- 1
    }
    formatted_df$well[i] <- current_well
    formatted_df$Type[i] <- "UNKN"
  }
  
  # Quantity stays empty ("")
  
  return(formatted_df)
}

# ---- List all .EDS files ----
# (Unchanged from original)
all_files <- fs::dir_ls(folder_path, recurse = FALSE, type = "file")
eds_files <- all_files[grepl("(?i)\\.eds$", fs::path_file(all_files), perl = TRUE)]
if(length(eds_files) == 0) stop("No .EDS files found in selected folder.")

# ---- Process each .EDS file (modified to track file type) ----
for(eds_file in eds_files){
  message("Processing: ", fs::path_file(eds_file))
  temp_dir <- tempfile(pattern = "eds_extract_")
  dir.create(temp_dir)
  
  # Extract the archive
  extract_ok <- FALSE
  try({
    utils::unzip(zipfile = eds_file, exdir = temp_dir)
    extract_ok <- TRUE
  }, silent = TRUE)
  if(!extract_ok){
    if(requireNamespace("archive", quietly = TRUE)){
      try({ archive::archive_extract(eds_file, dir = temp_dir); extract_ok <- TRUE }, silent = TRUE)
    }
  }
  if(!extract_ok){
    warning("Could not extract ", eds_file, " — skipping.")
    unlink(temp_dir, recursive = TRUE, force = TRUE)
    next
  }
  
  # Detect type by checking for specific files (recursively)
  all_extracted_files <- fs::dir_ls(temp_dir, recurse = TRUE, type = "file")
  has_txt <- any(grepl("(?i)analysis_result\\.txt$", fs::path_file(all_extracted_files), perl = TRUE))
  has_json <- any(grepl("(?i)relative_quantification_result\\.json$", fs::path_file(all_extracted_files), perl = TRUE))
  
  if(has_txt && has_json){
    warning("Both TXT and JSON files found in ", eds_file, " — ambiguous, skipping.")
    unlink(temp_dir, recursive = TRUE, force = TRUE)
    next
  } else if(!has_txt && !has_json){
    warning("Neither TXT nor JSON file found in ", eds_file, " — skipping.")
    unlink(temp_dir, recursive = TRUE, force = TRUE)
    next
  }
  
  # Process based on type and record file type
  df <- NULL
  file_type <- NULL
  if(has_txt){
    df <- process_txt_eds(eds_file, temp_dir)
    file_type <- "TXT (analysis_result.txt)"
  } else if(has_json){
    df <- process_json_eds(eds_file, temp_dir)
    file_type <- "JSON (relative_quantification_result.json)"
  }
  
  if(!is.null(df)){
    all_dfs <- append(all_dfs, list(df))
    eds_file_types[[basename(eds_file)]] <- file_type  # Store file type for this .EDS
  }
  
  # Cleanup temp
  unlink(temp_dir, recursive = TRUE, force = TRUE)
}

# ---- If no data processed, stop ----
# (Unchanged from original)
if(length(all_dfs) == 0){
  stop("No valid data extracted from any .EDS files.")
}

# ---- Stack all processed data frames ----
# (Unchanged from original)
stacked_data <- dplyr::bind_rows(all_dfs)

# ---- Save qLOOK_qBASE xlsx & txt (harmonize Quantity handling) ----
# (Unchanged from original)
stacked_for_txt <- stacked_data
if("Quantity" %in% names(stacked_for_txt)){
  stacked_for_txt$Quantity <- as.character(stacked_for_txt$Quantity)
  stacked_for_txt$Quantity[is.na(stacked_for_txt$Quantity)] <- ""
}
qbase_xlsx <- fs::path(folder_path, paste0("qLOOK_qBASE_", timestamp_global, ".xlsx"))
openxlsx::write.xlsx(stacked_data, qbase_xlsx, overwrite = TRUE)
created_files <- c(created_files, qbase_xlsx)

qbase_txt <- fs::path(folder_path, paste0("qLOOK_qBASE_", timestamp_global, ".txt"))
write.table(stacked_for_txt, file = qbase_txt, sep = "\t", row.names = FALSE, quote = FALSE)
created_files <- c(created_files, qbase_txt)

message("Saved stacked qLOOK_qBASE files: ", qbase_xlsx, " and ", qbase_txt)

# ---- Create wide qLOOK_Data (from stacked_data) ----
# (Unchanged from original)
stacked_data <- stacked_data %>%
  dplyr::mutate(Ct = suppressWarnings(as.numeric(Ct))) %>%
  dplyr::mutate(Ct = ifelse(is.na(Ct), 36, Ct))

agg <- stacked_data %>%
  dplyr::group_by(Sample, Gene) %>%
  dplyr::summarise(Ct = mean(Ct, na.rm = TRUE), .groups = "drop")

wide_data <- agg %>%
  tidyr::pivot_wider(names_from = Gene, values_from = Ct, names_repair = "unique") %>%
  dplyr::rename(Sample_Name = Sample)

qdata_file <- fs::path(folder_path, paste0("qLOOK_Data_", timestamp_global, ".xlsx"))
openxlsx::write.xlsx(wide_data, qdata_file, overwrite = TRUE)
created_files <- c(created_files, qdata_file)
message("Saved wide data file: ", qdata_file)

# ---- Run the reference gene analysis (unchanged from original) ----
# Read the data
ready_data <- openxlsx::read.xlsx(qdata_file, colNames = TRUE)

# Convert numeric-like columns (all except first column)
sample_col <- colnames(ready_data)[1]
gene_cols <- colnames(ready_data)[-1]

ready_data[, gene_cols] <- lapply(ready_data[, gene_cols, drop = FALSE], function(x) {
  as.numeric(as.character(x))
})

# Prepare ct_data
ct_data <- ready_data[, gene_cols, drop = FALSE]
ct_data <- ct_data[, colSums(!is.na(ct_data)) > 0, drop = FALSE]
rownames(ct_data) <- ready_data[[sample_col]]

# Delta Ct calculation
delta_ct_calc <- function(df){
  genes <- colnames(df)
  stability <- numeric(length(genes))
  for(i in seq_along(genes)){
    g1 <- df[[genes[i]]]
    others <- df[, -i, drop=FALSE]
    diff_matrix <- sweep(others, 1, g1, FUN = "-")
    sd_diff <- apply(diff_matrix, 2, sd, na.rm = TRUE)
    stability[i] <- mean(sd_diff, na.rm = TRUE)
  }
  result <- data.frame(Gene = genes, DeltaCtStability = stability)
  result <- result[order(result$DeltaCtStability), , drop = FALSE]
  return(result)
}
deltaCt_result <- delta_ct_calc(ct_data)

# geNorm calculation
geNorm_calc <- function(df){
  genes <- colnames(df)
  RQ <- 2^(-df)
  M_values <- numeric(length(genes)); names(M_values) <- genes
  remaining <- genes
  while(length(remaining) > 2){
    avg_sd <- sapply(remaining, function(g){
      others <- setdiff(remaining, g)
      log_ratios <- log2(sweep(RQ[, remaining, drop = FALSE], 1, RQ[, g], FUN = "/")[, others, drop = FALSE])
      mean(apply(log_ratios, 2, sd, na.rm = TRUE))
    })
    M_values[remaining] <- avg_sd
    remove_gene <- names(which.max(avg_sd))
    remaining <- setdiff(remaining, remove_gene)
  }
  final_sd <- sapply(remaining, function(g){
    others <- setdiff(remaining, g)
    log_ratios <- log2(sweep(RQ[, remaining, drop = FALSE], 1, RQ[, g], FUN = "/")[, others, drop = FALSE])
    mean(apply(log_ratios, 2, sd, na.rm = TRUE))
  })
  M_values[remaining] <- final_sd
  result <- data.frame(Gene = names(M_values), GeNormM = M_values)
  result <- result[order(result$GeNormM), , drop = FALSE]
  return(result)
}
geNorm_result <- geNorm_calc(ct_data)

# NormFinder calculation
normFinder_calc <- function(df){
  genes <- colnames(df)
  rho <- sapply(genes, function(g){
    sd(df[[g]], na.rm = TRUE)
  })
  result <- data.frame(Gene = genes, NormFinderRho = rho)
  result <- result[order(result$NormFinderRho), , drop = FALSE]
  return(result)
}
normFinder_result <- normFinder_calc(ct_data)

# Save results
output_file <- file.path(folder_path, paste0("qLOOK_RefGene_", timestamp_global, ".xlsx"))

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "DeltaCt")
openxlsx::writeData(wb, "DeltaCt", deltaCt_result)
openxlsx::addWorksheet(wb, "geNorm")
openxlsx::writeData(wb, "geNorm", geNorm_result)
openxlsx::addWorksheet(wb, "NormFinder")
openxlsx::writeData(wb, "NormFinder", normFinder_result)
openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)
created_files <- c(created_files, output_file)

message("✅ Gene stability results saved: ", output_file)

# ---- Logging (modified to include file type) ----
log_file <- file.path(folder_path, "qLOOK_Summary.txt")

# Original .EDS files with their types
eds_files_basenames <- basename(eds_files)
# Create a list of processed files with their types (only those successfully processed)
processed_eds <- names(eds_file_types)
eds_log_entries <- sapply(eds_files_basenames, function(f) {
  if(f %in% processed_eds) {
    paste(f, "processed as", eds_file_types[[f]])
  } else {
    paste(f, "skipped (no valid TXT or JSON found)")
  }
})

# Script run timestamp
script_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

# Prepare the log text
created_files_text <- paste("-", basename(created_files), collapse = "\n")
eds_files_text <- paste("-", eds_log_entries, collapse = "\n")

log_text <- paste0(
  "===== qLOOK - qPCR-LOg-boOK =====","\n",
  "===== EDS Processing Log =====", "\n",
  "\n--- ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ---\n",
  "Script: qLOOK_Module1_v1.0.R (v1.0)\n",
  "Script: qLOOK_RefGene_v1.R (v1)\n",
  "Original .EDS files processed:\n", eds_files_text, "\n", "\n",
  "If .TXT -> EDS was generated by ViiA7, StepOneplus, QuantStudio6, ABI7500, or ABI7500fast", "\n",
  "If .JSON -> EDS was generated by QuantStudio3", "\n", "\n",
  "Created files:\n", created_files_text, 
  "\n"
)

# Write or append to qLOOK_Summary.txt
if(file.exists(log_file)){
  existing_text <- readLines(log_file)
  new_text <- c(existing_text, "", "", log_text)
} else {
  new_text <- log_text
}

writeLines(new_text, log_file)
created_files <- c(created_files, log_file)

message("✅ Log updated: ", log_file)
message("\n🎉 All steps completed successfully!")