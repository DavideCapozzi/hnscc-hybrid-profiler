# src/01_data_processing.R
# ==============================================================================
# STEP 01: DATA PROCESSING
# Description: Data ingestion, Hybrid Quality Control, and Transform logic.
# ==============================================================================

source("R/utils_io.R")     
source("R/modules_coda.R") 
source("R/modules_qc.R")   

message("\n=== PIPELINE STEP 1: INGESTION + HYBRID QC + TRANSFORM ===")

# 1. Load Configuration & Data Initialization
# ------------------------------------------------------------------------------
if (!exists("config")) {
  args <- commandArgs(trailingOnly = TRUE)
  config_path <- if (length(args) > 0) args[1] else "config/global_params.yml"
  config <- load_config(config_path)
  if (is.null(config$run_mode)) config$run_mode <- "standard"
}

# Handle Longitudinal vs Cross-Sectional Data Ingestion
if (isTRUE(config$is_longitudinal)) {
  if (is.null(config$input_file_t0) || is.null(config$input_file_t1)) {
    stop("[FATAL] Configuration keys 'input_file_t0' and 'input_file_t1' are required for longitudinal mode.")
  }
  if (!file.exists(config$input_file_t0)) stop(sprintf("T0 Input file not found: %s", config$input_file_t0))
  if (!file.exists(config$input_file_t1)) stop(sprintf("T1 Input file not found: %s", config$input_file_t1))
  
  raw_t0 <- readxl::read_excel(config$input_file_t0) %>% dplyr::mutate(Timepoint = "T0")
  raw_t1 <- readxl::read_excel(config$input_file_t1) %>% dplyr::mutate(Timepoint = "T1")
  
  colnames(raw_t0)[1] <- "Patient_ID"
  colnames(raw_t1)[1] <- "Patient_ID"
  
  raw_data <- dplyr::bind_rows(raw_t0, raw_t1)
  raw_data$Sample_ID <- paste(as.character(raw_data$Patient_ID), raw_data$Timepoint, sep = "_")
  
} else {
  input_file <- config$input_file
  if (!file.exists(input_file)) stop(sprintf("Input file not found: %s", input_file))
  
  raw_data <- readxl::read_excel(input_file)
  colnames(raw_data)[1] <- "Patient_ID"
  raw_data$Timepoint <- "Baseline"
  
  if (any(duplicated(raw_data$Patient_ID))) {
    warning("[Data] Duplicated Patient_IDs detected. Applying make.unique() standardization.")
    raw_data$Patient_ID <- make.unique(as.character(raw_data$Patient_ID))
  }
  raw_data$Sample_ID <- as.character(raw_data$Patient_ID)
}

# 2. Setup Initial Matrices & Clinical Target Filtering
# ------------------------------------------------------------------------------
clin_col <- config$clinical$target_column
resp_lbl <- config$clinical$responder_label
nresp_lbl <- config$clinical$non_responder_label

if (!clin_col %in% colnames(raw_data)) {
  stop(sprintf("[FATAL] Target clinical column '%s' not found in the dataset.", clin_col))
}

# Apply mappings if defined in configuration
if (!is.null(config$clinical$mapping)) {
  mapped_resp <- as.character(unlist(config$clinical$mapping[[resp_lbl]]))
  mapped_nresp <- as.character(unlist(config$clinical$mapping[[nresp_lbl]]))
  
  raw_data <- raw_data %>%
    dplyr::mutate(!!clin_col := dplyr::case_when(
      as.character(.data[[clin_col]]) %in% mapped_resp ~ resp_lbl,
      as.character(.data[[clin_col]]) %in% mapped_nresp ~ nresp_lbl,
      TRUE ~ as.character(.data[[clin_col]])
    ))
}

raw_filtered <- raw_data %>%
  dplyr::filter(.data[[clin_col]] %in% c(resp_lbl, nresp_lbl)) %>%
  dplyr::mutate(Group = factor(.data[[clin_col]], levels = c(nresp_lbl, resp_lbl)))

if (nrow(raw_filtered) < 5) stop("[FATAL] Insufficient patient observations for specified clinical labels.")

# Dynamic Covariate Extraction
meta_cols <- c("Patient_ID", "Sample_ID", "Timepoint", "Group")
if (!is.null(config$clinical$covariates)) {
  cov_cols <- unlist(config$clinical$covariates)
  valid_covs <- intersect(cov_cols, colnames(raw_filtered))
  if (length(valid_covs) < length(cov_cols)) {
    warning("[Data] Discrepancy detected: Some configured covariates are missing from the raw dataset.")
  }
  meta_cols <- unique(c(meta_cols, valid_covs))
}

metadata_raw <- raw_filtered %>% dplyr::select(dplyr::all_of(meta_cols))

facs_feats <- unlist(config$features$facs)
sol_feats <- unlist(config$features$soluble)
target_markers <- unique(c(facs_feats, sol_feats))

missing_markers <- setdiff(target_markers, colnames(raw_filtered))
if (length(missing_markers) > 0) {
  warning(sprintf("[Data] %d configured markers not found in the input matrix. Proceeding with intersection.", length(missing_markers)))
  target_markers <- intersect(target_markers, colnames(raw_filtered))
  facs_feats <- intersect(facs_feats, colnames(raw_filtered))
  sol_feats <- intersect(sol_feats, colnames(raw_filtered))
}

mat_raw <- as.matrix(raw_filtered[, target_markers])
rownames(mat_raw) <- raw_filtered$Sample_ID

message(sprintf("[Data] Matrix initialized: %d Samples x %d Hybrid Markers", nrow(mat_raw), ncol(mat_raw)))

# 3. Hybrid Quality Control
# ------------------------------------------------------------------------------
message("[QC] Executing Split Quality Control Strategy...")

qc_config_basic <- config$qc
qc_config_basic$remove_outliers <- FALSE

message("   [QC-A] Enforcing Missingness Tolerances (Dual Domain Thresholds)...")
qc_result <- run_qc_pipeline(
  mat_raw = mat_raw, 
  metadata = metadata_raw, 
  qc_config = qc_config_basic,
  facs_features = facs_feats,
  soluble_features = sol_feats,
  stratification_col = "Group"
)

mat_clean_basic  <- qc_result$data
meta_clean_basic <- qc_result$metadata

if (isTRUE(config$qc$remove_outliers)) {
  message("   [QC-B] Generating hybrid PCA proxy for Multivariate Outlier Detection...")
  
  proxy_results <- perform_data_transformation(mat_clean_basic, config, mode = "fast")
  mat_proxy <- proxy_results$hybrid_data_z
  
  is_outlier <- detect_pca_outliers(
    mat = mat_proxy, 
    groups = meta_clean_basic$Group, 
    conf_level = config$qc$outlier_conf_level
  )
  
  if (any(is_outlier)) {
    out_pids <- rownames(mat_clean_basic)[is_outlier]
    message(sprintf("   [QC-B] Removing %d multivariate outliers based on robust Mahalanobis distance.", length(out_pids)))
    
    group_info <- as.character(meta_clean_basic$Group[is_outlier])
    na_pcts <- rowMeans(is.na(mat_clean_basic[is_outlier, , drop=FALSE])) * 100
    
    qc_result$report$dropped_rows_detail <- rbind(qc_result$report$dropped_rows_detail, data.frame(
      Patient_ID = out_pids, NA_Percent = round(na_pcts, 2),
      Reason = "Outlier", Original_Source = group_info, stringsAsFactors = FALSE
    ))
    qc_result$report$n_row_dropped <- nrow(qc_result$report$dropped_rows_detail)
    
    mat_clean_basic  <- mat_clean_basic[!is_outlier, , drop = FALSE]
    meta_clean_basic <- meta_clean_basic[!is_outlier, , drop = FALSE]
  }
}

mat_raw <- mat_clean_basic
metadata_raw <- meta_clean_basic
message(sprintf("[QC] Final Matrix Dimensions: %d Samples x %d Markers", nrow(mat_raw), ncol(mat_raw)))

# 4. Final Hybrid Transformation (Mode-Aware: Prevents BPCA Leakage)
# ------------------------------------------------------------------------------
transform_mode <- if (isTRUE(config$is_longitudinal)) "fast" else "complete"
message(sprintf("[Transform] Executing Final Hybrid Transformation (Mode: %s)...", transform_mode))

transform_results <- perform_data_transformation(mat_raw, config, mode = transform_mode)

mat_hybrid_raw <- transform_results$hybrid_data_raw
mat_hybrid_z   <- transform_results$hybrid_data_z
final_markers  <- transform_results$hybrid_markers

# 5. Export Final Artifacts
# ------------------------------------------------------------------------------
out_dir <- file.path(config$output_root, "01_data_processing")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

df_hybrid_raw <- cbind(metadata_raw, as.data.frame(mat_hybrid_raw))
df_hybrid_z   <- cbind(metadata_raw, as.data.frame(mat_hybrid_z))

processed_data <- list(
  metadata        = metadata_raw,
  markers         = colnames(mat_raw),      
  raw_matrix      = mat_raw, 
  hybrid_markers  = final_markers,
  hybrid_data_raw = df_hybrid_raw,
  hybrid_data_z   = df_hybrid_z,
  config          = config
)

out_file_data <- file.path(out_dir, sprintf("data_processed_%s_%s.rds", config$project_name, config$run_mode))
out_file_qc   <- file.path(out_dir, sprintf("QC_Filtering_Report_%s_%s.xlsx", config$project_name, config$run_mode))

saveRDS(processed_data, out_file_data)
save_qc_report(qc_result$report, out_file_qc)

message(sprintf("   [Output] Serialized data objects saved: %s", basename(out_file_data)))
message("=== STEP 1 COMPLETE ===\n")