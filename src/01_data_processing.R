# src/01_data_processing.R
# ==============================================================================
# STEP 01: DATA PROCESSING
# Description: Data ingestion, Hybrid QC, and Hybrid Transformation.
# ==============================================================================

source("R/utils_io.R")     
source("R/modules_coda.R") 
source("R/modules_qc.R")   

message("\n=== PIPELINE STEP 1: INGESTION + HYBRID QC + TRANSFORM ===")

# 1. Load Configuration & Data
# ------------------------------------------------------------------------------
config <- load_config("config/global_params.yml")

# Use basic readxl for single file (ignoring old multi-sheet logic)
input_file <- config$input_file
if (!file.exists(input_file)) stop(sprintf("Input file not found: %s", input_file))
raw_data <- readxl::read_excel(input_file)

# Standardize Patient ID column (assume it's the first column)
colnames(raw_data)[1] <- "Patient_ID"

# 2. Setup Initial Matrices & Clinical Filtering
# ------------------------------------------------------------------------------
clin_col <- config$clinical$target_column
resp_lbl <- config$clinical$responder_label
nresp_lbl <- config$clinical$non_responder_label

if (!clin_col %in% colnames(raw_data)) {
  stop(sprintf("[FATAL] Clinical target column '%s' not found in dataset.", clin_col))
}

# Filter to keep only target clinical labels and rename to 'Group'
raw_filtered <- raw_data %>%
  dplyr::filter(.data[[clin_col]] %in% c(resp_lbl, nresp_lbl)) %>%
  dplyr::mutate(Group = factor(.data[[clin_col]], levels = c(nresp_lbl, resp_lbl)))

if (nrow(raw_filtered) < 5) stop("[FATAL] Insufficient valid patients found for specified clinical labels.")

# Extract Features
facs_feats <- unlist(config$features$facs)
sol_feats <- unlist(config$features$soluble)
target_markers <- unique(c(facs_feats, sol_feats))

# Ensure markers exist in data
missing_markers <- setdiff(target_markers, colnames(raw_filtered))
if (length(missing_markers) > 0) {
  warning(sprintf("[Data] %d configured markers not found in Excel. They will be ignored.", length(missing_markers)))
  target_markers <- intersect(target_markers, colnames(raw_filtered))
  facs_feats <- intersect(facs_feats, colnames(raw_filtered))
  sol_feats <- intersect(sol_feats, colnames(raw_filtered))
}

# Build strictly numeric matrix
mat_raw <- as.matrix(raw_filtered[, target_markers])
rownames(mat_raw) <- raw_filtered$Patient_ID
metadata_raw <- raw_filtered %>% dplyr::select(Patient_ID, Group)

message(sprintf("[Data] Initial Matrix: %d Samples x %d Hybrid Markers", nrow(mat_raw), ncol(mat_raw)))

# 3. Hybrid Quality Control (Dual Thresholds)
# ------------------------------------------------------------------------------
message("[QC] Running Split QC Strategy...")

qc_config_basic <- config$qc
qc_config_basic$remove_outliers <- FALSE # Disable outlier detection for first pass

message("   [QC-A] Filtering Missingness with Dual Thresholds...")
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

# --- STEP B: Multivariate Outlier Detection ---
if (config$qc$remove_outliers) {
  message("   [QC-B] Generating hybrid proxy for Outlier Detection...")
  
  proxy_results <- perform_data_transformation(mat_clean_basic, config, mode = "fast")
  mat_proxy <- proxy_results$hybrid_data_z
  
  is_outlier <- detect_pca_outliers(
    mat = mat_proxy, 
    groups = meta_clean_basic$Group, 
    conf_level = config$qc$outlier_conf_level
  )
  
  if (any(is_outlier)) {
    out_pids <- rownames(mat_clean_basic)[is_outlier]
    message(sprintf("   [QC-B] Dropping %d multivariate outliers.", length(out_pids)))
    
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
message(sprintf("[QC] Final Clean Dimensions: %d Samples x %d Markers", nrow(mat_raw), ncol(mat_raw)))

# 4. Final Hybrid Transformation (Complete Mode with BPCA)
# ------------------------------------------------------------------------------
message("[Transform] Running Final Hybrid Data Transformation...")
transform_results <- perform_data_transformation(mat_raw, config, mode = "complete")

mat_hybrid_raw <- transform_results$hybrid_data_raw
mat_hybrid_z   <- transform_results$hybrid_data_z
final_markers  <- transform_results$hybrid_markers

# 5. Save Final Output
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

saveRDS(processed_data, file.path(out_dir, "data_processed.rds"))
save_qc_report(qc_result$report, file.path(out_dir, "QC_Filtering_Report.xlsx"))

message("=== STEP 1 COMPLETE ===\n")
