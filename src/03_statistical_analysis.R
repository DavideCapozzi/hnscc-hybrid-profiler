# src/03_statistical_analysis.R
# ==============================================================================
# STEP 03: STATISTICAL ANALYSIS & REPORTING
# Description: Extracts hybrid drivers of clinical response via sPLS-DA.
# ==============================================================================

source("R/utils_io.R")          
source("R/modules_multivariate.R")
source("R/modules_viz.R") 

message("\n=== PIPELINE STEP 3: STATISTICAL ANALYSIS & REPORTING ===")

# 1. INITIALIZATION & DATA SETUP
args <- commandArgs(trailingOnly = TRUE)
# 1. Load Config & Data
if (!exists("config")) {
  args <- commandArgs(trailingOnly = TRUE)
  config_path <- if (length(args) > 0) args[1] else "config/global_params.yml"
  config <- load_config(config_path)
}
input_file <- file.path(config$output_root, "01_data_processing", "data_processed.rds")

if (!file.exists(input_file)) stop("[FATAL] Step 01 output not found.")
DATA <- readRDS(input_file)

results_dir <- file.path(config$output_root, "03_statistical_analysis")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

df_global <- DATA$hybrid_data_z
safe_markers <- DATA$hybrid_markers
meta_stats <- DATA$metadata 

# Extract strict matrix with exact rownames to prevent downstream heatmap crashes
mat_z_global <- as.matrix(df_global[, safe_markers])
rownames(mat_z_global) <- meta_stats$Patient_ID

# Align colors safely (protecting against missing config keys)
colors_viz <- c()
if (!is.null(config$colors$groups[[config$clinical$responder_label]])) {
  colors_viz[config$clinical$responder_label] <- config$colors$groups[[config$clinical$responder_label]]
} else { colors_viz[config$clinical$responder_label] <- "blue" }

if (!is.null(config$colors$groups[[config$clinical$non_responder_label]])) {
  colors_viz[config$clinical$non_responder_label] <- config$colors$groups[[config$clinical$non_responder_label]]
} else { colors_viz[config$clinical$non_responder_label] <- "red" }

wb_master <- openxlsx::createWorkbook()

# 2. MULTIVARIATE HYBRID sPLS-DA
message("\n--- RUNNING HYBRID sPLS-DA (Responder vs Non-Responder) ---")

set.seed(config$stats$seed) 
tryCatch({
  pls_res <- run_splsda_model(
    data_z = mat_z_global, 
    metadata = meta_stats, 
    group_col = "Group", 
    n_comp = config$multivariate$n_comp,
    folds = config$multivariate$validation_folds, 
    n_repeat = config$multivariate$n_repeat_cv
  ) 
  
  top_drivers <- extract_plsda_loadings(pls_res)
  
  if (nrow(top_drivers) > 0) {
    # Flag marker origin (FACS vs Soluble) for interpretability
    facs_feats <- unlist(config$features$facs)
    top_drivers <- top_drivers %>%
      dplyr::mutate(Domain = ifelse(Marker %in% facs_feats, "FACS", "Soluble")) %>%
      dplyr::select(Domain, Marker, Importance, dplyr::everything())
    
    openxlsx::addWorksheet(wb_master, "sPLSDA_Drivers")
    openxlsx::writeData(wb_master, "sPLSDA_Drivers", top_drivers)
    
    # Safe extraction of n_top_boxplots to prevent NULL crash
    safe_top_n <- if(!is.null(config$multivariate$top_n_loadings)) config$multivariate$top_n_loadings else 15
    
    # Generate Plots
    viz_report_plsda(
      pls_res = pls_res, 
      drivers_df = top_drivers, 
      metadata_viz = meta_stats, 
      colors_viz = colors_viz, 
      out_path = file.path(results_dir, "Hybrid_sPLSDA_Results.pdf"),
      group_col = "Group", 
      n_top_boxplots = safe_top_n
    )
  }
}, error = function(e) message(paste("   [ERROR] sPLS-DA Failed:", e$message)))

# 3. EXPORT FINAL REPORT
report_path <- file.path(results_dir, "Hybrid_Statistical_Report.xlsx")
openxlsx::saveWorkbook(wb_master, report_path, overwrite = TRUE)
message(sprintf("\n   [Output] Final Statistical Report saved: %s", report_path))

message("=== STEP 3 COMPLETE ===\n")