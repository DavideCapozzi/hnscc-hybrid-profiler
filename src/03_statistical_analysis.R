# src/03_statistical_analysis.R
# ==============================================================================
# STEP 03: STATISTICAL ANALYSIS & REPORTING
# Description: Extracts hybrid drivers of clinical response via Sparse PLS-DA.
# ==============================================================================

source("R/utils_io.R")          
source("R/modules_multivariate.R")
source("R/modules_viz.R") 

message("\n=== PIPELINE STEP 3: STATISTICAL ANALYSIS & REPORTING ===")

# 1. Initialization & Data Loading
# ------------------------------------------------------------------------------
if (!exists("config")) {
  args <- commandArgs(trailingOnly = TRUE)
  config_path <- if (length(args) > 0) args[1] else "config/global_params.yml"
  config <- load_config(config_path)
}

input_file <- file.path(config$output_root, "01_data_processing", sprintf("data_processed_%s_standard.rds", config$project_name))
if (!file.exists(input_file)) stop("[FATAL] Standard data processing output not found.")

DATA <- readRDS(input_file)
results_dir <- file.path(config$output_root, "03_statistical_analysis")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

df_global <- DATA$hybrid_data_z
safe_markers <- DATA$hybrid_markers
meta_stats <- DATA$metadata 

# Isolate strict matrix to prevent downstream ComplexHeatmap dimensional mismatches
mat_z_global <- as.matrix(df_global[, safe_markers])

# Structural requirement: Enforce unique rownames for mixOmics C++ backend compatibility
rownames(mat_z_global) <- make.unique(as.character(meta_stats$Patient_ID))

colors_viz <- get_clinical_colors(config)
wb_master <- openxlsx::createWorkbook()

# 2. Multivariate Modeling (Hybrid sPLS-DA)
# ------------------------------------------------------------------------------
message(sprintf("\n--- RUNNING HYBRID sPLS-DA (%s vs %s) ---", config$clinical$responder_label, config$clinical$non_responder_label))

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
    # Flag marker structural domain (FACS vs Soluble) for biological interpretability
    facs_feats <- unlist(config$features$facs)
    top_drivers <- top_drivers %>%
      dplyr::mutate(Domain = ifelse(Marker %in% facs_feats, "FACS", "Soluble")) %>%
      dplyr::select(Domain, Marker, Importance, dplyr::everything())
    
    openxlsx::addWorksheet(wb_master, "sPLSDA_Drivers")
    openxlsx::writeData(wb_master, "sPLSDA_Drivers", top_drivers)
    
    df_perf <- extract_plsda_performance(pls_res)
    
    auc_list <- list()
    if (!is.null(pls_res$performance$auc)) {
      for (i in 1:pls_res$model$ncomp) {
        comp_name <- paste0("comp", i)
        if (comp_name %in% names(pls_res$performance$auc)) {
          auc_val <- tryCatch(pls_res$performance$auc[[comp_name]][1, "AUC"], error = function(e) NA)
          pval_val <- tryCatch(pls_res$performance$auc[[comp_name]][1, "p-value"], error = function(e) NA)
          auc_list[[i]] <- data.frame(Component = i, AUC = auc_val, P_Value = pval_val)
        }
      }
    }
    
    if (length(auc_list) > 0) {
      df_auc <- do.call(rbind, auc_list)
      df_perf <- merge(df_perf, df_auc, by = "Component", all.x = TRUE)
    }
    
    df_perf$keepX <- pls_res$model$keepX[1:pls_res$model$ncomp]
    
    openxlsx::addWorksheet(wb_master, "Model_Performance")
    openxlsx::writeData(wb_master, "Model_Performance", df_perf)
    
    # Generate machine-readable JSON object for downstream parsing
    machine_output <- list(
      project_name = config$project_name,
      clinical_target = config$clinical$target_column,
      model_type = "sPLS-DA",
      validation_method = df_perf$Validation_Method[1],
      cv_folds = config$multivariate$validation_folds,
      cv_repeats = config$multivariate$n_repeat_cv,
      performance_per_component = df_perf,
      top_drivers = top_drivers
    )
    
    json_path <- file.path(results_dir, sprintf("Machine_Metrics_%s.json", config$project_name))
    if (requireNamespace("jsonlite", quietly = TRUE)) {
      jsonlite::write_json(machine_output, json_path, pretty = TRUE, auto_unbox = TRUE)
      message(sprintf("   [Output] Serialized JSON metrics saved: %s", basename(json_path)))
    } else {
      saveRDS(machine_output, sub(".json$", ".rds", json_path))
      warning("[Dependencies] jsonlite package missing. Falling back to RDS format.")
    }
    
    safe_top_n <- if(!is.null(config$multivariate$top_n_loadings)) config$multivariate$top_n_loadings else 15
    
    # Render Analytical Report
    viz_report_plsda(
      pls_res = pls_res, 
      drivers_df = top_drivers, 
      metadata_viz = meta_stats, 
      colors_viz = colors_viz, 
      out_path = file.path(results_dir, sprintf("Hybrid_sPLSDA_Results_%s.pdf", config$project_name)),
      group_col = "Group", 
      n_top_boxplots = safe_top_n
    )
  }
}, error = function(e) message(paste("   [ERROR] sPLS-DA computation failed:", e$message)))

# 3. Export Global Statistics
# ------------------------------------------------------------------------------
report_path <- file.path(results_dir, sprintf("Hybrid_Statistical_Report_%s.xlsx", config$project_name))
openxlsx::saveWorkbook(wb_master, report_path, overwrite = TRUE)
message(sprintf("\n   [Output] Mathematical summary workbook saved: %s", report_path))

message("=== STEP 3 COMPLETE ===\n")