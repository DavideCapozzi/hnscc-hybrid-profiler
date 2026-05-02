# main.R
# ==============================================================================
# MAIN PIPELINE ORCHESTRATOR: clinical-onco-profiler
# Description: End-to-end execution for Hybrid Data (FACS + Solubles)
#              Goal: Identify drivers of clinical response via sPLS-DA & LMM.
#              Supports nested multi-experiment and flat single-cohort configs.
# ==============================================================================

# 1. Environment Setup
rm(list = ls())
graphics.off()

suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(here)
  library(openxlsx)
})

options(crayon.enabled = FALSE)

# Load BASE Configuration 
base_config_path <- here("config/global_params.yml")
if (!file.exists(base_config_path)) stop("[FATAL] Base Config file not found at expected path.")
base_config <- yaml::read_yaml(base_config_path)

# --- Configuration Router (Nested vs Flat Support) ---
experiments_list <- list()

if (!is.null(base_config$experiments) && length(base_config$experiments) > 0) {
  # Mode: Nested Multi-Experiment
  experiments_list <- base_config$experiments
  message("[System] Detected nested multi-experiment configuration.")
} else {
  # Mode: Flat Single-Cohort Configuration (Fallback)
  exp_name <- if (!is.null(base_config$project_name)) base_config$project_name else "Single_Cohort_Run"
  experiments_list[[exp_name]] <- base_config
  message("[System] Detected flat single-cohort configuration.")
}

# --- Module Loading ---
message("\n>>> LOADING MODULES <<<")
list.files(here("R"), pattern = "\\.R$", full.names = TRUE) %>% purrr::walk(source)
message("[System] Modules loaded successfully.")

# 2. Pipeline Execution Loop
# ------------------------------------------------------------------------------
for (exp_name in names(experiments_list)) {
  
  exp_cfg <- experiments_list[[exp_name]]
  
  # Inherit boolean flags with safety fallbacks
  run_standard     <- if (!is.null(exp_cfg$run_standard))    as.logical(exp_cfg$run_standard)    else TRUE
  run_longitudinal <- if (!is.null(exp_cfg$is_longitudinal)) as.logical(exp_cfg$is_longitudinal) else FALSE
  run_network      <- if (!is.null(exp_cfg$run_network))     as.logical(exp_cfg$run_network)     else FALSE
  
  message(sprintf("\n========================================================"))
  message(sprintf("STARTING EXPERIMENT BLOCK: %s", exp_name))
  message(sprintf("========================================================\n"))
  
  # ============================================================================
  # PASS 1: STANDARD CROSS-SECTIONAL PIPELINE (01, 02, 03[, 05])
  # Note: step 05 (network) runs here, before pass 2 step 04 (longitudinal).
  #       The two steps belong to independent analytical passes and different data modalities.
  # ============================================================================
  if (run_standard) {
    # Isolate configuration state for this pass
    config <- base_config
    config$project_name <- exp_name 
    config$run_mode <- "standard"   
    
    # Overwrite clinical logic if defined specifically in the experiment block
    if (!is.null(exp_cfg$clinical)) config$clinical <- exp_cfg$clinical
    
    config$is_longitudinal <- FALSE
    
    if (!is.null(exp_cfg$input_file)) config$input_file <- exp_cfg$input_file
    
    config$output_root <- file.path(base_config$output_root, config$project_name)
    if (!dir.exists(config$output_root)) dir.create(config$output_root, recursive = TRUE)
    
    log_file_std <- file.path(config$output_root, paste0("log_standard_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
    
    message(sprintf("\n--- [PASS 1] STANDARD PIPELINE: %s ---", config$project_name))
    cat(sprintf("=== STANDARD PIPELINE STARTED: %s ===\n", Sys.time()), file = log_file_std)
    
    tryCatch({
      withCallingHandlers({
        message(">>> RUNNING PHASE 1: DATA PROCESSING <<<")
        source(here("src/01_data_processing.R"), echo = FALSE, local = FALSE)

        message("\n>>> RUNNING PHASE 2: VISUALIZATION <<<")
        source(here("src/02_visualization.R"), echo = FALSE, local = FALSE)

        message("\n>>> RUNNING PHASE 3: STATISTICAL ANALYSIS & REPORTING <<<")
        source(here("src/03_statistical_analysis.R"), echo = FALSE, local = FALSE)

        if (run_network) {
          message("\n>>> RUNNING PHASE 5: DIFFERENTIAL NETWORK ANALYSIS <<<")
          source(here("src/05_network_analysis.R"), echo = FALSE, local = FALSE)
        }

        message(sprintf("\n[SUCCESS] Standard Pipeline completed for: %s", exp_name))
      }, message = function(m) cat(conditionMessage(m), file = log_file_std, append = TRUE, sep = "\n"),
      warning = function(w) cat(paste0("WARNING: ", conditionMessage(w)), file = log_file_std, append = TRUE, sep = "\n"))
    }, error = function(e) {
      err_msg <- paste0("\n[FATAL ERROR] Standard Pipeline stopped: ", e$message)
      cat(err_msg, file = log_file_std, append = TRUE, sep = "\n")
      message(err_msg)
    })
  }
  
  # ============================================================================
  # PASS 2: LONGITUDINAL PIPELINE (01 Joint, 04 LMM)
  # ============================================================================
  if (run_longitudinal) {
    # Isolate configuration state for this pass
    config <- base_config
    config$project_name <- exp_name     
    config$run_mode <- "longitudinal"   
    
    if (!is.null(exp_cfg$clinical)) config$clinical <- exp_cfg$clinical
    
    config$is_longitudinal <- TRUE
    
    # Critical Structural Rule: Disable Multivariate Outliers to prevent temporal censoring
    config$qc$remove_outliers <- FALSE
    
    if (!is.null(exp_cfg$input_file_t0)) config$input_file_t0 <- exp_cfg$input_file_t0
    if (!is.null(exp_cfg$input_file_t1)) config$input_file_t1 <- exp_cfg$input_file_t1
    
    config$output_root <- file.path(base_config$output_root, config$project_name)
    if (!dir.exists(config$output_root)) dir.create(config$output_root, recursive = TRUE)
    
    log_file_long <- file.path(config$output_root, paste0("log_longitudinal_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
    
    message(sprintf("\n--- [PASS 2] LONGITUDINAL PIPELINE: %s ---", config$project_name))
    cat(sprintf("=== LONGITUDINAL PIPELINE STARTED: %s ===\n", Sys.time()), file = log_file_long)
    
    tryCatch({
      withCallingHandlers({
        message(">>> RUNNING PHASE 1: JOINT DATA PROCESSING <<<")
        source(here("src/01_data_processing.R"), echo = FALSE, local = FALSE)
        
        message("\n>>> RUNNING PHASE 4: LONGITUDINAL LMM ANALYSIS <<<")
        source(here("src/04_longitudinal_lmm.R"), echo = FALSE, local = FALSE)
        
        message(sprintf("\n[SUCCESS] Longitudinal Pipeline completed for: %s", exp_name))
      }, message = function(m) cat(conditionMessage(m), file = log_file_long, append = TRUE, sep = "\n"),
      warning = function(w) cat(paste0("WARNING: ", conditionMessage(w)), file = log_file_long, append = TRUE, sep = "\n"))
    }, error = function(e) {
      err_msg <- paste0("\n[FATAL ERROR] Longitudinal Pipeline stopped: ", e$message)
      cat(err_msg, file = log_file_long, append = TRUE, sep = "\n")
      message(err_msg)
    })
  }
}

options(crayon.enabled = TRUE)
message("\n=== ALL EXPERIMENT BLOCKS COMPLETED ===")