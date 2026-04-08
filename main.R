# main.R
# ==============================================================================
# MAIN PIPELINE ORCHESTRATOR: onco-hybrid-profiler
# Description: End-to-end execution for Hybrid Data (FACS + Solubles)
#              Goal: Identify drivers of clinical response via sPLS-DA & LMM.
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
if (!file.exists(base_config_path)) stop("[FATAL] Base Config file not found!")
base_config <- yaml::read_yaml(base_config_path)

if (is.null(base_config$experiments) || length(base_config$experiments) == 0) {
  stop("[FATAL] No experiments defined in the base config file!")
}

# --- Module Loading ---
message("\n>>> LOADING MODULES <<<")
list.files(here("R"), pattern = "\\.R$", full.names = TRUE) %>% purrr::walk(source)
message("[System] Modules loaded successfully.")

# 2. Pipeline Execution Loop
# ------------------------------------------------------------------------------
for (exp_name in names(base_config$experiments)) {
  
  exp_cfg <- base_config$experiments[[exp_name]]
  
  run_standard <- if (!is.null(exp_cfg$run_standard)) as.logical(exp_cfg$run_standard) else TRUE
  run_longitudinal <- if (!is.null(exp_cfg$is_longitudinal)) as.logical(exp_cfg$is_longitudinal) else FALSE
  
  message(sprintf("\n========================================================"))
  message(sprintf("STARTING EXPERIMENT BLOCK: %s", exp_name))
  message(sprintf("========================================================\n"))
  
  # ============================================================================
  # PASS 1: STANDARD CROSS-SECTIONAL PIPELINE (01, 02, 03)
  # ============================================================================
  if (run_standard) {
    config <- base_config
    config$project_name <- exp_name # Shared root directory
    config$run_mode <- "standard"   # Filename isolation
    config$clinical <- exp_cfg$clinical
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
    config <- base_config
    config$project_name <- exp_name     # Shared root directory
    config$run_mode <- "longitudinal"   # Filename isolation
    config$clinical <- exp_cfg$clinical
    config$is_longitudinal <- TRUE
    
    # Critical: Disable Multivariate Outliers to prevent temporal censoring
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