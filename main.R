# main.R
# ==============================================================================
# MAIN PIPELINE ORCHESTRATOR: onco-hybrid-profiler
# Description: End-to-end execution for Hybrid Data (FACS + Solubles)
#              Goal: Identify drivers of clinical response via sPLS-DA.
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

# Load BASE Configuration (contains features, qc thresholds, colors, experiments, etc.)
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
  
  # A. Merge Configs (Base + Experiment specific)
  exp_cfg <- base_config$experiments[[exp_name]]
  config <- base_config # This creates the global config variable used by 01, 02, 03
  
  config$project_name <- exp_name
  config$clinical <- exp_cfg$clinical
  if (!is.null(exp_cfg$input_file)) config$input_file <- exp_cfg$input_file
  
  # Construct nested output directory (e.g., results/Toxicity_1v0)
  if (!is.null(config$project_name) && config$project_name != "") {
    config$output_root <- file.path(base_config$output_root, config$project_name)
  }
  
  # B. Logging Setup for this specific experiment
  out_root <- here(config$output_root)
  if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE)
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  log_file <- file.path(out_root, paste0("pipeline_log_", timestamp, ".txt"))
  
  message(sprintf("\n========================================================"))
  message(sprintf("STARTING EXPERIMENT: %s", config$project_name))
  message(sprintf("Log file: %s", log_file))
  message(sprintf("========================================================\n"))
  
  cat(sprintf("=== ONCO-HYBRID-PROFILER STARTED: %s ===\n", Sys.time()), file = log_file)
  cat(sprintf("Experiment: %s\n", config$project_name), file = log_file, append = TRUE)
  
  log_handler <- function(m) {
    cat(conditionMessage(m), file = log_file, append = TRUE, sep = "\n")
  }
  
  warn_handler <- function(w) {
    cat(paste0("WARNING: ", conditionMessage(w)), file = log_file, append = TRUE, sep = "\n")
  }
  
  # C. Run Scripts
  tryCatch({
    withCallingHandlers({
      
      message("\n>>> RUNNING PHASE 1: DATA PROCESSING <<<")
      source(here("src/01_data_processing.R"), echo = FALSE, local = FALSE)
      
      message("\n>>> RUNNING PHASE 2: VISUALIZATION <<<")
      source(here("src/02_visualization.R"), echo = FALSE, local = FALSE)
      
      message("\n>>> RUNNING PHASE 3: STATISTICAL ANALYSIS & REPORTING <<<")
      source(here("src/03_statistical_analysis.R"), echo = FALSE, local = FALSE)
      
      final_msg <- sprintf("\n=== EXPERIMENT FINISHED SUCCESSFULLY: %s ===", Sys.time())
      message(final_msg)
      
    }, message = log_handler, warning = warn_handler)
    
  }, error = function(e) {
    err_msg <- paste0("\n[FATAL ERROR] Experiment stopped: ", e$message)
    cat(err_msg, file = log_file, append = TRUE, sep = "\n")
    message(err_msg)
  })
}

options(crayon.enabled = TRUE)
message("\n=== ALL EXPERIMENTS COMPLETED ===")