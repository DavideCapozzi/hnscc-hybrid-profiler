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

# Load Configuration
config_path <- here("config/global_params.yml")
if (!file.exists(config_path)) stop("[FATAL] Config file not found!")
config <- yaml::read_yaml(config_path)

# 2. Logging Setup
# ------------------------------------------------------------------------------
out_root <- here(config$output_root)
if (!is.null(config$project_name) && config$project_name != "") {
  out_root <- paste0(out_root, "_", config$project_name)
}
if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(out_root, paste0("pipeline_log_", timestamp, ".txt"))

cat(sprintf("=== ONCO-HYBRID-PROFILER STARTED: %s ===\n", Sys.time()), file = log_file)
message(sprintf("[System] Saving full log to: %s", log_file))

log_handler <- function(m) {
  cat(conditionMessage(m), file = log_file, append = TRUE, sep = "\n")
}

warn_handler <- function(w) {
  cat(paste0("WARNING: ", conditionMessage(w)), file = log_file, append = TRUE, sep = "\n")
}

# 3. Pipeline Execution
# ------------------------------------------------------------------------------
tryCatch({
  withCallingHandlers({
    
    # --- Module Loading ---
    message("\n>>> LOADING MODULES <<<")
    list.files(here("R"), pattern = "\\.R$", full.names = TRUE) %>% purrr::walk(source)
    message("[System] Modules loaded successfully.")
    
    # --- PHASE 1: DATA INGESTION, QC & HYBRID TRANSFORM ---
    message("\n>>> RUNNING PHASE 1: DATA PROCESSING <<<")
    source(here("src/01_data_processing.R"), echo = FALSE)
    
    # --- PHASE 2: OVERVIEW VISUALIZATION ---
    message("\n>>> RUNNING PHASE 2: VISUALIZATION <<<")
    source(here("src/02_visualization.R"), echo = FALSE)
    
    # --- PHASE 3: MULTIVARIATE ANALYSIS (sPLS-DA) & REPORTING ---
    message("\n>>> RUNNING PHASE 3: STATISTICAL ANALYSIS & REPORTING <<<")
    source(here("src/03_statistical_analysis.R"), echo = FALSE)
    
    final_msg <- sprintf("\n=== PIPELINE FINISHED SUCCESSFULLY: %s ===", Sys.time())
    message(final_msg)
    
  }, message = log_handler, warning = warn_handler)
  
}, error = function(e) {
  err_msg <- paste0("\n[FATAL ERROR] Pipeline stopped: ", e$message)
  cat(err_msg, file = log_file, append = TRUE, sep = "\n")
  stop(e$message)
}, finally = {
  options(crayon.enabled = TRUE)
})
