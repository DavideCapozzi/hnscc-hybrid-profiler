# run_step04_only.R
# Minimal runner: re-executes Step 04 (LMM) for BestResponse_2v3_4 longitudinal pass only.
# Use when only config$clinical$sensitivity_covariates changed and Steps 01/02/03/05/06
# outputs are still valid. Replicates exactly the config setup that main.R performs.

library(here)
library(yaml)
library(purrr)
library(dplyr)

# Load all R modules (same as main.R)
list.files(here("R"), pattern = "\\.R$", full.names = TRUE) %>% purrr::walk(source)

base_config <- yaml::read_yaml(here("config/global_params.yml"))
exp_cfg     <- base_config$experiments$BestResponse_2v3_4

# Replicate main.R longitudinal pass config setup exactly
config <- base_config
config$project_name <- "BestResponse_2v3_4"
config$run_mode     <- "longitudinal"
if (!is.null(exp_cfg$clinical)) config$clinical <- exp_cfg$clinical
config$is_longitudinal   <- TRUE
config$qc$remove_outliers <- FALSE  # Critical structural rule: prevent temporal censoring
if (!is.null(exp_cfg$input_file_t0)) config$input_file_t0 <- exp_cfg$input_file_t0
if (!is.null(exp_cfg$input_file_t1)) config$input_file_t1 <- exp_cfg$input_file_t1
config$output_root <- file.path(base_config$output_root, config$project_name)

message("=== STEP 04 ONLY — BestResponse_2v3_4 longitudinal pass ===")
message(sprintf("Sensitivity covariates: %s",
                paste(config$clinical$sensitivity_covariates, collapse = ", ")))

source(here("src/04_longitudinal_lmm.R"), echo = FALSE, local = FALSE)

message("\n=== Step 04 complete. Check results/BestResponse_2v3_4/04_longitudinal_analysis/ ===")
