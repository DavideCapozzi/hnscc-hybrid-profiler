# run_step06_only.R
# Minimal runner: re-executes Step 06 (ML) for BestResponse_2v3_4 standard pass only.
# Use when Step 01/03/04 outputs are still valid (the LMM JSON gate is the only
# upstream dependency). Replicates exactly the config setup that main.R performs.

library(here)
library(yaml)
library(purrr)
library(dplyr)

list.files(here("R"), pattern = "\\.R$", full.names = TRUE) %>% purrr::walk(source)

base_config <- yaml::read_yaml(here("config/global_params.yml"))
exp_cfg     <- base_config$experiments$BestResponse_2v3_4

config <- base_config
config$project_name <- "BestResponse_2v3_4"
config$run_mode     <- "standard"
if (!is.null(exp_cfg$clinical)) config$clinical <- exp_cfg$clinical
config$is_longitudinal <- FALSE
if (!is.null(exp_cfg$input_file_t0)) config$input_file_t0 <- exp_cfg$input_file_t0
if (!is.null(exp_cfg$input_file_t1)) config$input_file_t1 <- exp_cfg$input_file_t1
config$output_root <- file.path(base_config$output_root, config$project_name)

message("=== STEP 06 ONLY — BestResponse_2v3_4 standard pass ===")
message(sprintf("LMM gate: %s",
                file.path(config$output_root, "04_longitudinal_analysis",
                          sprintf("Machine_Metrics_LMM_%s.json", config$project_name))))

source(here("src/06_machine_learning.R"), echo = FALSE, local = FALSE)

message("\n=== Step 06 complete. Check results/BestResponse_2v3_4/06_machine_learning/ ===")
