# src/04_longitudinal_lmm.R
# ==============================================================================
# STEP 04: LONGITUDINAL LMM ANALYSIS
# Description: Loads joint processed data to evaluate Temporal x Clinical interaction.
# ==============================================================================

source("R/utils_io.R")
source("R/modules_longitudinal.R")

message("\n=== PIPELINE STEP 4: LONGITUDINAL ANALYSIS (LMM) ===")

# 1. Configuration & Data Loading
if (!exists("config")) stop("[FATAL] Configuration not loaded.")

input_rds <- file.path(config$output_root, "01_data_processing", sprintf("data_processed_%s_longitudinal.rds", config$project_name))
if (!file.exists(input_rds)) stop(sprintf("[FATAL] Step 01 longitudinal output not found: %s", input_rds))

out_dir <- file.path(config$output_root, "04_longitudinal_analysis")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(sprintf("[Data] Loading Joint Processed Dataset: %s", basename(input_rds)))
DATA <- readRDS(input_rds)

df_long <- DATA$hybrid_data_raw
meta_long <- DATA$metadata

df_model <- cbind(meta_long, as.data.frame(df_long[, DATA$hybrid_markers]))

df_model$Timepoint <- factor(df_model$Timepoint, levels = c("T0", "T1"))
df_model$Patient_ID <- as.factor(df_model$Patient_ID)

ref_level <- config$clinical$non_responder_label
avail_levels <- unique(as.character(df_model$Group))

if (ref_level %in% avail_levels) {
  df_model$Group <- factor(df_model$Group, levels = c(ref_level, setdiff(avail_levels, ref_level)))
} else {
  df_model$Group <- factor(df_model$Group)
}

message(sprintf("   [Data] Longitudinal Matrix formed: %d total observations (T0: %d, T1: %d)", 
                nrow(df_model), sum(df_model$Timepoint == "T0"), sum(df_model$Timepoint == "T1")))

# Extract covariates from config if defined
covariates_list <- if (!is.null(config$clinical$covariates)) unlist(config$clinical$covariates) else NULL

# 2. Execution of Linear Mixed Models
message("\n[Stats] Running Linear Mixed Models (LMM) for all markers...")

results_list <- list()
pb <- txtProgressBar(min = 0, max = length(DATA$hybrid_markers), style = 3)

for (i in seq_along(DATA$hybrid_markers)) {
  mk <- DATA$hybrid_markers[i]
  res <- fit_feature_lmm(data_long = df_model, feature = mk, group_col = "Group", 
                         time_col = "Timepoint", id_col = "Patient_ID", 
                         covariates = covariates_list)
  results_list[[mk]] <- res
  setTxtProgressBar(pb, i)
}
close(pb)

df_results <- do.call(rbind, results_list)
rownames(df_results) <- NULL

# 3. FDR Correction 
df_results <- df_results %>%
  dplyr::filter(!is.na(P_Value_Interaction)) %>%
  dplyr::mutate(FDR_Interaction = p.adjust(P_Value_Interaction, method = "BH")) %>%
  dplyr::arrange(P_Value_Interaction)

n_sig_raw <- sum(df_results$P_Value_Interaction < 0.05, na.rm = TRUE)
n_sig_fdr <- sum(df_results$FDR_Interaction < 0.05, na.rm = TRUE)

message(sprintf("\n[Stats] LMM Analysis complete. %d markers significant (p < 0.05), %d survive FDR correction.", 
                n_sig_raw, n_sig_fdr))

# 4. Exporting Results
json_path <- file.path(out_dir, sprintf("Machine_Metrics_LMM_%s.json", config$project_name))
machine_output <- list(
  project_name = config$project_name,
  clinical_target = config$clinical$target_column,
  model_type = "LMM_Interaction",
  n_observations = nrow(df_model),
  n_patients = length(unique(df_model$Patient_ID)),
  significant_features_fdr = n_sig_fdr,
  full_results = df_results
)

if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(machine_output, json_path, pretty = TRUE, auto_unbox = TRUE)
}

excel_path <- file.path(out_dir, sprintf("Longitudinal_LMM_Report_%s.xlsx", config$project_name))
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "LMM_Interaction_Results")
openxlsx::writeData(wb, "LMM_Interaction_Results", df_results)

sig_style <- openxlsx::createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
p_col_idx <- which(names(df_results) == "P_Value_Interaction")

if (length(p_col_idx) > 0) {
  openxlsx::conditionalFormatting(wb, "LMM_Interaction_Results", cols = p_col_idx, 
                                  rows = 2:(nrow(df_results)+1), rule = "< 0.05", style = sig_style)
}
openxlsx::saveWorkbook(wb, excel_path, overwrite = TRUE)
message(sprintf("   [Output] Full statistics saved: %s", basename(excel_path)))

# 5. Volcano Plot Generation
plot_path <- file.path(out_dir, sprintf("Volcano_LMM_%s.pdf", config$project_name))
pdf(plot_path, width = 9, height = 7)
tryCatch({
  p_volcano <- plot_lmm_volcano(df_results, title = sprintf("LMM Interaction: Time x %s", config$clinical$target_column))
  if (!is.null(p_volcano)) print(p_volcano)
}, error = function(e) warning(paste("Volcano plot failed:", e$message)))
dev.off()

message("=== STEP 4 COMPLETE ===\n")