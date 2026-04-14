# src/04_longitudinal_lmm.R
# ==============================================================================
# STEP 04: LONGITUDINAL LMM ANALYSIS
# Description: Fits Linear Mixed Models to evaluate Time x Clinical Interaction.
#              Includes LOO Sensitivity Analysis to protect against outlier bias.
# ==============================================================================

source("R/utils_io.R")
source("R/modules_longitudinal.R")
source("R/modules_viz.R")

message("\n=== PIPELINE STEP 4: LONGITUDINAL ANALYSIS (LMM) ===")

# 1. Configuration & Data Parsing
# ------------------------------------------------------------------------------
if (!exists("config")) stop("[FATAL] Global configuration object not detected in environment.")

input_rds <- file.path(config$output_root, "01_data_processing", sprintf("data_processed_%s_longitudinal.rds", config$project_name))
if (!file.exists(input_rds)) stop(sprintf("[FATAL] Step 01 longitudinal processed dataset not found at: %s", input_rds))

out_dir <- file.path(config$output_root, "04_longitudinal_analysis")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(sprintf("[Data] Loading Joint Pharmacodynamic Dataset: %s", basename(input_rds)))
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

message(sprintf("   [Data] Longitudinal dimensions established: %d matrix observations (T0: %d, T1: %d)", 
                nrow(df_model), sum(df_model$Timepoint == "T0"), sum(df_model$Timepoint == "T1")))

covariates_list <- if (!is.null(config$clinical$covariates)) unlist(config$clinical$covariates) else NULL

# --- SAFETY CHECK: Prevent Silent Cohort Decimation via Covariates ---
safe_max_na_covariates <- if(!is.null(config$qc$max_na_covariates)) config$qc$max_na_covariates else 0.10

if (!is.null(covariates_list) && length(intersect(covariates_list, colnames(df_model))) > 0) {
  valid_covs <- intersect(covariates_list, colnames(df_model))
  na_rows <- sum(!complete.cases(df_model[, valid_covs, drop = FALSE]))
  pct_na <- na_rows / nrow(df_model)
  
  if (pct_na > safe_max_na_covariates) {
    stop(sprintf("[FATAL] Security Gate Triggered: Covariates induce %.1f%% missingness, exceeding the safety threshold (%.1f%%). Aborting to prevent silent cohort decimation.", pct_na * 100, safe_max_na_covariates * 100))
  } else if (pct_na > 0) {
    message(sprintf("   [Data] Notice: Covariate inclusion will drop %d observations (%.1f%% of cohort).", na_rows, pct_na * 100))
  }
}

colors_viz <- get_clinical_colors(config)

# 2. Linear Mixed Models Computation
# ------------------------------------------------------------------------------
message("\n[Stats] Initiating Linear Mixed Models (LMM) for longitudinal variance extraction...")

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

# 3. Multi-Testing Correction (FDR)
# ------------------------------------------------------------------------------
df_results <- df_results %>%
  dplyr::filter(!is.na(P_Value_Interaction)) %>%
  dplyr::mutate(FDR_Interaction = p.adjust(P_Value_Interaction, method = "BH")) %>%
  dplyr::arrange(P_Value_Interaction)

n_sig_raw <- sum(df_results$P_Value_Interaction < 0.05, na.rm = TRUE)
n_sig_fdr <- sum(df_results$FDR_Interaction < 0.05, na.rm = TRUE)

message(sprintf("\n[Stats] LMM integration successful. Identified %d raw significant markers, %d survive FDR adjustment.", 
                n_sig_raw, n_sig_fdr))

# 4. LOO Sensitivity Security Check (FDR subset only)
# ------------------------------------------------------------------------------
df_results$Max_P_Value_LOO <- NA

if (n_sig_fdr > 0) {
  message("   [Stats] Executing Leave-One-Out (LOO) Sensitivity protocol on topological drivers...")
  sig_markers <- df_results$Marker[which(df_results$FDR_Interaction < 0.05)]
  
  for (mk in sig_markers) {
    max_p <- run_loo_sensitivity(data_long = df_model, feature = mk, group_col = "Group", 
                                 time_col = "Timepoint", id_col = "Patient_ID", 
                                 covariates = covariates_list)
    
    df_results$Max_P_Value_LOO[df_results$Marker == mk] <- max_p
    
    if (!is.na(max_p) && max_p < 0.05) {
      message(sprintf("      -> %s: Structurally Robust (Max P-Value = %.4f)", mk, max_p))
    } else {
      warning(sprintf("      -> %s: OUTLIER BIAS DETECTED (Max P-Value spikes to %.4f upon LOO)", mk, max_p))
    }
  }
}

# 5. Output Serialization
# ------------------------------------------------------------------------------
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
message(sprintf("   [Output] Mathematical models matrix saved: %s", basename(excel_path)))

# 6. Volcano Plot & Trajectory Rendering
# ------------------------------------------------------------------------------
plot_path <- file.path(out_dir, sprintf("Volcano_LMM_%s.pdf", config$project_name))
pdf(plot_path, width = 9, height = 7)
tryCatch({
  p_volcano <- plot_lmm_volcano(df_results, title = sprintf("LMM Pharmacodynamics: Time x %s", config$clinical$target_column))
  if (!is.null(p_volcano)) print(p_volcano)
}, error = function(e) warning(paste("Volcano plot rendering failed:", e$message)))
dev.off()

message("\n[Viz] Generating Patient Trajectory vectors...")

top_df <- df_results %>% dplyr::filter(FDR_Interaction < 0.05)
sig_type <- "FDR"

if (nrow(top_df) == 0) {
  top_df <- df_results %>% dplyr::arrange(P_Value_Interaction) %>% head(4)
  sig_type <- "RAW"
  message("   [Viz] Null set for FDR boundaries. Defaulting plot limits to top 4 absolute P-values.")
}

if (nrow(top_df) > 0) {
  traj_path <- file.path(out_dir, sprintf("Trajectories_LMM_%s.pdf", config$project_name))
  pdf(traj_path, width = 8, height = 6)
  
  for (i in 1:nrow(top_df)) {
    mk <- top_df$Marker[i]
    pval_disp <- if (sig_type == "FDR") top_df$FDR_Interaction[i] else top_df$P_Value_Interaction[i]
    
    tryCatch({
      p_traj <- plot_lmm_trajectories(
        data_long = df_model,
        feature = mk,
        group_col = "Group",
        time_col = "Timepoint",
        id_col = "Patient_ID",
        colors = colors_viz,
        p_val = pval_disp
      )
      print(p_traj)
    }, error = function(e) warning(sprintf("Trajectory matrix error on marker %s: %s", mk, e$message)))
  }
  
  dev.off()
  message(sprintf("   [Output] Temporal trajectory visualizations exported: %s", basename(traj_path)))
}

message("=== STEP 4 COMPLETE ===\n")