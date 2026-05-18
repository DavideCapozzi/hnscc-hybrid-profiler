# src/06_machine_learning.R
# ==============================================================================
# STEP 06: HYBRID TWO-STEP MACHINE LEARNING CLASSIFIER
# Description: Nested-LOOCV classification on cross-sectional T0 data using
#              only markers that survived the LMM LOO robustness gate (Step 04).
#              Feature selection is frozen before any CV iteration — no leakage.
#              Primary: Elastic Net Logistic Regression. Secondary: SVM-RBF.
#              Degrades gracefully when no LOO-robust markers exist for an
#              experiment (e.g., Comorbidity, Toxicity with FDR > 0.05).
# ==============================================================================

source("R/utils_io.R")
source("R/modules_ml.R")
source("R/modules_viz.R")

message("\n=== PIPELINE STEP 6: MACHINE LEARNING CLASSIFICATION ===")

# 1. Configuration Guard
# ------------------------------------------------------------------------------
if (!exists("config")) stop("[FATAL] Global configuration object not detected in environment.")

ml_cfg <- if (!is.null(config$machine_learning)) config$machine_learning else list()

# 2. Load LMM LOO-Robust Features (auto-detect, graceful degradation)
# ------------------------------------------------------------------------------
lmm_json_path <- file.path(
  config$output_root, "04_longitudinal_analysis",
  sprintf("Machine_Metrics_LMM_%s.json", config$project_name)
)

fdr_thresh <- if (!is.null(ml_cfg$fdr_threshold)) as.numeric(ml_cfg$fdr_threshold) else 0.05
loo_thresh <- if (!is.null(ml_cfg$loo_threshold)) as.numeric(ml_cfg$loo_threshold) else 0.05

lmm_robust <- load_lmm_robust_features(lmm_json_path, fdr_thresh, loo_thresh)

if (lmm_robust$n_robust == 0) {
  message(sprintf(
    "[ML] Experiment '%s': no LOO-robust markers found. Step 06 complete with no output.",
    config$project_name
  ))
} else {

  # 3. Load Step 01 Standard RDS (cross-sectional T0 data)
  # Note: this is intentionally the STANDARD (not longitudinal) output. The LMM
  # feature gate was derived from longitudinal T0/T1 data in Step 04, keeping the
  # feature selection and classification data modalities strictly separate.
  # ------------------------------------------------------------------------------
  input_rds <- file.path(
    config$output_root, "01_data_processing",
    sprintf("data_processed_%s_standard.rds", config$project_name)
  )
  if (!file.exists(input_rds)) {
    stop(sprintf("[FATAL] Step 01 standard RDS not found at: %s", input_rds))
  }

  DATA <- readRDS(input_rds)

  out_dir <- file.path(config$output_root, "06_machine_learning")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  colors_viz <- get_clinical_colors(config)

  # 4. Build Feature Matrix
  # ------------------------------------------------------------------------------
  include_interactions <- if (!is.null(ml_cfg$include_interactions)) {
    isTRUE(ml_cfg$include_interactions)
  } else TRUE

  cor_threshold <- if (!is.null(ml_cfg$collinearity_threshold)) {
    as.numeric(ml_cfg$collinearity_threshold)
  } else 0.85

  ml_matrix <- build_ml_matrix(DATA, lmm_robust$markers,
                                include_interactions = include_interactions,
                                cor_threshold        = cor_threshold)
  X <- ml_matrix$X
  y <- ml_matrix$y

  message(sprintf("[ML] Input: %d samples x %d features (main: %d) | Classes: %s",
                  nrow(X), ncol(X), ml_matrix$n_main,
                  paste(names(table(y)), table(y), sep = "=", collapse = ", ")))

  # 5. Nested-LOOCV — Elastic Net Logistic Regression (primary)
  # ------------------------------------------------------------------------------
  glmnet_cfg  <- if (!is.null(ml_cfg$glmnet)) ml_cfg$glmnet else list()
  alpha_grid  <- if (!is.null(glmnet_cfg$alpha_grid))  as.numeric(unlist(glmnet_cfg$alpha_grid)) else c(0, 0.5, 1)
  n_lambda    <- if (!is.null(glmnet_cfg$n_lambda))    as.integer(glmnet_cfg$n_lambda)            else 100L
  inner_k_en  <- if (!is.null(glmnet_cfg$inner_folds)) as.integer(glmnet_cfg$inner_folds)         else 5L

  message("\n[ML] Running Nested-LOOCV: Elastic Net Logistic Regression...")
  set.seed(config$stats$seed)
  res_glmnet <- run_nested_loocv_glmnet(
    X           = X,
    y           = y,
    alpha_grid  = alpha_grid,
    n_lambda    = n_lambda,
    inner_folds = inner_k_en,
    seed        = config$stats$seed
  )

  metrics_glmnet <- compute_classification_metrics(
    res_glmnet$y_true, res_glmnet$predicted_probs, res_glmnet$positive_label
  )
  message(sprintf(
    "   [Elastic Net] AUC=%.3f [%.3f–%.3f] | BalAcc=%.3f | BER=%.3f | Sens=%.3f | Spec=%.3f",
    metrics_glmnet$auc,
    metrics_glmnet$auc_ci[1], metrics_glmnet$auc_ci[3],
    metrics_glmnet$balanced_accuracy,
    metrics_glmnet$ber,
    metrics_glmnet$sensitivity,
    metrics_glmnet$specificity
  ))

  # 6. Nested-LOOCV — SVM with RBF Kernel (secondary comparison)
  # ------------------------------------------------------------------------------
  svm_cfg    <- if (!is.null(ml_cfg$svm)) ml_cfg$svm else list()
  C_grid     <- if (!is.null(svm_cfg$C_grid))     as.numeric(unlist(svm_cfg$C_grid))     else c(0.01, 0.1, 1, 10, 100)
  gamma_grid <- if (!is.null(svm_cfg$gamma_grid)) as.numeric(unlist(svm_cfg$gamma_grid)) else c(0.01, 0.1, 1, 10)
  inner_k_sv <- if (!is.null(svm_cfg$inner_folds)) as.integer(svm_cfg$inner_folds)       else 5L

  message("\n[ML] Running Nested-LOOCV: SVM with RBF Kernel...")
  set.seed(config$stats$seed)
  res_svm <- run_nested_loocv_svm(
    X           = X,
    y           = y,
    C_grid      = C_grid,
    gamma_grid  = gamma_grid,
    inner_folds = inner_k_sv,
    seed        = config$stats$seed
  )

  metrics_svm <- compute_classification_metrics(
    res_svm$y_true, res_svm$predicted_probs, res_svm$positive_label
  )
  message(sprintf(
    "   [SVM-RBF]     AUC=%.3f [%.3f–%.3f] | BalAcc=%.3f | BER=%.3f | Sens=%.3f | Spec=%.3f",
    metrics_svm$auc,
    metrics_svm$auc_ci[1], metrics_svm$auc_ci[3],
    metrics_svm$balanced_accuracy,
    metrics_svm$ber,
    metrics_svm$sensitivity,
    metrics_svm$specificity
  ))

  # 7. Determine Primary Method
  # Elastic Net is default primary. Swap only if SVM-RBF exceeds it by > 0.05 AUC.
  # ------------------------------------------------------------------------------
  primary_method <- if (
    !is.na(metrics_glmnet$auc) && !is.na(metrics_svm$auc) &&
    metrics_svm$auc > metrics_glmnet$auc + 0.05
  ) "SVM-RBF" else "Elastic-Net"

  message(sprintf("\n[ML] Primary method designated: %s", primary_method))

  # 8. Coefficient Stability Analysis (Elastic Net only)
  # ------------------------------------------------------------------------------
  df_coef_stability <- NULL
  if (!is.null(res_glmnet$coef_matrix)) {
    coef_mat <- res_glmnet$coef_matrix

    df_coef_stability <- data.frame(
      Feature   = colnames(coef_mat),
      Mean_Coef = round(colMeans(coef_mat, na.rm = TRUE), 5),
      SD_Coef   = round(apply(coef_mat, 2, sd, na.rm = TRUE), 5),
      Sign_Flip = apply(coef_mat, 2, function(x) {
        x_nz <- x[!is.na(x) & x != 0]
        length(x_nz) >= 2 && any(x_nz > 0) && any(x_nz < 0)
      }),
      stringsAsFactors = FALSE
    )

    unstable <- df_coef_stability$Feature[df_coef_stability$Sign_Flip]
    unstable <- unstable[unstable != "(Intercept)"]
    if (length(unstable) > 0) {
      warning(sprintf(
        "[ML] Sign flip detected in Elastic Net coefficients for: %s — interpret with caution.",
        paste(unstable, collapse = ", ")
      ))
    }
  }

  # 9. Permutation AUC Test
  # Uses the already-computed out-of-fold probabilities — no model re-fitting.
  # ------------------------------------------------------------------------------
  n_perm <- if (!is.null(ml_cfg$n_perm)) as.integer(ml_cfg$n_perm) else 2000L

  message(sprintf("\n[ML] Running permutation AUC test (n_perm=%d)...", n_perm))
  set.seed(config$stats$seed)
  perm_glmnet <- run_permutation_auc_test(
    res_glmnet$y_true, res_glmnet$predicted_probs, res_glmnet$positive_label,
    n_perm = n_perm, seed = config$stats$seed
  )
  perm_svm <- run_permutation_auc_test(
    res_svm$y_true, res_svm$predicted_probs, res_svm$positive_label,
    n_perm = n_perm, seed = config$stats$seed
  )
  message(sprintf(
    "   [Perm] Elastic Net p=%.4f | SVM-RBF p=%.4f",
    perm_glmnet$p_value, perm_svm$p_value
  ))

  # 10. Univariate AUC and LOO-Threshold Classification
  # X_main = main effects only (no interaction columns) — avoids any selection bias
  # from the pairwise terms constructed in build_ml_matrix().
  # ------------------------------------------------------------------------------
  X_main <- X[, seq_len(ml_matrix$n_main), drop = FALSE]

  message("\n[ML] Running univariate AUC analysis (parameter-free)...")
  uni_auc_df <- run_univariate_auc(X_main, y, res_glmnet$positive_label)
  message(sprintf("   [Univariate AUC] %s",
                  paste(sprintf("%s=%.3f", uni_auc_df$Marker, uni_auc_df$AUC),
                        collapse = " | ")))

  message("\n[ML] Running univariate LOO-threshold classification...")
  set.seed(config$stats$seed)
  uni_thresh_df <- run_univariate_loo_threshold(
    X_main, y, res_glmnet$positive_label, seed = config$stats$seed
  )
  message(sprintf("   [LOO Threshold] %s",
                  paste(sprintf("%s BalAcc=%.3f", uni_thresh_df$Marker,
                                uni_thresh_df$Balanced_Accuracy), collapse = " | ")))

  # 10b. Clinical Benchmark Comparison (optional, config-driven)
  # Compares the primary model against a clinical standard biomarker (e.g. PD-L1)
  # read from config$clinical$benchmark_column in the raw input Excel.
  # The DeLong test is run on the available-cases subset only; it may be underpowered
  # and should be reported as descriptive when p > 0.05.
  # ------------------------------------------------------------------------------
  benchmark_result <- NULL
  bench_col   <- config$clinical$benchmark_column
  bench_label <- if (!is.null(config$clinical$benchmark_label)) {
    config$clinical$benchmark_label
  } else bench_col

  if (!is.null(bench_col) && !is.null(config$input_file_t0)) {
    message(sprintf("\n[ML] Clinical benchmark: '%s'...", bench_label))
    primary_res <- if (primary_method == "SVM-RBF") res_svm else res_glmnet
    benchmark_result <- tryCatch(
      run_clinical_benchmark(
        DATA            = DATA,
        primary_probs   = primary_res$predicted_probs,
        y_primary       = primary_res$y_true,
        positive_label  = res_glmnet$positive_label,
        input_file      = config$input_file_t0,
        benchmark_col   = bench_col,
        benchmark_label = bench_label
      ),
      error = function(e) {
        warning(sprintf("[ML] Clinical benchmark failed: %s", e$message))
        NULL
      }
    )
  }

  # 10c. Benchmark-Stratified + Combined Model Information Gain
  # Activated when a clinical benchmark column is configured. Adds two
  # complementary perspectives beyond the basic DeLong comparison:
  #   - Stratified analysis: clinical-strata bins (e.g. PD-L1 TPS) + per-bin
  #     subgroup AUC of the model (most clinically actionable: benchmark-low)
  #   - Combined model: logistic LOOCV with benchmark + gate features, with
  #     IDI / continuous NRI and patient-level bootstrap CIs
  # ------------------------------------------------------------------------------
  stratified_result   <- NULL
  combined_result     <- NULL

  if (!is.null(bench_col) && !is.null(config$input_file_t0)) {
    message(sprintf("\n[ML] Benchmark-stratified subgroup analysis: '%s'...", bench_label))
    stratified_result <- tryCatch(
      run_pdl1_stratified(
        DATA            = DATA,
        X_main          = X_main,
        y               = y,
        positive_label  = res_glmnet$positive_label,
        input_file      = config$input_file_t0,
        benchmark_col   = bench_col,
        benchmark_label = bench_label,
        C_grid          = C_grid,
        gamma_grid      = gamma_grid,
        inner_folds     = inner_k_sv,
        seed            = config$stats$seed
      ),
      error = function(e) {
        warning(sprintf("[ML] Stratified analysis failed: %s", e$message))
        NULL
      }
    )

    if (!is.null(stratified_result)) {
      bc <- stratified_result$binary_cut
      message(sprintf("   [Stratified] Binary cut at %g: BalAcc=%.3f Sens=%.3f Spec=%.3f Fisher p=%.4f",
                      bc$threshold, bc$balanced_accuracy, bc$sensitivity, bc$specificity, bc$fisher_p))
      sl <- stratified_result$subgroup_low; sh <- stratified_result$subgroup_high
      if (!is.na(sl$auc)) message(sprintf("   [Stratified] Subgroup low  (n=%d): SVM AUC=%.3f", sl$n, sl$auc))
      if (!is.na(sh$auc)) message(sprintf("   [Stratified] Subgroup high (n=%d): SVM AUC=%.3f", sh$n, sh$auc))
    }

    message(sprintf("\n[ML] Combined model (information gain over '%s')...", bench_label))
    combined_result <- tryCatch(
      run_combined_benchmark_model(
        DATA            = DATA,
        X_main          = X_main,
        y               = y,
        positive_label  = res_glmnet$positive_label,
        input_file      = config$input_file_t0,
        benchmark_col   = bench_col,
        benchmark_label = bench_label,
        n_boot          = 1000L,
        C_grid          = C_grid,
        gamma_grid      = gamma_grid,
        inner_folds     = inner_k_sv,
        seed            = config$stats$seed
      ),
      error = function(e) {
        warning(sprintf("[ML] Combined model analysis failed: %s", e$message))
        NULL
      }
    )

    if (!is.null(combined_result)) {
      ig <- combined_result$information_gain
      lg <- combined_result$logistic
      sv <- combined_result$svm_comparison
      message(sprintf("   [Combined] Logistic %s vs %s+gate: AUC %.3f -> %.3f (ΔAUC=%+.3f, DeLong p=%.4f)",
                      bench_label, bench_label, lg$auc_benchmark, lg$auc_combined,
                      lg$delta_auc, lg$delong_p))
      message(sprintf("   [Combined] IDI=%+.4f [%+.4f, %+.4f]  bootstrap p=%.4f",
                      ig$IDI, ig$IDI_95CI[[1]], ig$IDI_95CI[[2]], ig$IDI_bootstrap_p))
      message(sprintf("   [Combined] cNRI=%+.4f [%+.4f, %+.4f]  bootstrap p=%.4f",
                      ig$cNRI, ig$cNRI_95CI[[1]], ig$cNRI_95CI[[2]], ig$cNRI_bootstrap_p))
      message(sprintf("   [Combined] SVM features-only=%.3f vs features+benchmark=%.3f (ΔAUC=%+.3f, DeLong p=%.4f)",
                      sv$auc_features_only, sv$auc_features_with_bench,
                      sv$delta_auc, sv$delong_p))
    }
  }

  # 11. Nested LOO Validation (optional, config-driven)
  # ------------------------------------------------------------------------------
  nested_validation <- NULL
  run_nested <- isTRUE(config$run_nested_loocv_validation)

  if (run_nested) {
    message("\n[ML] Running fully-nested LOO validation (LMM gate inside outer fold)...")
    lmm_json_for_nested <- file.path(config$output_root, "04_longitudinal_analysis",
                                     sprintf("Machine_Metrics_LMM_%s.json", config$project_name))

    DATA_LONG_nested <- tryCatch({
      rds_long <- file.path(config$output_root, "01_data_processing",
                            sprintf("data_processed_%s_longitudinal.rds", config$project_name))
      if (file.exists(rds_long)) readRDS(rds_long) else NULL
    }, error = function(e) NULL)

    if (!is.null(DATA_LONG_nested)) {
      nested_validation <- tryCatch(
        run_nested_loocv_svm_validated(
          DATA_T0       = DATA,
          DATA_LONG     = DATA_LONG_nested,
          fdr_thr       = fdr_thresh,
          loo_thr       = loo_thresh,
          cor_threshold = cor_threshold,
          C_grid        = C_grid,
          gamma_grid    = gamma_grid,
          inner_folds   = inner_k_sv,
          seed          = config$stats$seed
        ),
        error = function(e) {
          warning(sprintf("[ML] Nested LOO validation failed: %s", e$message))
          NULL
        }
      )

      if (!is.null(nested_validation)) {
        m <- nested_validation$metrics
        message(sprintf(
          "   [ML] Nested LOO validation complete — AUC=%.4f [%.3f–%.3f]  BalAcc=%.4f",
          m$auc, m$auc_ci[1], m$auc_ci[3], m$balanced_accuracy
        ))
        message(sprintf("   [ML] Gate stability: %d / %d folds with non-empty gate",
                        nrow(DATA$metadata) - nested_validation$n_empty_folds,
                        nrow(DATA$metadata)))
      }
    } else {
      warning("[ML] Nested LOO validation skipped: longitudinal data not found.")
    }
  }

  # 12. Export Machine-Readable JSON
  # ------------------------------------------------------------------------------
  machine_output <- list(
    project_name         = config$project_name,
    clinical_target      = config$clinical$target_column,
    model_type           = "Hybrid Two-Step Nested-LOOCV",
    n_samples            = nrow(X),
    n_features_main      = ml_matrix$n_main,
    n_features_total     = ncol(X),
    include_interactions = include_interactions,
    robust_markers       = as.list(lmm_robust$markers),
    primary_method       = primary_method,
    elastic_net          = list(
      method  = res_glmnet$method,
      metrics = metrics_glmnet
    ),
    svm_rbf              = list(
      method  = res_svm$method,
      metrics = metrics_svm
    ),
    scaling_note = paste(
      "Z-scores computed globally in Step 01 (deterministic centering, not a learned parameter).",
      "For linear classifiers, global z-scoring is an affine transformation that cancels out;",
      "for SVM-RBF, leakage potential was ruled out empirically (diag_11): refitting with strict",
      "per-fold z-scoring from raw transformed input yields identical AUC (DeLong p=1.0).",
      "Feature gate derived from longitudinal data (Step 04) — disjoint from classifier training data.",
      if (!is.null(nested_validation))
        sprintf("Leakage further quantified: fully-nested LOO validation yielded AUC=%.3f [%.3f-%.3f] with %d%% gate stability for primary marker (see nested_loocv_validation).",
                nested_validation$metrics$auc,
                nested_validation$metrics$auc_ci[1],
                nested_validation$metrics$auc_ci[3],
                if (nrow(nested_validation$gate_stability) > 0) nested_validation$gate_stability$Pct[1] else 0L)
      else ""
    ),
    permutation_test = list(
      elastic_net = perm_glmnet,
      svm_rbf     = perm_svm,
      n_perm      = n_perm
    ),
    univariate_auc            = uni_auc_df,
    univariate_loo_threshold  = uni_thresh_df,
    clinical_benchmark        = if (!is.null(benchmark_result)) {
      list(
        label                 = benchmark_result$label,
        column                = benchmark_result$column,
        n_valid               = benchmark_result$n_valid,
        n_na                  = benchmark_result$n_na,
        auc                   = benchmark_result$auc,
        auc_ci                = as.list(benchmark_result$auc_ci),
        primary_auc_on_subset = benchmark_result$primary_auc_on_subset,
        delong_p              = benchmark_result$delong_p,
        note                  = sprintf(
          "Primary model AUC=%.3f vs %s AUC=%.3f on %d/%d cases with non-NA values (DeLong p=%.4f). Test may be underpowered.",
          benchmark_result$primary_auc_on_subset,
          benchmark_result$label,
          benchmark_result$auc,
          benchmark_result$n_valid,
          nrow(X),
          if (is.na(benchmark_result$delong_p)) 0 else benchmark_result$delong_p
        )
      )
    } else NULL,
    benchmark_stratified = if (!is.null(stratified_result)) {
      # Strip per-patient raw vectors before JSON serialization (kept in-memory for plotting)
      sr_json <- stratified_result
      for (sub in c("subgroup_low", "subgroup_high")) {
        sr_json[[sub]]$predicted_probs <- NULL
        sr_json[[sub]]$y_true          <- NULL
        sr_json[[sub]]$benchmark_vals  <- NULL
      }
      sr_json
    } else NULL,
    benchmark_combined   = if (!is.null(combined_result)) {
      cr_json <- combined_result
      cr_json$plot_data <- NULL
      cr_json
    } else NULL,
    nested_loocv_validation = if (!is.null(nested_validation)) {
      nv_m <- nested_validation$metrics
      list(
        method            = "SVM-RBF Fully-Nested LOO (LMM gate re-selected inside each outer fold)",
        auc               = nv_m$auc,
        auc_ci            = as.list(nv_m$auc_ci),
        balanced_accuracy = nv_m$balanced_accuracy,
        ber               = nv_m$ber,
        sensitivity       = nv_m$sensitivity,
        specificity       = nv_m$specificity,
        gate_stability    = nested_validation$gate_stability,
        n_empty_folds     = nested_validation$n_empty_folds,
        note              = sprintf(
          "Fully-nested validation: LMM feature selection repeated within each outer LOO fold on n-1 patients. AUC=%.3f vs fixed-gate AUC=%.3f (delta=%.3f). Primary marker gate stability: %d%%.",
          nv_m$auc,
          if (primary_method == "SVM-RBF") metrics_svm$auc else metrics_glmnet$auc,
          nv_m$auc - (if (primary_method == "SVM-RBF") metrics_svm$auc else metrics_glmnet$auc),
          if (nrow(nested_validation$gate_stability) > 0) nested_validation$gate_stability$Pct[1] else 0L
        )
      )
    } else NULL
  )

  json_path <- file.path(out_dir, sprintf("Machine_Metrics_ML_%s.json", config$project_name))
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(machine_output, json_path, pretty = TRUE, auto_unbox = TRUE)
    message(sprintf("   [Output] ML metrics JSON saved: %s", basename(json_path)))
  } else {
    saveRDS(machine_output, sub("\\.json$", ".rds", json_path))
    warning("[ML] jsonlite unavailable — metrics saved as RDS.")
  }

  # 12. Export Excel Report
  # ------------------------------------------------------------------------------
  pos_lbl <- res_glmnet$positive_label
  neg_lbl <- setdiff(levels(res_glmnet$y_true), pos_lbl)

  df_preds <- data.frame(
    Patient_ID      = as.character(DATA$metadata$Patient_ID),
    True_Group      = as.character(y),
    Prob_ElasticNet = round(res_glmnet$predicted_probs, 4),
    Prob_SVM_RBF    = round(res_svm$predicted_probs,    4),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      Pred_ElasticNet = ifelse(Prob_ElasticNet >= metrics_glmnet$threshold, pos_lbl, neg_lbl),
      Pred_SVM_RBF    = ifelse(Prob_SVM_RBF    >= metrics_svm$threshold,    pos_lbl, neg_lbl),
      Correct_ElasticNet = (Pred_ElasticNet == True_Group),
      Correct_SVM_RBF    = (Pred_SVM_RBF    == True_Group)
    )

  df_perf_summary <- data.frame(
    Method            = c(res_glmnet$method, res_svm$method),
    Is_Primary        = c(primary_method == "Elastic-Net", primary_method == "SVM-RBF"),
    AUC               = c(metrics_glmnet$auc,               metrics_svm$auc),
    AUC_CI_Lower      = c(metrics_glmnet$auc_ci[1],         metrics_svm$auc_ci[1]),
    AUC_CI_Upper      = c(metrics_glmnet$auc_ci[3],         metrics_svm$auc_ci[3]),
    Balanced_Accuracy = c(metrics_glmnet$balanced_accuracy, metrics_svm$balanced_accuracy),
    BER               = c(metrics_glmnet$ber,               metrics_svm$ber),
    Sensitivity       = c(metrics_glmnet$sensitivity,       metrics_svm$sensitivity),
    Specificity       = c(metrics_glmnet$specificity,       metrics_svm$specificity),
    Youden_Threshold  = c(metrics_glmnet$threshold,         metrics_svm$threshold),
    stringsAsFactors  = FALSE
  )

  wb <- openxlsx::createWorkbook()

  openxlsx::addWorksheet(wb, "Performance_Summary")
  openxlsx::writeData(wb, "Performance_Summary", df_perf_summary)

  openxlsx::addWorksheet(wb, "Patient_Predictions")
  openxlsx::writeData(wb, "Patient_Predictions", df_preds)

  openxlsx::addWorksheet(wb, "Robust_Markers_Gate")
  openxlsx::writeData(wb, "Robust_Markers_Gate", lmm_robust$df_robust)

  if (!is.null(df_coef_stability)) {
    openxlsx::addWorksheet(wb, "ElasticNet_Coef_Stability")
    openxlsx::writeData(wb, "ElasticNet_Coef_Stability", df_coef_stability)
  }

  openxlsx::addWorksheet(wb, "Permutation_Test")
  openxlsx::writeData(wb, "Permutation_Test", data.frame(
    Method       = c("Elastic Net", "SVM-RBF"),
    Observed_AUC = c(perm_glmnet$observed_auc, perm_svm$observed_auc),
    P_Value      = c(perm_glmnet$p_value,      perm_svm$p_value),
    N_Perm       = n_perm,
    stringsAsFactors = FALSE
  ))

  openxlsx::addWorksheet(wb, "Univariate_AUC")
  openxlsx::writeData(wb, "Univariate_AUC", uni_auc_df)

  openxlsx::addWorksheet(wb, "Univariate_LOO_Threshold")
  openxlsx::writeData(wb, "Univariate_LOO_Threshold", uni_thresh_df)

  if (!is.null(benchmark_result)) {
    df_bench_report <- data.frame(
      Metric = c(
        "Benchmark biomarker", "Source column",
        "N valid (non-NA)", "N missing (NA)",
        "Benchmark AUC", "Benchmark AUC CI Lower", "Benchmark AUC CI Upper",
        "Primary model AUC (same subset)", "DeLong p-value (primary vs benchmark)",
        "Methodological note"
      ),
      Value = c(
        benchmark_result$label,
        benchmark_result$column,
        as.character(benchmark_result$n_valid),
        as.character(benchmark_result$n_na),
        as.character(round(benchmark_result$auc,                   4)),
        as.character(round(benchmark_result$auc_ci[1],             4)),
        as.character(round(benchmark_result$auc_ci[3],             4)),
        as.character(round(benchmark_result$primary_auc_on_subset, 4)),
        if (is.na(benchmark_result$delong_p)) "NA" else
          as.character(round(benchmark_result$delong_p, 5)),
        sprintf(
          "DeLong test run on available-cases subset (n=%d). Benchmark is univariate; primary model is multivariate nested-LOOCV. AUC difference is descriptive when p > 0.05.",
          benchmark_result$n_valid
        )
      ),
      stringsAsFactors = FALSE
    )
    openxlsx::addWorksheet(wb, "Clinical_Benchmark")
    openxlsx::writeData(wb, "Clinical_Benchmark", df_bench_report)
  }

  if (!is.null(stratified_result)) {
    bc <- stratified_result$binary_cut
    sl <- stratified_result$subgroup_low
    sh <- stratified_result$subgroup_high
    df_strat_summary <- data.frame(
      Metric = c(
        "Benchmark biomarker", "N valid", "Bin breaks", "Bin labels",
        "Fisher exact p (3 bins)", "Cochran-Armitage trend p",
        sprintf("Binary cut threshold (>=%s)", bc$threshold),
        "Binary cut: N", "Binary cut: N high", "Binary cut: N low",
        "Binary cut: Balanced Accuracy", "Binary cut: Sensitivity",
        "Binary cut: Specificity", "Binary cut: Fisher p", "Binary cut: Odds Ratio",
        "Subgroup low: N", "Subgroup low: AUC",
        "Subgroup low: AUC CI Lower", "Subgroup low: AUC CI Upper",
        "Subgroup high: N", "Subgroup high: AUC",
        "Subgroup high: AUC CI Lower", "Subgroup high: AUC CI Upper",
        "Note"
      ),
      Value = c(
        stratified_result$label, as.character(stratified_result$n_valid),
        paste(stratified_result$bin_breaks, collapse = " / "),
        paste(stratified_result$bin_labels, collapse = " | "),
        as.character(stratified_result$fisher_3bins_p),
        as.character(stratified_result$ca_trend_p),
        as.character(bc$threshold),
        as.character(bc$n), as.character(bc$n_high), as.character(bc$n_low),
        as.character(bc$balanced_accuracy), as.character(bc$sensitivity),
        as.character(bc$specificity), as.character(bc$fisher_p),
        as.character(bc$odds_ratio),
        as.character(sl$n),
        if (is.na(sl$auc)) "NA" else as.character(sl$auc),
        if (any(is.na(sl$auc_ci))) "NA" else as.character(sl$auc_ci[1]),
        if (any(is.na(sl$auc_ci))) "NA" else as.character(sl$auc_ci[3]),
        as.character(sh$n),
        if (is.na(sh$auc)) "NA" else as.character(sh$auc),
        if (any(is.na(sh$auc_ci))) "NA" else as.character(sh$auc_ci[1]),
        if (any(is.na(sh$auc_ci))) "NA" else as.character(sh$auc_ci[3]),
        stratified_result$note
      ),
      stringsAsFactors = FALSE
    )
    openxlsx::addWorksheet(wb, "Benchmark_Stratified")
    openxlsx::writeData(wb, "Benchmark_Stratified", df_strat_summary)
    openxlsx::addWorksheet(wb, "Benchmark_Stratified_Bins")
    openxlsx::writeData(wb, "Benchmark_Stratified_Bins", stratified_result$bin_crosstab)
  }

  if (!is.null(combined_result)) {
    lg <- combined_result$logistic
    ig <- combined_result$information_gain
    sv <- combined_result$svm_comparison
    df_comb_summary <- data.frame(
      Metric = c(
        "Benchmark biomarker", "N valid (non-NA)",
        sprintf("Logistic AUC %s alone",          combined_result$label),
        sprintf("Logistic AUC %s+gate combined",  combined_result$label),
        "Logistic delta AUC", "Logistic DeLong p",
        "IDI", "IDI 95% CI Lower", "IDI 95% CI Upper", "IDI bootstrap p",
        "cNRI", "cNRI 95% CI Lower", "cNRI 95% CI Upper", "cNRI bootstrap p",
        "Bootstrap iterations (IDI/cNRI)",
        "SVM AUC features only (same subset)",
        sprintf("SVM AUC features + %s",   combined_result$label),
        "SVM delta AUC", "SVM DeLong p",
        "Note"
      ),
      Value = c(
        combined_result$label, as.character(combined_result$n_valid),
        as.character(lg$auc_benchmark), as.character(lg$auc_combined),
        as.character(lg$delta_auc), as.character(lg$delong_p),
        as.character(ig$IDI), as.character(ig$IDI_95CI[[1]]),
        as.character(ig$IDI_95CI[[2]]), as.character(ig$IDI_bootstrap_p),
        as.character(ig$cNRI), as.character(ig$cNRI_95CI[[1]]),
        as.character(ig$cNRI_95CI[[2]]), as.character(ig$cNRI_bootstrap_p),
        as.character(ig$n_boot),
        as.character(sv$auc_features_only),
        as.character(sv$auc_features_with_bench),
        as.character(sv$delta_auc), as.character(sv$delong_p),
        combined_result$note
      ),
      stringsAsFactors = FALSE
    )
    openxlsx::addWorksheet(wb, "Benchmark_Combined")
    openxlsx::writeData(wb, "Benchmark_Combined", df_comb_summary)
  }

  if (!is.null(nested_validation)) {
    nv_m      <- nested_validation$metrics
    df_nested <- data.frame(
      Metric = c("Method", "AUC", "AUC CI Lower", "AUC CI Upper",
                 "Balanced Accuracy", "BER", "Sensitivity", "Specificity",
                 "Empty Gate Folds", "Total Folds",
                 "Note"),
      Value  = c(
        "SVM-RBF Fully-Nested LOO",
        round(nv_m$auc, 4), round(nv_m$auc_ci[1], 4), round(nv_m$auc_ci[3], 4),
        round(nv_m$balanced_accuracy, 4), round(nv_m$ber, 4),
        round(nv_m$sensitivity, 4), round(nv_m$specificity, 4),
        nested_validation$n_empty_folds, nrow(DATA$metadata),
        sprintf("LMM gate re-selected inside each fold. Fixed-gate AUC=%.4f. Delta=%.4f.",
                if (primary_method == "SVM-RBF") metrics_svm$auc else metrics_glmnet$auc,
                nv_m$auc - (if (primary_method == "SVM-RBF") metrics_svm$auc else metrics_glmnet$auc))
      ),
      stringsAsFactors = FALSE
    )
    openxlsx::addWorksheet(wb, "Nested_LOO_Validation")
    openxlsx::writeData(wb, "Nested_LOO_Validation", df_nested)

    openxlsx::addWorksheet(wb, "Gate_Stability")
    openxlsx::writeData(wb, "Gate_Stability", nested_validation$gate_stability)
  }

  excel_path <- file.path(out_dir, sprintf("ML_Classification_Report_%s.xlsx", config$project_name))
  openxlsx::saveWorkbook(wb, excel_path, overwrite = TRUE)
  message(sprintf("   [Output] Classification report saved: %s", basename(excel_path)))

  # 11. PDF Visualizations
  # ------------------------------------------------------------------------------
  results_list_plot <- list(
    `Elastic Net` = res_glmnet,
    `SVM-RBF`     = res_svm
  )

  pdf(file.path(out_dir, sprintf("ROC_ML_%s.pdf", config$project_name)), width = 8, height = 7)
  tryCatch({
    p_roc <- plot_ml_roc(
      results_list   = results_list_plot,
      colors_viz     = colors_viz,
      title          = sprintf("Nested-LOOCV ROC: %s vs %s\n(%s)",
                               config$clinical$responder_label,
                               config$clinical$non_responder_label,
                               config$project_name),
      benchmark_list = if (!is.null(benchmark_result)) list(benchmark_result) else NULL
    )
    if (!is.null(p_roc)) print(p_roc)
  }, error = function(e) warning(paste("ROC plot failed:", e$message)))
  dev.off()

  pdf(file.path(out_dir, sprintf("Predictions_ML_%s.pdf", config$project_name)), width = 10, height = 6)
  tryCatch({
    df_preds_plot <- df_preds %>%
      dplyr::mutate(Group = factor(True_Group, levels = levels(y)))
    p_pred <- plot_ml_predictions(
      predictions_df = df_preds_plot,
      colors_viz     = colors_viz,
      title          = sprintf("Out-of-Fold Predicted Probabilities: %s", config$project_name)
    )
    if (!is.null(p_pred)) print(p_pred)
  }, error = function(e) warning(paste("Predictions plot failed:", e$message)))
  dev.off()

  pdf(file.path(out_dir, sprintf("Univariate_ML_%s.pdf", config$project_name)), width = 8, height = 7)
  tryCatch({
    p_uni <- plot_univariate_roc(
      X_main         = X_main,
      y              = y,
      positive_label = res_glmnet$positive_label,
      univariate_df  = uni_auc_df,
      colors_viz     = colors_viz,
      title          = sprintf("Univariate ROC: %s vs %s\n(%s)",
                               config$clinical$responder_label,
                               config$clinical$non_responder_label,
                               config$project_name)
    )
    if (!is.null(p_uni)) print(p_uni)
  }, error = function(e) warning(paste("Univariate ROC plot failed:", e$message)))
  dev.off()

  # Benchmark-stratified figure (clinical-strata + subgroup ROC)
  if (!is.null(stratified_result)) {
    pdf(file.path(out_dir, sprintf("Benchmark_Stratified_%s.pdf", config$project_name)),
        width = 13, height = 6)
    tryCatch({
      p_strat <- plot_benchmark_stratified(
        stratified_result = stratified_result,
        colors_viz        = colors_viz,
        positive_label    = res_glmnet$positive_label,
        title             = sprintf("%s strata and gate-model subgroup performance (%s)",
                                    stratified_result$label, config$project_name)
      )
      if (!is.null(p_strat)) print(p_strat)
    }, error = function(e) warning(paste("Stratified figure failed:", e$message)))
    dev.off()
  }

  # Combined-model information-gain figure (3-ROC overlay + IDI/cNRI forest)
  if (!is.null(combined_result)) {
    pdf(file.path(out_dir, sprintf("InformationGain_%s.pdf", config$project_name)),
        width = 13, height = 6)
    tryCatch({
      p_ig <- plot_combined_information_gain(
        combined_result = combined_result,
        title           = sprintf("Information gain over %s (%s)",
                                  combined_result$label, config$project_name)
      )
      if (!is.null(p_ig)) print(p_ig)
    }, error = function(e) warning(paste("InformationGain figure failed:", e$message)))
    dev.off()
  }

  bench_line <- if (!is.null(benchmark_result)) {
    sprintf(" | %s AUC=%.3f (n=%d, DeLong p=%.4f)",
            benchmark_result$label,
            benchmark_result$auc,
            benchmark_result$n_valid,
            if (is.na(benchmark_result$delong_p)) 0 else benchmark_result$delong_p)
  } else ""

  message(sprintf(
    "\n[ML] Summary — %s | Primary: %s | AUC(EN)=%.3f | AUC(SVM)=%.3f | Perm p(SVM)=%.4f%s",
    config$project_name, primary_method,
    metrics_glmnet$auc, metrics_svm$auc, perm_svm$p_value,
    bench_line
  ))
}

message("=== STEP 6 COMPLETE ===\n")
