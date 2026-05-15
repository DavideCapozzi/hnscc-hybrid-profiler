# diagnostics/diag_03_nested_lmm_loocv.R
# ==============================================================================
# DIAGNOSTIC TEST 3: Fully-nested LOO — LMM feature selection inside the loop
#
# For each outer LOO fold (leave out patient i):
#   1. Run LMM on the remaining n-1 patients (longitudinal data)
#   2. BH-correct and select FDR < fdr_thr markers
#   3. Run LOO sensitivity check on those markers (still on n-1 patients)
#   4. Apply collinearity filter (|r| > 0.85 → drop higher mean-|r| marker)
#   5. If >= 1 feature survives: train SVM-RBF on T0 of n-1 patients → predict patient i
#   6. If 0 features: predict 0.5 (uninformative)
#
# Reports: AUC distribution, how often each marker appears as a gate marker,
#          comparison with current AUC=0.716 (fixed gate).
#
# NOTE: This is slow (75 outer folds × 74 inner LOO folds × n_markers LMMs).
#       Estimated runtime: 30–60 min. Progress is reported every 5 folds.
#
# Usage: Rscript diagnostics/diag_03_nested_lmm_loocv.R
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(lmerTest)
  library(lme4)
  library(dplyr)
  library(e1071)
  library(pROC)
})

# ── Parameters ────────────────────────────────────────────────────────────────
fdr_thr  <- 0.05
loo_thr  <- 0.05
cor_thr  <- 0.85
seed     <- 2026

# ── Load data ─────────────────────────────────────────────────────────────────
cat("Loading data...\n")
obj_std  <- readRDS(here("results/BestResponse_2v3_4/01_data_processing/data_processed_BestResponse_2v3_4_standard.rds"))
obj_long <- readRDS(here("results/BestResponse_2v3_4/01_data_processing/data_processed_BestResponse_2v3_4_longitudinal.rds"))

# Cross-sectional T0 data (for classifier)
# hybrid_data_z contains meta columns (Patient_ID, Sample_ID, Timepoint, Group)
# plus marker columns — extract only markers; use intersection of both datasets
META_COLS <- c("Patient_ID", "Sample_ID", "Timepoint", "Group")

dz_std_full <- obj_std$hybrid_data_z
dz_lng_full <- obj_long$hybrid_data_z
meta_t0     <- obj_std$metadata
meta_lng    <- obj_long$metadata

std_markers <- setdiff(colnames(dz_std_full), META_COLS)
lng_markers <- setdiff(colnames(dz_lng_full), META_COLS)
all_markers <- intersect(std_markers, lng_markers)   # features present in both

data_z_t0   <- as.matrix(dz_std_full[, all_markers, drop = FALSE])
y_all       <- ifelse(meta_t0$Group == "RP", 1, 0)
pid_t0      <- meta_t0$Patient_ID
n_t0        <- nrow(meta_t0)

# Longitudinal data (for LMM feature selection)
data_z_lng  <- as.data.frame(dz_lng_full[, all_markers, drop = FALSE])

cat(sprintf("  T0 classifier: %d patients x %d markers\n", n_t0, length(all_markers)))
cat(sprintf("  Longitudinal:  %d obs x %d markers (%d unique patients)\n\n",
            nrow(meta_lng), ncol(data_z_lng), length(unique(meta_lng$Patient_ID))))

# ── Helper: collinearity filter ───────────────────────────────────────────────
filter_collinear <- function(X, markers, threshold = 0.85) {
  if (length(markers) <= 1) return(markers)
  sub <- X[, markers, drop = FALSE]
  if (any(is.na(sub))) sub <- sub[complete.cases(sub), ]
  if (nrow(sub) < 5) return(markers)
  cm <- cor(sub, use = "pairwise.complete.obs")
  diag(cm) <- 0
  mean_abs_r <- rowMeans(abs(cm))
  kept <- markers
  repeat {
    sub_cm <- abs(cor(X[, kept, drop = FALSE], use = "pairwise.complete.obs"))
    diag(sub_cm) <- 0
    if (max(sub_cm, na.rm = TRUE) <= threshold) break
    pair_idx <- which(sub_cm == max(sub_cm, na.rm = TRUE), arr.ind = TRUE)[1, ]
    drop_idx <- which.max(mean_abs_r[kept])
    kept <- kept[-drop_idx]
    if (length(kept) == 0) break
  }
  kept
}

# ── Helper: LMM feature selection on a training subset ───────────────────────
select_lmm_features <- function(df_lng_sub, markers, fdr_thr, loo_thr) {
  results <- lapply(markers, function(mk) {
    vals <- as.numeric(df_lng_sub[[mk]])
    val_sd <- sd(vals, na.rm = TRUE)
    if (is.na(val_sd) || val_sd == 0) val_sd <- 1

    df_m <- data.frame(
      Value = vals / val_sd,
      Time  = as.factor(df_lng_sub$Timepoint),
      Group = as.factor(df_lng_sub$Group),
      ID    = as.factor(df_lng_sub$Patient_ID)
    )
    df_m <- df_m[complete.cases(df_m), ]
    if (nrow(df_m) < 10 || length(unique(df_m$Group)) < 2) return(data.frame(Marker=mk, P=NA))

    tryCatch({
      mod <- suppressMessages(suppressWarnings(
        lmerTest::lmer(Value ~ Time * Group + (1 | ID), data = df_m,
                       REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE))
      ))
      ct  <- summary(mod)$coefficients
      idx <- grep("Time.*:Group", rownames(ct))
      if (length(idx) != 1) return(data.frame(Marker=mk, P=NA))
      p_col <- grep("Pr\\(>\\|t\\|\\)", colnames(ct))
      data.frame(Marker = mk, P = ct[idx, p_col])
    }, error = function(e) data.frame(Marker = mk, P = NA))
  })

  res_df <- bind_rows(results)
  res_df$FDR <- p.adjust(res_df$P, method = "BH")
  fdr_pass <- res_df$Marker[!is.na(res_df$FDR) & res_df$FDR < fdr_thr]
  if (length(fdr_pass) == 0) return(character(0))

  # LOO sensitivity on FDR-passing markers
  loo_pass <- c()
  for (mk in fdr_pass) {
    unique_ids <- unique(df_lng_sub$Patient_ID)
    max_p <- 0
    for (drop_id in unique_ids) {
      sub2 <- df_lng_sub[df_lng_sub$Patient_ID != drop_id, ]
      vals2 <- as.numeric(sub2[[mk]])
      val_sd2 <- sd(vals2, na.rm = TRUE)
      if (is.na(val_sd2) || val_sd2 == 0) val_sd2 <- 1
      df_m2 <- data.frame(
        Value = vals2 / val_sd2,
        Time  = as.factor(sub2$Timepoint),
        Group = as.factor(sub2$Group),
        ID    = as.factor(sub2$Patient_ID)
      )
      df_m2 <- df_m2[complete.cases(df_m2), ]
      if (nrow(df_m2) < 8 || length(unique(df_m2$Group)) < 2) next
      res_loo <- tryCatch({
        mod2 <- suppressMessages(suppressWarnings(
          lmerTest::lmer(Value ~ Time * Group + (1 | ID), data = df_m2,
                         REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE))
        ))
        ct2  <- summary(mod2)$coefficients
        idx2 <- grep("Time.*:Group", rownames(ct2))
        if (length(idx2) != 1) NA else ct2[idx2, grep("Pr\\(>\\|t\\|\\)", colnames(ct2))]
      }, error = function(e) NA)
      if (!is.na(res_loo) && res_loo > max_p) max_p <- res_loo
    }
    if (max_p < loo_thr) loo_pass <- c(loo_pass, mk)
  }
  loo_pass
}

# ── Helper: train SVM-RBF inner CV ───────────────────────────────────────────
svm_inner_cv <- function(X_train, y_train) {
  C_grid     <- c(0.01, 0.1, 1, 10, 100)
  gamma_grid <- c(0.01, 0.1, 1, 10)
  inner_folds <- 5
  if (length(unique(y_train)) < 2) return(list(C = 1, gamma = 0.1))

  set.seed(seed)
  fold_ids <- sample(rep(1:inner_folds, length.out = nrow(X_train)))
  best_C <- 1; best_gamma <- 0.1; best_acc <- -Inf

  for (C in C_grid) {
    for (gam in gamma_grid) {
      acc_inner <- 0; n_inner <- 0
      for (f in 1:inner_folds) {
        xtr <- X_train[fold_ids != f, , drop = FALSE]
        xv  <- X_train[fold_ids == f, , drop = FALSE]
        ytr <- y_train[fold_ids != f]
        yv  <- y_train[fold_ids == f]
        if (length(unique(ytr)) < 2) next
        fit_i <- tryCatch(
          svm(xtr, as.factor(ytr), kernel = "radial", cost = C, gamma = gam, probability = TRUE),
          error = function(e) NULL
        )
        if (is.null(fit_i)) next
        pv <- predict(fit_i, xv)
        acc_inner <- acc_inner + mean(pv == as.factor(yv))
        n_inner <- n_inner + 1
      }
      if (n_inner == 0) next
      acc_avg <- acc_inner / n_inner
      if (acc_avg > best_acc) { best_acc <- acc_avg; best_C <- C; best_gamma <- gam }
    }
  }
  list(C = best_C, gamma = best_gamma)
}

# ── Outer LOO loop ────────────────────────────────────────────────────────────
cat("Starting fully-nested LOO (this may take 30-60 min)...\n")
cat(sprintf("  Outer folds: %d patients\n", n_t0))
cat(sprintf("  FDR threshold: %.2f  |  LOO threshold: %.2f  |  Collinearity: %.2f\n\n",
            fdr_thr, loo_thr, cor_thr))

probs_nested  <- numeric(n_t0)
gate_markers_per_fold <- vector("list", n_t0)

t_start <- proc.time()["elapsed"]

for (i in seq_len(n_t0)) {
  pid_test <- pid_t0[i]

  # T0 training data (n-1 patients)
  idx_train_t0 <- setdiff(seq_len(n_t0), i)
  X_train_t0   <- as.matrix(data_z_t0[idx_train_t0, , drop = FALSE])
  y_train       <- y_all[idx_train_t0]
  X_test_t0     <- as.matrix(data_z_t0[i, , drop = FALSE])

  # Longitudinal training data (same n-1 patients, both timepoints)
  # data_z_lng contains only marker columns; combine with meta_lng for LMM
  lng_mask     <- meta_lng$Patient_ID != pid_test
  df_lng_train <- cbind(
    as.data.frame(meta_lng[lng_mask, ]),
    as.data.frame(data_z_lng[lng_mask, , drop = FALSE])
  )

  # Step A: LMM feature selection on longitudinal n-1
  gate <- select_lmm_features(df_lng_train, all_markers, fdr_thr, loo_thr)

  # Step B: Collinearity filter
  if (length(gate) > 1) {
    gate <- filter_collinear(X_train_t0, gate, threshold = cor_thr)
  }

  gate_markers_per_fold[[i]] <- gate

  if (length(gate) == 0) {
    probs_nested[i] <- 0.5
    if (i %% 5 == 0 || i == n_t0) {
      elapsed <- proc.time()["elapsed"] - t_start
      cat(sprintf("  Fold %2d/%d  gate=EMPTY  [%.0fs elapsed]\n", i, n_t0, elapsed))
    }
    next
  }

  # Step C: Train SVM-RBF on T0 training data with gate features
  X_tr  <- X_train_t0[, gate, drop = FALSE]
  X_te  <- X_test_t0[,  gate, drop = FALSE]

  # Z-score on training fold
  mu   <- colMeans(X_tr, na.rm = TRUE)
  sds  <- apply(X_tr, 2, sd, na.rm = TRUE)
  sds[sds == 0] <- 1
  X_tr_s <- scale(X_tr, center = mu, scale = sds)
  X_te_s <- scale(X_te, center = mu, scale = sds)

  # Inner CV
  best_params <- svm_inner_cv(X_tr_s, y_train)

  fit <- tryCatch(
    svm(X_tr_s, as.factor(y_train), kernel = "radial",
        cost = best_params$C, gamma = best_params$gamma, probability = TRUE),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    probs_nested[i] <- 0.5
  } else {
    pred_obj <- predict(fit, X_te_s, probability = TRUE)
    pm <- attr(pred_obj, "probabilities")
    probs_nested[i] <- if ("1" %in% colnames(pm)) pm[, "1"] else 0.5
  }

  if (i %% 5 == 0 || i == n_t0) {
    elapsed <- proc.time()["elapsed"] - t_start
    cat(sprintf("  Fold %2d/%d  gate=[%s]  [%.0fs elapsed]\n",
                i, n_t0, paste(gate, collapse = ","), elapsed))
  }
}

# ── Results ───────────────────────────────────────────────────────────────────
cat("\n", strrep("=", 70), "\n")
cat("RESULTS — Fully-nested LOO (LMM feature selection inside loop)\n")
cat(strrep("=", 70), "\n\n")

roc_nested <- pROC::roc(y_all, probs_nested, direction = "<", quiet = TRUE)
ci_nested  <- pROC::ci.auc(roc_nested, method = "bootstrap", boot.n = 2000,
                            conf.level = 0.95, progress = "none")
thr_nested <- pROC::coords(roc_nested, "best", ret = c("threshold","sensitivity","specificity"))
pred_bin   <- as.integer(probs_nested >= thr_nested$threshold[1])
bal_acc    <- mean(c(mean(pred_bin[y_all == 1] == 1), mean(pred_bin[y_all == 0] == 0)))

cat(sprintf("  AUC        : %.4f [%.3f–%.3f]\n",
            as.numeric(roc_nested$auc), as.numeric(ci_nested[1]), as.numeric(ci_nested[3])))
cat(sprintf("  Bal. Acc   : %.4f\n", bal_acc))
cat(sprintf("  Sensitivity: %.4f\n", thr_nested$sensitivity[1]))
cat(sprintf("  Specificity: %.4f\n\n", thr_nested$specificity[1]))

# Gate marker frequency
cat("Gate marker selection frequency across folds:\n")
all_gates <- unlist(gate_markers_per_fold)
freq_tab  <- sort(table(all_gates), decreasing = TRUE)
for (nm in names(freq_tab)) {
  cat(sprintf("  %-16s  %3d / %d folds (%.0f%%)\n",
              nm, freq_tab[[nm]], n_t0, 100 * freq_tab[[nm]] / n_t0))
}

n_empty <- sum(sapply(gate_markers_per_fold, length) == 0)
cat(sprintf("\n  Folds with empty gate: %d / %d (%.0f%%)\n", n_empty, n_t0, 100*n_empty/n_t0))

cat(sprintf("\n  Reference (fixed gate, current pipeline): AUC=0.7158 [0.592–0.839]\n"))
cat(sprintf("  Difference: %.4f\n",
            as.numeric(roc_nested$auc) - 0.7158))

# DeLong test: nested vs reference (we approximate with fixed-gate probs from JSON)
cat("\nNote: DeLong test vs fixed-gate requires re-running pipeline with saved probs.\n")
cat("      The AUC difference above gives the unbiased leakage correction estimate.\n")

# Save probs for potential follow-up DeLong
saveRDS(list(probs = probs_nested, y = y_all, gate_per_fold = gate_markers_per_fold),
        here("diagnostics/diag_03_nested_probs.rds"))
cat("\nPredicted probabilities saved to diagnostics/diag_03_nested_probs.rds\n")
cat("Done.\n")
