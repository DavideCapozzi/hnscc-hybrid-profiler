# diagnostics/diag_09_gate_composite.R
# ==============================================================================
# DIAGNOSTIC 9: Simple gate-averaged KI67 composite score
#
# Hypothesis: a 1-feature classifier built as the mean z-score of the LMM gate
# markers (KI67NAIVE + KI67CD4) may match the 2-feature SVM-RBF baseline
# (AUC=0.716) while simplifying the narrative ("1 composite biomarker").
#
# diag_02 already showed that PC1 of the *full* KI67 family is worse (AUC=0.28)
# because it absorbs noise from CD137KI67/PD1KI67 (no LMM signal). This test
# restricts the composite to the FDR+LOO-robust gate markers only.
#
# Output: diagnostics/diag_09_output.txt
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(pROC)
  library(e1071)
})

OUT_FILE <- here("diagnostics/diag_09_output.txt")
con <- file(OUT_FILE, open = "wt")
sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)

# ── Load cross-sectional T0 data ─────────────────────────────────────────────
obj_std <- readRDS(here(
  "results/BestResponse_2v3_4/01_data_processing/data_processed_BestResponse_2v3_4_standard.rds"
))
data_z <- obj_std$hybrid_data_z
meta   <- obj_std$metadata

cat(sprintf("Cross-sectional T0: %d patients x %d markers\n", nrow(meta), ncol(data_z)))

GATE <- c("KI67NAIVE", "KI67CD4")
stopifnot(all(GATE %in% colnames(data_z)))

y <- ifelse(meta$Group == "RP", 1, 0)
cat(sprintf("Class distribution — RP: %d  |  SD_PD: %d\n\n", sum(y == 1), sum(y == 0)))

# ── Composite: mean of gate marker z-scores (already z-scored at Step 01) ────
composite_score <- rowMeans(as.matrix(data_z[, GATE]), na.rm = FALSE)

cat("=== Composite definition ===\n")
cat("  score_i = mean(z(KI67NAIVE_i), z(KI67CD4_i))\n")
cat(sprintf("  range = [%.2f, %.2f]\n", min(composite_score), max(composite_score)))
cat(sprintf("  cor(KI67NAIVE, KI67CD4) = %.3f\n\n",
            cor(data_z[, "KI67NAIVE"], data_z[, "KI67CD4"], use = "complete.obs")))

# ── Model A: ROC of composite alone (direct, no model fit needed) ────────────
cat("=== Model A: Composite score, direct ROC (no model) ===\n")
roc_comp <- pROC::roc(y, composite_score, direction = ">", quiet = TRUE)
ci_comp  <- pROC::ci.auc(roc_comp, method = "bootstrap", boot.n = 2000,
                         conf.level = 0.95, progress = "none")
thr_comp <- pROC::coords(roc_comp, "best", ret = c("threshold", "sensitivity", "specificity"))
pred_comp <- as.integer(composite_score <= thr_comp$threshold[1])  # lower score → RP
bal_acc_comp <- mean(c(mean(pred_comp[y == 1] == 1), mean(pred_comp[y == 0] == 0)))
cat(sprintf("  AUC=%.4f [%.3f-%.3f]  BalAcc=%.4f  Sens=%.3f  Spec=%.3f\n\n",
            as.numeric(roc_comp$auc), ci_comp[1], ci_comp[3],
            bal_acc_comp, thr_comp$sensitivity[1], thr_comp$specificity[1]))

# ── Helper: nested-LOOCV SVM-RBF ────────────────────────────────────────────
run_svm_loocv <- function(X, y, label = "model") {
  n <- nrow(X)
  probs <- numeric(n)
  C_grid <- c(0.01, 0.1, 1, 10, 100)
  gamma_grid <- c(0.01, 0.1, 1, 10)

  for (i in seq_len(n)) {
    X_train <- X[-i, , drop = FALSE]; X_test <- X[i, , drop = FALSE]
    y_train <- y[-i]; y_test <- y[i]

    mu  <- colMeans(X_train, na.rm = TRUE)
    sds <- apply(X_train, 2, sd, na.rm = TRUE); sds[sds == 0] <- 1
    X_train_s <- scale(X_train, center = mu, scale = sds)
    X_test_s  <- scale(X_test,  center = mu, scale = sds)

    set.seed(2026 + i)
    fold_ids <- sample(rep(1:5, length.out = nrow(X_train_s)))
    best_C <- 1; best_gamma <- 0.1; best_acc <- -Inf
    for (C in C_grid) for (gam in gamma_grid) {
      acc_inner <- 0
      for (f in 1:5) {
        xtr <- X_train_s[fold_ids != f, , drop = FALSE]
        xv  <- X_train_s[fold_ids == f, , drop = FALSE]
        ytr <- y_train[fold_ids != f]; yv <- y_train[fold_ids == f]
        if (length(unique(ytr)) < 2) next
        fit_in <- tryCatch(
          svm(xtr, as.factor(ytr), kernel = "radial", cost = C, gamma = gam, probability = TRUE),
          error = function(e) NULL)
        if (is.null(fit_in)) next
        pv <- predict(fit_in, xv)
        acc_inner <- acc_inner + mean(pv == as.factor(yv))
      }
      acc_inner <- acc_inner / 5
      if (acc_inner > best_acc) { best_acc <- acc_inner; best_C <- C; best_gamma <- gam }
    }
    fit <- tryCatch(
      svm(X_train_s, as.factor(y_train), kernel = "radial",
          cost = best_C, gamma = best_gamma, probability = TRUE),
      error = function(e) NULL)
    if (is.null(fit)) { probs[i] <- 0.5; next }
    pred_obj <- predict(fit, X_test_s, probability = TRUE)
    prob_mat <- attr(pred_obj, "probabilities")
    probs[i] <- if ("1" %in% colnames(prob_mat)) prob_mat[, "1"] else 0.5
  }

  roc_obj <- pROC::roc(y, probs, direction = "<", quiet = TRUE)
  ci_obj  <- pROC::ci.auc(roc_obj, method = "bootstrap", boot.n = 2000,
                          conf.level = 0.95, progress = "none")
  thr <- pROC::coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"))
  pred_bin <- as.integer(probs >= thr$threshold[1])
  bal_acc <- mean(c(mean(pred_bin[y == 1] == 1), mean(pred_bin[y == 0] == 0)))
  cat(sprintf("  %-30s AUC=%.4f [%.3f-%.3f]  BalAcc=%.4f  Sens=%.3f  Spec=%.3f\n",
              label, as.numeric(roc_obj$auc), ci_obj[1], ci_obj[3],
              bal_acc, thr$sensitivity[1], thr$specificity[1]))
  invisible(list(auc = as.numeric(roc_obj$auc), ci = as.numeric(ci_obj),
                 bal_acc = bal_acc, probs = probs, roc = roc_obj))
}

# ── Model B: SVM-RBF on composite (1 feature) ────────────────────────────────
cat("=== Model B: SVM-RBF on composite (1 feature, nested LOOCV) ===\n")
res_comp_svm <- run_svm_loocv(matrix(composite_score, ncol = 1,
                                     dimnames = list(NULL, "KI67_gate_avg")),
                              y, "Composite (1 feat)")

# ── Model C: Baseline — 2 separate features ──────────────────────────────────
cat("\n=== Model C: SVM-RBF on KI67NAIVE + KI67CD4 (current baseline) ===\n")
X_base <- as.matrix(data_z[, GATE])
res_base <- run_svm_loocv(X_base, y, "Baseline (2 feats)")

# ── Model D: Logistic regression on composite (linear, simpler) ──────────────
cat("\n=== Model D: Logistic regression on composite (LOOCV) ===\n")
df_lr <- data.frame(y = y, score = composite_score)
probs_lr <- numeric(length(y))
for (i in seq_along(y)) {
  fit_lr <- glm(y ~ score, data = df_lr[-i, ], family = binomial())
  probs_lr[i] <- predict(fit_lr, newdata = df_lr[i, , drop = FALSE], type = "response")
}
roc_lr <- pROC::roc(y, probs_lr, direction = "<", quiet = TRUE)
ci_lr  <- pROC::ci.auc(roc_lr, method = "bootstrap", boot.n = 2000,
                       conf.level = 0.95, progress = "none")
cat(sprintf("  Logistic (composite, LOOCV)    AUC=%.4f [%.3f-%.3f]\n",
            as.numeric(roc_lr$auc), ci_lr[1], ci_lr[3]))

# ── DeLong tests ─────────────────────────────────────────────────────────────
cat("\n=== DeLong tests: composite vs current baseline ===\n")
for (rcomp in list(list(r = res_comp_svm, l = "SVM (composite)"),
                   list(r = list(roc = roc_comp), l = "ROC composite (raw)"),
                   list(r = list(roc = roc_lr),  l = "Logistic (composite)"))) {
  tryCatch({
    dt <- pROC::roc.test(res_base$roc, rcomp$r$roc, method = "delong", paired = TRUE)
    cat(sprintf("  Baseline (2 feats) vs %-22s  DeLong p=%.4f\n", rcomp$l, dt$p.value))
  }, error = function(e) cat(sprintf("  DeLong failed: %s\n", rcomp$l)))
}

# ── Interpretation ──────────────────────────────────────────────────────────
cat("\n=== Interpretation ===\n")
cat("If composite-1feat AUC >= baseline-2feats AUC - 0.02, the composite is a\n")
cat("narratively cleaner choice (single biomarker, simpler clinical adoption).\n")
cat("If AUC drops > 0.05, keep the 2-marker baseline (SVM-RBF non-linear boundary).\n")
cat("\nDone.\n")
