# diagnostics/diag_11_perfold_zscore.R
# ==============================================================================
# DIAGNOSTIC 11: Per-fold z-scoring vs global z-scoring (leakage quantification)
#
# The current pipeline computes z-scores once in Step 01 over ALL patients
# (including the held-out fold during nested LOOCV). This is a known minor
# leakage source documented in the Step 06 JSON.
#
# This test refits the SVM-RBF classifier with z-scoring done strictly INSIDE
# each outer LOO fold (training-only) and compares AUC to the current pipeline
# value (0.716). If ΔAUC < 0.02, the global-z caveat can be defensibly retained;
# if > 0.02, per-fold scaling should be integrated in Step 06.
#
# Note: the diag_03 helper already z-scores per-fold inside the SVM-RBF inner
# loop using TRAINING fold statistics. The remaining question is whether the
# Step 01 *raw input* (logit/log2 transformed but z-scored globally) carries
# residual information from the held-out fold.
#
# Implementation: load the RAW (non-z) hybrid matrix, z-score on training fold
# only, then run nested LOOCV identical to Step 06.
#
# Output: diagnostics/diag_11_output.txt
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(pROC)
  library(e1071)
})

OUT_FILE <- here("diagnostics/diag_11_output.txt")
con <- file(OUT_FILE, open = "wt")
sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)

# ── Load cross-sectional T0 data (RAW hybrid scale, pre-z-score) ─────────────
obj_std <- readRDS(here(
  "results/BestResponse_2v3_4/01_data_processing/data_processed_BestResponse_2v3_4_standard.rds"
))
data_raw <- obj_std$hybrid_data_raw   # logit/log2 transformed, NOT z-scored
data_z   <- obj_std$hybrid_data_z     # globally z-scored (current pipeline)
meta     <- obj_std$metadata

cat(sprintf("Cross-sectional T0: %d patients\n", nrow(meta)))
cat(sprintf("hybrid_data_raw cols: %d  |  hybrid_data_z cols: %d\n",
            ncol(data_raw), ncol(data_z)))

GATE <- c("KI67NAIVE", "KI67CD4")
y <- ifelse(meta$Group == "RP", 1, 0)

stopifnot(all(GATE %in% colnames(data_raw)))
stopifnot(all(GATE %in% colnames(data_z)))

X_raw <- as.matrix(data_raw[, GATE])
X_z   <- as.matrix(data_z[,   GATE])

cat(sprintf("Class distribution: RP=%d  SD_PD=%d\n\n", sum(y == 1), sum(y == 0)))

# ── Helper: nested LOOCV with configurable scaling source ───────────────────
# scale_source: "raw" (recompute z inside each fold) or "z" (use pre-z input,
# still re-center on training fold to match the current diag_03 behavior).
run_svm_loocv <- function(X, y, scale_source = "raw", label = "model") {
  n <- nrow(X)
  probs <- numeric(n)
  C_grid <- c(0.01, 0.1, 1, 10, 100)
  gamma_grid <- c(0.01, 0.1, 1, 10)

  for (i in seq_len(n)) {
    X_train <- X[-i, , drop = FALSE]; X_test <- X[i, , drop = FALSE]
    y_train <- y[-i]

    # Z-score using TRAINING-FOLD statistics only (both branches do this;
    # only the input differs: raw transformed vs globally-z input).
    mu  <- colMeans(X_train, na.rm = TRUE)
    sds <- apply(X_train, 2, sd, na.rm = TRUE); sds[sds == 0] <- 1
    X_train_s <- scale(X_train, center = mu, scale = sds)
    X_test_s  <- scale(X_test,  center = mu, scale = sds)

    set.seed(2026 + i)
    fold_ids <- sample(rep(1:5, length.out = nrow(X_train_s)))
    best_C <- 1; best_gamma <- 0.1; best_acc <- -Inf
    for (C in C_grid) for (gam in gamma_grid) {
      acc <- 0
      for (f in 1:5) {
        xtr <- X_train_s[fold_ids != f, , drop = FALSE]
        xv  <- X_train_s[fold_ids == f, , drop = FALSE]
        ytr <- y_train[fold_ids != f]; yv <- y_train[fold_ids == f]
        if (length(unique(ytr)) < 2) next
        fit <- tryCatch(svm(xtr, as.factor(ytr), kernel = "radial",
                            cost = C, gamma = gam, probability = TRUE),
                        error = function(e) NULL)
        if (is.null(fit)) next
        acc <- acc + mean(predict(fit, xv) == as.factor(yv))
      }
      acc <- acc / 5
      if (acc > best_acc) { best_acc <- acc; best_C <- C; best_gamma <- gam }
    }
    fit <- tryCatch(svm(X_train_s, as.factor(y_train), kernel = "radial",
                        cost = best_C, gamma = best_gamma, probability = TRUE),
                    error = function(e) NULL)
    if (is.null(fit)) { probs[i] <- 0.5; next }
    pr_mat <- attr(predict(fit, X_test_s, probability = TRUE), "probabilities")
    probs[i] <- if ("1" %in% colnames(pr_mat)) pr_mat[, "1"] else 0.5
  }

  roc_obj <- pROC::roc(y, probs, direction = "<", quiet = TRUE)
  ci_obj  <- pROC::ci.auc(roc_obj, method = "bootstrap", boot.n = 2000,
                          conf.level = 0.95, progress = "none")
  thr <- pROC::coords(roc_obj, "best", ret = c("threshold","sensitivity","specificity"))
  pred_bin <- as.integer(probs >= thr$threshold[1])
  bal_acc <- mean(c(mean(pred_bin[y == 1] == 1), mean(pred_bin[y == 0] == 0)))
  cat(sprintf("  %-35s AUC=%.4f [%.3f-%.3f]  BalAcc=%.4f  Sens=%.3f  Spec=%.3f\n",
              label, as.numeric(roc_obj$auc), ci_obj[1], ci_obj[3],
              bal_acc, thr$sensitivity[1], thr$specificity[1]))
  invisible(list(roc = roc_obj, auc = as.numeric(roc_obj$auc),
                 ci = as.numeric(ci_obj), probs = probs))
}

# ── Model A: Per-fold z-scoring from RAW transformed scale ───────────────────
cat("=== Model A: Per-fold z-scoring from RAW hybrid (no global z leakage) ===\n")
res_raw <- run_svm_loocv(X_raw, y, scale_source = "raw",
                         label = "Per-fold-z from RAW")

# ── Model B: Re-center pre-z input on training fold (matches diag_03/Step 06)
cat("\n=== Model B: Globally-z input + per-fold re-center (current pipeline) ===\n")
res_z <- run_svm_loocv(X_z, y, scale_source = "z",
                       label = "Globally-z + per-fold")

# ── DeLong test ─────────────────────────────────────────────────────────────
cat("\n=== DeLong test: leakage impact ===\n")
dt <- pROC::roc.test(res_z$roc, res_raw$roc, method = "delong", paired = TRUE)
cat(sprintf("  Current pipeline vs Per-fold-raw: ΔAUC=%+.4f  DeLong p=%.4f\n",
            res_raw$auc - res_z$auc, dt$p.value))

cat("\n=== Interpretation ===\n")
cat("- Model A: STRICT (no global z) — refits z on training fold from raw input\n")
cat("- Model B: Current pipeline behavior on cached z-scored input\n")
cat("- |ΔAUC| < 0.02: global-z leakage is negligible; current pipeline is defensible\n")
cat("- |ΔAUC| > 0.02 or DeLong p<0.05: integrate per-fold scaling in Step 06\n")

cat("\nDone.\n")
