# diagnostics/diag_02_ki67_composite.R
# ==============================================================================
# DIAGNOSTIC TEST 2: KI67 composite score (PC1) as classifier feature
#
# Computes PC1 of the KI67 subfamily at T0, then runs nested-LOOCV with:
#   (a) PC1 alone
#   (b) PC1 + KI67NAIVE (best univariate)
# and compares AUC to the current 2-feature SVM-RBF (0.716).
#
# The PC1 is computed on training folds only to avoid leakage.
#
# Usage: Rscript diagnostics/diag_02_ki67_composite.R
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(e1071)
  library(pROC)
})

setwd("/home/laboratorio/projects/clinical-onco-profiler")

# ── Load cross-sectional T0 data ─────────────────────────────────────────────
obj_std <- readRDS("results/BestResponse_2v3_4/01_data_processing/data_processed_BestResponse_2v3_4_standard.rds")
data_z  <- obj_std$hybrid_data_z
meta    <- obj_std$metadata

cat(sprintf("Cross-sectional T0: %d patients x %d markers\n\n", nrow(meta), ncol(data_z)))

# KI67 subfamily — markers with KI67 in name or known KI67 co-expressed
ki67_markers <- c("KI67", "KI67CD4", "KI67CD8", "KI67CM", "KI67EFF",
                  "KI657EMRA", "KI67NAIVE", "CD28KI67", "PD1KI67", "CD137KI67")
ki67_markers <- intersect(ki67_markers, colnames(data_z))
cat("KI67 subfamily markers used for PCA:\n")
cat(" ", paste(ki67_markers, collapse = ", "), "\n\n")

# Binary outcome
y <- ifelse(meta$Group == "RP", 1, 0)
cat(sprintf("Class distribution — RP: %d  |  SD_PD: %d\n\n", sum(y==1), sum(y==0)))

# ── Helper: SVM-RBF nested LOOCV with inner CV ───────────────────────────────
run_svm_loocv <- function(X, y, label = "model") {
  n <- nrow(X)
  preds <- numeric(n)
  probs <- numeric(n)

  C_grid     <- c(0.01, 0.1, 1, 10, 100)
  gamma_grid <- c(0.01, 0.1, 1, 10)

  for (i in seq_len(n)) {
    X_train <- X[-i, , drop = FALSE]
    X_test  <- X[i,  , drop = FALSE]
    y_train <- y[-i]
    y_test  <- y[i]

    # Z-score on training fold
    mu  <- colMeans(X_train, na.rm = TRUE)
    sds <- apply(X_train, 2, sd, na.rm = TRUE)
    sds[sds == 0] <- 1
    X_train_s <- scale(X_train, center = mu, scale = sds)
    X_test_s  <- scale(X_test,  center = mu, scale = sds)

    # Inner 5-fold CV for C and gamma
    set.seed(2026 + i)
    inner_folds <- 5
    fold_ids <- sample(rep(1:inner_folds, length.out = nrow(X_train_s)))
    best_C <- 1; best_gamma <- 0.1; best_acc <- -Inf
    for (C in C_grid) {
      for (gam in gamma_grid) {
        acc_inner <- 0
        for (f in 1:inner_folds) {
          xtr <- X_train_s[fold_ids != f, , drop = FALSE]
          xv  <- X_train_s[fold_ids == f, , drop = FALSE]
          ytr <- y_train[fold_ids != f]
          yv  <- y_train[fold_ids == f]
          if (length(unique(ytr)) < 2) next
          fit_inner <- tryCatch(
            svm(xtr, as.factor(ytr), kernel = "radial", cost = C, gamma = gam,
                probability = TRUE),
            error = function(e) NULL
          )
          if (is.null(fit_inner)) next
          pv <- predict(fit_inner, xv)
          acc_inner <- acc_inner + mean(pv == as.factor(yv))
        }
        acc_inner <- acc_inner / inner_folds
        if (acc_inner > best_acc) {
          best_acc <- acc_inner; best_C <- C; best_gamma <- gam
        }
      }
    }

    # Train final model on full training fold
    fit <- tryCatch(
      svm(X_train_s, as.factor(y_train), kernel = "radial",
          cost = best_C, gamma = best_gamma, probability = TRUE),
      error = function(e) NULL
    )
    if (is.null(fit)) { probs[i] <- 0.5; next }

    pred_obj <- predict(fit, X_test_s, probability = TRUE)
    prob_mat <- attr(pred_obj, "probabilities")
    probs[i] <- if ("1" %in% colnames(prob_mat)) prob_mat[, "1"] else 0.5
  }

  roc_obj <- pROC::roc(y, probs, direction = "<", quiet = TRUE)
  ci_obj  <- pROC::ci.auc(roc_obj, method = "bootstrap", boot.n = 2000,
                           conf.level = 0.95, progress = "none")

  # Optimal threshold (Youden)
  thr_obj  <- pROC::coords(roc_obj, "best", ret = c("threshold","sensitivity","specificity"))
  pred_bin <- as.integer(probs >= thr_obj$threshold[1])
  bal_acc  <- mean(c(
    mean(pred_bin[y == 1] == 1),
    mean(pred_bin[y == 0] == 0)
  ))

  cat(sprintf("\n%-35s  AUC=%.4f [%.3f–%.3f]  BalAcc=%.4f  Sens=%.3f  Spec=%.3f\n",
              label,
              as.numeric(roc_obj$auc),
              as.numeric(ci_obj[1]),
              as.numeric(ci_obj[3]),
              bal_acc,
              thr_obj$sensitivity[1],
              thr_obj$specificity[1]))
  invisible(list(auc = as.numeric(roc_obj$auc), ci = as.numeric(ci_obj),
                 bal_acc = bal_acc, probs = probs, roc = roc_obj))
}

# ── Model A: Current baseline (KI67NAIVE + KI67CD4) ─────────────────────────
cat("=== Model A: Current baseline (KI67NAIVE + KI67CD4) ===\n")
X_baseline <- as.matrix(data_z[, c("KI67NAIVE", "KI67CD4")])
res_baseline <- run_svm_loocv(X_baseline, y, label = "KI67NAIVE + KI67CD4 (current)")

# ── Model B: PC1 of KI67 family (LOO-clean PCA) ─────────────────────────────
cat("\n=== Model B: PC1 of full KI67 family (LOO-clean PCA) ===\n")
cat("Note: PCA fitted on training fold at each LOO iteration\n")

n <- nrow(data_z)
ki67_mat <- as.matrix(data_z[, ki67_markers])
pc1_scores <- numeric(n)

for (i in seq_len(n)) {
  train_mat <- ki67_mat[-i, , drop = FALSE]
  test_vec  <- ki67_mat[i, , drop = FALSE]
  # Remove columns with zero variance in training
  col_sds <- apply(train_mat, 2, sd, na.rm = TRUE)
  keep <- col_sds > 0 & !is.na(col_sds)
  if (sum(keep) < 2) { pc1_scores[i] <- 0; next }
  pca_fit <- prcomp(train_mat[, keep, drop = FALSE], scale. = TRUE, center = TRUE)
  pc1_scores[i] <- predict(pca_fit, test_vec[, keep, drop = FALSE])[1, 1]
}

cat(sprintf("\nPC1 variance explained (full data PCA for reference):\n"))
pca_full <- prcomp(ki67_mat, scale. = TRUE, center = TRUE)
cat(sprintf("  PC1: %.1f%%\n", 100 * pca_full$sdev[1]^2 / sum(pca_full$sdev^2)))
cat("  PC1 loadings (top contributors):\n")
loadings_sorted <- sort(abs(pca_full$rotation[, 1]), decreasing = TRUE)
for (nm in names(loadings_sorted)[1:min(6, length(loadings_sorted))]) {
  cat(sprintf("    %-14s  %.3f\n", nm, pca_full$rotation[nm, 1]))
}

X_pc1 <- matrix(pc1_scores, ncol = 1)
colnames(X_pc1) <- "KI67_PC1"
res_pc1 <- run_svm_loocv(X_pc1, y, label = "KI67_PC1 only")

# ── Model C: PC1 + KI67NAIVE (best univariate) ───────────────────────────────
cat("\n=== Model C: PC1 + KI67NAIVE ===\n")
X_pc1_naive <- cbind(pc1_scores, as.numeric(data_z[, "KI67NAIVE"]))
colnames(X_pc1_naive) <- c("KI67_PC1", "KI67NAIVE")
res_pc1_naive <- run_svm_loocv(X_pc1_naive, y, label = "KI67_PC1 + KI67NAIVE")

# ── Model D: All KI67 markers directly (no compression) ──────────────────────
cat("\n=== Model D: All KI67 markers directly ===\n")
X_all_ki67 <- as.matrix(data_z[, ki67_markers])
res_all <- run_svm_loocv(X_all_ki67, y, label = paste("All", length(ki67_markers), "KI67 markers"))

# ── DeLong test: PC1 vs baseline ─────────────────────────────────────────────
cat("\n\n=== DeLong tests vs current baseline ===\n")
for (res in list(list(r=res_pc1, l="PC1 only"),
                 list(r=res_pc1_naive, l="PC1+KI67NAIVE"),
                 list(r=res_all, l="All KI67"))) {
  tryCatch({
    dt <- pROC::roc.test(res_baseline$roc, res$r$roc, method = "delong", paired = TRUE)
    cat(sprintf("  Baseline vs %-20s  DeLong p=%.4f\n", res$l, dt$p.value))
  }, error = function(e) cat(sprintf("  DeLong test failed for %s\n", res$l)))
}

cat("\nDone.\n")
