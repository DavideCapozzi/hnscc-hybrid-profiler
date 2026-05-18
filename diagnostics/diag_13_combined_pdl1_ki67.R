# diagnostics/diag_13_combined_pdl1_ki67.R
# ==============================================================================
# DIAGNOSTIC 13: PD-L1 alone vs PD-L1 + KI67 score (combined model)
#
# The most common reviewer question for a candidate ICI biomarker is:
#   "Does it add independent information ON TOP of PD-L1?"
#
# Answers this by:
#   1. Fitting a logistic regression with PD-L1 alone (n=60 with PD-L1 data)
#   2. Fitting logistic regression with PD-L1 + KI67NAIVE + KI67CD4
#   3. Computing LOOCV AUC for both, DeLong test, IDI, continuous NRI
#   4. Fitting SVM-RBF on the combined feature set (3 features) as the
#      non-linear alternative — directly comparable to the primary model
#
# Output: diagnostics/diag_13_output.txt
# ==============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readxl); library(pROC); library(e1071)
})

OUT_FILE <- here("diagnostics/diag_13_output.txt")
con <- file(OUT_FILE, open = "wt")
sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)

# ── Load ─────────────────────────────────────────────────────────────────────
obj_std <- readRDS(here(
  "results/BestResponse_2v3_4/01_data_processing/data_processed_BestResponse_2v3_4_standard.rds"
))
data_z <- obj_std$hybrid_data_z
meta   <- obj_std$metadata
raw_t0 <- read_excel(here("data/Dati_NSCLC_standardizzati_anonimi_T0.xlsx"))

df <- meta %>% select(Patient_ID, Group) %>%
  left_join(raw_t0 %>% select(Patient_ID, PD_L1), by = "Patient_ID") %>%
  mutate(y = ifelse(Group == "RP", 1, 0))

# Restrict to PD-L1-available subset (apples-to-apples comparison)
keep <- !is.na(df$PD_L1)
df_k <- df[keep, ]
y    <- df_k$y
pdl1 <- df_k$PD_L1
ki67naive <- as.numeric(data_z[keep, "KI67NAIVE"])
ki67cd4   <- as.numeric(data_z[keep, "KI67CD4"])

cat(sprintf("Cohort with PD-L1 data: n=%d  (RP=%d  SD_PD=%d)\n",
            nrow(df_k), sum(y == 1), sum(y == 0)))

# ── Helper: LOOCV logistic regression ────────────────────────────────────────
loocv_glm <- function(formula_str, data, y) {
  n <- nrow(data); probs <- numeric(n)
  for (i in seq_len(n)) {
    fit <- glm(as.formula(formula_str), data = data[-i, ], family = binomial())
    probs[i] <- predict(fit, newdata = data[i, , drop = FALSE], type = "response")
  }
  probs
}

df_glm <- data.frame(y = y, PD_L1 = pdl1,
                     KI67NAIVE = ki67naive, KI67CD4 = ki67cd4)

cat("\n=== Model 1: Logistic — PD-L1 alone ===\n")
p1 <- loocv_glm("y ~ PD_L1", df_glm, y)
roc1 <- pROC::roc(y, p1, direction = "<", quiet = TRUE)
ci1  <- pROC::ci.auc(roc1, method = "bootstrap", boot.n = 2000,
                     conf.level = 0.95, progress = "none")
cat(sprintf("  AUC=%.4f [%.3f-%.3f]\n", as.numeric(roc1$auc), ci1[1], ci1[3]))

cat("\n=== Model 2: Logistic — PD-L1 + KI67NAIVE + KI67CD4 ===\n")
p2 <- loocv_glm("y ~ PD_L1 + KI67NAIVE + KI67CD4", df_glm, y)
roc2 <- pROC::roc(y, p2, direction = "<", quiet = TRUE)
ci2  <- pROC::ci.auc(roc2, method = "bootstrap", boot.n = 2000,
                     conf.level = 0.95, progress = "none")
cat(sprintf("  AUC=%.4f [%.3f-%.3f]\n", as.numeric(roc2$auc), ci2[1], ci2[3]))

# DeLong (paired since same patients)
dt <- pROC::roc.test(roc1, roc2, method = "delong", paired = TRUE)
cat(sprintf("\n  DeLong (PD-L1 vs PD-L1+KI67): ΔAUC=%+.4f  p=%.4f\n",
            as.numeric(roc2$auc) - as.numeric(roc1$auc), dt$p.value))

# IDI / continuous NRI (Pencina 2008)
idi_components <- function(p_old, p_new, y) {
  is_evt <- y == 1
  imp_evt <- mean(p_new[is_evt]) - mean(p_old[is_evt])
  red_nev <- mean(p_old[!is_evt]) - mean(p_new[!is_evt])
  list(IDI = imp_evt + red_nev, IDI_event = imp_evt, IDI_nonevent = red_nev)
}
cont_nri <- function(p_old, p_new, y) {
  # Continuous NRI: proportion of events with increased predicted prob minus
  # proportion of non-events with increased predicted prob
  is_evt <- y == 1
  up_evt <- mean((p_new - p_old)[is_evt] > 0) -
            mean((p_new - p_old)[is_evt] < 0)
  dn_nev <- mean((p_new - p_old)[!is_evt] < 0) -
            mean((p_new - p_old)[!is_evt] > 0)
  list(NRI = up_evt + dn_nev, NRI_event = up_evt, NRI_nonevent = dn_nev)
}

idi <- idi_components(p1, p2, y)
nri <- cont_nri(p1, p2, y)
cat(sprintf("\n  IDI = %+.4f  (event:%+.4f  non-event:%+.4f)\n",
            idi$IDI, idi$IDI_event, idi$IDI_nonevent))
cat(sprintf("  cNRI = %+.4f  (event:%+.4f  non-event:%+.4f)\n",
            nri$NRI, nri$NRI_event, nri$NRI_nonevent))

# Bootstrap IDI/NRI 95% CI (1000 resamples, patient-level)
set.seed(2026)
B <- 1000
idi_boot <- numeric(B); nri_boot <- numeric(B)
for (b in seq_len(B)) {
  idx <- sample(seq_along(y), length(y), replace = TRUE)
  yb <- y[idx]; if (length(unique(yb)) < 2) { idi_boot[b] <- NA; next }
  idi_boot[b] <- idi_components(p1[idx], p2[idx], yb)$IDI
  nri_boot[b] <- cont_nri(p1[idx], p2[idx], yb)$NRI
}
cat(sprintf("\n  IDI 95%% bootstrap CI: [%+.4f, %+.4f]  one-sided p>0: %.4f\n",
            quantile(idi_boot, 0.025, na.rm = TRUE),
            quantile(idi_boot, 0.975, na.rm = TRUE),
            mean(idi_boot <= 0, na.rm = TRUE)))
cat(sprintf("  cNRI 95%% bootstrap CI: [%+.4f, %+.4f]  one-sided p>0: %.4f\n",
            quantile(nri_boot, 0.025, na.rm = TRUE),
            quantile(nri_boot, 0.975, na.rm = TRUE),
            mean(nri_boot <= 0, na.rm = TRUE)))

# ── Model 3: SVM-RBF on combined features ───────────────────────────────────
cat("\n=== Model 3: SVM-RBF — PD-L1 + KI67NAIVE + KI67CD4 (nested LOOCV) ===\n")
# Need PD-L1 z-scored to be on comparable scale with the gate
pdl1_z <- as.numeric(scale(pdl1))  # standardized PD-L1
X_full <- cbind(PDL1 = pdl1_z, KI67NAIVE = ki67naive, KI67CD4 = ki67cd4)

run_svm <- function(X, y) {
  n <- nrow(X); probs <- numeric(n)
  for (i in seq_len(n)) {
    X_train <- X[-i, , drop = FALSE]; X_test <- X[i, , drop = FALSE]
    y_train <- y[-i]
    mu <- colMeans(X_train); sds <- apply(X_train, 2, sd); sds[sds == 0] <- 1
    X_train_s <- scale(X_train, center = mu, scale = sds)
    X_test_s  <- scale(X_test,  center = mu, scale = sds)
    set.seed(2026 + i)
    fold_ids <- sample(rep(1:5, length.out = nrow(X_train_s)))
    best_C <- 1; best_g <- 0.1; best_acc <- -Inf
    for (C in c(0.01, 0.1, 1, 10, 100)) for (g in c(0.01, 0.1, 1, 10)) {
      acc <- 0
      for (f in 1:5) {
        xtr <- X_train_s[fold_ids != f, , drop = FALSE]
        xv  <- X_train_s[fold_ids == f, , drop = FALSE]
        ytr <- y_train[fold_ids != f]; yv <- y_train[fold_ids == f]
        if (length(unique(ytr)) < 2) next
        fit <- tryCatch(svm(xtr, as.factor(ytr), kernel = "radial",
                            cost = C, gamma = g, probability = TRUE),
                        error = function(e) NULL)
        if (is.null(fit)) next
        acc <- acc + mean(predict(fit, xv) == as.factor(yv))
      }
      if (acc / 5 > best_acc) { best_acc <- acc / 5; best_C <- C; best_g <- g }
    }
    fit <- tryCatch(svm(X_train_s, as.factor(y_train), kernel = "radial",
                        cost = best_C, gamma = best_g, probability = TRUE),
                    error = function(e) NULL)
    if (is.null(fit)) { probs[i] <- 0.5; next }
    pr <- attr(predict(fit, X_test_s, probability = TRUE), "probabilities")
    probs[i] <- if ("1" %in% colnames(pr)) pr[, "1"] else 0.5
  }
  probs
}

p3 <- run_svm(X_full, y)
roc3 <- pROC::roc(y, p3, direction = "<", quiet = TRUE)
ci3  <- pROC::ci.auc(roc3, method = "bootstrap", boot.n = 2000,
                     conf.level = 0.95, progress = "none")
cat(sprintf("  AUC=%.4f [%.3f-%.3f]\n", as.numeric(roc3$auc), ci3[1], ci3[3]))

# DeLong vs KI67-only baseline (on same n=60 subset)
X_ki <- cbind(KI67NAIVE = ki67naive, KI67CD4 = ki67cd4)
p_ki <- run_svm(X_ki, y)
roc_ki <- pROC::roc(y, p_ki, direction = "<", quiet = TRUE)
ci_ki  <- pROC::ci.auc(roc_ki, method = "bootstrap", boot.n = 2000,
                       conf.level = 0.95, progress = "none")
cat(sprintf("  KI67-only on n=60 subset: AUC=%.4f [%.3f-%.3f]\n",
            as.numeric(roc_ki$auc), ci_ki[1], ci_ki[3]))

dt_ki_full <- pROC::roc.test(roc_ki, roc3, method = "delong", paired = TRUE)
dt_pdl1_full <- pROC::roc.test(roc1, roc3, method = "delong", paired = TRUE)
cat(sprintf("\n  DeLong KI67-only vs PDL1+KI67 (SVM): ΔAUC=%+.4f  p=%.4f\n",
            as.numeric(roc3$auc) - as.numeric(roc_ki$auc), dt_ki_full$p.value))
cat(sprintf("  DeLong PDL1-only vs PDL1+KI67 (SVM): ΔAUC=%+.4f  p=%.4f\n",
            as.numeric(roc3$auc) - as.numeric(roc1$auc), dt_pdl1_full$p.value))

cat("\n=== Interpretation ===\n")
cat("- Logistic IDI > 0 and cNRI > 0: KI67 reclassifies patients beyond PD-L1\n")
cat("- IDI bootstrap CI excluding 0: rigorous evidence of added information\n")
cat("- SVM combined > SVM KI67-only: PD-L1 contributes orthogonal info\n")
cat("- SVM combined ~ SVM KI67-only: KI67 already captures the PD-L1 signal —\n")
cat("  this is the strongest narrative (KI67 dynamics dominate baseline phenotype)\n")

cat("\nDone.\n")
