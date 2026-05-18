# diagnostics/diag_12_pdl1_stratified.R
# ==============================================================================
# DIAGNOSTIC 12: PD-L1 categorical (clinical TPS bins) vs continuous
#
# The current Step 06 benchmark treats PD-L1 as continuous (AUC=0.594, n=60,
# DeLong p=0.151 vs primary). Clinical decision-making in NSCLC uses TPS bins:
#   - high   : TPS ≥ 50%  (pembrolizumab monotherapy indication)
#   - low    : TPS 1-49%  (combo therapy)
#   - neg    : TPS < 1%   (chemo-immuno or chemo)
#
# Test:
#   1. Fisher / Cochran-Armitage trend tests for PD-L1 bins vs Best Response
#   2. KI67 SVM AUC restricted to PD-L1 < 50% subgroup (where PD-L1 is a weak
#      signal and a complementary biomarker is most valuable)
#   3. PD-L1 binary (≥50% vs <50%) accuracy as direct clinical comparator
#
# Output: diagnostics/diag_12_output.txt
# ==============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readxl); library(pROC); library(e1071)
})

OUT_FILE <- here("diagnostics/diag_12_output.txt")
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
pdl1   <- raw_t0 %>% select(Patient_ID, PD_L1)

df <- meta %>%
  select(Patient_ID, Group) %>%
  left_join(pdl1, by = "Patient_ID") %>%
  mutate(
    y = ifelse(Group == "RP", 1, 0),
    PDL1_bin3 = cut(PD_L1, breaks = c(-Inf, 0.99, 49.99, Inf),
                    labels = c("neg(<1%)", "low(1-49%)", "high(>=50%)")),
    PDL1_bin2 = ifelse(PD_L1 >= 50, "high(>=50%)", "low(<50%)")
  )

cat(sprintf("Cohort: %d total | %d with PD-L1 data | %d NA\n",
            nrow(df), sum(!is.na(df$PD_L1)), sum(is.na(df$PD_L1))))

cat("\n=== PD-L1 distribution ===\n")
print(summary(df$PD_L1))

# ── 1. Frequency table & Fisher test on 3 bins ───────────────────────────────
cat("\n=== Cross-tab: PD-L1 3 bins x Response ===\n")
tab3 <- table(df$PDL1_bin3, df$Group, useNA = "no")
print(tab3)

cat("\n  Response rate (RP %) by PD-L1 bin:\n")
rr <- prop.table(tab3, margin = 1)[, "RP"] * 100
for (k in seq_along(rr)) {
  cat(sprintf("    %-12s  RP rate = %5.1f%%  (n=%d)\n",
              names(rr)[k], rr[k], sum(tab3[k, ])))
}

fish3 <- fisher.test(tab3)
cat(sprintf("\n  Fisher exact test (3 bins): p=%.4f\n", fish3$p.value))

# Cochran-Armitage trend test (PD-L1 ordinal vs binary outcome)
ca_score <- c("neg(<1%)" = 0, "low(1-49%)" = 1, "high(>=50%)" = 2)
ca_x <- ca_score[as.character(df$PDL1_bin3)]
ca_df <- data.frame(x = ca_x, y = df$y) %>% filter(!is.na(x))
ca_fit <- prop.trend.test(x = tapply(ca_df$y, ca_df$x, sum),
                          n = tapply(ca_df$y, ca_df$x, length))
cat(sprintf("  Cochran-Armitage trend (0/1/2 scoring): p=%.4f\n", ca_fit$p.value))

# ── 2. Binary cut at 50% (pembro-mono indication threshold) ──────────────────
cat("\n=== Binary cut at TPS=50% (clinical pembro-mono threshold) ===\n")
tab2 <- table(df$PDL1_bin2, df$Group, useNA = "no")
print(tab2)
fish2 <- fisher.test(tab2)
cat(sprintf("  Fisher exact (TPS>=50 vs <50): p=%.4f  OR=%.2f\n",
            fish2$p.value, fish2$estimate))

# PD-L1 binary as "classifier" — sens/spec
df_b <- df %>% filter(!is.na(PDL1_bin2))
pred_pdl1 <- ifelse(df_b$PDL1_bin2 == "high(>=50%)", 1, 0)
sens_pdl1 <- mean(pred_pdl1[df_b$y == 1] == 1)
spec_pdl1 <- mean(pred_pdl1[df_b$y == 0] == 0)
bal_pdl1  <- mean(c(sens_pdl1, spec_pdl1))
cat(sprintf("  As binary classifier: BalAcc=%.3f  Sens=%.3f  Spec=%.3f  (n=%d)\n",
            bal_pdl1, sens_pdl1, spec_pdl1, nrow(df_b)))

# ── 3. KI67 SVM-RBF restricted to PD-L1 < 50% subgroup ───────────────────────
cat("\n=== KI67 model performance in PD-L1 <50% subgroup ===\n")
cat("(The subgroup where PD-L1 alone is insufficient for treatment decision)\n")
keep_sub <- !is.na(df$PD_L1) & df$PD_L1 < 50
n_sub <- sum(keep_sub)
y_sub <- df$y[keep_sub]
cat(sprintf("  Subgroup n=%d  (RP=%d, SD_PD=%d)\n",
            n_sub, sum(y_sub == 1), sum(y_sub == 0)))

X_sub <- as.matrix(data_z[keep_sub, c("KI67NAIVE", "KI67CD4")])

run_svm_loocv <- function(X, y) {
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
        ytr_f <- factor(ytr, levels = c(0, 1))
        yv_f  <- factor(yv,  levels = c(0, 1))
        if (length(unique(ytr_f)) < 2) next
        fit <- tryCatch(svm(xtr, ytr_f, kernel = "radial",
                            cost = C, gamma = g, probability = TRUE),
                        error = function(e) NULL)
        if (is.null(fit)) next
        pv <- predict(fit, xv)
        acc <- acc + mean(as.character(pv) == as.character(yv_f))
      }
      if (acc / 5 > best_acc) { best_acc <- acc / 5; best_C <- C; best_g <- g }
    }
    fit <- tryCatch(svm(X_train_s, factor(y_train, levels = c(0, 1)),
                        kernel = "radial",
                        cost = best_C, gamma = best_g, probability = TRUE),
                    error = function(e) NULL)
    if (is.null(fit)) { probs[i] <- 0.5; next }
    pr <- attr(predict(fit, X_test_s, probability = TRUE), "probabilities")
    probs[i] <- if ("1" %in% colnames(pr)) pr[, "1"] else 0.5
  }
  probs
}

probs_sub <- run_svm_loocv(X_sub, y_sub)
roc_sub <- pROC::roc(y_sub, probs_sub, direction = "<", quiet = TRUE)
ci_sub  <- pROC::ci.auc(roc_sub, method = "bootstrap", boot.n = 2000,
                         conf.level = 0.95, progress = "none")
cat(sprintf("  SVM-RBF on KI67 gate, PD-L1<50%% subgroup: AUC=%.4f [%.3f-%.3f]\n",
            as.numeric(roc_sub$auc), ci_sub[1], ci_sub[3]))

# Same for PD-L1 >= 50%
keep_high <- !is.na(df$PD_L1) & df$PD_L1 >= 50
if (sum(keep_high) >= 10) {
  y_h <- df$y[keep_high]
  X_h <- as.matrix(data_z[keep_high, c("KI67NAIVE", "KI67CD4")])
  cat(sprintf("\n  Subgroup PD-L1>=50%% n=%d  (RP=%d, SD_PD=%d)\n",
              sum(keep_high), sum(y_h == 1), sum(y_h == 0)))
  if (length(unique(y_h)) == 2 && min(table(y_h)) >= 3) {
    probs_h <- run_svm_loocv(X_h, y_h)
    roc_h <- pROC::roc(y_h, probs_h, direction = "<", quiet = TRUE)
    cat(sprintf("  SVM-RBF on KI67 gate, PD-L1>=50%% subgroup: AUC=%.4f\n",
                as.numeric(roc_h$auc)))
  } else {
    cat("  (Skipping: degenerate class distribution)\n")
  }
}

cat("\n=== Interpretation ===\n")
cat("- 3-bin Fisher / CA trend test: clinically standard PD-L1 strata vs response\n")
cat("- KI67 model AUC in PD-L1<50% subgroup is the most clinically actionable\n")
cat("  number: PD-L1 alone is weak in this group, so a complementary marker\n")
cat("  changing treatment selection has direct clinical implication.\n")
cat("- If KI67 subgroup AUC > 0.65, the biomarker is most useful in TPS<50%\n")
cat("  patients — a strong differentiator for the manuscript discussion.\n")

cat("\nDone.\n")
