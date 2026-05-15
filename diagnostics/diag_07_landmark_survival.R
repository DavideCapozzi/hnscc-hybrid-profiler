# diagnostics/diag_07_landmark_survival.R
# Diagnostic: landmark Cox analysis using KI67 delta (T1 - T0) as predictor of PFS/OS.
#
# Purpose:
#   Tests whether the treatment-induced KI67 change (pharmacodynamic Δ) predicts
#   survival beyond what the binary response category (RP vs SD/PD) already explains.
#   Uses landmark analysis anchored at T1 to avoid immortal-time bias.
#
#   Key question: does a deeper KI67 reset at T1 independently predict longer PFS,
#   even after adjusting for response group? If yes: Δ is a pharmacodynamic
#   predictor of survival, not just a correlate of binary response.
#
# Analyses:
#   1. Compute KI67NAIVE_delta and KI67CD4_delta for paired patients (n=50)
#   2. KM curves by KI67NAIVE_delta above/below median
#   3. Univariate Cox: PFS ~ KI67NAIVE_delta (landmark from T1)
#   4. Adjusted Cox: PFS ~ KI67NAIVE_delta + Group
#   5. C-statistic comparison: Group alone vs Group + KI67NAIVE_delta
#   6. Repeat key models for OS
#
# Output: diagnostics/diag_07_output.txt

library(here)
library(readxl)
library(survival)
library(dplyr)

OUT_FILE <- here("diagnostics/diag_07_output.txt")
con <- file(OUT_FILE, open = "wt")
sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)

# ── 1. Load data ──────────────────────────────────────────────────────────────
DATA_L <- readRDS(here(
  "results/BestResponse_2v3_4/01_data_processing/data_processed_BestResponse_2v3_4_longitudinal.rds"
))
hz   <- DATA_L$hybrid_data_z
meta <- DATA_L$metadata

df_raw_t0 <- read_excel(here("data/Dati_NSCLC_standardizzati_anonimi_T0.xlsx"))

surv_map <- df_raw_t0 %>%
  select(Patient_ID, PFS, OS, Dead = `Alive_0/Dead_1`) %>%
  mutate(
    Dead     = as.integer(Dead),
    event_os = Dead,
    event_pfs = as.integer(Dead == 1 | PFS < OS)
  )

# ── 2. Compute Δ for paired patients ──────────────────────────────────────────
df_t0 <- hz %>% filter(Timepoint == "T0") %>%
  select(Patient_ID, Group, KI67NAIVE_T0 = KI67NAIVE,
         KI67CD4_T0 = KI67CD4, CD28KI67_T0 = CD28KI67)

df_t1 <- hz %>% filter(Timepoint == "T1") %>%
  select(Patient_ID, KI67NAIVE_T1 = KI67NAIVE,
         KI67CD4_T1 = KI67CD4, CD28KI67_T1 = CD28KI67)

df_delta <- inner_join(df_t0, df_t1, by = "Patient_ID") %>%
  mutate(
    KI67NAIVE_delta = KI67NAIVE_T1 - KI67NAIVE_T0,
    KI67CD4_delta   = KI67CD4_T1   - KI67CD4_T0,
    CD28KI67_delta  = CD28KI67_T1  - CD28KI67_T0
  ) %>%
  left_join(surv_map, by = "Patient_ID") %>%
  filter(!is.na(PFS)) %>%
  mutate(Group = factor(Group, levels = c("SD_PD", "RP")))

cat("=== Landmark Survival Analysis (Cox) ===\n")
cat("Paired patients (T0+T1):", nrow(df_delta), "\n")
cat("Deaths:", sum(df_delta$Dead == 1, na.rm = TRUE), "\n")
cat("Median PFS (months):", median(df_delta$PFS), "\n")
cat("Groups:", table(df_delta$Group), "\n\n")

cat("KI67NAIVE_delta summary:\n")
print(summary(df_delta$KI67NAIVE_delta))

# ── 3. KM strata by delta above/below median ──────────────────────────────────
cat("\n=== KM log-rank: PFS by KI67NAIVE_delta (above/below median) ===\n")
med_delta <- median(df_delta$KI67NAIVE_delta, na.rm = TRUE)
df_delta$delta_group <- ifelse(df_delta$KI67NAIVE_delta < med_delta, "Reset_high", "Reset_low")
cat("Median delta:", round(med_delta, 3), "\n")
cat("N per delta group:", table(df_delta$delta_group), "\n")

km_delta <- survfit(Surv(PFS, event_pfs) ~ delta_group, data = df_delta)
lr_delta  <- survdiff(Surv(PFS, event_pfs) ~ delta_group, data = df_delta)
cat("KM median PFS by delta group:\n")
print(summary(km_delta)$table[, c("records", "events", "median")])
cat("Log-rank p:", round(1 - pchisq(lr_delta$chisq, df = 1), 4), "\n")

# ── 4. Cox models ─────────────────────────────────────────────────────────────
fit_models <- function(time_col, event_col, label) {
  cat("\n===", label, "===\n")

  # Use a single complete-case dataset for all models to allow LRT comparison
  df_cox <- df_delta %>%
    filter(!is.na(.data[[time_col]]), !is.na(.data[[event_col]]),
           !is.na(KI67NAIVE_delta), !is.na(KI67CD4_delta))

  cat("n (complete cases):", nrow(df_cox), "| events:", sum(df_cox[[event_col]]), "\n")

  # Model A: Group only (baseline)
  mA <- coxph(Surv(df_cox[[time_col]], df_cox[[event_col]]) ~ Group,
               data = df_cox, ties = "efron")

  # Model B: KI67NAIVE_delta only
  mB <- coxph(Surv(df_cox[[time_col]], df_cox[[event_col]]) ~ KI67NAIVE_delta,
               data = df_cox, ties = "efron")

  # Model C: Group + KI67NAIVE_delta (key model)
  mC <- coxph(Surv(df_cox[[time_col]], df_cox[[event_col]]) ~ Group + KI67NAIVE_delta,
               data = df_cox, ties = "efron")

  # Model D: KI67 composite (KI67NAIVE_delta + KI67CD4_delta)
  mD <- coxph(Surv(df_cox[[time_col]], df_cox[[event_col]]) ~
                Group + KI67NAIVE_delta + KI67CD4_delta,
               data = df_cox, ties = "efron")

  for (nm in c("A: Group only", "B: KI67NAIVE_delta only",
                "C: Group + KI67NAIVE_delta", "D: Group + KI67NAIVE + KI67CD4")) {
    m <- list(mA, mB, mC, mD)[[which(c("A", "B", "C", "D") == substr(nm, 1, 1))]]
    cstat <- concordance(m)$concordance
    cat("\nModel", nm, "\n")
    coefs <- summary(m)$coefficients
    print(round(coefs[, c("coef", "exp(coef)", "se(coef)", "Pr(>|z|)")], 4))
    cat("C-statistic:", round(cstat, 3), "\n")
  }

  # LRT: does adding delta improve on Group alone?
  lrt <- anova(mA, mC)
  cat("\nLRT Group-only vs Group+delta: p =", round(lrt[2, "Pr(>|Chi|)"], 4), "\n")
}

fit_models("PFS", "event_pfs", "PFS landmark analysis")
fit_models("OS",  "event_os",  "OS landmark analysis")

# ── 5. Correlation: delta vs PFS (continuous) ─────────────────────────────────
cat("\n=== Spearman correlation: KI67NAIVE_delta vs PFS ===\n")
sp <- cor.test(df_delta$KI67NAIVE_delta, df_delta$PFS,
               method = "spearman", exact = FALSE)
cat("rho =", round(sp$estimate, 3), ", p =", round(sp$p.value, 4), "\n")

cat("\n=== Spearman correlation: KI67NAIVE_delta vs OS ===\n")
sp2 <- cor.test(df_delta$KI67NAIVE_delta, df_delta$OS,
                method = "spearman", exact = FALSE)
cat("rho =", round(sp2$estimate, 3), ", p =", round(sp2$p.value, 4), "\n")

# ── 6. Delta by group: sanity check ───────────────────────────────────────────
cat("\n=== KI67NAIVE_delta by response group ===\n")
print(tapply(df_delta$KI67NAIVE_delta, df_delta$Group,
             function(x) c(mean = round(mean(x, na.rm=TRUE), 3),
                           median = round(median(x, na.rm=TRUE), 3),
                           n = sum(!is.na(x)))))
cat("\nWilcoxon test RP vs SD_PD:\n")
wt <- wilcox.test(KI67NAIVE_delta ~ Group, data = df_delta, exact = FALSE)
cat("W =", wt$statistic, ", p =", round(wt$p.value, 4), "\n")

# ── 7. Dropout T1 bias check ──────────────────────────────────────────────────
cat("\n=== T1 dropout bias check ===\n")
all_t0 <- hz %>% filter(Timepoint == "T0") %>%
  left_join(surv_map, by = "Patient_ID") %>%
  mutate(has_T1 = Patient_ID %in% df_delta$Patient_ID,
         Group = meta$Group[match(Patient_ID, meta$Patient_ID)])

cat("PFS by T1 availability:\n")
print(tapply(all_t0$PFS, all_t0$has_T1, function(x)
  c(median = median(x, na.rm=TRUE), n = length(x))))

cat("\nGroup distribution — with T1 vs without T1:\n")
print(table(all_t0$has_T1, all_t0$Group))
