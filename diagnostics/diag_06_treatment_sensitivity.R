# diagnostics/diag_06_treatment_sensitivity.R
# Diagnostic: LMM with treatment type as covariate.
#
# Purpose:
#   Tests whether the Timepoint x Group interaction for KI67NAIVE, KI67CD4,
#   and CD28KI67 survives after adjusting for treatment regime
#   (1=Pembro mono, 2=Pembro+CHT, 3=Ipi/Nivo+CHT).
#   Addresses the key reviewer concern: is the KI67 reset immune-specific
#   or confounded by cytotoxic chemotherapy?
#
# Output: diagnostics/diag_06_output.txt

library(here)
library(readxl)
library(lmerTest)
library(dplyr)

OUT_FILE <- here("diagnostics/diag_06_output.txt")
con <- file(OUT_FILE, open = "wt")
sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)

# ── 1. Load processed longitudinal RDS ────────────────────────────────────────
DATA_L <- readRDS(here(
  "results/BestResponse_2v3_4/01_data_processing/data_processed_BestResponse_2v3_4_longitudinal.rds"
))
meta <- DATA_L$metadata
hz   <- DATA_L$hybrid_data_z

df_raw_t0 <- read_excel(here("data/Dati_NSCLC_standardizzati_anonimi_T0.xlsx"))

trt_col <- "Terapia_Pembro_1_Pembro_CHT_2_ipi/nivo_CHT_3;_Atezo_CHT_4"

trt_map <- df_raw_t0 %>%
  select(Patient_ID, Terapia = all_of(trt_col)) %>%
  mutate(Terapia = factor(Terapia, levels = c(1, 2, 3),
                          labels = c("Pembro_mono", "Pembro_CHT", "IpiNivo_CHT")))

# ── 2. Build long format data frame ───────────────────────────────────────────
df_long <- hz %>%
  select(Patient_ID, Timepoint, Group, KI67NAIVE, KI67CD4, CD28KI67) %>%
  left_join(trt_map, by = "Patient_ID") %>%
  filter(!is.na(Terapia)) %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("T0", "T1")),
    Group     = factor(Group,     levels = c("SD_PD", "RP")),
    Patient_ID = factor(Patient_ID)
  )

n_with_trt <- length(unique(df_long$Patient_ID))
cat("=== LMM with Treatment Covariate ===\n")
cat("Patients with treatment info:", n_with_trt, "\n")
cat("Treatment distribution:\n")
print(table(unique(df_long[, c("Patient_ID", "Terapia")])$Terapia))
cat("\n")

# ── 3. Fit LMMs with and without treatment covariate ─────────────────────────
MARKERS <- c("KI67NAIVE", "KI67CD4", "CD28KI67")

results <- lapply(MARKERS, function(marker) {
  df_m <- df_long %>% filter(!is.na(.data[[marker]]))

  # Model without treatment (baseline)
  f_base <- as.formula(paste0(marker, " ~ Timepoint * Group + (1 | Patient_ID)"))
  m_base <- tryCatch(lmer(f_base, data = df_m, REML = FALSE), error = function(e) NULL)

  # Model with treatment as covariate
  f_trt  <- as.formula(paste0(marker, " ~ Timepoint * Group + Terapia + (1 | Patient_ID)"))
  m_trt  <- tryCatch(lmer(f_trt,  data = df_m, REML = FALSE), error = function(e) NULL)

  if (is.null(m_base) || is.null(m_trt)) {
    cat(marker, ": model failed\n")
    return(NULL)
  }

  coef_base <- coef(summary(m_base))
  coef_trt  <- coef(summary(m_trt))

  int_row_base <- coef_base["TimepointT1:GroupRP", , drop = FALSE]
  int_row_trt  <- coef_trt["TimepointT1:GroupRP",  , drop = FALSE]

  trt_rows <- coef_trt[grep("Terapia", rownames(coef_trt)), , drop = FALSE]

  list(marker = marker, base = int_row_base, trt = int_row_trt,
       trt_coefs = trt_rows, lrt_p = anova(m_base, m_trt)$`Pr(>Chisq)`[2])
})

# ── 4. Print comparison table ─────────────────────────────────────────────────
cat("=== Timepoint x Group interaction: baseline vs treatment-adjusted ===\n\n")
cat(sprintf("%-12s  %-8s %-8s %-8s   %-8s %-8s %-8s   LRT_p(treatment)\n",
            "Marker", "beta_base", "SE_base", "p_base",
            "beta_adj", "SE_adj", "p_adj"))
cat(strrep("-", 85), "\n")

for (r in results) {
  if (is.null(r)) next
  b <- r$base[1, ]
  a <- r$trt[1, ]
  cat(sprintf("%-12s  %8.3f %8.3f %8.4f   %8.3f %8.3f %8.4f   %.4f\n",
              r$marker,
              b["Estimate"], b["Std. Error"], b["Pr(>|t|)"],
              a["Estimate"], a["Std. Error"], a["Pr(>|t|)"],
              r$lrt_p))
}

cat("\n=== Treatment regime coefficients (adjusted model) ===\n")
for (r in results) {
  if (is.null(r)) next
  cat("\n", r$marker, ":\n", sep = "")
  print(round(r$trt_coefs[, c("Estimate", "Std. Error", "Pr(>|t|)")], 4))
}

cat("\n=== Interpretation ===\n")
cat("- beta_base:  Timepoint x Group estimate in standard LMM (no treatment covariate)\n")
cat("- beta_adj:   Same estimate after adding Terapia as fixed effect\n")
cat("- LRT_p:      Likelihood ratio test p-value for adding treatment to the model\n")
cat("- If beta_adj ~ beta_base and p_adj remains < 0.05, the KI67 signal is\n")
cat("  robust to treatment confounding.\n")
