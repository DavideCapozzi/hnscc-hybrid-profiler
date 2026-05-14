# diagnostics/diag_01_lmm_covariates.R
# ==============================================================================
# DIAGNOSTIC TEST 1: LMM with clinical covariates + marginal R²
#
# Tests whether KI67 markers (KI67NAIVE, CD28KI67, KI67CD4) retain significance
# after adjusting for clinical covariates (PD-L1, PS, treatment, sex, smoking).
# Also computes marginal R² (Nakagawa & Schielzeth 2013) for each model.
#
# Usage: Rscript diagnostics/diag_01_lmm_covariates.R
# ==============================================================================

suppressPackageStartupMessages({
  library(lmerTest)
  library(lme4)
  library(dplyr)
  library(readxl)
  library(tidyr)
})

setwd("/home/laboratorio/projects/clinical-onco-profiler")

# ── Marginal R² (manual implementation, Nakagawa & Schielzeth 2013) ─────────
r2_nakagawa <- function(mod) {
  tryCatch({
    # Fixed effects variance
    sigma2_f <- var(predict(mod, re.form = NA))
    # Random effects variance (Patient intercept)
    vc <- as.data.frame(lme4::VarCorr(mod))
    sigma2_r <- sum(vc$vcov[vc$grp != "Residual"])
    sigma2_e <- sigma(mod)^2
    total <- sigma2_f + sigma2_r + sigma2_e
    list(
      R2m = round(sigma2_f / total, 4),
      R2c = round((sigma2_f + sigma2_r) / total, 4)
    )
  }, error = function(e) list(R2m = NA, R2c = NA))
}

# ── Load clinical metadata ───────────────────────────────────────────────────
cat("Loading raw clinical data...\n")
df_raw_t0 <- read_excel("data/Dati_NSCLC_standardizzati_T0.xlsx")
df_raw_t1 <- read_excel("data/Dati_NSCLC_standardizzati_T1.xlsx")

# Standardise clinical covariate column names
df_clin <- df_raw_t0 %>%
  select(
    Patient_ID,
    PD_L1,
    PS          = PS_prima_immuno,
    Therapy     = `Terapia_Pembro_1_Pembro_CHT_2_ipi/nivo_CHT_3;_Atezo_CHT_4`,
    Sex         = `Sex_M_1;_F_2`,
    Smoking     = `Fumatori_0_no;_1_ex;2_si`
  ) %>%
  mutate(
    PS       = as.numeric(PS),
    PD_L1    = as.numeric(PD_L1),
    Therapy  = as.factor(Therapy),
    Sex      = as.factor(Sex),
    Smoking  = as.factor(Smoking)
  )

cat(sprintf("  Clinical data: %d patients\n", nrow(df_clin)))
cat(sprintf("  PD_L1 available: %d / %d\n", sum(!is.na(df_clin$PD_L1)), nrow(df_clin)))
cat(sprintf("  PS available:    %d / %d\n", sum(!is.na(df_clin$PS)), nrow(df_clin)))
cat(sprintf("  Therapy avail:   %d / %d\n", sum(!is.na(df_clin$Therapy)), nrow(df_clin)))
cat(sprintf("  Sex available:   %d / %d\n", sum(!is.na(df_clin$Sex)), nrow(df_clin)))
cat(sprintf("  Smoking avail:   %d / %d\n\n", sum(!is.na(df_clin$Smoking)), nrow(df_clin)))

# ── Load longitudinal processed data ────────────────────────────────────────
cat("Loading longitudinal processed data...\n")
obj_long <- readRDS("results/BestResponse_2v3_4/01_data_processing/data_processed_BestResponse_2v3_4_longitudinal.rds")

data_z  <- obj_long$hybrid_data_z    # z-scored feature matrix
meta    <- obj_long$metadata         # Patient_ID, Timepoint, Group

cat(sprintf("  Longitudinal data: %d obs, %d markers\n\n",
            nrow(meta), ncol(data_z)))

# data_z already includes Patient_ID, Timepoint, Group columns
df_long <- dplyr::left_join(as.data.frame(data_z), df_clin, by = "Patient_ID")

# ── Define markers and covariate sets ───────────────────────────────────────
# Focus on the 3 FDR-significant + LOO-robust markers + the KI67 family trend markers
focal_markers <- c("KI67NAIVE", "CD28KI67", "KI67CD4",
                   "KI67CD8", "KI67CM", "KI657EMRA", "KI67EFF")

covariate_sets <- list(
  "Baseline (no covariates)" = NULL,
  "+ PD_L1"                  = c("PD_L1"),
  "+ PD_L1 + PS"             = c("PD_L1", "PS"),
  "+ PD_L1 + PS + Therapy"   = c("PD_L1", "PS", "Therapy"),
  "+ All (PD_L1+PS+Therapy+Sex+Smoking)" = c("PD_L1", "PS", "Therapy", "Sex", "Smoking")
)

# ── Helper: fit LMM, extract interaction term + R² ──────────────────────────
fit_lmm_with_r2 <- function(df, feature, covariates = NULL) {
  vals <- as.numeric(df[[feature]])
  val_sd <- sd(vals, na.rm = TRUE)
  if (is.na(val_sd) || val_sd == 0) val_sd <- 1

  df_model <- data.frame(
    Value    = vals / val_sd,
    Time     = as.factor(df$Timepoint),
    Group    = as.factor(df$Group),
    ID       = as.factor(df$Patient_ID)
  )
  valid_covs <- character(0)
  if (!is.null(covariates)) {
    valid_covs <- intersect(covariates, names(df))
    for (cv in valid_covs) df_model[[cv]] <- df[[cv]]
  }
  df_model <- df_model[complete.cases(df_model), ]

  if (nrow(df_model) < 10) return(NULL)

  formula_str <- "Value ~ Time * Group"
  if (length(valid_covs) > 0)
    formula_str <- paste(formula_str, "+", paste(valid_covs, collapse = " + "))
  formula_str <- paste(formula_str, "+ (1 | ID)")

  tryCatch({
    mod <- suppressMessages(suppressWarnings(
      lmerTest::lmer(as.formula(formula_str), data = df_model,
                     REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE))
    ))
    ct   <- summary(mod)$coefficients
    idx  <- grep("Time.*:Group", rownames(ct))
    if (length(idx) != 1) return(NULL)

    r2  <- r2_nakagawa(mod)
    list(
      beta  = ct[idx, "Estimate"] * val_sd,
      se    = ct[idx, "Std. Error"] * val_sd,
      tval  = ct[idx, grep("t value", colnames(ct))],
      pval  = ct[idx, grep("Pr\\(>\\|t\\|\\)", colnames(ct))],
      n     = nrow(df_model),
      R2m   = r2$R2m,
      R2c   = r2$R2c,
      singular = lme4::isSingular(mod)
    )
  }, error = function(e) NULL)
}

# ── Run all combinations ─────────────────────────────────────────────────────
cat("Running LMM across covariate sets...\n")
cat(strrep("=", 90), "\n")

results_all <- list()

for (cov_name in names(covariate_sets)) {
  covs <- covariate_sets[[cov_name]]
  cat(sprintf("\n[%s]\n", cov_name))
  cat(sprintf("  %-14s %8s %8s %8s %8s %6s %6s %s\n",
              "Marker", "Beta", "SE", "t", "p", "R2m", "R2c", "Singular"))
  cat(sprintf("  %s\n", strrep("-", 72)))

  for (mk in focal_markers) {
    res <- fit_lmm_with_r2(df_long, mk, covs)
    if (is.null(res)) {
      cat(sprintf("  %-14s  FAILED\n", mk))
      next
    }
    sig_flag <- if (!is.na(res$pval) && res$pval < 0.05) " *" else ""
    cat(sprintf("  %-14s %8.4f %8.4f %8.3f %8.4f %6.3f %6.3f  %s%s\n",
                mk, res$beta, res$se, res$tval, res$pval,
                res$R2m, res$R2c,
                if (res$singular) "[sing]" else "",
                sig_flag))
    results_all[[paste(cov_name, mk, sep = "|")]] <- c(
      marker = mk, model = cov_name,
      beta = res$beta, pval = res$pval, R2m = res$R2m, R2c = res$R2c
    )
  }
}

# ── FDR after adjustment ──────────────────────────────────────────────────────
cat("\n\n", strrep("=", 90), "\n")
cat("FDR correction within each covariate model (focal 7 markers)\n")
cat(strrep("=", 90), "\n")

for (cov_name in names(covariate_sets)) {
  covs <- covariate_sets[[cov_name]]
  pvals <- sapply(focal_markers, function(mk) {
    res <- fit_lmm_with_r2(df_long, mk, covs)
    if (is.null(res)) NA else res$pval
  })
  fdrs <- p.adjust(pvals, method = "BH")
  n_sig_fdr <- sum(fdrs < 0.05, na.rm = TRUE)
  cat(sprintf("\n[%s] → %d markers survive FDR<0.05\n", cov_name, n_sig_fdr))
  for (i in seq_along(focal_markers)) {
    star <- if (!is.na(fdrs[i]) && fdrs[i] < 0.05) " ***" else
            if (!is.na(pvals[i]) && pvals[i] < 0.05) " (nom. sig)" else ""
    cat(sprintf("  %-14s  raw p=%6.4f  FDR=%6.4f%s\n",
                focal_markers[i], pvals[i], fdrs[i], star))
  }
}

cat("\n\nDone.\n")
