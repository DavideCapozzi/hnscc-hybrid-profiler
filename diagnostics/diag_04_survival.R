# diagnostics/diag_04_survival.R
# Diagnostic: PFS/OS survival analysis with KI67NAIVE and KI67CD4
# Uses processed T0 data (n=70 post-QC) matched to raw clinical survival columns.
# Output mirrored to diagnostics/diag_04_output.txt

library(readxl)
library(survival)
library(dplyr)

OUT_FILE <- "diagnostics/diag_04_output.txt"
con <- file(OUT_FILE, open = "wt")
sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)

DATA    <- readRDS("results/BestResponse_2v3_4/01_data_processing/data_processed_BestResponse_2v3_4_standard.rds")
df_raw  <- read_excel("data/Dati_NSCLC_standardizzati_T0.xlsx")

# --- Prepare survival data ---
df_surv <- df_raw %>%
  select(Patient_ID,
         PFS, OS,
         event_os = `Alive_0/Dead_1`) %>%
  mutate(
    event_os  = as.integer(event_os),
    event_pfs = as.integer(event_os == 1 | (event_os == 0 & PFS < OS))
  )

# --- Merge with processed T0 markers ---
meta  <- DATA$metadata
hz    <- DATA$hybrid_data_z

df_ml <- data.frame(
  Patient_ID = meta$Patient_ID,
  Group      = as.character(meta$Group),
  KI67NAIVE  = hz$KI67NAIVE,
  KI67CD4    = hz$KI67CD4,
  stringsAsFactors = FALSE
)

df <- inner_join(df_ml, df_surv, by = "Patient_ID")

# Remove patients with zero follow-up (survival models require time > 0)
zero_time <- df$Patient_ID[df$OS <= 0 | df$PFS <= 0]
if (length(zero_time) > 0) {
  cat(sprintf("[NOTE] Excluding %d patient(s) with zero follow-up time: %s\n",
              length(zero_time), paste(zero_time, collapse = ", ")))
  df <- df %>% filter(OS > 0, PFS > 0)
}

cat(sprintf("Matched patients: %d | OS events: %d | PFS events: %d\n",
            nrow(df), sum(df$event_os), sum(df$event_pfs)))
cat(sprintf("Median OS: %.1f months | Median PFS: %.1f months\n",
            median(df$OS), median(df$PFS)))

# --- Cox models ---
run_cox <- function(time, event, predictor, marker_name, endpoint) {
  fit <- coxph(Surv(time, event) ~ predictor)
  s   <- summary(fit)
  hr  <- s$conf.int[1, c("exp(coef)", "lower .95", "upper .95")]
  p   <- s$coefficients[1, "Pr(>|z|)"]
  cat(sprintf("  Cox %s ~ %s: HR=%.3f [%.3f-%.3f] p=%.4f\n",
              endpoint, marker_name, hr[1], hr[2], hr[3], p))
  invisible(list(hr = hr, p = p, fit = fit))
}

cat("\n=== Cox PH — KI67NAIVE (z-score) ===\n")
cox_os_ki67n  <- run_cox(df$OS,  df$event_os,  df$KI67NAIVE, "KI67NAIVE", "OS")
cox_pfs_ki67n <- run_cox(df$PFS, df$event_pfs, df$KI67NAIVE, "KI67NAIVE", "PFS")

cat("\n=== Cox PH — KI67CD4 (z-score) ===\n")
cox_os_ki67c  <- run_cox(df$OS,  df$event_os,  df$KI67CD4, "KI67CD4", "OS")
cox_pfs_ki67c <- run_cox(df$PFS, df$event_pfs, df$KI67CD4, "KI67CD4", "PFS")

# --- Optimal cutpoint (Youden on Kaplan-Meier log-rank) ---
find_best_cut <- function(time, event, marker, n_cuts = 50) {
  cuts  <- quantile(marker, probs = seq(0.2, 0.8, length.out = n_cuts))
  chisqs <- sapply(cuts, function(c) {
    grp <- ifelse(marker >= c, "High", "Low")
    if (length(unique(grp)) < 2) return(NA_real_)
    tryCatch(survdiff(Surv(time, event) ~ grp)$chisq, error = function(e) NA_real_)
  })
  chisqs <- as.numeric(chisqs)
  valid  <- which(!is.na(chisqs))
  if (length(valid) == 0) return(list(cut = NA, chisq = NA, p = NA))
  best_i <- valid[which.max(chisqs[valid])]
  list(cut = cuts[best_i], chisq = chisqs[best_i],
       p = pchisq(chisqs[best_i], df = 1, lower.tail = FALSE))
}

cat("\n=== Optimal cutpoint — KI67NAIVE ===\n")
cut_os_n  <- find_best_cut(df$OS,  df$event_os,  df$KI67NAIVE)
cut_pfs_n <- find_best_cut(df$PFS, df$event_pfs, df$KI67NAIVE)
cat(sprintf("  OS  best cut=%.3f (z-score) | log-rank chi2=%.2f | p=%.4f\n",
            cut_os_n$cut, cut_os_n$chisq, cut_os_n$p))
cat(sprintf("  PFS best cut=%.3f (z-score) | log-rank chi2=%.2f | p=%.4f\n",
            cut_pfs_n$cut, cut_pfs_n$chisq, cut_pfs_n$p))

cat("\n=== Optimal cutpoint — KI67CD4 ===\n")
cut_os_c  <- find_best_cut(df$OS,  df$event_os,  df$KI67CD4)
cut_pfs_c <- find_best_cut(df$PFS, df$event_pfs, df$KI67CD4)
cat(sprintf("  OS  best cut=%.3f (z-score) | log-rank chi2=%.2f | p=%.4f\n",
            cut_os_c$cut, cut_os_c$chisq, cut_os_c$p))
cat(sprintf("  PFS best cut=%.3f (z-score) | log-rank chi2=%.2f | p=%.4f\n",
            cut_pfs_c$cut, cut_pfs_c$chisq, cut_pfs_c$p))

# --- Concordance index ---
cat("\n=== Concordance index (C-stat) ===\n")
cat(sprintf("  OS  KI67NAIVE: C=%.3f\n", concordance(Surv(df$OS, df$event_os) ~ df$KI67NAIVE)$concordance))
cat(sprintf("  PFS KI67NAIVE: C=%.3f\n", concordance(Surv(df$PFS, df$event_pfs) ~ df$KI67NAIVE)$concordance))
cat(sprintf("  OS  KI67CD4:   C=%.3f\n", concordance(Surv(df$OS, df$event_os) ~ df$KI67CD4)$concordance))
cat(sprintf("  PFS KI67CD4:   C=%.3f\n", concordance(Surv(df$PFS, df$event_pfs) ~ df$KI67CD4)$concordance))

# --- KM curves with optimal cut ---
km_summary <- function(time, event, marker, cut_v, marker_name, endpoint) {
  grp    <- factor(ifelse(marker >= cut_v,
                          paste0("High_", marker_name),
                          paste0("Low_",  marker_name)))
  km     <- survfit(Surv(time, event) ~ grp)
  med    <- summary(km)$table[, "median"]
  med[is.nan(med) | is.na(med)] <- Inf   # median not reached
  sd_res <- survdiff(Surv(time, event) ~ grp)
  p_lr   <- pchisq(sd_res$chisq, df = 1, lower.tail = FALSE)
  # survfit prefixes strata names with "grp="; match by grep
  low_key  <- names(med)[grepl("Low_",  names(med))]
  high_key <- names(med)[grepl("High_", names(med))]
  med_low  <- if (length(low_key)  && is.finite(med[low_key]))  sprintf("%.1f", med[low_key])  else "NR"
  med_high <- if (length(high_key) && is.finite(med[high_key])) sprintf("%.1f", med[high_key]) else "NR"
  low_lbl  <- paste0("Low_",  marker_name)
  high_lbl <- paste0("High_", marker_name)
  cat(sprintf("  %s: Low median=%s | High median=%s | log-rank p=%.4f\n",
              endpoint, med_low, med_high, p_lr))
  cat(sprintf("       n(Low)=%d n(High)=%d\n",
              sum(grp == low_lbl), sum(grp == high_lbl)))
}

cat("\n=== Kaplan-Meier with optimal cut (KI67NAIVE) ===\n")
km_summary(df$OS,  df$event_os,  df$KI67NAIVE, cut_os_n$cut,  "KI67NAIVE", "OS")
km_summary(df$PFS, df$event_pfs, df$KI67NAIVE, cut_pfs_n$cut, "KI67NAIVE", "PFS")

cat("\n=== Kaplan-Meier with optimal cut (KI67CD4) ===\n")
km_summary(df$OS,  df$event_os,  df$KI67CD4, cut_os_c$cut,  "KI67CD4", "OS")
km_summary(df$PFS, df$event_pfs, df$KI67CD4, cut_pfs_c$cut, "KI67CD4", "PFS")

# --- Multivariate Cox (KI67NAIVE + PD-L1 + PS) ---
cat("\n=== Multivariate Cox OS: KI67NAIVE + PD_L1 + PS ===\n")
df_mv <- df %>%
  left_join(df_raw %>% select(Patient_ID, PD_L1, PS = PS_prima_immuno), by = "Patient_ID") %>%
  filter(!is.na(PD_L1), !is.na(PS))
cat(sprintf("  n for multivariate: %d\n", nrow(df_mv)))
if (nrow(df_mv) >= 20) {
  fit_mv <- coxph(Surv(OS, event_os) ~ KI67NAIVE + PD_L1 + PS, data = df_mv)
  s_mv   <- summary(fit_mv)
  out_mv <- cbind(
    s_mv$conf.int[,  c("exp(coef)", "lower .95", "upper .95"), drop = FALSE],
    "Pr(>|z|)" = s_mv$coefficients[, "Pr(>|z|)"]
  )
  print(round(out_mv, 4))
}
