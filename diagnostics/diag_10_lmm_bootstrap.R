# diagnostics/diag_10_lmm_bootstrap.R
# ==============================================================================
# DIAGNOSTIC 10: Patient-level cluster bootstrap of LMM betas
#
# For each primary marker (KI67NAIVE, CD28KI67, KI67CD4) draws B=1000 bootstrap
# samples by *resampling patient IDs with replacement* (cluster bootstrap —
# preserves the paired T0/T1 structure inside each patient). Reports:
#   - median, 2.5%/97.5% bootstrap CI of β_interaction
#   - bootstrap p-value (proportion of resamples with β ≥ 0)
#   - % of bootstrap resamples in which the marker still passes FDR<0.05
#
# This is the standard sensitivity check reviewers ask for to confirm the LMM
# point estimates are not driven by a few influential patients.
#
# Output: diagnostics/diag_10_output.txt
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(lmerTest)
})

OUT_FILE <- here("diagnostics/diag_10_output.txt")
con <- file(OUT_FILE, open = "wt")
sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)

# Source the LMM helper (reuse production fit_feature_lmm to avoid drift)
source(here("R/modules_longitudinal.R"))

DATA_L <- readRDS(here(
  "results/BestResponse_2v3_4/01_data_processing/data_processed_BestResponse_2v3_4_longitudinal.rds"
))
hz <- DATA_L$hybrid_data_z

META_COLS   <- c("Patient_ID", "Sample_ID", "Timepoint", "Group")
ALL_MARKERS <- setdiff(colnames(hz), META_COLS)
PRIMARY <- c("KI67NAIVE", "CD28KI67", "KI67CD4")
stopifnot(all(PRIMARY %in% ALL_MARKERS))

cat(sprintf("Longitudinal data: %d obs, %d patients, %d markers\n",
            nrow(hz), length(unique(hz$Patient_ID)), length(ALL_MARKERS)))
cat(sprintf("Primary markers: %s\n", paste(PRIMARY, collapse = ", ")))

B <- 500  # bootstrap iterations (~13 min runtime)
set.seed(2026)

unique_pids <- unique(hz$Patient_ID)
n_pat <- length(unique_pids)

cat(sprintf("\nRunning %d cluster bootstrap iterations (patient-level resampling)...\n", B))
cat("Each iteration: refit LMM on all markers + recompute BH-FDR (matches pipeline)\n\n")

# Pre-allocate: per-iteration table of (marker, beta, p-raw, FDR)
boot_beta <- matrix(NA_real_, nrow = B, ncol = length(PRIMARY),
                    dimnames = list(NULL, PRIMARY))
boot_fdr  <- matrix(NA_real_, nrow = B, ncol = length(PRIMARY),
                    dimnames = list(NULL, PRIMARY))

t0 <- Sys.time()
for (b in seq_len(B)) {
  # Cluster bootstrap: resample patient IDs with replacement
  sampled_pids <- sample(unique_pids, n_pat, replace = TRUE)

  # Re-build data preserving the within-patient T0/T1 structure
  # (each sampled patient brings all their rows; duplicates get new pseudo-IDs)
  rows_list <- lapply(seq_along(sampled_pids), function(j) {
    pid <- sampled_pids[j]
    sub <- hz[hz$Patient_ID == pid, , drop = FALSE]
    sub$Patient_ID <- paste0(pid, "_b", j)  # unique pseudo-ID
    sub
  })
  hz_boot <- do.call(rbind, rows_list)

  # Fit LMM on all markers (mirror Step 04 pipeline so FDR is comparable)
  res_list <- lapply(ALL_MARKERS, function(m) {
    fit_feature_lmm(hz_boot, m,
                    group_col = "Group", time_col = "Timepoint",
                    id_col = "Patient_ID")
  })
  res_df <- do.call(rbind, res_list)
  # BH-FDR over the full marker set (same as pipeline)
  res_df$FDR <- p.adjust(res_df$P_Value_Interaction, method = "BH")

  for (m in PRIMARY) {
    row <- res_df[res_df$Marker == m, ]
    if (nrow(row) == 1 && !is.na(row$Estimate_Interaction)) {
      boot_beta[b, m] <- row$Estimate_Interaction
      boot_fdr[b, m]  <- row$FDR
    }
  }

  if (b %% 100 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    eta <- elapsed / b * (B - b)
    cat(sprintf("  iter %4d/%d  [%.0fs elapsed, ~%.0fs remaining]\n", b, B, elapsed, eta))
  }
}

cat(sprintf("\nDone in %.0fs.\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# ── Summary table ───────────────────────────────────────────────────────────
cat("\n=== Bootstrap β summary (cluster bootstrap, B=", B, ") ===\n", sep = "")
cat(sprintf("%-12s  %-10s  %-10s  %-22s  %-12s  %-14s\n",
            "Marker", "obs_beta", "boot_med", "boot_95CI", "boot_p", "%FDR<0.05"))
cat(strrep("-", 95), "\n")

# Observed (full-data) betas — use pipeline JSON for accuracy
obs <- list(KI67NAIVE = -1.3946, CD28KI67 = -0.9101, KI67CD4 = -1.0848)

for (m in PRIMARY) {
  b_vec <- boot_beta[, m]; b_vec <- b_vec[!is.na(b_vec)]
  f_vec <- boot_fdr[, m];  f_vec <- f_vec[!is.na(f_vec)]
  if (length(b_vec) == 0) { cat(sprintf("%-12s  NO CONVERGENT MODELS\n", m)); next }

  ci <- quantile(b_vec, c(0.025, 0.975))
  boot_p <- mean(b_vec >= 0)  # one-sided: how often does the effect flip?
  boot_p <- 2 * min(boot_p, 1 - boot_p)  # two-sided bootstrap p
  pct_fdr <- mean(f_vec < 0.05) * 100

  cat(sprintf("%-12s  %+8.3f    %+8.3f    [%+6.3f, %+6.3f]    %.4f       %5.1f%% (%d/%d)\n",
              m, obs[[m]], median(b_vec), ci[1], ci[2], boot_p,
              pct_fdr, sum(f_vec < 0.05), length(f_vec)))
}

# ── Save bootstrap distributions for downstream forest plot ─────────────────
boot_out <- list(
  primary  = PRIMARY,
  obs_beta = obs,
  beta_mat = boot_beta,
  fdr_mat  = boot_fdr,
  B        = B,
  seed     = 2026
)
saveRDS(boot_out, here("diagnostics/diag_10_lmm_bootstrap.rds"))
cat("\nBootstrap distributions saved to diagnostics/diag_10_lmm_bootstrap.rds\n")

# ── Interpretation ──────────────────────────────────────────────────────────
cat("\n=== Interpretation ===\n")
cat("- boot_med: bootstrap median estimate of beta (should agree with obs_beta)\n")
cat("- boot_95CI: 95% bootstrap CI — should exclude zero for robustness\n")
cat("- boot_p: two-sided proportion-based bootstrap p-value\n")
cat("- %FDR<0.05: fraction of resamples in which the marker is still 'significant'\n")
cat("  after BH-correction across the full marker panel — a stricter criterion\n")
cat("  than the standard FDR check on the original data\n")

cat("\nDone.\n")
