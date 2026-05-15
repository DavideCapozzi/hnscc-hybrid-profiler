# diagnostics/diag_05_threeway_lmm.R
# Diagnostic: three-group LMM (RP vs SD vs PD) for KI67 and CD137 marker families.
#
# Purpose:
#   KI67 block  — tests whether the Timepoint x Group interaction is RP-specific
#                 or graded (RP > SD > PD); reconciles with Gelibter 2024 DCR contrast.
#   CD137 block — tests whether T cell activation (4-1BB) dynamics show
#                 group-specific patterns independent of proliferation.
#
# NOTE: CD137 markers are loaded from raw Excel files, logit-transformed and
#   z-scored inline because the current pipeline RDS predates their addition.
#   Re-run main.R to get pipeline-integrated CD137 z-scores.
#
# Output: diagnostics/diag_05_output.txt

library(readxl)
library(lmerTest)
library(dplyr)
library(tidyr)

OUT_FILE <- "diagnostics/diag_05_output.txt"
con <- file(OUT_FILE, open = "wt")
sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)

EPSILON <- 1e-4
logit_pct <- function(x) {
  p <- as.numeric(x) / 100
  p <- pmax(pmin(p, 1 - EPSILON), EPSILON)
  log(p / (1 - p))
}

# ── 1. Load processed longitudinal data (KI67 z-scores + QC'd patient list) ───
DATA_L <- readRDS(
  "results/BestResponse_2v3_4/01_data_processing/data_processed_BestResponse_2v3_4_longitudinal.rds"
)
df_raw_t0 <- read_excel("data/Dati_NSCLC_standardizzati_anonimi_T0.xlsx")
df_raw_t1 <- read_excel("data/Dati_NSCLC_standardizzati_anonimi_T1.xlsx")

meta <- DATA_L$metadata
hz   <- DATA_L$hybrid_data_z

# 3-group response coding from raw T0 (PD = reference level)
resp_map <- df_raw_t0 %>%
  select(Patient_ID,
         resp_code = `Best_response_rate_RC_1_RP_2_SD_3_PD_4`) %>%
  mutate(Group3 = factor(
    case_when(resp_code == 2 ~ "RP",
              resp_code == 3 ~ "SD",
              resp_code == 4 ~ "PD"),
    levels = c("PD", "SD", "RP")
  ))

# ── 2. KI67 long data (pipeline z-scores from RDS) ────────────────────────────
markers_ki67 <- c("KI67NAIVE", "KI67CD4", "CD28KI67",
                  "KI67CD8", "KI67CM", "KI657EMRA", "KI67EFF", "KI67")

df_ki67 <- hz %>%
  select(Patient_ID, Timepoint, any_of(markers_ki67)) %>%
  left_join(resp_map %>% select(Patient_ID, Group3), by = "Patient_ID") %>%
  mutate(Timepoint = factor(Timepoint, levels = c("T0", "T1"))) %>%
  filter(!is.na(Group3))

# ── 3. CD137 long data (raw Excel → logit → global z-score) ──────────────────
markers_cd137 <- c("CD137tot", "CD137CD4", "CD137CD8",
                   "CD137KI67", "CD137PD1",
                   "CD137CM", "CD137EFF", "CD137EMRA", "CD137NAIVE")

valid_patients <- unique(hz$Patient_ID)

build_cd137_long <- function(df_t0, df_t1, markers, valid_pts) {
  avail <- intersect(markers, intersect(colnames(df_t0), colnames(df_t1)))
  if (length(avail) == 0) return(NULL)

  t0 <- df_t0 %>%
    filter(Patient_ID %in% valid_pts) %>%
    select(Patient_ID, all_of(avail)) %>%
    mutate(Timepoint = "T0")
  t1 <- df_t1 %>%
    filter(Patient_ID %in% valid_pts) %>%
    select(Patient_ID, all_of(avail)) %>%
    mutate(Timepoint = "T1")

  long <- bind_rows(t0, t1) %>%
    pivot_longer(all_of(avail), names_to = "marker", values_to = "raw_pct") %>%
    mutate(logit_val = logit_pct(raw_pct)) %>%
    group_by(marker) %>%
    mutate(z = (logit_val - mean(logit_val, na.rm = TRUE)) /
                sd(logit_val,   na.rm = TRUE)) %>%
    ungroup() %>%
    select(Patient_ID, Timepoint, marker, z) %>%
    pivot_wider(names_from = marker, values_from = z)

  long
}

df_cd137_wide <- build_cd137_long(df_raw_t0, df_raw_t1, markers_cd137, valid_patients)

if (!is.null(df_cd137_wide)) {
  df_cd137 <- df_cd137_wide %>%
    left_join(resp_map %>% select(Patient_ID, Group3), by = "Patient_ID") %>%
    mutate(Timepoint = factor(Timepoint, levels = c("T0", "T1"))) %>%
    filter(!is.na(Group3))
  cat(sprintf("[INFO] CD137 markers available: %s\n",
              paste(intersect(markers_cd137, colnames(df_cd137)), collapse = ", ")))
} else {
  df_cd137 <- NULL
  cat("[WARN] No CD137 markers found in raw Excel — check column names.\n")
}

# ── 4. Sample sizes ───────────────────────────────────────────────────────────
cat("\n=== Sample sizes (from pipeline QC) ===\n")
df_ki67 %>% distinct(Patient_ID, Group3) %>% count(Group3) %>% print()

# ── 5. LMM helpers ───────────────────────────────────────────────────────────
run_lmm3 <- function(marker, data) {
  d <- data %>%
    select(Patient_ID, Timepoint, Group3, y = all_of(marker)) %>%
    drop_na(y)
  fit <- tryCatch(
    lmer(y ~ Timepoint * Group3 + (1 | Patient_ID), data = d, REML = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  coefs   <- summary(fit)$coefficients
  rp_row  <- grep("TimepointT1:Group3RP", rownames(coefs), value = TRUE)
  sd_row  <- grep("TimepointT1:Group3SD", rownames(coefs), value = TRUE)
  extract <- function(row) {
    if (length(row) == 0 || !row %in% rownames(coefs))
      return(c(beta = NA, se = NA, t = NA, p = NA))
    r <- coefs[row, ]
    c(beta = unname(r["Estimate"]), se = unname(r["Std. Error"]),
      t    = unname(r["t value"]),  p  = unname(r["Pr(>|t|)"]))
  }
  list(marker = marker, n_obs = nrow(d),
       rp_vs_pd = extract(rp_row), sd_vs_pd = extract(sd_row),
       singular = isSingular(fit))
}

run_lmm_binary <- function(marker, data, group_var) {
  d <- data %>%
    select(Patient_ID, Timepoint, Group = all_of(group_var), y = all_of(marker)) %>%
    drop_na(y)
  fit <- tryCatch(
    lmer(y ~ Timepoint * Group + (1 | Patient_ID), data = d, REML = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  coefs <- summary(fit)$coefficients
  row   <- grep("TimepointT1:Group", rownames(coefs), value = TRUE)[1]
  if (is.na(row)) return(NULL)
  r <- coefs[row, ]
  list(marker = marker, beta = r["Estimate"],
       t = r["t value"], p = r["Pr(>|t|)"],
       singular = isSingular(fit))
}

# ── 6. Generic analysis runner ────────────────────────────────────────────────
run_marker_analysis <- function(markers, label, df_base) {
  avail <- intersect(markers, colnames(df_base))
  if (length(avail) == 0) {
    cat(sprintf("\n[SKIP] %s — no markers in data.\n", label))
    return(invisible(NULL))
  }

  cat(sprintf("\n\n══════════════════════════════════════════════════════\n"))
  cat(sprintf("  %s\n", label))
  cat(sprintf("══════════════════════════════════════════════════════\n"))

  # --- 3-group LMM ---
  cat("\n--- Three-group LMM: Timepoint x Group (reference = PD) ---\n")
  cat(sprintf("%-14s  %7s %6s %6s  |  %7s %6s %6s  | Singular\n",
              "Marker", "RP_β", "RP_t", "RP_p", "SD_β", "SD_t", "SD_p"))
  cat(strrep("-", 74), "\n")
  for (mk in avail) {
    r <- run_lmm3(mk, df_base)
    if (is.null(r)) next
    cat(sprintf("%-14s  %+7.3f %6.2f %6.4f  |  %+7.3f %6.2f %6.4f  | %s\n",
                r$marker,
                r$rp_vs_pd["beta"], r$rp_vs_pd["t"], r$rp_vs_pd["p"],
                r$sd_vs_pd["beta"], r$sd_vs_pd["t"], r$sd_vs_pd["p"],
                if (r$singular) "YES" else "no"))
  }

  # --- DCR replication ---
  cat("\n--- DCR replication: RP+SD vs PD (Gelibter 2024 contrast) ---\n")
  df_dcr <- df_base %>%
    mutate(Group_DCR = factor(
      ifelse(Group3 %in% c("RP", "SD"), "DCR", "PD"),
      levels = c("PD", "DCR")
    ))
  cat(sprintf("%-14s  %7s %6s %6s\n", "Marker", "DCR_β", "DCR_t", "DCR_p"))
  cat(strrep("-", 40), "\n")
  for (mk in avail) {
    r <- run_lmm_binary(mk, df_dcr, "Group_DCR")
    if (is.null(r)) next
    cat(sprintf("%-14s  %+7.3f %6.2f %6.4f\n", r$marker, r$beta, r$t, r$p))
  }

  # --- Marginal means ---
  cat("\n--- Marginal mean change T0→T1 by group (mean ± SE of delta) ---\n")
  cat(sprintf("%-14s  %10s %10s %10s\n", "Marker", "RP", "SD", "PD"))
  cat(strrep("-", 50), "\n")
  for (mk in avail) {
    deltas <- df_base %>%
      select(Patient_ID, Timepoint, Group3, y = all_of(mk)) %>%
      drop_na(y) %>%
      pivot_wider(names_from = Timepoint, values_from = y) %>%
      drop_na(T0, T1) %>%
      mutate(delta = T1 - T0) %>%
      group_by(Group3) %>%
      summarise(m = mean(delta), se = sd(delta) / sqrt(n()), .groups = "drop")
    fmt <- function(g) {
      row <- deltas[deltas$Group3 == g, ]
      if (nrow(row) == 0) return("  NA")
      sprintf("%+.2f(%.2f)", row$m, row$se)
    }
    cat(sprintf("%-14s  %10s %10s %10s\n", mk, fmt("RP"), fmt("SD"), fmt("PD")))
  }
}

# ── 7. Run both families ──────────────────────────────────────────────────────
run_marker_analysis(markers_ki67, "KI67 proliferation subfamily", df_ki67)

if (!is.null(df_cd137)) {
  run_marker_analysis(markers_cd137, "CD137 (4-1BB) activation subfamily", df_cd137)
} else {
  cat("\n[SKIP] CD137 block — raw Excel columns not found.\n")
}
