# diagnostics/diag_08_monotherapy_subgroup.R
# Diagnostic: LMM in ICI monotherapy subgroup (Pembrolizumab only).
#
# Purpose:
#   Tests whether the KI67 proliferative reset survives in the "pure ICI"
#   subgroup (Pembro monotherapy, n≈24), where cytotoxic chemotherapy is absent.
#   This is the cleanest biological test of an immune-mediated mechanism.
#
#   Expected result (if signal is real): same direction and comparable magnitude
#   of Timepoint x Group interaction, even if underpowered for significance.
#   A null result here could indicate the reset is chemo-dependent.
#
# Output: diagnostics/diag_08_output.txt

library(here)
library(readxl)
library(lmerTest)
library(dplyr)

OUT_FILE <- here("diagnostics/diag_08_output.txt")
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
trt_col <- "Terapia_Pembro_1_Pembro_CHT_2_ipi/nivo_CHT_3;_Atezo_CHT_4"

trt_map <- df_raw_t0 %>%
  select(Patient_ID, Terapia = all_of(trt_col))

df_all <- hz %>%
  left_join(trt_map, by = "Patient_ID") %>%
  mutate(
    Timepoint  = factor(Timepoint, levels = c("T0", "T1")),
    Group      = factor(Group, levels = c("SD_PD", "RP")),
    Patient_ID = factor(Patient_ID)
  )

# ── 2. Subgroup tables ────────────────────────────────────────────────────────
df_mono <- df_all %>% filter(Terapia == 1)

cat("=== ICI Monotherapy Subgroup (Pembro only) ===\n")
cat("Total rows (T0+T1):", nrow(df_mono), "\n")
ids_mono <- unique(df_mono$Patient_ID)
cat("Unique patients:", length(ids_mono), "\n")

meta_mono <- df_mono %>%
  filter(Timepoint == "T0") %>%
  select(Patient_ID, Group) %>%
  distinct()
cat("Group distribution (T0 only):\n")
print(table(meta_mono$Group))
cat("\n")

# ── 3. Fit LMMs ───────────────────────────────────────────────────────────────
MARKERS <- c("KI67NAIVE", "KI67CD4", "CD28KI67", "KI67CM", "KI67EFF",
             "KI67CD8", "KI67")

fit_lmm_both <- function(marker, df_full, df_sub) {
  f <- as.formula(paste0(marker, " ~ Timepoint * Group + (1 | Patient_ID)"))

  m_full <- tryCatch(lmer(f, data = df_full %>% filter(!is.na(.data[[marker]])),
                           REML = FALSE), error = function(e) NULL)
  m_mono <- tryCatch(lmer(f, data = df_sub  %>% filter(!is.na(.data[[marker]])),
                           REML = FALSE), error = function(e) NULL)

  extract <- function(m) {
    if (is.null(m)) return(c(beta = NA, se = NA, p = NA))
    co <- coef(summary(m))
    if (!"TimepointT1:GroupRP" %in% rownames(co)) return(c(beta = NA, se = NA, p = NA))
    row <- co["TimepointT1:GroupRP", ]
    c(beta = unname(row["Estimate"]), se = unname(row["Std. Error"]), p = unname(row["Pr(>|t|)"]))
  }

  list(marker = marker,
       full = extract(m_full),
       mono = extract(m_mono),
       n_mono = sum(!is.na(df_sub[[marker]]) & df_sub$Timepoint == "T0"))
}

results <- lapply(MARKERS, fit_lmm_both,
                  df_full = df_all %>% filter(!is.na(Terapia)),
                  df_sub  = df_mono)

# ── 4. Print comparison table ─────────────────────────────────────────────────
cat("=== Timepoint x Group interaction: full cohort vs Pembro-mono subgroup ===\n")
cat(sprintf("%-12s  %6s  %-8s %-8s %-8s   %-8s %-8s %-8s\n",
            "Marker", "n_mono",
            "beta_full", "SE_full", "p_full",
            "beta_mono", "SE_mono", "p_mono"))
cat(strrep("-", 80), "\n")

for (r in results) {
  sig_full <- ifelse(!is.na(r$full["p"]) & r$full["p"] < 0.05, "*", " ")
  sig_mono <- ifelse(!is.na(r$mono["p"]) & r$mono["p"] < 0.05, "*", " ")
  cat(sprintf("%-12s  %6d  %8.3f %8.3f %7.4f%s  %8.3f %8.3f %7.4f%s\n",
              r$marker, r$n_mono,
              r$full["beta"], r$full["se"], r$full["p"], sig_full,
              r$mono["beta"], r$mono["se"], r$mono["p"], sig_mono))
}

cat("\n* = p < 0.05\n")

# ── 5. Direction consistency check ────────────────────────────────────────────
cat("\n=== Direction consistency (full vs mono) ===\n")
for (r in results) {
  if (any(is.na(c(r$full["beta"], r$mono["beta"])))) next
  same_dir <- sign(r$full["beta"]) == sign(r$mono["beta"])
  cat(sprintf("  %-12s  full=%.3f  mono=%.3f  same_direction=%s\n",
              r$marker, r$full["beta"], r$mono["beta"],
              ifelse(same_dir, "YES", "NO *** DISCORDANT ***")))
}

cat("\n=== Notes ===\n")
cat("- n_mono = patients with T0 observation in monotherapy subgroup\n")
cat("- Full cohort uses only patients with non-NA Terapia (excludes NA=9 patients)\n")
cat("- Underpowered for significance in subgroup; direction is the key metric\n")
cat("- Discordant direction would challenge ICI-specificity of the KI67 signal\n")
