# clinical-onco-profiler

An R-based, config-driven pipeline for immunological profiling of oncology cohorts. Designed for hybrid immunophenotyping data (FACS flow cytometry + soluble biomarkers), it performs:

- **Cross-sectional analysis** (Steps 01–03): QC → visualization → sPLS-DA multivariate modeling
- **Longitudinal analysis** (Steps 01 + 04): joint T0/T1 processing → Linear Mixed Models (LMM) with Leave-One-Out sensitivity and covariate adjustment
- **Differential network analysis** (Step 05, optional): bootstrap partial correlation baselines per clinical group → permutation-tested differential overlay → Cytoscape-ready CSV exports
- **Machine learning classification** (Step 06): nested-LOOCV Elastic Net + SVM-RBF classifier gated on LMM LOO-robust features

Applied to NSCLC patients treated with immune checkpoint inhibitors (ICI) to identify response biomarkers. Converging evidence across three analytical layers (sPLS-DA, LMM, ML) points to a KI67 proliferative reset signature discriminating responders from non-responders.

---

## Requirements

- [conda](https://docs.conda.io) or [micromamba](https://mamba.readthedocs.io)
- R 4.3.x (managed via the conda environment)

---

## Installation

```bash
conda env create -f env/environment.yml
conda activate clinical-onco-profiler
```

Or with micromamba:

```bash
micromamba env create -f env/environment.yml
micromamba activate clinical-onco-profiler
```

The environment pins R 4.3.x with Bioconductor packages (`mixOmics`, `pcaMethods`, `speckle`, `ALDEx2`) and key CRAN packages (`lmerTest`, `lme4`, `glmnet`, `e1071`, `pROC`, `igraph`).

---

## Usage

```bash
# Full pipeline (reads config/global_params.yml)
Rscript main.R

# With a custom config
Rscript main.R config/global_params_hnscc.yml
```

The pipeline auto-detects whether the config is **nested multi-experiment** (has an `experiments:` key) or **flat single-cohort**, then runs each experiment through two independent passes (standard cross-sectional + longitudinal), followed by an optional ML pass.

---

## Project structure

```
config/                    Pipeline configuration (YAML)
  global_params.yml        Primary multi-experiment config
  global_params_hnscc.yml  Flat single-cohort config (legacy)

data/                      Input Excel files (not tracked in git — add your own)

diagnostics/               Supplementary stand-alone analyses
  diag_01_lmm_covariates.R LMM with clinical covariates (PD-L1, PS, etc.)
  diag_02_ki67_composite.R KI67 subfamily PCA as classifier feature
  diag_03_nested_lmm_loocv.R Fully-nested LOO (LMM gate inside loop)
  diag_04_survival.R       OS/PFS survival analysis with KI67 markers
  diag_05_threeway_lmm.R   Three-group LMM: RP vs SD vs PD

env/
  environment.yml          Conda environment specification

R/                         Reusable module library (sourced by main.R at startup)
  modules_qc.R             PCA outlier detection, QC pipeline
  modules_coda.R           Compositional data transformation, BPCA imputation
  modules_multivariate.R   sPLS-DA fitting and extraction
  modules_longitudinal.R   LMM fitting, LOO sensitivity, R² (Nakagawa 2013)
  modules_network.R        Bootstrap partial correlations, differential overlay
  modules_ml.R             Nested-LOOCV classifiers, collinearity filter, metrics
  modules_viz.R            ggplot rendering functions
  utils_io.R               Config loading, QC report writing

src/                       Pipeline step scripts (sourced in order by main.R)
  01_data_processing.R
  02_visualization.R
  03_statistical_analysis.R
  04_longitudinal_lmm.R
  05_network_analysis.R
  06_machine_learning.R

results/                   All pipeline outputs (not tracked in git)
  <experiment_name>/
    01_data_processing/    Processed .rds objects
    02_visualization/      Distribution and PCA PDFs
    03_statistical_analysis/ sPLS-DA results, JSON drivers
    04_longitudinal_lmm/   LMM tables (Excel + JSON), volcano/trajectory PDFs
    05_network_analysis/   Cytoscape CSVs, network PDFs
    06_machine_learning/   Classifier results (Excel, PDFs, JSON)
```

---

## Configuration

All pipeline parameters are controlled by `config/global_params.yml`. Key sections:

| Section | Description |
|---------|-------------|
| `experiments` | Named experiment blocks (override clinical mapping, input files, flags) |
| `features.facs` / `features.soluble` | Marker panel definitions |
| `clinical.target_column` | Column encoding the response group |
| `clinical.mapping` | Maps raw integer codes to responder/non-responder labels |
| `qc` | Missingness thresholds, outlier detection settings |
| `network` | Bootstrap/permutation counts, edge thresholds, Cytoscape export filters |
| `machine_learning` | FDR/LOO thresholds, collinearity filter, CV grids |

See `config/global_params.yml` for inline documentation of each parameter.

---

## Input data format

Input files are Excel workbooks (`.xlsx`) with one row per patient, one column per marker. Required columns: `Patient_ID` plus all marker columns listed under `features.facs` and `features.soluble` in the config. Longitudinal analyses expect separate T0 and T1 files with matching `Patient_ID` values.

Input files are **not tracked in git** (the `data/` directory is gitignored) to protect patient data.

---

## Supplementary diagnostics

Scripts in `diagnostics/` are stand-alone analyses that extend or validate the main pipeline results. They read processed `.rds` outputs from `results/` and can be run independently:

```bash
Rscript diagnostics/diag_01_lmm_covariates.R
Rscript diagnostics/diag_04_survival.R
```

These scripts require the main pipeline to have been run first (so that `results/` outputs exist).

---

## Citation

If you use this pipeline, please cite:

> [Manuscript in preparation] — Capozzi D. et al., *Multi-layer immunological profiling identifies a KI67 proliferative reset signature associated with ICI response in NSCLC*. 2025.

---

## License

This project is released for academic use. Contact the authors for other uses.
