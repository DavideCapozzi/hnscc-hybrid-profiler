# R/modules_longitudinal.R
# ==============================================================================
# LONGITUDINAL LMM MODULE
# Description: Wrapper functions for Linear Mixed Models handling time & group.
# Dependencies: lmerTest, dplyr, ggplot2, ggrepel
# ==============================================================================

library(dplyr)
library(ggplot2)
library(ggrepel)

#' @title Marginal and Conditional R² for LMM (Nakagawa & Schielzeth 2013)
#' @description Computes variance-partition R² without requiring external packages.
#' @param mod A fitted lmerMod object.
#' @return Named list with R2m (marginal, fixed effects only) and R2c (conditional).
r2_nakagawa <- function(mod) {
  tryCatch({
    sigma2_f <- var(predict(mod, re.form = NA))
    vc       <- as.data.frame(lme4::VarCorr(mod))
    sigma2_r <- sum(vc$vcov[vc$grp != "Residual"])
    sigma2_e <- sigma(mod)^2
    total    <- sigma2_f + sigma2_r + sigma2_e
    list(R2m = round(sigma2_f / total, 4),
         R2c = round((sigma2_f + sigma2_r) / total, 4))
  }, error = function(e) list(R2m = NA_real_, R2c = NA_real_))
}

#' @title Run Linear Mixed Model on a single feature
#' @description Fits LMM with an interaction term between Time and Clinical Group.
#' Includes internal standard deviation scaling and Marginal Time Effect extraction.
#' @param data_long Dataframe in long format.
#' @param feature String. Name of the column containing feature values.
#' @param group_col String. Name of the column containing clinical labels.
#' @param time_col String. Name of the column containing timepoints.
#' @param id_col String. Name of the column containing Patient IDs.
#' @param covariates Vector of strings. Names of clinical covariates to include.
#' @return A one-row dataframe with model statistics.
fit_feature_lmm <- function(data_long, feature, group_col = "Group", 
                            time_col = "Timepoint", id_col = "Patient_ID",
                            covariates = NULL) {
  
  if (!requireNamespace("lmerTest", quietly = TRUE)) {
    stop("Package 'lmerTest' is required for p-value calculation in LMM.")
  }
  
  raw_vals <- as.numeric(data_long[[feature]])
  val_sd <- sd(raw_vals, na.rm = TRUE)
  if (is.na(val_sd) || val_sd == 0) val_sd <- 1 # Safety fallback for zero variance
  
  df_list <- list(
    Value = raw_vals / val_sd, # Internal scaling for optimizer convergence
    Time = as.factor(data_long[[time_col]]),
    Group = as.factor(data_long[[group_col]]),
    ID = as.factor(data_long[[id_col]])
  )
  
  valid_covs <- c()
  if (!is.null(covariates) && length(covariates) > 0) {
    valid_covs <- intersect(covariates, colnames(data_long))
    for (cov in valid_covs) {
      df_list[[cov]] <- data_long[[cov]]
    }
  }
  
  df_model <- as.data.frame(df_list, check.names = FALSE)
  df_model <- df_model[complete.cases(df_model), ]
  
  result <- data.frame(
    Marker = feature,
    Estimate_Interaction = NA,
    T_Value_Interaction = NA,
    Std_Error = NA,
    P_Value_Interaction = NA,
    Estimate_Time_Main = NA,
    P_Value_Time_Main = NA,
    Model_Converged = FALSE,
    Is_Singular = NA,
    N_Observations = nrow(df_model),
    R2m = NA_real_,
    R2c = NA_real_
  )
  
  if (nrow(df_model) < 10) return(result)
  
  tryCatch({
    # --- MODEL 1: Interaction Model (Primary Objective) ---
    formula_int <- "Value ~ Time * Group"
    if (length(valid_covs) > 0) formula_int <- paste(formula_int, "+", paste(sprintf("`%s`", valid_covs), collapse = " + "))
    formula_int <- paste(formula_int, "+ (1 | ID)")
    
    mod_int <- suppressMessages(suppressWarnings(
      lmerTest::lmer(as.formula(formula_int), data = df_model, 
                     REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE))
    ))
    
    result$Is_Singular <- lme4::isSingular(mod_int)
    
    coef_table_int <- summary(mod_int)$coefficients
    interaction_idx <- grep("Time.*:Group", rownames(coef_table_int))
    
    if (length(interaction_idx) == 1) {
      # Reverse scaling to restore native hybrid logit/log2 effect sizes
      result$Estimate_Interaction <- coef_table_int[interaction_idx, "Estimate"] * val_sd
      result$Std_Error <- coef_table_int[interaction_idx, "Std. Error"] * val_sd
      
      t_col <- grep("t value", colnames(coef_table_int))
      if (length(t_col) == 1) result$T_Value_Interaction <- coef_table_int[interaction_idx, t_col]
      
      p_col <- grep("Pr\\(>\\|t\\|\\)", colnames(coef_table_int))
      if (length(p_col) == 1) {
        result$P_Value_Interaction <- coef_table_int[interaction_idx, p_col]
        result$Model_Converged <- TRUE
      }
    }

    r2 <- r2_nakagawa(mod_int)
    result$R2m <- r2$R2m
    result$R2c <- r2$R2c
    
    # --- MODEL 2: Marginal Time Model (Positive Control) ---
    formula_time <- "Value ~ Time"
    if (length(valid_covs) > 0) formula_time <- paste(formula_time, "+", paste(sprintf("`%s`", valid_covs), collapse = " + "))
    formula_time <- paste(formula_time, "+ (1 | ID)")
    
    mod_time <- suppressMessages(suppressWarnings(
      lmerTest::lmer(as.formula(formula_time), data = df_model, 
                     REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE))
    ))
    
    coef_table_time <- summary(mod_time)$coefficients
    time_idx <- grep("Time", rownames(coef_table_time))
    
    if (length(time_idx) == 1) {
      result$Estimate_Time_Main <- coef_table_time[time_idx, "Estimate"]
      p_col <- grep("Pr\\(>\\|t\\|\\)", colnames(coef_table_time))
      if (length(p_col) == 1) result$P_Value_Time_Main <- coef_table_time[time_idx, p_col]
    }
    
  }, error = function(e) {
    # Silent fail on calculation error, returning NA defaults
  })
  
  return(result)
}

#' @title Run Leave-One-Out (LOO) Sensitivity Analysis
#' @description Iteratively drops one patient at a time to test interaction robustness.
#' @return Numeric. The maximum (worst-case) interaction P-value found across all iterations.
run_loo_sensitivity <- function(data_long, feature, group_col = "Group", 
                                time_col = "Timepoint", id_col = "Patient_ID",
                                covariates = NULL) {
  
  unique_ids <- unique(data_long[[id_col]])
  max_p_val <- 0
  
  for (drop_id in unique_ids) {
    # Create LOO dataset
    df_subset <- data_long[data_long[[id_col]] != drop_id, ]
    
    # Fit model silently
    res <- fit_feature_lmm(data_long = df_subset, feature = feature, 
                           group_col = group_col, time_col = time_col, 
                           id_col = id_col, covariates = covariates)
    
    # Track the worst p-value
    current_p <- res$P_Value_Interaction
    if (!is.na(current_p) && current_p > max_p_val) {
      max_p_val <- current_p
    }
  }
  
  # max_p_val stays 0 only when every LOO fold returned NA (all models failed).
  # Return NA so the gate filter (!is.na) correctly excludes the marker.
  return(if (max_p_val == 0) NA else max_p_val)
}

#' @title Cluster Bootstrap CI for LMM Interaction Betas
#' @description
#' Patient-level cluster bootstrap (resample patient IDs with replacement,
#' preserving the within-patient T0/T1 pairing). For each bootstrap iteration
#' refits the full LMM panel and re-applies BH-FDR, then reports for the
#' supplied target markers:
#'   - bootstrap median estimate
#'   - 2.5% / 97.5% bootstrap quantiles (95% CI)
#'   - two-sided proportion-based bootstrap p-value
#'   - fraction of resamples in which the marker still passes FDR < fdr_threshold
#'
#' @param data_long Dataframe in long format (one row per (Patient_ID, Timepoint)).
#' @param all_markers Character vector of all markers fitted in the primary LMM
#'   panel (needed for BH-FDR re-computation inside each bootstrap iteration).
#' @param target_markers Character vector of markers to report results for
#'   (typically the FDR-significant subset from the primary run).
#' @param group_col,time_col,id_col Column names in data_long.
#' @param covariates Optional clinical covariates passed to fit_feature_lmm.
#' @param n_boot Integer. Number of bootstrap iterations (default 500).
#' @param fdr_threshold Numeric. FDR cutoff for the per-iteration significance
#'   tally (default 0.05).
#' @param seed Integer. RNG seed for reproducibility (default 2026).
#' @param progress_message Logical. Print every 100 iterations (default TRUE).
#' @return A list with:
#'   - summary_df: one row per target marker (Marker, Median_Beta_Boot,
#'                 CI_Lower_2.5, CI_Upper_97.5, Bootstrap_P, Pct_FDR_Significant,
#'                 N_Valid_Iterations)
#'   - beta_matrix: B x length(target_markers) matrix of bootstrap betas
#'   - fdr_matrix:  B x length(target_markers) matrix of bootstrap FDR values
#'   - n_boot, seed, fdr_threshold
run_lmm_bootstrap_ci <- function(data_long,
                                 all_markers,
                                 target_markers,
                                 group_col       = "Group",
                                 time_col        = "Timepoint",
                                 id_col          = "Patient_ID",
                                 covariates      = NULL,
                                 n_boot          = 500L,
                                 fdr_threshold   = 0.05,
                                 seed            = 2026L,
                                 progress_message = TRUE) {

  target_markers <- intersect(target_markers, all_markers)
  if (length(target_markers) == 0) {
    return(list(summary_df = data.frame(), beta_matrix = NULL, fdr_matrix = NULL,
                n_boot = 0L, seed = seed, fdr_threshold = fdr_threshold))
  }

  unique_pids <- unique(data_long[[id_col]])
  n_pat       <- length(unique_pids)

  beta_mat <- matrix(NA_real_, nrow = n_boot, ncol = length(target_markers),
                     dimnames = list(NULL, target_markers))
  fdr_mat  <- matrix(NA_real_, nrow = n_boot, ncol = length(target_markers),
                     dimnames = list(NULL, target_markers))

  set.seed(seed)
  t0 <- Sys.time()

  for (b in seq_len(n_boot)) {
    sampled_pids <- sample(unique_pids, n_pat, replace = TRUE)

    rows_list <- lapply(seq_along(sampled_pids), function(j) {
      pid <- sampled_pids[j]
      sub <- data_long[data_long[[id_col]] == pid, , drop = FALSE]
      sub[[id_col]] <- paste0(pid, "_b", j)  # unique pseudo-ID per draw
      sub
    })
    data_boot <- do.call(rbind, rows_list)

    res_list <- lapply(all_markers, function(mk) {
      fit_feature_lmm(data_long = data_boot, feature = mk,
                      group_col = group_col, time_col = time_col,
                      id_col = id_col, covariates = covariates)
    })
    res_df <- do.call(rbind, res_list)
    res_df$FDR <- p.adjust(res_df$P_Value_Interaction, method = "BH")

    for (m in target_markers) {
      row <- res_df[res_df$Marker == m, , drop = FALSE]
      if (nrow(row) == 1 && !is.na(row$Estimate_Interaction)) {
        beta_mat[b, m] <- row$Estimate_Interaction
        fdr_mat[b, m]  <- row$FDR
      }
    }

    if (progress_message && b %% 100 == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      eta     <- elapsed / b * (n_boot - b)
      message(sprintf("      [Bootstrap] iter %d/%d  [%.0fs elapsed, ~%.0fs remaining]",
                      b, n_boot, elapsed, eta))
    }
  }

  summary_rows <- lapply(target_markers, function(m) {
    b_vec <- beta_mat[, m]; b_vec <- b_vec[!is.na(b_vec)]
    f_vec <- fdr_mat[, m];  f_vec <- f_vec[!is.na(f_vec)]
    if (length(b_vec) == 0) {
      return(data.frame(Marker = m, Median_Beta_Boot = NA_real_,
                        CI_Lower_2.5 = NA_real_, CI_Upper_97.5 = NA_real_,
                        Bootstrap_P = NA_real_, Pct_FDR_Significant = NA_real_,
                        N_Valid_Iterations = 0L, stringsAsFactors = FALSE))
    }
    ci      <- quantile(b_vec, c(0.025, 0.975))
    one_sided <- mean(b_vec >= 0)
    two_sided <- 2 * min(one_sided, 1 - one_sided)
    pct_fdr <- mean(f_vec < fdr_threshold) * 100
    data.frame(
      Marker             = m,
      Median_Beta_Boot   = round(median(b_vec), 4),
      CI_Lower_2.5       = round(ci[1], 4),
      CI_Upper_97.5      = round(ci[2], 4),
      Bootstrap_P        = round(two_sided, 4),
      Pct_FDR_Significant = round(pct_fdr, 1),
      N_Valid_Iterations = length(b_vec),
      stringsAsFactors   = FALSE
    )
  })
  summary_df <- do.call(rbind, summary_rows)

  list(
    summary_df    = summary_df,
    beta_matrix   = beta_mat,
    fdr_matrix    = fdr_mat,
    n_boot        = n_boot,
    seed          = seed,
    fdr_threshold = fdr_threshold
  )
}


#' @title Plot Longitudinal Trajectories (Spaghetti + Boxplot)
#' @description Visualizes patient trajectories over time, split by group.
#' @param data_long Dataframe in long format.
#' @param feature String. Marker name.
#' @param group_col String.
#' @param time_col String.
#' @param id_col String.
#' @param colors Named vector of colors.
#' @param p_val Optional numeric. P-value or FDR to display in subtitle.
#' @return A ggplot object.
plot_lmm_trajectories <- function(data_long, feature, group_col = "Group", 
                                  time_col = "Timepoint", id_col = "Patient_ID",
                                  colors = NULL, p_val = NULL) {
  
  require(ggplot2)
  
  # Ensure NAs are dropped for pure visualization
  plot_df <- data_long[!is.na(data_long[[feature]]), ]
  
  sub_title <- "Patient trajectories over time"
  if (!is.null(p_val)) {
    metric_name <- if (p_val < 0.05) "FDR" else "Interaction P-Value"
    sub_title <- sprintf("%s: %.4f", metric_name, p_val)
  }
  
  p <- ggplot(plot_df, aes(x = .data[[time_col]], y = .data[[feature]], fill = .data[[group_col]])) +
    
    # Base Boxplot distribution
    geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.4) +
    
    # Spaghetti Lines (Connecting individual patients)
    geom_line(aes(group = .data[[id_col]], color = .data[[group_col]]), alpha = 0.3, linewidth = 0.6) +
    
    # Individual Points
    geom_point(aes(fill = .data[[group_col]]), shape = 21, size = 2.5, color = "white", stroke = 0.3) +
    
    # Split by clinical group
    facet_wrap(as.formula(paste("~", group_col))) +
    
    # Aesthetics scaling
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    
    labs(
      title = paste("Trajectory:", feature),
      subtitle = sub_title,
      x = "Timepoint", 
      y = "Expression Level (Hybrid Scale)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "gray95"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40")
    )
  
  return(p)
}

#' @title Volcano Plot for LMM Results
#' @description Plots the t-statistic vs -log10(FDR) of the interaction terms.
#' @param results_df Dataframe output from LMM loop.
#' @param title Plot title.
#' @return A ggplot object.
plot_lmm_volcano <- function(results_df, title = "Longitudinal Interaction (LMM)") {
  
  plot_df <- results_df %>% filter(!is.na(P_Value_Interaction), !is.na(T_Value_Interaction))
  if (nrow(plot_df) == 0) return(NULL)
  
  plot_df$logP <- -log10(plot_df$P_Value_Interaction)
  plot_df$logFDR <- -log10(plot_df$FDR_Interaction)
  
  plot_df$Significance <- "Not Significant"
  plot_df$Significance[plot_df$FDR_Interaction < 0.05 & plot_df$T_Value_Interaction > 0] <- "Positive Interaction"
  plot_df$Significance[plot_df$FDR_Interaction < 0.05 & plot_df$T_Value_Interaction < 0] <- "Negative Interaction"
  
  if (all(plot_df$Significance == "Not Significant")) {
    plot_df$Significance[plot_df$P_Value_Interaction < 0.05 & plot_df$T_Value_Interaction > 0] <- "Positive (Raw P < 0.05)"
    plot_df$Significance[plot_df$P_Value_Interaction < 0.05 & plot_df$T_Value_Interaction < 0] <- "Negative (Raw P < 0.05)"
  }
  
  colors <- c(
    "Not Significant" = "gray70",
    "Positive Interaction" = "#B2182B",
    "Negative Interaction" = "#2166AC",
    "Positive (Raw P < 0.05)" = "#F4A582",
    "Negative (Raw P < 0.05)" = "#92C5DE"
  )
  
  p <- ggplot(plot_df, aes(x = T_Value_Interaction, y = logFDR, color = Significance)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred") +
    geom_point(size = 3, alpha = 0.8) +
    geom_text_repel(
      data = subset(plot_df, Significance != "Not Significant" | P_Value_Interaction < 0.01),
      aes(label = Marker), size = 3, show.legend = FALSE, max.overlaps = 20
    ) +
    scale_color_manual(values = colors) +
    labs(
      title = title,
      subtitle = "X: t-statistic of Time:Group Interaction | Y: -log10(FDR)",
      x = "LMM Interaction (t-statistic)",
      y = "-log10(FDR)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  return(p)
}


#' @title Forest Plot of Bootstrap LMM Beta Estimates
#' @description
#' Renders a horizontal forest plot of the bootstrap-derived interaction beta
#' estimates with 95% bootstrap CIs. Each row is one target marker; the
#' observed (point) estimate from the primary LMM is plotted as a coloured
#' point, the bootstrap CI as a horizontal segment, and the bootstrap median
#' as a small dashed vertical tick. A reference line at beta=0 marks the null.
#'
#' @param boot_summary Dataframe from run_lmm_bootstrap_ci()$summary_df. Must
#'   contain Marker, Median_Beta_Boot, CI_Lower_2.5, CI_Upper_97.5,
#'   Pct_FDR_Significant.
#' @param observed_df Optional dataframe with columns Marker and
#'   Estimate_Interaction (primary LMM point estimate). If provided, plotted
#'   alongside the bootstrap median.
#' @param title Plot title.
#' @return A ggplot object, or NULL if no rows.
plot_lmm_forest <- function(boot_summary, observed_df = NULL,
                            title = "Bootstrap 95% CI of LMM Interaction Betas") {
  if (is.null(boot_summary) || nrow(boot_summary) == 0) return(NULL)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)

  plot_df <- boot_summary
  if (!is.null(observed_df) && all(c("Marker", "Estimate_Interaction") %in% colnames(observed_df))) {
    plot_df <- merge(plot_df,
                     observed_df[, c("Marker", "Estimate_Interaction")],
                     by = "Marker", all.x = TRUE)
  } else {
    plot_df$Estimate_Interaction <- plot_df$Median_Beta_Boot
  }

  # Order by observed effect size (ascending — most negative at top)
  plot_df$Marker <- factor(plot_df$Marker,
                           levels = plot_df$Marker[order(plot_df$Estimate_Interaction)])

  plot_df$Label <- sprintf("%s  (%.0f%% FDR<0.05)", plot_df$Marker, plot_df$Pct_FDR_Significant)

  p <- ggplot(plot_df, aes(y = Marker)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = CI_Lower_2.5, xmax = CI_Upper_97.5),
                   height = 0.25, color = "gray30") +
    geom_point(aes(x = Median_Beta_Boot), shape = 4, size = 2.5, color = "gray40") +
    geom_point(aes(x = Estimate_Interaction), shape = 19, size = 3.5, color = "#B2182B") +
    geom_text(aes(x = CI_Upper_97.5,
                  label = sprintf("  %.0f%% FDR", Pct_FDR_Significant)),
              hjust = 0, size = 3.2, color = "gray30") +
    labs(
      title    = title,
      subtitle = "Red dot: observed beta | X: bootstrap median | Bar: 95% bootstrap CI",
      x        = "Time x Group interaction (beta, hybrid scale)",
      y        = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      panel.grid.minor = element_blank(),
      axis.text.y   = element_text(size = 11)
    ) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(8, 60, 8, 8))

  return(p)
}