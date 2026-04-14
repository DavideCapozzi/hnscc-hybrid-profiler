# R/modules_longitudinal.R
# ==============================================================================
# LONGITUDINAL LMM MODULE
# Description: Wrapper functions for Linear Mixed Models handling time & group.
# Dependencies: lmerTest, dplyr, ggplot2, ggrepel
# ==============================================================================

library(dplyr)
library(ggplot2)
library(ggrepel)

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
  
  df_model <- as.data.frame(df_list)
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
    N_Observations = nrow(df_model)
  )
  
  if (nrow(df_model) < 10) return(result)
  
  tryCatch({
    # --- MODEL 1: Interaction Model (Primary Objective) ---
    formula_int <- "Value ~ Time * Group"
    if (length(valid_covs) > 0) formula_int <- paste(formula_int, "+", paste(valid_covs, collapse = " + "))
    formula_int <- paste(formula_int, "+ (1 | ID)")
    
    mod_int <- suppressMessages(suppressWarnings(
      lmerTest::lmer(as.formula(formula_int), data = df_model, 
                     REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE))
    ))
    
    result$Is_Singular <- lme4::isSingular(mod_int)
    
    coef_table_int <- summary(mod_int)$coefficients
    interaction_idx <- grep("Time.*:Group", rownames(coef_table_int))
    
    if (length(interaction_idx) == 1) {
      result$Estimate_Interaction <- coef_table_int[interaction_idx, "Estimate"]
      result$Std_Error <- coef_table_int[interaction_idx, "Std. Error"]
      
      t_col <- grep("t value", colnames(coef_table_int))
      if (length(t_col) == 1) result$T_Value_Interaction <- coef_table_int[interaction_idx, t_col]
      
      p_col <- grep("Pr\\(>\\|t\\|\\)", colnames(coef_table_int))
      if (length(p_col) == 1) {
        result$P_Value_Interaction <- coef_table_int[interaction_idx, p_col]
        result$Model_Converged <- TRUE
      }
    }
    
    # --- MODEL 2: Marginal Time Model (Positive Control) ---
    formula_time <- "Value ~ Time"
    if (length(valid_covs) > 0) formula_time <- paste(formula_time, "+", paste(valid_covs, collapse = " + "))
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
  
  return(if (max_p_val == 0) NA else max_p_val)
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
    metric_name <- if (p_val < 0.05 && p_val > 0.0000) "FDR/P-Value" else "Interaction P-Value"
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