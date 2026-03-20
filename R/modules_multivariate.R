# R/modules_multivariate.R
# ==============================================================================
# MULTIVARIATE ANALYSIS MODULE (sPLS-DA)
# Description: Wrappers for Sparse PLS-DA with automatic feature tuning.
# Dependencies: mixOmics, dplyr
# ==============================================================================

library(dplyr)
library(mixOmics)

#' @title Run Sparse PLS-DA (sPLS-DA) with Tuning
#' @description 
#' Fits a supervised sparse PLS-DA model.
#' Performs automatic tuning to select the optimal number of variables (keepX)
#' for each component using cross-validation.
#' 
#' @param data_z Matrix of Z-scored data (Samples x Features).
#' @param metadata Dataframe containing metadata (must match data_z rows).
#' @param group_col Name of the column in metadata to use as the class/Y.
#' @param n_comp Number of components to compute (default: 2).
#' @param validation_method CV method ("Mfold" or "loo"). 
#' @param folds Number of folds for Mfold CV (default: 5).
#' @param n_repeat Number of repetitions for Cross-Validation (default: 50).
#' @return A list containing the final optimized model, performance metrics, and tuning results.
run_splsda_model <- function(data_z, metadata, group_col = "Group", n_comp = 2, 
                             validation_method = "Mfold", folds = 5, n_repeat = 50) {
  
  # Ensure mixOmics is loaded
  requireNamespace("mixOmics", quietly = TRUE)
  
  Y <- as.factor(metadata[[group_col]])
  X <- as.matrix(data_z)
  
  # Validation: Dimensions
  if (nrow(X) != length(Y)) stop("Dimension mismatch between Data and Group labels.")
  
  # Validation: Minimum samples for Mfold
  min_samples_per_class <- min(table(Y))
  if (validation_method == "Mfold" && min_samples_per_class < folds) {
    message(sprintf("   [sPLS-DA] Warning: Smallest class has %d samples. Adjusting folds to %d.", 
                    min_samples_per_class, min_samples_per_class))
    folds <- min_samples_per_class
  }
  
  # Enforce absolute minimum folds required by mixOmics for Mfold CV
  if (validation_method == "Mfold" && folds < 3) {
    stop("Minimum 3 samples per class required for M-fold cross-validation. Analysis aborted.")
  }
  
  # 1. Tuning Step: Determine optimal keepX
  message(sprintf("   [sPLS-DA] Tuning optimal features (keepX) for %d components...", n_comp))
  
  # Define grid of keepX to test: from 5 markers up to all markers
  n_vars <- ncol(X)
  
  if (n_vars < 5) {
    list_keepX <- c(n_vars)
  } else {
    # Sequence of 5s
    list_keepX <- seq(5, min(30, n_vars), by = 5)
    
    # Always include the full set of variables if small enough, or at least the max cap
    max_val <- min(30, n_vars)
    if (!max_val %in% list_keepX) {
      list_keepX <- c(list_keepX, max_val)
    }
    
    # Strictly include n_vars if it's <= 30 and not in list
    if (n_vars <= 30 && !n_vars %in% list_keepX) {
      list_keepX <- c(list_keepX, n_vars)
    }
  }
  list_keepX <- sort(unique(list_keepX))
  
  # Run Tuning (tune.splsda)
  tune_splsda <- mixOmics::tune.splsda(
    X = X, 
    Y = Y, 
    ncomp = n_comp, 
    test.keepX = list_keepX, 
    validation = validation_method, 
    folds = folds, 
    dist = "max.dist", 
    progressBar = FALSE,
    nrepeat = n_repeat  
  )
  
  choice_keepX <- tune_splsda$choice.keepX
  message(sprintf("      -> Optimal keepX selected: Comp1=%d, Comp2=%d", 
                  choice_keepX[1], choice_keepX[2]))
  
  # 2. Final Model Fitting
  message("   [sPLS-DA] Fitting final sparse model with selected parameters...")
  final_model <- mixOmics::splsda(X, Y, ncomp = n_comp, keepX = choice_keepX)
  
  # 3. Final Performance Evaluation (Error Rate)
  # We wrap this in tryCatch as perf() can sometimes fail on edge cases
  perf_splsda <- tryCatch({
    mixOmics::perf(final_model, validation = validation_method, folds = folds, 
                   progressBar = FALSE, nrepeat = n_repeat, auc = TRUE) 
  }, error = function(e) {
    message(paste("      [WARN] Performance evaluation failed:", e$message))
    return(NULL) 
  })
  
  return(list(model = final_model, performance = perf_splsda, tuning = tune_splsda))
}

#' @title Extract PLS-DA Performance Metrics
#' @description Returns BER (Balanced Error Rate) from the performance object.
extract_plsda_performance <- function(pls_result) {
  perf <- pls_result$performance
  n_comp <- pls_result$model$ncomp
  
  # Handle failure case
  if (is.null(perf)) {
    return(data.frame(
      Component = 1:n_comp,
      Overall_BER = NA,
      Validation_Method = "FAILED"
    ))
  }
  
  # Extract BER safely
  ber_vals <- rep(NA, n_comp)
  
  tryCatch({
    if (!is.null(perf$error.rate$BER)) {
      # Check if matrix or vector
      if (is.matrix(perf$error.rate$BER)) {
        # Usually extract 'max.dist' distance metric for sPLS-DA
        raw_ber <- perf$error.rate$BER[, "max.dist"]
        if(length(raw_ber) >= n_comp) ber_vals <- raw_ber[1:n_comp]
      } else if (is.vector(perf$error.rate$BER)) {
        ber_vals[1] <- perf$error.rate$BER["max.dist"]
      }
    }
  }, error = function(e) {
    message("      [WARN] BER extraction encountered an issue. Using NAs.")
  })
  
  # Extract validation method string
  val_method <- "Unknown"
  if (!is.null(perf$validation)) {
    val_method <- as.character(perf$validation)[1]
  }
  
  df_perf <- data.frame(
    Component = 1:n_comp,
    Overall_BER = ber_vals,
    Validation_Method = val_method 
  )
  
  return(df_perf)
}

#' @title Extract Discriminant Features (Sparse Loadings)
#' @description 
#' Extracts loading weights from the sPLS-DA model. 
#' Filters out variables with zero weight (not selected by the model).
#' 
#' @param pls_result The list returned by run_splsda_model.
#' @return A dataframe of non-zero loadings.
extract_plsda_loadings <- function(pls_result) {
  model <- pls_result$model
  
  # Check if loadings exist
  if(is.null(model$loadings$X)) return(data.frame())
  
  # Convert loadings matrix to dataframe
  loadings_raw <- as.data.frame(model$loadings$X)
  loadings_raw$Marker <- rownames(loadings_raw)
  
  # Reshape to long format or keep wide? 
  # Strategy: Create a clean table with Comp1 and Comp2 weights
  
  df_clean <- loadings_raw %>%
    dplyr::select(Marker, everything()) %>%
    dplyr::rename(Comp1_Weight = comp1)
  
  if(model$ncomp > 1) {
    df_clean <- df_clean %>% dplyr::rename(Comp2_Weight = comp2)
  }
  
  # Calculate Importance (Absolute weight on Comp1)
  df_clean$Importance <- abs(df_clean$Comp1_Weight)
  
  # Filter: In sPLS-DA, we only care about non-zero weights
  # We keep rows where at least one component has a non-zero weight
  df_clean <- df_clean %>%
    filter(Importance > 0 | (if("Comp2_Weight" %in% names(.)) abs(Comp2_Weight) > 0 else FALSE)) %>%
    arrange(desc(Importance)) %>%
    mutate(
      Direction = ifelse(Comp1_Weight > 0, "Positive_Assoc", "Negative_Assoc")
    )
  
  return(df_clean)
}