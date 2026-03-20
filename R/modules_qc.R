# R/modules_qc.R
# ==============================================================================
# QUALITY CONTROL MODULE
# Description: Functions for data filtering (Variance, Missingness, Outliers).
# ==============================================================================

#' @title Detect Multivariate Outliers (PCA-based Mahalanobis)
#' @description 
#' Identifies outliers within specific groups using Robust Mahalanobis distance 
#' calculated on the first Principal Components (PCs).
#' 
#' @param mat Numeric matrix (Samples x Markers).
#' @param groups Vector of group labels corresponding to rows of mat.
#' @param conf_level Confidence level for Chi-squared cutoff (default 0.99).
#' @return A logical vector (TRUE = Outlier, FALSE = Keep).
detect_pca_outliers <- function(mat, groups, conf_level = 0.99) {
  
  is_outlier <- rep(FALSE, nrow(mat))
  
  # Ensure groups are factors
  groups <- as.factor(groups)
  
  for (g in levels(groups)) {
    # Get indices for current group
    idx <- which(groups == g)
    
    # Skip if too few samples (need at least 5 for PCA/Covariance stability)
    if (length(idx) < 5) {
      next
    }
    
    sub_mat <- mat[idx, , drop = FALSE]
    
    # 1. Handle NAs: First remove columns that are 100% NA in this subgroup
    na_counts <- colSums(is.na(sub_mat))
    valid_cols <- na_counts < nrow(sub_mat)
    sub_mat <- sub_mat[, valid_cols, drop = FALSE]
    
    if (ncol(sub_mat) < 2) next
    
    # 2. Impute remaining sporadic NAs (Median Imputation)
    if (any(is.na(sub_mat))) {
      sub_mat <- apply(sub_mat, 2, function(x) {
        if (all(is.na(x))) return(x) 
        x[is.na(x)] <- median(x, na.rm = TRUE)
        return(x)
      })
    }
    
    # 3. Remove zero-variance columns locally (Robust check)
    vars <- apply(sub_mat, 2, var, na.rm = TRUE)
    keep_vars <- !is.na(vars) & vars > 1e-12
    
    sub_mat <- sub_mat[, keep_vars, drop = FALSE]
    
    if (ncol(sub_mat) < 2) next 
    
    # Run PCA 
    tryCatch({
      pca_res <- prcomp(sub_mat, scale. = TRUE, center = TRUE)
      
      max_available_pcs <- ncol(pca_res$x)
      n_pc <- min(ncol(sub_mat), length(idx) - 2, 5) 
      
      if (n_pc > max_available_pcs) n_pc <- max_available_pcs
      if (n_pc < 1) n_pc <- 1
      
      scores <- pca_res$x[, 1:n_pc, drop = FALSE]
      
      # Calculate Robust Mahalanobis Distance
      center_vec <- colMeans(scores)
      cov_mat <- cov(scores)
      
      d2 <- mahalanobis(scores, center_vec, cov_mat)
      
      # Cutoff based on Chi-Squared distribution
      cutoff <- qchisq(conf_level, df = n_pc)
      
      # Flag local outliers
      local_outliers <- d2 > cutoff
      is_outlier[idx[local_outliers]] <- TRUE
      
    }, error = function(e) {
      warning(sprintf("[QC] Outlier detection failed for group '%s': %s", g, e$message))
    })
  }
  
  return(is_outlier)
}

#' @title Run Quality Control Pipeline (Hybrid Version)
#' @description 
#' Filters the data matrix based on variance, missingness, and multivariate outliers.
#' Applies dual missingness thresholds for FACS and Soluble features.
#' 
#' @param mat_raw Raw numeric matrix (Samples x Markers).
#' @param metadata Dataframe containing Group information.
#' @param qc_config List containing thresholds.
#' @param facs_features Character vector of FACS marker names.
#' @param soluble_features Character vector of Soluble marker names.
#' @param stratification_col String. Column name for detailed group breakdown.
#' @param dropped_markers_apriori Dataframe of markers dropped by config.
#' @param dropped_samples_apriori Dataframe of samples dropped by config.
#' @return A list containing the filtered 'data' and the 'report'.
run_qc_pipeline <- function(mat_raw, metadata, qc_config, 
                            facs_features = NULL, soluble_features = NULL,
                            stratification_col = "Group",
                            dropped_markers_apriori = NULL,
                            dropped_samples_apriori = NULL) {
  
  message("[QC] Running Hybrid Quality Control...")
  
  get_group_counts <- function(meta, col_name) {
    if (is.null(col_name) || !col_name %in% colnames(meta)) return(NULL)
    return(table(meta[[col_name]]))
  }
  
  init_counts <- get_group_counts(metadata, stratification_col)
  
  qc_summary <- list(
    n_row_init = nrow(mat_raw),
    n_col_init = ncol(mat_raw),
    breakdown_init = init_counts, 
    dropped_markers_apriori = dropped_markers_apriori,
    dropped_samples_apriori = dropped_samples_apriori, 
    n_col_zerovar = 0,
    dropped_rows_detail = data.frame(Patient_ID = character(), NA_Percent = numeric(), 
                                     Reason = character(), Original_Source = character(),
                                     stringsAsFactors = FALSE),
    dropped_cols_detail = data.frame(Marker = character(), NA_Percent = numeric(), Category = character()),
    n_row_dropped = 0,
    n_col_dropped = 0
  )
  
  curr_mat <- mat_raw
  curr_meta <- metadata
  
  # 1. Remove Constant Columns (Zero Variance)
  col_vars <- apply(curr_mat, 2, var, na.rm = TRUE)
  const_cols <- which(col_vars == 0 | is.na(col_vars))
  
  if (length(const_cols) > 0) {
    message(sprintf("   [QC] Dropping %d markers with zero variance.", length(const_cols)))
    qc_summary$n_col_zerovar <- length(const_cols)
    curr_mat <- curr_mat[, -const_cols, drop = FALSE]
  }
  
  # 2. Filter by Missingness (Rows - Patients)
  row_na_freq <- rowMeans(is.na(curr_mat))
  drop_row_na <- which(row_na_freq > qc_config$max_na_row_pct)
  
  if (length(drop_row_na) > 0) {
    pids <- rownames(curr_mat)[drop_row_na]
    message(sprintf("   [QC] Dropping %d patients with >%.0f%% total missingness.", 
                    length(pids), qc_config$max_na_row_pct * 100))
    
    group_info <- if(!is.null(stratification_col) && stratification_col %in% colnames(curr_meta)) {
      as.character(curr_meta[[stratification_col]][drop_row_na])
    } else { rep("N/A", length(pids)) }
    
    qc_summary$dropped_rows_detail <- rbind(qc_summary$dropped_rows_detail, data.frame(
      Patient_ID = pids,
      NA_Percent = round(row_na_freq[drop_row_na] * 100, 2),
      Reason = "High Missingness",
      Original_Source = group_info,
      stringsAsFactors = FALSE
    ))
    
    curr_mat <- curr_mat[-drop_row_na, , drop = FALSE]
    curr_meta <- curr_meta[-drop_row_na, , drop = FALSE]
  }
  
  # 3. Filter by Missingness (Cols - Markers) with DUAL THRESHOLDS
  col_na_freq <- colMeans(is.na(curr_mat))
  drop_col_idx <- c()
  drop_col_details <- list()
  
  current_cols <- colnames(curr_mat)
  
  # A. FACS Columns
  if (!is.null(facs_features) && length(facs_features) > 0) {
    facs_in_mat <- intersect(current_cols, facs_features)
    if (length(facs_in_mat) > 0) {
      facs_na <- col_na_freq[facs_in_mat]
      drop_facs <- names(which(facs_na > qc_config$max_na_col_facs_pct))
      if (length(drop_facs) > 0) {
        drop_col_idx <- c(drop_col_idx, drop_facs)
        drop_col_details[[1]] <- data.frame(
          Marker = drop_facs, 
          NA_Percent = round(facs_na[drop_facs] * 100, 2),
          Category = "FACS"
        )
      }
    }
  }
  
  # B. Soluble Columns
  if (!is.null(soluble_features) && length(soluble_features) > 0) {
    sol_in_mat <- intersect(current_cols, soluble_features)
    if (length(sol_in_mat) > 0) {
      sol_na <- col_na_freq[sol_in_mat]
      drop_sol <- names(which(sol_na > qc_config$max_na_col_sol_pct))
      if (length(drop_sol) > 0) {
        drop_col_idx <- c(drop_col_idx, drop_sol)
        drop_col_details[[2]] <- data.frame(
          Marker = drop_sol, 
          NA_Percent = round(sol_na[drop_sol] * 100, 2),
          Category = "Soluble"
        )
      }
    }
  }
  
  if (length(drop_col_idx) > 0) {
    message(sprintf("   [QC] Dropping %d markers due to domain-specific missingness thresholds.", length(drop_col_idx)))
    qc_summary$dropped_cols_detail <- do.call(rbind, drop_col_details)
    qc_summary$n_col_dropped <- length(drop_col_idx)
    
    keep_cols <- setdiff(colnames(curr_mat), drop_col_idx)
    curr_mat <- curr_mat[, keep_cols, drop = FALSE]
  }
  
  # 4. Filter Multivariate Outliers (Per Group)
  if (!is.null(qc_config$remove_outliers) && qc_config$remove_outliers) {
    message("   [QC] Checking for multivariate outliers (PCA-based)...")
    
    if ("Group" %in% names(curr_meta)) {
      is_outlier <- detect_pca_outliers(curr_mat, curr_meta$Group, conf_level = qc_config$outlier_conf_level)
      
      if (any(is_outlier)) {
        out_pids <- rownames(curr_mat)[is_outlier]
        message(sprintf("   [QC] Dropping %d outliers detected based on group distribution.", length(out_pids)))
        
        target_info_col <- if(!is.null(stratification_col) && stratification_col %in% colnames(curr_meta)) stratification_col else "Group"
        group_info <- as.character(curr_meta[[target_info_col]][is_outlier])
        
        qc_summary$dropped_rows_detail <- rbind(qc_summary$dropped_rows_detail, data.frame(
          Patient_ID = out_pids,
          NA_Percent = round(rowMeans(is.na(curr_mat[is_outlier,,drop=FALSE])) * 100, 2),
          Reason = "Outlier",
          Original_Source = group_info,
          stringsAsFactors = FALSE
        ))
        
        curr_mat <- curr_mat[!is_outlier, , drop = FALSE]
        curr_meta <- curr_meta[!is_outlier, , drop = FALSE]
      }
    } else {
      warning("[QC] 'Group' column not found in metadata. Skipping outlier detection.")
    }
  }
  
  group_map <- unique(data.frame(
    Subgroup = if(!is.null(stratification_col) && stratification_col %in% colnames(curr_meta)) as.character(curr_meta[[stratification_col]]) else "All",
    Group = if("Group" %in% colnames(curr_meta)) as.character(curr_meta$Group) else "All",
    stringsAsFactors = FALSE
  ))
  
  qc_summary$n_row_dropped <- nrow(qc_summary$dropped_rows_detail)
  qc_summary$n_row_final <- nrow(curr_mat)
  qc_summary$n_col_final <- ncol(curr_mat)
  qc_summary$final_markers_names <- colnames(curr_mat)
  qc_summary$breakdown_final <- get_group_counts(curr_meta, stratification_col)
  qc_summary$group_mapping <- group_map
  
  message(sprintf("   [QC] Final Dimensions: %d Samples x %d Markers", nrow(curr_mat), ncol(curr_mat)))
  
  return(list(
    data = curr_mat,
    metadata = curr_meta, 
    report = qc_summary
  ))
}