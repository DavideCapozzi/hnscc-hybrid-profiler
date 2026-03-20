# R/modules_coda.R
# ==============================================================================
# DATA TRANSFORMATION MODULE
# Description: Handling imputation (BPCA) and transformation (Logit & Log2).
# Dependencies: dplyr, pcaMethods
# ==============================================================================

library(dplyr)

# --- HELPER FUNCTIONS ---

#' @title Internal Median Fallback Helper
#' @description Fills NAs with column medians. Used when BPCA fails.
#' @param mat Numeric matrix.
#' @return Numeric matrix with NAs filled.
.impute_fallback_median <- function(mat) {
  apply(mat, 2, function(x) {
    if(all(is.na(x))) return(rep(0, length(x))) 
    x[is.na(x)] <- median(x, na.rm = TRUE)
    return(x)
  })
}

# --- MAIN FUNCTIONS ---

#' @title Bayesian PCA Imputation (Robust & Fail-Safe)
#' @description Imputes missing values using Bayesian PCA.
#' @param mat Numeric matrix (samples x markers).
#' @param nPcs Number of components. "auto" uses kEstimate (CV).
#' @param seed Integer. Seed for reproducibility.
#' @return A numeric matrix with NAs filled.
#' @export
impute_matrix_bpca <- function(mat, nPcs = "auto", seed = 123) {
  
  if (!any(is.na(mat)) && !any(is.infinite(mat))) return(mat)
  
  if (nrow(mat) < 5 || ncol(mat) < 2) {
    warning("[Impute] Matrix too small (<5 samples or <2 cols). Falling back to Median.")
    return(.impute_fallback_median(mat)) 
  }
  
  mat[is.infinite(mat)] <- NA
  
  if (!requireNamespace("pcaMethods", quietly = TRUE)) stop("Package 'pcaMethods' required.")
  
  vars <- apply(mat, 2, var, na.rm = TRUE)
  is_const <- !is.na(vars) & (vars < 1e-12)
  const_cols <- names(which(is_const))
  var_cols   <- names(which(!is_const))
  
  mat_bpca_in <- mat[, var_cols, drop = FALSE]
  
  if (ncol(mat_bpca_in) < 2) {
    warning("[Impute] Too few variable columns after filtering. Using Median fallback.")
    return(.impute_fallback_median(mat))
  }
  
  hard_limit_k <- min(nrow(mat_bpca_in), ncol(mat_bpca_in)) - 1
  if (hard_limit_k < 1) hard_limit_k <- 1
  
  k_used <- 2 
  
  if (is.null(nPcs) || nPcs == "auto") {
    k_used <- tryCatch({
      set.seed(seed)
      max_k <- min(5, hard_limit_k)
      
      utils::capture.output({
        k_est <- pcaMethods::kEstimate(mat_bpca_in, method = "bpca", evalPcs = 1:max_k, verbose = FALSE, scale = "none") 
      })
      suggested_k <- k_est$bestNPcs
      
      if (length(suggested_k) > 1) {
        warning(sprintf("   [Impute] kEstimate ambiguous. Using conservative K=%d.", min(suggested_k)))
        min(suggested_k)
      } else if (length(suggested_k) == 0 || is.na(suggested_k)) {
        2
      } else if (suggested_k > hard_limit_k) {
        hard_limit_k
      } else {
        suggested_k
      }
    }, error = function(e) {
      warning("   [Impute] kEstimate failed. Falling back to K=2.")
      2 
    })
  } else {
    k_used <- as.numeric(nPcs)
  }
  
  if (k_used > hard_limit_k) k_used <- hard_limit_k
  
  mat_bpca_out <- mat_bpca_in 
  
  tryCatch({
    set.seed(seed)
    res_bpca <- pcaMethods::pca(mat_bpca_in, method = "bpca", nPcs = k_used, verbose = FALSE)
    mat_bpca_out <- as.matrix(pcaMethods::completeObs(res_bpca))
    
  }, error = function(e) {
    warning(sprintf("   [Impute] BPCA execution failed. Applying Median Fallback."))
    mat_bpca_out <<- .impute_fallback_median(mat_bpca_in)
  })
  
  mat_final <- mat 
  mat_final[, var_cols] <- mat_bpca_out
  
  if (length(const_cols) > 0) {
    for (col in const_cols) {
      val <- mean(mat[, col], na.rm = TRUE)
      if (is.na(val)) val <- 0 
      mat_final[is.na(mat_final[, col]), col] <- val
    }
  }
  
  return(mat_final)
}

#' @title Robust Logit Transformation
#' @description 
#' Transforms values into Logit space (Log-Odds).
#' Handles scale normalization and mathematically sound epsilon bounding.
#' 
#' @param mat Numeric matrix.
#' @param epsilon Boundary protection value.
#' @param input_type String. "percentage" (divides by 100) or "proportion" (keeps as is).
#' @return Numeric matrix in logit space.
coda_transform_logit <- function(mat, epsilon = 1e-6, input_type = "percentage") {
  
  p <- mat
  
  # 1. Normalize Scale and Epsilon alignment
  if (input_type == "percentage") {
    if (any(mat > 100, na.rm = TRUE)) {
      stop("[Transform] Critical Error: Values > 100 detected in percentage mode.")
    }
    p <- mat / 100
    epsilon <- epsilon / 100 
  } else {
    if (any(mat > 1.0 + epsilon, na.rm = TRUE)) {
      stop("[Transform] Critical Error: Values > 1.0 detected in proportion mode.")
    }
  }
  
  # 2. Boundary Protection (Clamping)
  n_low  <- sum(p < epsilon, na.rm = TRUE)
  n_high <- sum(p > (1 - epsilon), na.rm = TRUE)
  
  if (n_low > 0 || n_high > 0) {
    message(sprintf("   [Transform-FACS] Clamped %d boundary values.", n_low + n_high))
  }
  
  p[p < epsilon] <- epsilon
  p[p > (1 - epsilon)] <- 1 - epsilon
  
  # 3. Transform
  logit_mat <- log(p / (1 - p))
  return(logit_mat)
}

#' @title Perform Hybrid Data Transformation Pipeline
#' @description 
#' Splits the matrix into FACS and Soluble blocks, applies independent 
#' transformations (Logit vs Log2), merges them, imputes via BPCA, 
#' and standardizes globally (Z-score).
#' @export
perform_data_transformation <- function(mat_raw, config, mode = "complete") {
  
  if (!mode %in% c("complete", "fast")) stop("Mode must be 'complete' or 'fast'.")
  
  # --- Setup Parameters ---
  facs_feats <- unlist(config$features$facs)
  sol_feats <- unlist(config$features$soluble)
  
  facs_method <- if(!is.null(config$transformation$facs_method)) config$transformation$facs_method else "logit"
  facs_fmt <- if(!is.null(config$transformation$facs_format)) config$transformation$facs_format else "percentage"
  sol_method <- if(!is.null(config$transformation$soluble_method)) config$transformation$soluble_method else "log2"
  
  # Dynamic Epsilon Calculation (Calculated strictly on FACS strictly positive data)
  facs_cols_present <- intersect(colnames(mat_raw), facs_feats)
  raw_eps_conf <- if(!is.null(config$transformation$epsilon)) config$transformation$epsilon else 1e-6
  
  if (is.character(raw_eps_conf) && tolower(raw_eps_conf) == "auto" && length(facs_cols_present) > 0) {
    facs_data <- mat_raw[, facs_cols_present, drop = FALSE]
    pos_vals <- facs_data[facs_data > 0 & !is.na(facs_data)]
    eps_val <- if (length(pos_vals) > 0) as.numeric(quantile(pos_vals, 0.01) / 2) else 1e-6
    message(sprintf("   [Transform] Epsilon set to 'auto'. Calculated data-driven boundary for FACS: %s", format(eps_val, scientific = TRUE)))
  } else {
    eps_val <- as.numeric(raw_eps_conf)
  }
  
  message(sprintf("\n[Transform] Starting Hybrid Strategy (%s mode)...", mode))
  rn_safe <- rownames(mat_raw)
  mat_transformed_parts <- list()
  
  # --- 1. FACS BLOCK PROCESSING ---
  if (length(facs_cols_present) > 0) {
    mat_facs <- mat_raw[, facs_cols_present, drop = FALSE]
    
    if (facs_method == "logit") {
      if (mode == "complete") {
        mins <- apply(mat_facs, 2, function(x) {
          pos <- x[x > 0 & !is.na(x)]
          if(length(pos) == 0) return(eps_val)
          return(min(pos) / 2)
        })
        for (j in 1:ncol(mat_facs)) {
          zeros <- which(mat_facs[, j] == 0 & !is.na(mat_facs[, j]))
          if (length(zeros) > 0) mat_facs[zeros, j] <- mins[j]
        }
      } else {
        mat_facs[mat_facs <= 0] <- eps_val 
      }
      
      mat_facs_trans <- coda_transform_logit(mat_facs, epsilon = eps_val, input_type = facs_fmt)
      mat_transformed_parts[["FACS"]] <- mat_facs_trans
      message(sprintf("   [Transform] Processed %d FACS markers (Method: %s).", ncol(mat_facs), facs_method))
    }
  }
  
  # --- 2. SOLUBLE BLOCK PROCESSING ---
  sol_cols_present <- intersect(colnames(mat_raw), sol_feats)
  if (length(sol_cols_present) > 0) {
    mat_sol <- mat_raw[, sol_cols_present, drop = FALSE]
    
    if (sol_method == "log2") {
      mat_sol_trans <- log2(mat_sol + 1)
      message(sprintf("   [Transform] Processed %d Soluble markers (Method: log2(x+1)).", ncol(mat_sol)))
    } else if (sol_method == "none") {
      mat_sol_trans <- mat_sol
      message(sprintf("   [Transform] Processed %d Soluble markers (Method: None/Raw).", ncol(mat_sol)))
    } else {
      stop(sprintf("[Transform] Unknown soluble transformation method: %s", sol_method))
    }
    mat_transformed_parts[["Soluble"]] <- mat_sol_trans
  }
  
  # --- 3. RE-MERGE ---
  if (length(mat_transformed_parts) == 0) stop("[FATAL] No valid features processed in transformation step.")
  mat_merged <- do.call(cbind, mat_transformed_parts)
  
  expected_cols <- intersect(colnames(mat_raw), colnames(mat_merged))
  mat_merged <- mat_merged[, expected_cols, drop = FALSE]
  
  # --- 4. JOINT IMPUTATION (BPCA) ---
  if (mode == "complete") {
    seed_val <- if(!is.null(config$stats$seed)) config$stats$seed else 123
    message("   [Transform] Running Joint BPCA Imputation on combined data...")
    mat_imputed <- impute_matrix_bpca(mat_merged, nPcs = "auto", seed = seed_val)
  } else {
    mat_imputed <- mat_merged
  }
  
  rownames(mat_imputed) <- rn_safe
  
  # --- 5. GLOBAL Z-SCORE ---
  mat_final_raw <- mat_imputed[, sort(colnames(mat_imputed)), drop = FALSE]
  mat_final_z <- NULL
  
  if (mode == "complete") {
    col_vars_post <- apply(mat_final_raw, 2, var, na.rm = TRUE)
    valid_cols_post <- !is.na(col_vars_post) & (col_vars_post > 1e-12)
    
    if (sum(!valid_cols_post) > 0) {
      n_dropped <- sum(!valid_cols_post)
      message(sprintf("   [Transform] Dropping %d columns with zero variance post-imputation.", n_dropped))
      mat_final_raw <- mat_final_raw[, valid_cols_post, drop = FALSE]
    }
    
    if (ncol(mat_final_raw) == 0) stop("[FATAL] All columns dropped due to zero variance.")
    
    mat_final_z <- scale(mat_final_raw)
    attr(mat_final_z, "scaled:center") <- NULL
    attr(mat_final_z, "scaled:scale") <- NULL
    mat_final_z <- as.matrix(mat_final_z)
    rownames(mat_final_z) <- rn_safe
  } else {
    mat_final_z <- mat_final_raw
  }
  
  return(list(
    hybrid_data_raw = mat_final_raw,
    hybrid_data_z = mat_final_z,
    hybrid_markers = colnames(mat_final_z)
  ))
}