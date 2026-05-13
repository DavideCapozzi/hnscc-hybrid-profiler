# R/modules_ml.R
# ==============================================================================
# MACHINE LEARNING MODULE (Hybrid Two-Step Classifier)
# Description: Nested-LOOCV classification on LOO-robust biomarkers from LMM.
#              Primary method: Elastic Net Logistic Regression (glmnet).
#              Secondary method: SVM with RBF Kernel (e1071).
# Dependencies: dplyr, ggplot2, tidyr, glmnet, e1071, pROC
# ==============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)


#' @title Load LMM LOO-Robust Features
#' @description
#' Reads the Step 04 JSON and extracts markers that survive both FDR correction
#' and Leave-One-Out sensitivity testing. Returns an empty list when no markers
#' qualify, enabling graceful step-level skip without pipeline interruption.
#'
#' @param lmm_json_path String. Path to Machine_Metrics_LMM_<name>.json from Step 04.
#' @param fdr_threshold Numeric. Maximum FDR_Interaction to pass the gate (default 0.05).
#' @param loo_threshold Numeric. Maximum Max_P_Value_LOO to pass the gate (default 0.05).
#' @return A named list: markers (character), n_robust (integer), df_robust (data.frame).
#' @export
load_lmm_robust_features <- function(lmm_json_path,
                                     fdr_threshold = 0.05,
                                     loo_threshold = 0.05) {
  empty <- list(markers = character(0), n_robust = 0L, df_robust = data.frame())

  if (!file.exists(lmm_json_path)) {
    message(sprintf("   [ML] LMM JSON not found: %s. Skipping ML analysis.", lmm_json_path))
    return(empty)
  }

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    message("   [ML] Package 'jsonlite' required to read LMM JSON. Skipping ML analysis.")
    return(empty)
  }

  metrics <- jsonlite::read_json(lmm_json_path, simplifyVector = TRUE)
  df      <- as.data.frame(metrics$full_results)

  if (nrow(df) == 0 || !"FDR_Interaction" %in% names(df)) {
    message("   [ML] LMM JSON contains no results. Skipping ML analysis.")
    return(empty)
  }

  # jsonlite drops all-NULL columns when simplifyVector = TRUE; restore as NA
  if (!"Max_P_Value_LOO" %in% names(df)) df$Max_P_Value_LOO <- NA_real_

  df_robust <- df %>%
    dplyr::filter(
      !is.na(FDR_Interaction),
      !is.na(Max_P_Value_LOO),
      FDR_Interaction < fdr_threshold,
      Max_P_Value_LOO < loo_threshold
    ) %>%
    dplyr::arrange(FDR_Interaction)

  if (nrow(df_robust) == 0) {
    message(sprintf(
      "   [ML] No markers survive LOO gate (FDR < %.2f AND LOO p < %.2f). Skipping ML analysis.",
      fdr_threshold, loo_threshold
    ))
    return(empty)
  }

  message(sprintf("   [ML] %d LOO-robust marker(s) selected for ML feature gate:", nrow(df_robust)))
  for (i in seq_len(nrow(df_robust))) {
    message(sprintf("         -> %s  (FDR=%.4f | LOO_max_p=%.4f)",
                    df_robust$Marker[i],
                    df_robust$FDR_Interaction[i],
                    df_robust$Max_P_Value_LOO[i]))
  }

  list(markers = df_robust$Marker, n_robust = nrow(df_robust), df_robust = df_robust)
}


#' @title Filter Collinear Features
#' @description
#' Iteratively removes near-collinear markers from a feature matrix. At each
#' step, the pair with the highest absolute Pearson correlation is identified;
#' the marker with the higher MEAN absolute correlation to all other remaining
#' markers is dropped (the more globally redundant one). Repeats until no pair
#' exceeds the threshold. This correctly handles chains of collinearity without
#' relying on univariate AUC, which can be noisy in small samples.
#'
#' @param X_main Numeric matrix. Main-effect feature matrix (samples x features).
#' @param cor_threshold Numeric. Drop one of a pair when |r| exceeds this value (default 0.85).
#' @return A named list: X (filtered matrix), kept (character), dropped (character).
#' @export
filter_collinear_features <- function(X_main, cor_threshold = 0.85) {
  markers <- colnames(X_main)
  dropped <- character(0)

  repeat {
    if (length(markers) <= 1) break

    cm <- abs(cor(X_main[, markers, drop = FALSE], method = "pearson"))
    diag(cm) <- 0
    max_r <- max(cm, na.rm = TRUE)
    if (!is.finite(max_r) || max_r < cor_threshold) break

    # Identify the most correlated pair
    idx    <- which(cm == max_r, arr.ind = TRUE)[1, ]
    m1     <- markers[idx[1]]
    m2     <- markers[idx[2]]
    others <- markers[!markers %in% c(m1, m2)]

    if (length(others) == 0) {
      # Only two markers left and they're collinear — keep the one with lower
      # correlation to itself (identical: keep m1 by convention, drop m2)
      dropped <- c(dropped, m2)
      markers <- m1
      break
    }

    # Drop whichever has higher mean |r| to all other markers (more redundant)
    mean_r_m1 <- mean(abs(cm[m1, others]), na.rm = TRUE)
    mean_r_m2 <- mean(abs(cm[m2, others]), na.rm = TRUE)
    to_drop   <- if (mean_r_m1 >= mean_r_m2) m1 else m2

    dropped <- c(dropped, to_drop)
    markers <- markers[markers != to_drop]
  }

  list(X = X_main[, markers, drop = FALSE], kept = markers, dropped = dropped)
}


#' @title Build ML Feature Matrix
#' @description
#' Extracts the Z-scored feature columns for the LOO-robust markers from the
#' Step 01 processed data object. Applies a collinearity filter before adding
#' interaction terms to prevent near-duplicate features from degrading the
#' classifier. Optionally appends pairwise interaction terms guided by the
#' network co-activation topology.
#'
#' @param DATA Named list. Processed data object from Step 01 (standard mode).
#' @param robust_markers Character vector. Markers selected by the LMM LOO gate.
#' @param include_interactions Logical. Append all pairwise products (default TRUE).
#' @param cor_threshold Numeric. Collinearity filter threshold (default 0.85).
#'   Set to 1 to disable filtering.
#' @return A named list: X (matrix), y (factor), feature_names (character), n_main (integer).
#' @export
build_ml_matrix <- function(DATA, robust_markers, include_interactions = TRUE,
                            cor_threshold = 0.85) {
  # Resolve marker availability against the transformed Z-score matrix
  all_z_cols  <- DATA$hybrid_markers
  available   <- intersect(robust_markers, all_z_cols)

  if (length(available) < length(robust_markers)) {
    missing <- setdiff(robust_markers, all_z_cols)
    warning(sprintf("[ML] %d marker(s) not found in Step 01 Z-score matrix: %s",
                    length(missing), paste(missing, collapse = ", ")))
  }
  if (length(available) == 0) stop("[ML][FATAL] No robust markers present in data matrix.")

  X <- as.matrix(DATA$hybrid_data_z[, available, drop = FALSE])
  mode(X) <- "numeric"
  rownames(X) <- make.unique(as.character(DATA$metadata$Patient_ID))

  # Collinearity filter: drop the more redundant of any near-duplicate pair
  if (length(available) >= 2 && is.finite(cor_threshold) && cor_threshold < 1) {
    filt <- filter_collinear_features(X, cor_threshold = cor_threshold)
    if (length(filt$dropped) > 0) {
      message(sprintf(
        "   [ML] Collinearity filter (|r| > %.2f): removed %s (kept %s)",
        cor_threshold,
        paste(filt$dropped, collapse = ", "),
        paste(filt$kept,    collapse = ", ")
      ))
    }
    X         <- filt$X
    available <- filt$kept
  }

  n_main <- length(available)

  if (include_interactions && n_main >= 2) {
    combos <- utils::combn(available, 2, simplify = FALSE)
    for (pair in combos) {
      col_name        <- paste(pair[1], pair[2], sep = "_x_")
      X               <- cbind(X, as.numeric(X[, pair[1]]) * as.numeric(X[, pair[2]]))
      colnames(X)[ncol(X)] <- col_name
    }
    message(sprintf("   [ML] Feature matrix: %d main effects + %d interaction term(s) = %d total",
                    n_main, length(combos), ncol(X)))
  } else {
    message(sprintf("   [ML] Feature matrix: %d main effect(s) (interactions disabled or p < 2)", n_main))
  }

  y <- DATA$metadata$Group  # Factor: levels = c(non_responder, responder)

  list(X = X, y = y, feature_names = colnames(X), n_main = n_main)
}


#' @title Run Nested-LOOCV with Elastic Net Logistic Regression
#' @description
#' Outer loop: Leave-One-Out (n folds). Inner loop: k-fold cross-validation via
#' glmnet::cv.glmnet to jointly select alpha (Ridge/Elastic Net/LASSO) and lambda.
#' No feature re-selection occurs inside the loop; the feature set is frozen prior
#' to LOOCV, ensuring no data leakage from the classification data.
#'
#' @param X Numeric matrix. Feature matrix (samples x features).
#' @param y Factor. Binary outcome (two levels; positive class = second level).
#' @param alpha_grid Numeric vector. Alpha values to try in inner CV (default c(0, 0.5, 1)).
#' @param n_lambda Integer. Number of lambda values on the regularization path (default 100).
#' @param inner_folds Integer. k for inner k-fold CV (default 5).
#' @param seed Integer. Random seed for reproducibility.
#' @return A named list: method, predicted_probs, y_true, positive_label, coef_matrix, feature_names.
#' @export
run_nested_loocv_glmnet <- function(X, y,
                                    alpha_grid  = c(0, 0.5, 1),
                                    n_lambda    = 100L,
                                    inner_folds = 5L,
                                    seed        = 2026) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("[ML] Package 'glmnet' required. Add r-glmnet to env/environment.yml and rebuild.")
  }

  n           <- nrow(X)
  pos_label   <- levels(y)[2]
  y_bin       <- as.integer(y) - 1L  # 0 = reference, 1 = positive class

  pred_probs  <- numeric(n)
  coef_list   <- vector("list", n)
  coef_names  <- c("(Intercept)", colnames(X))

  message(sprintf("   [Elastic Net] Outer LOO: %d folds | Inner CV: %d folds | alpha grid: {%s}",
                  n, inner_folds, paste(alpha_grid, collapse = ", ")))

  for (i in seq_len(n)) {
    X_train <- X[-i, , drop = FALSE]
    y_train <- y_bin[-i]
    X_test  <- X[i,  , drop = FALSE]

    n_min_class  <- min(sum(y_train == 0L), sum(y_train == 1L))
    safe_folds   <- max(2L, min(as.integer(inner_folds), n_min_class))

    # Inner CV: select best (alpha, lambda) by maximum inner AUC
    best_alpha      <- alpha_grid[1]
    best_lambda     <- NULL
    best_inner_auc  <- -Inf

    for (a in alpha_grid) {
      set.seed(seed + i)
      cv_fit <- tryCatch(
        glmnet::cv.glmnet(
          x            = X_train,
          y            = as.numeric(y_train),
          family       = "binomial",
          alpha        = a,
          nfolds       = safe_folds,
          type.measure = "auc",
          nlambda      = as.integer(n_lambda)
        ),
        error = function(e) NULL
      )

      if (!is.null(cv_fit)) {
        inner_auc <- suppressWarnings(max(cv_fit$cvm, na.rm = TRUE))
        if (is.finite(inner_auc) && inner_auc > best_inner_auc) {
          best_inner_auc <- inner_auc
          best_alpha     <- a
          best_lambda    <- cv_fit$lambda.min
        }
      }
    }

    # Fit final model on full training fold with best hyperparameters
    set.seed(seed + i)
    final_fit <- tryCatch(
      glmnet::glmnet(X_train, as.numeric(y_train), family = "binomial",
                     alpha = best_alpha, lambda = best_lambda),
      error = function(e) NULL
    )

    if (!is.null(final_fit) && !is.null(best_lambda)) {
      pred_probs[i] <- tryCatch(
        as.numeric(predict(final_fit, newx = X_test, s = best_lambda, type = "response")),
        error = function(e) 0.5
      )
      coef_list[[i]] <- tryCatch(
        as.numeric(coef(final_fit, s = best_lambda)),
        error = function(e) rep(NA_real_, length(coef_names))
      )
    } else {
      pred_probs[i]  <- 0.5
      coef_list[[i]] <- rep(NA_real_, length(coef_names))
    }
  }

  coef_matrix           <- do.call(rbind, coef_list)
  colnames(coef_matrix) <- coef_names

  list(
    method        = "Elastic Net Logistic Regression (Nested-LOOCV)",
    predicted_probs = pred_probs,
    y_true        = y,
    positive_label = pos_label,
    coef_matrix   = coef_matrix,
    feature_names = colnames(X)
  )
}


#' @title Run Nested-LOOCV with SVM (RBF Kernel)
#' @description
#' Outer loop: Leave-One-Out (n folds). Inner loop: grid search over (C, gamma)
#' via e1071::tune() with k-fold CV. Probability estimates obtained via Platt
#' scaling (probability = TRUE in e1071::svm). Serves as a non-linear comparison
#' arm to validate whether a linear decision boundary is sufficient.
#'
#' @param X Numeric matrix. Feature matrix (samples x features).
#' @param y Factor. Binary outcome (two levels; positive class = second level).
#' @param C_grid Numeric vector. Cost values for inner grid search.
#' @param gamma_grid Numeric vector. Gamma values for inner grid search.
#' @param inner_folds Integer. k for inner k-fold CV (default 5).
#' @param seed Integer. Random seed for reproducibility.
#' @return A named list: method, predicted_probs, y_true, positive_label, coef_matrix (NULL), feature_names.
#' @export
run_nested_loocv_svm <- function(X, y,
                                 C_grid      = c(0.01, 0.1, 1, 10, 100),
                                 gamma_grid  = c(0.01, 0.1, 1, 10),
                                 inner_folds = 5L,
                                 seed        = 2026) {
  if (!requireNamespace("e1071", quietly = TRUE)) {
    stop("[ML] Package 'e1071' required. Add r-e1071 to env/environment.yml and rebuild.")
  }

  n          <- nrow(X)
  pos_label  <- levels(y)[2]
  pred_probs <- numeric(n)

  message(sprintf("   [SVM-RBF] Outer LOO: %d folds | Inner CV: %d folds | Grid: C(%d) x gamma(%d)",
                  n, inner_folds, length(C_grid), length(gamma_grid)))

  for (i in seq_len(n)) {
    X_train <- X[-i, , drop = FALSE]
    y_train <- y[-i]
    X_test  <- X[i,  , drop = FALSE]

    n_min_class <- min(table(y_train))
    safe_folds  <- max(2L, min(as.integer(inner_folds), n_min_class))

    set.seed(seed + i)
    tune_res <- tryCatch(
      suppressWarnings(e1071::tune(
        e1071::svm,
        train.x    = X_train,
        train.y    = y_train,
        kernel     = "radial",
        ranges     = list(cost = C_grid, gamma = gamma_grid),
        tunecontrol = e1071::tune.control(sampling = "cross", cross = safe_folds)
      )),
      error = function(e) NULL
    )

    if (!is.null(tune_res)) {
      bp <- tune_res$best.parameters

      set.seed(seed + i)
      final_svm <- tryCatch(
        e1071::svm(x = X_train, y = y_train, kernel = "radial",
                   cost = bp$cost, gamma = bp$gamma, probability = TRUE),
        error = function(e) NULL
      )

      if (!is.null(final_svm)) {
        pred_obj <- tryCatch(
          predict(final_svm, newdata = X_test, probability = TRUE),
          error = function(e) NULL
        )

        if (!is.null(pred_obj)) {
          prob_mat       <- attr(pred_obj, "probabilities")
          pred_probs[i]  <- if (!is.null(prob_mat) && pos_label %in% colnames(prob_mat)) {
            as.numeric(prob_mat[1, pos_label])
          } else 0.5
        } else {
          pred_probs[i] <- 0.5
        }
      } else {
        pred_probs[i] <- 0.5
      }
    } else {
      pred_probs[i] <- 0.5
    }
  }

  list(
    method         = "SVM with RBF Kernel (Nested-LOOCV)",
    predicted_probs = pred_probs,
    y_true         = y,
    positive_label = pos_label,
    coef_matrix    = NULL,
    feature_names  = colnames(X)
  )
}


#' @title Compute Classification Performance Metrics
#' @description
#' Computes AUC with 95% CI (DeLong method), balanced accuracy, sensitivity,
#' and specificity at the Youden-J optimal threshold from out-of-fold predictions.
#'
#' @param y_true Factor. True class labels.
#' @param y_pred_prob Numeric vector. Predicted probabilities for the positive class.
#' @param positive_label String. The positive class label (second factor level).
#' @return A named list: auc, auc_ci (length-3 numeric: lower, estimate, upper),
#'   balanced_accuracy, ber, sensitivity, specificity, threshold.
#' @export
compute_classification_metrics <- function(y_true, y_pred_prob, positive_label) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("[ML] Package 'pROC' required. Add r-proc to env/environment.yml and rebuild.")
  }

  y_bin   <- as.integer(y_true == positive_label)
  na_mask <- !is.na(y_pred_prob)

  if (sum(na_mask) < length(y_bin)) {
    warning(sprintf("[ML] %d NA predictions replaced with 0.5 before metric computation.",
                    sum(!na_mask)))
    y_pred_prob[!na_mask] <- 0.5
  }

  roc_obj <- tryCatch(
    pROC::roc(response = y_bin, predictor = y_pred_prob, direction = "<", quiet = TRUE),
    error = function(e) {
      warning(paste("[ML] pROC::roc() failed:", e$message))
      return(NULL)
    }
  )

  if (is.null(roc_obj)) {
    return(list(auc = NA_real_, auc_ci = c(NA, NA, NA),
                balanced_accuracy = NA_real_, ber = NA_real_,
                sensitivity = NA_real_, specificity = NA_real_, threshold = 0.5))
  }

  auc_val <- as.numeric(pROC::auc(roc_obj))
  auc_ci  <- tryCatch(
    as.numeric(pROC::ci.auc(roc_obj, method = "delong")),
    error = function(e) c(NA_real_, auc_val, NA_real_)
  )

  # Youden-J optimal threshold
  coords_df <- tryCatch(
    pROC::coords(roc_obj, x = "best", best.method = "youden",
                 ret = c("threshold", "sensitivity", "specificity"),
                 transpose = FALSE),
    error = function(e) NULL
  )

  if (!is.null(coords_df) && nrow(coords_df) > 0) {
    opt_thresh <- as.numeric(coords_df$threshold[1])
    opt_sens   <- as.numeric(coords_df$sensitivity[1])
    opt_spec   <- as.numeric(coords_df$specificity[1])
  } else {
    opt_thresh <- 0.5
    opt_sens   <- NA_real_
    opt_spec   <- NA_real_
  }

  bal_acc <- if (!is.na(opt_sens) && !is.na(opt_spec)) (opt_sens + opt_spec) / 2 else NA_real_

  list(
    auc               = auc_val,
    auc_ci            = auc_ci,
    balanced_accuracy = bal_acc,
    ber               = if (!is.na(bal_acc)) 1 - bal_acc else NA_real_,
    sensitivity       = opt_sens,
    specificity       = opt_spec,
    threshold         = opt_thresh
  )
}


#' @title Plot Nested-LOOCV ROC Curves
#' @description
#' Overlays ROC curves for multiple methods on a single ggplot with AUC
#' annotated in the legend. Designed to compare Elastic Net and SVM-RBF outputs.
#'
#' @param results_list Named list of ML result objects (each with predicted_probs,
#'   y_true, positive_label, method).
#' @param colors_viz Named character vector of clinical group colours (optional).
#' @param title String. Plot title.
#' @return A ggplot object, or NULL on failure.
#' @export
plot_ml_roc <- function(results_list,
                        colors_viz = NULL,
                        title      = "Nested-LOOCV ROC Curves") {
  if (!requireNamespace("pROC", quietly = TRUE)) return(NULL)

  method_colors <- c(
    "Elastic Net" = "#2171B5",
    "SVM-RBF"     = "#CB181D"
  )

  roc_data <- purrr::map_dfr(names(results_list), function(nm) {
    res   <- results_list[[nm]]
    y_bin <- as.integer(res$y_true == res$positive_label)

    tryCatch({
      roc_obj <- pROC::roc(response = y_bin, predictor = res$predicted_probs,
                           direction = "<", quiet = TRUE)
      auc_val <- round(as.numeric(pROC::auc(roc_obj)), 3)

      data.frame(
        FPR    = 1 - roc_obj$specificities,
        TPR    = roc_obj$sensitivities,
        Method = sprintf("%s\n(AUC = %.3f)", nm, auc_val),
        Key    = nm,
        stringsAsFactors = FALSE
      )
    }, error = function(e) data.frame())
  })

  if (nrow(roc_data) == 0) return(NULL)

  p <- ggplot2::ggplot(roc_data, ggplot2::aes(x = FPR, y = TPR, color = Method)) +
    ggplot2::geom_line(linewidth = 1.2, na.rm = TRUE) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "grey55", linewidth = 0.6) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0),
                                 breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0),
                                 breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(title = title, x = "1 - Specificity (FPR)",
                  y = "Sensitivity (TPR)", color = NULL) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.direction = "vertical",
      plot.title       = ggplot2::element_text(size = 11, face = "bold")
    )

  p
}


#' @title Plot ML Predicted Probabilities by Clinical Group
#' @description
#' Jitter strip plot of out-of-fold predicted probabilities faceted by method.
#' A dashed reference line at 0.5 indicates the neutral decision boundary.
#'
#' @param predictions_df Data.frame with columns: Group, Prob_ElasticNet, Prob_SVM_RBF.
#' @param colors_viz Named character vector mapping group labels to hex colours.
#' @param title String. Plot title.
#' @return A ggplot object, or NULL on failure.
#' @export
plot_ml_predictions <- function(predictions_df,
                                colors_viz = NULL,
                                title      = "Predicted Probabilities by Clinical Group") {
  tryCatch({
    df_long <- predictions_df %>%
      dplyr::select(Group, dplyr::starts_with("Prob_")) %>%
      tidyr::pivot_longer(
        cols      = dplyr::starts_with("Prob_"),
        names_to  = "Method",
        values_to = "Probability"
      ) %>%
      dplyr::mutate(
        Method = dplyr::recode(Method,
          "Prob_ElasticNet" = "Elastic Net",
          "Prob_SVM_RBF"    = "SVM-RBF"
        ),
        Group = factor(Group)
      )

    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Group, y = Probability, color = Group)) +
      ggplot2::geom_jitter(width = 0.18, alpha = 0.80, size = 2.2, na.rm = TRUE) +
      ggplot2::stat_summary(fun = median, geom = "crossbar", width = 0.4,
                            color = "grey25", linewidth = 0.5, fatten = 1.5) +
      ggplot2::facet_wrap(~Method) +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed",
                          color = "grey45", linewidth = 0.6) +
      ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      ggplot2::labs(title = title, x = NULL,
                    y = sprintf("Predicted P(Positive Class)")) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(legend.position = "none",
                     plot.title = ggplot2::element_text(size = 11, face = "bold"))

    if (!is.null(colors_viz) && length(colors_viz) > 0) {
      avail <- intersect(levels(df_long$Group), names(colors_viz))
      if (length(avail) == length(levels(df_long$Group))) {
        p <- p + ggplot2::scale_color_manual(values = colors_viz)
      }
    }

    p
  }, error = function(e) {
    warning(paste("[ML] plot_ml_predictions failed:", e$message))
    return(NULL)
  })
}


#' @title Permutation AUC Test
#' @description
#' Empirical p-value for the observed AUC by permuting class labels against
#' fixed out-of-fold predicted probabilities. No ML models are re-fitted.
#' Uses the inclusive formula p = (B_extreme + 1) / (n_perm + 1).
#'
#' @param y_true Factor. True class labels (same order as predicted_probs).
#' @param predicted_probs Numeric vector. Out-of-fold predicted probabilities.
#' @param positive_label String. The positive class label.
#' @param n_perm Integer. Number of label permutations (default 2000).
#' @param seed Integer. Random seed for reproducibility.
#' @return A named list: observed_auc, p_value, n_perm, method_name.
#' @export
run_permutation_auc_test <- function(y_true, predicted_probs, positive_label,
                                     n_perm = 2000L, seed = 2026) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("[ML] Package 'pROC' required for permutation AUC test.")
  }

  y_bin    <- as.integer(y_true == positive_label)
  roc_obs  <- tryCatch(
    pROC::roc(response = y_bin, predictor = predicted_probs, direction = "<", quiet = TRUE),
    error = function(e) NULL
  )
  if (is.null(roc_obs)) return(list(observed_auc = NA_real_, p_value = NA_real_, n_perm = n_perm))

  obs_auc <- as.numeric(pROC::auc(roc_obs))
  n       <- length(y_bin)

  set.seed(seed)
  perm_aucs <- vapply(seq_len(n_perm), function(b) {
    y_perm <- sample(y_bin, size = n, replace = FALSE)
    roc_p  <- tryCatch(
      pROC::roc(response = y_perm, predictor = predicted_probs, direction = "<", quiet = TRUE),
      error = function(e) NULL
    )
    if (is.null(roc_p)) return(NA_real_)
    as.numeric(pROC::auc(roc_p))
  }, numeric(1))

  perm_aucs <- perm_aucs[!is.na(perm_aucs)]
  p_val     <- (sum(perm_aucs >= obs_auc) + 1L) / (length(perm_aucs) + 1L)

  list(observed_auc = obs_auc, p_value = round(p_val, 5), n_perm = n_perm)
}


#' @title Univariate AUC (Parameter-Free)
#' @description
#' Computes AUC with 95% CI (DeLong) for each marker independently using the
#' raw Z-score as predictor on all n samples. No hyperparameter tuning occurs,
#' so no LOOCV is needed — the estimate is unbiased by construction.
#'
#' @param X_main Numeric matrix. Main-effect columns only (no interaction terms).
#' @param y Factor. Binary outcome.
#' @param positive_label String. The positive class label.
#' @return A data.frame with columns: Marker, AUC, CI_Lower, CI_Upper.
#' @export
run_univariate_loo_auc <- function(X_main, y, positive_label) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("[ML] Package 'pROC' required for univariate AUC.")
  }

  y_bin   <- as.integer(y == positive_label)
  markers <- colnames(X_main)

  results <- purrr::map_dfr(markers, function(mk) {
    x_vec <- as.numeric(X_main[, mk])

    roc_fwd <- tryCatch(
      pROC::roc(response = y_bin, predictor = x_vec, direction = "<", quiet = TRUE),
      error = function(e) NULL
    )
    roc_rev <- tryCatch(
      pROC::roc(response = y_bin, predictor = x_vec, direction = ">", quiet = TRUE),
      error = function(e) NULL
    )

    auc_fwd <- if (!is.null(roc_fwd)) as.numeric(pROC::auc(roc_fwd)) else 0
    auc_rev <- if (!is.null(roc_rev)) as.numeric(pROC::auc(roc_rev)) else 0
    roc_obj <- if (auc_fwd >= auc_rev) roc_fwd else roc_rev
    auc_val <- max(auc_fwd, auc_rev)

    ci <- tryCatch(
      as.numeric(pROC::ci.auc(roc_obj, method = "delong")),
      error = function(e) c(NA_real_, auc_val, NA_real_)
    )

    data.frame(
      Marker   = mk,
      AUC      = round(auc_val, 4),
      CI_Lower = round(ci[1],   4),
      CI_Upper = round(ci[3],   4),
      stringsAsFactors = FALSE
    )
  })

  results
}


#' @title Univariate LOO Threshold Classification
#' @description
#' For each marker, runs outer LOO: selects the Youden-J optimal threshold on
#' the training fold ROC and applies it blindly to the held-out sample.
#' Threshold is the only learned parameter; LOO guarantees unbiased evaluation.
#'
#' @param X_main Numeric matrix. Main-effect columns only (no interaction terms).
#' @param y Factor. Binary outcome.
#' @param positive_label String. The positive class label.
#' @param seed Integer. Random seed (for reproducibility of any ties).
#' @return A data.frame with columns: Marker, AUC, Balanced_Accuracy, BER,
#'   Sensitivity, Specificity.
#' @export
run_univariate_loo_threshold <- function(X_main, y, positive_label, seed = 2026) {
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("[ML] Package 'pROC' required for LOO threshold analysis.")
  }

  y_bin   <- as.integer(y == positive_label)
  markers <- colnames(X_main)
  n       <- nrow(X_main)

  results <- purrr::map_dfr(markers, function(mk) {
    x_vec    <- as.numeric(X_main[, mk])
    preds    <- numeric(n)

    for (i in seq_len(n)) {
      x_train <- x_vec[-i]
      y_train <- y_bin[-i]

      roc_train <- tryCatch(
        pROC::roc(response = y_train, predictor = x_train, quiet = TRUE),
        error = function(e) NULL
      )

      if (is.null(roc_train)) { preds[i] <- 0.5; next }

      coords_df <- tryCatch(
        pROC::coords(roc_train, x = "best", best.method = "youden",
                     ret = "threshold", transpose = FALSE),
        error = function(e) NULL
      )
      thresh <- if (!is.null(coords_df) && nrow(coords_df) > 0) {
        as.numeric(coords_df$threshold[1])
      } else {
        median(x_train)
      }

      preds[i] <- if (x_vec[i] >= thresh) 1L else 0L
    }

    metrics <- compute_classification_metrics(y, as.numeric(preds), positive_label)

    data.frame(
      Marker            = mk,
      AUC               = round(metrics$auc,               4),
      Balanced_Accuracy = round(metrics$balanced_accuracy, 4),
      BER               = round(metrics$ber,               4),
      Sensitivity       = round(metrics$sensitivity,       4),
      Specificity       = round(metrics$specificity,       4),
      stringsAsFactors  = FALSE
    )
  })

  results
}


#' @title Plot Univariate ROC Curves
#' @description
#' Overlaid ROC curves for each marker from the univariate AUC analysis.
#' Uses the same visual style as plot_ml_roc() (bw theme, diagonal reference,
#' AUC in legend). Reconstructs ROC objects from the full dataset Z-scores.
#'
#' @param X_main Numeric matrix. Main-effect columns (same used in run_univariate_loo_auc).
#' @param y Factor. Binary outcome.
#' @param positive_label String. The positive class label.
#' @param univariate_df Data.frame returned by run_univariate_loo_auc() (for AUC labels).
#' @param colors_viz Named character vector (optional, for marker colours).
#' @param title String. Plot title.
#' @return A ggplot object, or NULL on failure.
#' @export
plot_univariate_roc <- function(X_main, y, positive_label, univariate_df,
                                colors_viz = NULL,
                                title      = "Univariate ROC Curves") {
  if (!requireNamespace("pROC", quietly = TRUE)) return(NULL)

  y_bin   <- as.integer(y == positive_label)
  markers <- colnames(X_main)

  roc_data <- purrr::map_dfr(markers, function(mk) {
    x_vec <- as.numeric(X_main[, mk])

    roc_fwd <- tryCatch(pROC::roc(response = y_bin, predictor = x_vec, direction = "<", quiet = TRUE), error = function(e) NULL)
    roc_rev <- tryCatch(pROC::roc(response = y_bin, predictor = x_vec, direction = ">", quiet = TRUE), error = function(e) NULL)
    auc_fwd <- if (!is.null(roc_fwd)) as.numeric(pROC::auc(roc_fwd)) else 0
    auc_rev <- if (!is.null(roc_rev)) as.numeric(pROC::auc(roc_rev)) else 0
    roc_obj <- if (auc_fwd >= auc_rev) roc_fwd else roc_rev
    auc_val <- max(auc_fwd, auc_rev)

    if (is.null(roc_obj)) return(data.frame())

    data.frame(
      FPR    = 1 - roc_obj$specificities,
      TPR    = roc_obj$sensitivities,
      Marker = sprintf("%s (AUC = %.3f)", mk, auc_val),
      stringsAsFactors = FALSE
    )
  })

  if (nrow(roc_data) == 0) return(NULL)

  p <- ggplot2::ggplot(roc_data, ggplot2::aes(x = FPR, y = TPR, color = Marker)) +
    ggplot2::geom_line(linewidth = 1.1, na.rm = TRUE) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "grey55", linewidth = 0.6) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0),
                                 breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0),
                                 breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(title = title, x = "1 - Specificity (FPR)",
                  y = "Sensitivity (TPR)", color = NULL) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.direction = "vertical",
      plot.title       = ggplot2::element_text(size = 11, face = "bold")
    )

  p
}
