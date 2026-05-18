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
      # Only two markers left and both are collinear with each other.
      # With no other markers to compute mean |r| against, break the tie
      # by convention (drop m2, keep m1).
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
#' Overlays ROC curves for multiple methods on a single ggplot with AUC annotated
#' in the legend. ML model curves are drawn as solid lines; optional clinical
#' benchmark curves are drawn as dashed lines at a distinct colour to signal their
#' different nature (univariate biomarker, potentially on a smaller subset).
#'
#' @param results_list Named list of ML result objects (each with predicted_probs,
#'   y_true, positive_label, method).
#' @param colors_viz Named character vector of clinical group colours (optional).
#' @param title String. Plot title.
#' @param benchmark_list Optional list of benchmark result objects from
#'   run_clinical_benchmark(). Each element must contain predicted_probs, y_true,
#'   positive_label, label, n_valid, and direction fields.
#' @return A ggplot object, or NULL on failure.
#' @export
plot_ml_roc <- function(results_list,
                        colors_viz     = NULL,
                        title          = "Nested-LOOCV ROC Curves",
                        benchmark_list = NULL) {
  if (!requireNamespace("pROC", quietly = TRUE)) return(NULL)

  method_colors    <- c("Elastic Net" = "#2171B5", "SVM-RBF" = "#CB181D")
  benchmark_colors <- c("#6A3D9A", "#FF7F00", "#33A02C")

  # ── ML model ROC data (solid curves) ────────────────────────────────────────
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
        stringsAsFactors = FALSE
      )
    }, error = function(e) data.frame())
  })

  if (nrow(roc_data) == 0) return(NULL)

  # ── Benchmark ROC data (dashed curves, possibly on a subset) ────────────────
  bench_data <- data.frame()
  if (!is.null(benchmark_list)) {
    bench_data <- purrr::map_dfr(seq_along(benchmark_list), function(bi) {
      br    <- benchmark_list[[bi]]
      y_bin <- as.integer(br$y_true == br$positive_label)
      x_vec <- as.numeric(br$predicted_probs)
      dir   <- if (!is.null(br$direction)) br$direction else "<"
      tryCatch({
        roc_obj <- pROC::roc(response = y_bin, predictor = x_vec,
                             direction = dir, quiet = TRUE)
        auc_val <- round(as.numeric(pROC::auc(roc_obj)), 3)
        data.frame(
          FPR    = 1 - roc_obj$specificities,
          TPR    = roc_obj$sensitivities,
          Method = sprintf("%s\n(AUC = %.3f, n=%d)", br$label, auc_val, br$n_valid),
          stringsAsFactors = FALSE
        )
      }, error = function(e) data.frame())
    })
  }

  # ── Colour map ───────────────────────────────────────────────────────────────
  # ML methods: keyed by short name prefix; benchmarks: cycle through palette.
  color_map <- stats::setNames(
    sapply(unique(roc_data$Method), function(m) {
      key <- sub("\n.*$", "", m)
      if (key %in% names(method_colors)) method_colors[[key]] else "#636363"
    }),
    unique(roc_data$Method)
  )
  if (nrow(bench_data) > 0) {
    bench_methods <- unique(bench_data$Method)
    color_map[bench_methods] <- benchmark_colors[
      ((seq_along(bench_methods) - 1L) %% length(benchmark_colors)) + 1L
    ]
  }

  # ── Build plot: two separate geom_line layers preserve linetype clarity ──────
  p <- ggplot2::ggplot(mapping = ggplot2::aes(x = FPR, y = TPR, color = Method)) +
    ggplot2::geom_line(data = roc_data, linewidth = 1.2, na.rm = TRUE)

  if (nrow(bench_data) > 0) {
    p <- p + ggplot2::geom_line(data = bench_data, linewidth = 1.0,
                                linetype = "dashed", na.rm = TRUE)
  }

  p <- p +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                         color = "grey55", linewidth = 0.5) +
    ggplot2::scale_color_manual(values = color_map, name = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0),
                                 breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0),
                                 breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(
      title    = title,
      x        = "1 - Specificity (FPR)",
      y        = "Sensitivity (TPR)",
      color    = NULL,
      caption  = if (nrow(bench_data) > 0) "Dashed: clinical benchmark (univariate, available-cases subset)" else NULL
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.direction = "vertical",
      plot.title       = ggplot2::element_text(size = 11, face = "bold"),
      plot.caption     = ggplot2::element_text(size = 8, color = "grey45")
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


#' @title Run Clinical Benchmark Comparison
#' @description
#' Loads a clinical biomarker column from the raw input Excel, merges by Patient_ID,
#' computes its univariate AUC on available cases, and performs a DeLong test against
#' the primary model's out-of-fold predicted probabilities. Returns NULL gracefully when
#' the column is absent, has too many NAs, or the input file cannot be read.
#'
#' @param DATA Named list. Processed data object from Step 01 (for Patient_ID and Group).
#' @param primary_probs Numeric vector. Out-of-fold predicted probabilities from the primary model.
#' @param y_primary Factor. True class labels (same order as primary_probs).
#' @param positive_label String. The positive class label.
#' @param input_file String. Path to the raw input Excel (config$input_file_t0).
#' @param benchmark_col String. Column name for the clinical biomarker.
#' @param benchmark_label String. Human-readable label for reporting.
#' @param min_n Integer. Minimum valid cases required to proceed (default 10).
#' @return A named list: label, auc, auc_ci, n_valid, n_na, delong_p, predicted_probs,
#'   y_true, positive_label. Returns NULL on failure.
#' @export
run_clinical_benchmark <- function(DATA, primary_probs, y_primary, positive_label,
                                   input_file, benchmark_col, benchmark_label = NULL,
                                   min_n = 10L) {
  if (!requireNamespace("pROC",    quietly = TRUE)) return(NULL)
  if (!requireNamespace("readxl",  quietly = TRUE)) return(NULL)
  if (!requireNamespace("dplyr",   quietly = TRUE)) return(NULL)

  if (is.null(benchmark_label)) benchmark_label <- benchmark_col

  if (!file.exists(input_file)) {
    message(sprintf("   [ML Benchmark] Input file not found: %s", input_file))
    return(NULL)
  }

  df_raw <- tryCatch(
    readxl::read_excel(input_file, sheet = 1),
    error = function(e) { message(sprintf("   [ML Benchmark] Cannot read %s: %s", input_file, e$message)); NULL }
  )
  if (is.null(df_raw)) return(NULL)

  if (!benchmark_col %in% names(df_raw)) {
    message(sprintf("   [ML Benchmark] Column '%s' not found in %s", benchmark_col, basename(input_file)))
    return(NULL)
  }

  # Merge: raw file → processed metadata by Patient_ID
  df_bench <- df_raw %>%
    dplyr::select(Patient_ID, dplyr::all_of(benchmark_col)) %>%
    dplyr::mutate(Patient_ID = as.character(Patient_ID))

  df_meta <- DATA$metadata %>%
    dplyr::mutate(Patient_ID = as.character(Patient_ID)) %>%
    dplyr::left_join(df_bench, by = "Patient_ID")

  bench_vals <- as.numeric(df_meta[[benchmark_col]])
  n_na    <- sum(is.na(bench_vals))
  n_valid <- sum(!is.na(bench_vals))

  message(sprintf("   [ML Benchmark] '%s': %d valid / %d total (NA=%d)",
                  benchmark_label, n_valid, nrow(df_meta), n_na))

  if (n_valid < min_n) {
    message(sprintf("   [ML Benchmark] Insufficient valid cases (%d < %d). Skipping.", n_valid, min_n))
    return(NULL)
  }

  valid_mask <- !is.na(bench_vals)
  y_bench    <- y_primary[valid_mask]
  y_bin      <- as.integer(y_bench == positive_label)
  x_bench    <- bench_vals[valid_mask]

  # Univariate AUC (direction auto-selected for maximum AUC)
  roc_fwd <- tryCatch(pROC::roc(y_bin, x_bench, direction = "<", quiet = TRUE), error = function(e) NULL)
  roc_rev <- tryCatch(pROC::roc(y_bin, x_bench, direction = ">", quiet = TRUE), error = function(e) NULL)
  auc_fwd <- if (!is.null(roc_fwd)) as.numeric(pROC::auc(roc_fwd)) else 0
  auc_rev <- if (!is.null(roc_rev)) as.numeric(pROC::auc(roc_rev)) else 0
  roc_best  <- if (auc_fwd >= auc_rev) roc_fwd else roc_rev
  bench_auc <- max(auc_fwd, auc_rev)

  ci_bench <- tryCatch(
    as.numeric(pROC::ci.auc(roc_best, method = "delong")),
    error = function(e) c(NA_real_, bench_auc, NA_real_)
  )

  # DeLong test: primary model (subset) vs benchmark
  probs_subset <- primary_probs[valid_mask]
  y_bin_sub    <- as.integer(y_primary[valid_mask] == positive_label)
  roc_primary  <- tryCatch(
    pROC::roc(y_bin_sub, probs_subset, direction = "<", quiet = TRUE),
    error = function(e) NULL
  )
  primary_auc_subset <- if (!is.null(roc_primary)) as.numeric(pROC::auc(roc_primary)) else NA_real_

  delong_result <- if (!is.null(roc_primary) && !is.null(roc_best)) {
    tryCatch(pROC::roc.test(roc_primary, roc_best, method = "delong"), error = function(e) NULL)
  } else NULL

  delong_p <- if (!is.null(delong_result)) round(delong_result$p.value, 5) else NA_real_

  message(sprintf(
    "   [ML Benchmark] %s: AUC=%.3f [%.3f-%.3f] | Primary(n=%d)=%.3f | DeLong p=%.4f",
    benchmark_label, bench_auc, ci_bench[1], ci_bench[3],
    n_valid, primary_auc_subset, if (is.na(delong_p)) 0 else delong_p
  ))

  list(
    label                 = benchmark_label,
    column                = benchmark_col,
    auc                   = bench_auc,
    auc_ci                = ci_bench,
    n_valid               = n_valid,
    n_na                  = n_na,
    primary_auc_on_subset = primary_auc_subset,
    delong_p              = delong_p,
    direction             = if (auc_fwd >= auc_rev) "<" else ">",
    predicted_probs       = x_bench,
    y_true                = y_bench,
    positive_label        = positive_label
  )
}


#' @title Benchmark-Stratified Subgroup Analysis
#' @description
#' Stratifies the cohort by clinical-standard categorical bins of the benchmark
#' biomarker (e.g. PD-L1 TPS at <1%, 1-49%, >=50%) and quantifies:
#'   - Frequency table of bins vs outcome (Fisher exact + Cochran-Armitage trend)
#'   - Binary cut performance ("benchmark high" vs "benchmark low") as classifier
#'   - Primary-model AUC restricted to the "benchmark low" subgroup, which is
#'     the clinically actionable population where the standard biomarker is
#'     insufficient for treatment decision
#'
#' @param DATA Step 01 standard DATA object (with metadata$Patient_ID).
#' @param X_main Numeric matrix of main-effect features (n_patients x n_main).
#' @param y Factor outcome aligned with X_main.
#' @param positive_label Positive class label.
#' @param input_file Path to raw clinical Excel (T0) with benchmark column.
#' @param benchmark_col Column name in raw Excel.
#' @param benchmark_label Optional human-readable label.
#' @param bin_breaks Numeric vector of bin breakpoints (default c(0.99, 49.99)
#'   for PD-L1 TPS clinical strata <1%, 1-49%, >=50%).
#' @param bin_labels Character vector of bin labels (default
#'   c("neg(<1%)", "low(1-49%)", "high(>=50%)")).
#' @param high_threshold Numeric. Cutoff defining "benchmark high" for the
#'   binary clinical-cut classifier (default 50).
#' @param C_grid,gamma_grid,inner_folds,seed SVM tuning hyperparameters.
#' @param min_subgroup_n Minimum patients per subgroup to attempt SVM fit
#'   (default 20). Smaller subgroups are reported with raw counts only.
#' @return A list with crosstab, fisher_p, ca_trend_p, binary_cut metrics, and
#'   per-subgroup nested-LOOCV SVM AUC (when subgroup_n >= min_subgroup_n).
#'   Returns NULL if benchmark column is unavailable.
#' @export
run_pdl1_stratified <- function(DATA, X_main, y, positive_label,
                                input_file, benchmark_col, benchmark_label = NULL,
                                bin_breaks     = c(0.99, 49.99),
                                bin_labels     = c("neg(<1%)", "low(1-49%)", "high(>=50%)"),
                                high_threshold = 50,
                                C_grid         = c(0.01, 0.1, 1, 10, 100),
                                gamma_grid     = c(0.01, 0.1, 1, 10),
                                inner_folds    = 5L,
                                seed           = 2026,
                                min_subgroup_n = 20L) {
  if (!requireNamespace("readxl", quietly = TRUE)) return(NULL)
  if (!requireNamespace("dplyr",  quietly = TRUE)) return(NULL)
  if (!requireNamespace("pROC",   quietly = TRUE)) return(NULL)

  if (is.null(benchmark_label)) benchmark_label <- benchmark_col
  if (!file.exists(input_file)) {
    message(sprintf("   [ML Stratified] Input file not found: %s", input_file))
    return(NULL)
  }

  df_raw <- tryCatch(readxl::read_excel(input_file, sheet = 1),
                     error = function(e) NULL)
  if (is.null(df_raw) || !benchmark_col %in% names(df_raw)) {
    message(sprintf("   [ML Stratified] Column '%s' not available — skipping.", benchmark_col))
    return(NULL)
  }

  df_bench <- df_raw %>%
    dplyr::select(Patient_ID, dplyr::all_of(benchmark_col)) %>%
    dplyr::mutate(Patient_ID = as.character(Patient_ID))

  df_meta <- DATA$metadata %>%
    dplyr::mutate(Patient_ID = as.character(Patient_ID)) %>%
    dplyr::left_join(df_bench, by = "Patient_ID")

  bench_vals <- as.numeric(df_meta[[benchmark_col]])
  n_valid <- sum(!is.na(bench_vals))
  if (n_valid < min_subgroup_n) {
    message(sprintf("   [ML Stratified] Insufficient valid cases (%d). Skipping.", n_valid))
    return(NULL)
  }

  # ---- Categorical bins (3 levels) ----
  bins <- cut(bench_vals,
              breaks = c(-Inf, bin_breaks, Inf),
              labels = bin_labels, include.lowest = TRUE)

  y_chr  <- as.character(y)
  pos    <- positive_label
  neg    <- setdiff(unique(y_chr), pos)
  tab3   <- table(bins, factor(y_chr, levels = c(neg, pos)), useNA = "no")
  fish3  <- tryCatch(fisher.test(tab3)$p.value, error = function(e) NA_real_)

  ca_p <- NA_real_
  if (nrow(tab3) >= 2) {
    rs <- rowSums(tab3)
    if (all(rs > 0)) {
      ca_x <- tab3[, pos]
      ca_n <- rs
      ca_p <- tryCatch(prop.trend.test(x = ca_x, n = ca_n)$p.value,
                       error = function(e) NA_real_)
    }
  }

  rr_per_bin <- as.numeric(prop.table(tab3, margin = 1)[, pos]) * 100
  bin_summary <- data.frame(
    Bin            = rownames(tab3),
    N              = as.numeric(rowSums(tab3)),
    N_Responder    = as.numeric(tab3[, pos]),
    Response_Rate  = round(rr_per_bin, 1),
    stringsAsFactors = FALSE
  )

  # ---- Binary cut (clinical threshold) as a classifier ----
  high_mask <- !is.na(bench_vals) & bench_vals >= high_threshold
  low_mask  <- !is.na(bench_vals) & bench_vals <  high_threshold
  y_b <- y_chr[high_mask | low_mask]
  pred_high <- ifelse(bench_vals[high_mask | low_mask] >= high_threshold, pos, neg)
  sens_b <- mean(pred_high[y_b == pos] == pos)
  spec_b <- mean(pred_high[y_b == neg] == neg)
  bal_b  <- (sens_b + spec_b) / 2
  tab2   <- table(factor(pred_high, levels = c(neg, pos)),
                  factor(y_b, levels = c(neg, pos)))
  fish2  <- tryCatch(fisher.test(tab2), error = function(e) NULL)
  binary_cut <- list(
    threshold         = high_threshold,
    n                 = length(y_b),
    n_high            = sum(high_mask),
    n_low             = sum(low_mask),
    balanced_accuracy = round(bal_b, 4),
    sensitivity       = round(sens_b, 4),
    specificity       = round(spec_b, 4),
    fisher_p          = if (!is.null(fish2)) round(fish2$p.value, 4) else NA_real_,
    odds_ratio        = if (!is.null(fish2)) round(as.numeric(fish2$estimate), 4) else NA_real_
  )

  # ---- Per-subgroup nested-LOOCV SVM on the LMM-gate features ----
  run_subgroup_svm <- function(mask, label) {
    n_s <- sum(mask)
    if (n_s < min_subgroup_n) {
      return(list(n = n_s,
                  n_responder = sum(y[mask] == pos),
                  auc = NA_real_, auc_ci = c(NA_real_, NA_real_, NA_real_),
                  predicted_probs = NULL, y_true = NULL, benchmark_vals = NULL,
                  note = sprintf("n=%d below min_subgroup_n=%d — not fitted.",
                                 n_s, min_subgroup_n)))
    }
    y_sub <- y[mask]
    if (length(unique(y_sub)) < 2 || min(table(y_sub)) < 3) {
      return(list(n = n_s,
                  n_responder = sum(y_sub == pos),
                  auc = NA_real_, auc_ci = c(NA_real_, NA_real_, NA_real_),
                  predicted_probs = NULL, y_true = NULL, benchmark_vals = NULL,
                  note = "Degenerate class distribution — not fitted."))
    }
    X_sub <- X_main[mask, , drop = FALSE]
    res <- run_nested_loocv_svm(
      X = X_sub, y = factor(as.character(y_sub), levels = c(neg, pos)),
      C_grid = C_grid, gamma_grid = gamma_grid,
      inner_folds = inner_folds, seed = seed
    )
    metrics <- compute_classification_metrics(res$y_true, res$predicted_probs, pos)
    list(
      n           = n_s,
      n_responder = sum(y_sub == pos),
      auc         = round(metrics$auc, 4),
      auc_ci      = round(metrics$auc_ci, 4),
      balanced_accuracy = round(metrics$balanced_accuracy, 4),
      sensitivity = round(metrics$sensitivity, 4),
      specificity = round(metrics$specificity, 4),
      predicted_probs = res$predicted_probs,
      y_true          = as.character(res$y_true),
      benchmark_vals  = as.numeric(bench_vals[mask]),
      note            = sprintf("SVM-RBF nested-LOOCV on KI67 gate, %s subgroup.", label)
    )
  }

  sub_low  <- run_subgroup_svm(low_mask,  sprintf("%s < %g", benchmark_label, high_threshold))
  sub_high <- run_subgroup_svm(high_mask, sprintf("%s >= %g", benchmark_label, high_threshold))

  list(
    label              = benchmark_label,
    column             = benchmark_col,
    n_valid            = n_valid,
    bin_breaks         = bin_breaks,
    bin_labels         = bin_labels,
    bin_crosstab       = bin_summary,
    fisher_3bins_p     = round(fish3, 4),
    ca_trend_p         = round(ca_p, 4),
    binary_cut         = binary_cut,
    subgroup_low       = sub_low,
    subgroup_high      = sub_high,
    note               = sprintf(
      "Stratified analysis on n=%d valid %s values. Categorical bins at %s. Binary clinical cut at %g. Subgroup SVM uses LMM gate features.",
      n_valid, benchmark_label,
      paste(bin_breaks, collapse = "/"), high_threshold
    )
  )
}


#' @title Combined Benchmark + Model Information Gain (IDI/NRI)
#' @description
#' Quantifies how much the LMM-gate features add over a continuous clinical
#' benchmark biomarker (e.g. PD-L1). Two complementary perspectives:
#'   1. Logistic LOOCV: benchmark alone vs benchmark + main features. Computes
#'      IDI (Pencina 2008) and continuous NRI with patient-level bootstrap CIs.
#'   2. SVM-RBF nested-LOOCV: same-subset comparison of model-with-benchmark
#'      vs model-without (informs whether the benchmark adds noise to the
#'      non-linear classifier).
#'
#' @param DATA Step 01 DATA object (for metadata).
#' @param X_main Numeric matrix of main-effect features (already z-scored at
#'   Step 01) aligned with DATA$metadata$Patient_ID.
#' @param y Factor outcome aligned with X_main.
#' @param positive_label Positive class label.
#' @param input_file Path to raw clinical Excel (T0).
#' @param benchmark_col Column name in raw Excel.
#' @param benchmark_label Optional label.
#' @param n_boot Number of patient-level bootstrap iterations for IDI/cNRI CI.
#' @param C_grid,gamma_grid,inner_folds SVM tuning.
#' @param seed RNG seed.
#' @return List with logistic, idi/nri, and SVM comparison. NULL if benchmark
#'   unavailable.
#' @export
run_combined_benchmark_model <- function(DATA, X_main, y, positive_label,
                                         input_file, benchmark_col,
                                         benchmark_label = NULL,
                                         n_boot          = 1000L,
                                         C_grid          = c(0.01, 0.1, 1, 10, 100),
                                         gamma_grid      = c(0.01, 0.1, 1, 10),
                                         inner_folds     = 5L,
                                         seed            = 2026) {
  if (!requireNamespace("readxl", quietly = TRUE)) return(NULL)
  if (!requireNamespace("dplyr",  quietly = TRUE)) return(NULL)
  if (!requireNamespace("pROC",   quietly = TRUE)) return(NULL)

  if (is.null(benchmark_label)) benchmark_label <- benchmark_col
  if (!file.exists(input_file)) return(NULL)

  df_raw <- tryCatch(readxl::read_excel(input_file, sheet = 1),
                     error = function(e) NULL)
  if (is.null(df_raw) || !benchmark_col %in% names(df_raw)) return(NULL)

  df_bench <- df_raw %>%
    dplyr::select(Patient_ID, dplyr::all_of(benchmark_col)) %>%
    dplyr::mutate(Patient_ID = as.character(Patient_ID))
  df_meta <- DATA$metadata %>%
    dplyr::mutate(Patient_ID = as.character(Patient_ID)) %>%
    dplyr::left_join(df_bench, by = "Patient_ID")

  bench_vals <- as.numeric(df_meta[[benchmark_col]])
  valid_mask <- !is.na(bench_vals)
  if (sum(valid_mask) < 20) return(NULL)

  pos <- positive_label
  neg <- setdiff(levels(y), pos)
  y_sub <- y[valid_mask]
  y_bin <- as.integer(y_sub == pos)
  x_bench <- bench_vals[valid_mask]
  X_sub  <- X_main[valid_mask, , drop = FALSE]

  # ---- Logistic LOOCV: benchmark alone vs benchmark + main features ----
  df_lr <- data.frame(y = y_bin, bench = x_bench, X_sub, check.names = FALSE)
  feat_names_lr <- make.names(colnames(X_sub), unique = TRUE)
  colnames(df_lr)[3:ncol(df_lr)] <- feat_names_lr

  formula_bench <- as.formula("y ~ bench")
  formula_comb  <- as.formula(paste("y ~ bench +", paste(feat_names_lr, collapse = " + ")))

  loocv_glm <- function(formula_obj, df) {
    n <- nrow(df); probs <- numeric(n)
    for (i in seq_len(n)) {
      fit <- suppressWarnings(
        glm(formula_obj, data = df[-i, , drop = FALSE], family = binomial())
      )
      probs[i] <- predict(fit, newdata = df[i, , drop = FALSE], type = "response")
    }
    probs
  }

  p_bench <- loocv_glm(formula_bench, df_lr)
  p_comb  <- loocv_glm(formula_comb,  df_lr)

  roc_bench <- pROC::roc(y_bin, p_bench, direction = "<", quiet = TRUE)
  roc_comb  <- pROC::roc(y_bin, p_comb,  direction = "<", quiet = TRUE)
  auc_bench <- as.numeric(pROC::auc(roc_bench))
  auc_comb  <- as.numeric(pROC::auc(roc_comb))
  ci_bench  <- tryCatch(as.numeric(pROC::ci.auc(roc_bench, method = "delong")),
                        error = function(e) c(NA_real_, auc_bench, NA_real_))
  ci_comb   <- tryCatch(as.numeric(pROC::ci.auc(roc_comb,  method = "delong")),
                        error = function(e) c(NA_real_, auc_comb,  NA_real_))
  delong_p  <- tryCatch(pROC::roc.test(roc_bench, roc_comb,
                                       method = "delong", paired = TRUE)$p.value,
                        error = function(e) NA_real_)

  # ---- IDI / Continuous NRI ----
  idi_components <- function(p_old, p_new, y_bin) {
    evt <- y_bin == 1
    imp_evt <- mean(p_new[evt]) - mean(p_old[evt])
    red_nev <- mean(p_old[!evt]) - mean(p_new[!evt])
    list(IDI = imp_evt + red_nev, IDI_event = imp_evt, IDI_nonevent = red_nev)
  }
  cnri_components <- function(p_old, p_new, y_bin) {
    evt <- y_bin == 1
    up_evt <- mean((p_new - p_old)[evt] > 0) - mean((p_new - p_old)[evt] < 0)
    dn_nev <- mean((p_new - p_old)[!evt] < 0) - mean((p_new - p_old)[!evt] > 0)
    list(NRI = up_evt + dn_nev, NRI_event = up_evt, NRI_nonevent = dn_nev)
  }

  idi <- idi_components(p_bench, p_comb, y_bin)
  nri <- cnri_components(p_bench, p_comb, y_bin)

  # Patient-level bootstrap CI for IDI / cNRI
  set.seed(seed)
  n_y <- length(y_bin)
  idi_boot <- numeric(n_boot); nri_boot <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample(n_y, n_y, replace = TRUE)
    if (length(unique(y_bin[idx])) < 2) { idi_boot[b] <- NA_real_; nri_boot[b] <- NA_real_; next }
    idi_boot[b] <- idi_components(p_bench[idx], p_comb[idx], y_bin[idx])$IDI
    nri_boot[b] <- cnri_components(p_bench[idx], p_comb[idx], y_bin[idx])$NRI
  }
  idi_ci <- as.numeric(quantile(idi_boot, c(0.025, 0.975), na.rm = TRUE))
  nri_ci <- as.numeric(quantile(nri_boot, c(0.025, 0.975), na.rm = TRUE))
  idi_p  <- 2 * min(mean(idi_boot <= 0, na.rm = TRUE), mean(idi_boot > 0, na.rm = TRUE))
  nri_p  <- 2 * min(mean(nri_boot <= 0, na.rm = TRUE), mean(nri_boot > 0, na.rm = TRUE))

  # ---- SVM-RBF on the same subset: benchmark + features vs features only ----
  bench_z <- as.numeric(scale(x_bench))  # standardize benchmark to feature scale
  X_with  <- cbind(BENCHMARK = bench_z, X_sub)
  res_with <- run_nested_loocv_svm(X = X_with, y = y_sub,
                                   C_grid = C_grid, gamma_grid = gamma_grid,
                                   inner_folds = inner_folds, seed = seed)
  res_only <- run_nested_loocv_svm(X = X_sub, y = y_sub,
                                   C_grid = C_grid, gamma_grid = gamma_grid,
                                   inner_folds = inner_folds, seed = seed)
  m_with <- compute_classification_metrics(res_with$y_true, res_with$predicted_probs, pos)
  m_only <- compute_classification_metrics(res_only$y_true, res_only$predicted_probs, pos)
  roc_with <- tryCatch(pROC::roc(as.integer(res_with$y_true == pos),
                                 res_with$predicted_probs, direction = "<", quiet = TRUE),
                       error = function(e) NULL)
  roc_only <- tryCatch(pROC::roc(as.integer(res_only$y_true == pos),
                                 res_only$predicted_probs, direction = "<", quiet = TRUE),
                       error = function(e) NULL)
  svm_delong_p <- if (!is.null(roc_with) && !is.null(roc_only)) {
    tryCatch(pROC::roc.test(roc_only, roc_with, method = "delong", paired = TRUE)$p.value,
             error = function(e) NA_real_)
  } else NA_real_

  list(
    label   = benchmark_label,
    column  = benchmark_col,
    n_valid = sum(valid_mask),
    logistic = list(
      auc_benchmark           = round(auc_bench, 4),
      auc_benchmark_ci        = as.list(round(ci_bench, 4)),
      auc_combined            = round(auc_comb, 4),
      auc_combined_ci         = as.list(round(ci_comb, 4)),
      delta_auc               = round(auc_comb - auc_bench, 4),
      delong_p                = round(delong_p, 5)
    ),
    information_gain = list(
      IDI                = round(idi$IDI, 4),
      IDI_event          = round(idi$IDI_event, 4),
      IDI_nonevent       = round(idi$IDI_nonevent, 4),
      IDI_95CI           = as.list(round(idi_ci, 4)),
      IDI_bootstrap_p    = round(idi_p, 4),
      cNRI               = round(nri$NRI, 4),
      cNRI_event         = round(nri$NRI_event, 4),
      cNRI_nonevent      = round(nri$NRI_nonevent, 4),
      cNRI_95CI          = as.list(round(nri_ci, 4)),
      cNRI_bootstrap_p   = round(nri_p, 4),
      n_boot             = n_boot
    ),
    svm_comparison = list(
      auc_features_only        = round(m_only$auc, 4),
      auc_features_with_bench  = round(m_with$auc, 4),
      delta_auc                = round(m_with$auc - m_only$auc, 4),
      delong_p                 = round(svm_delong_p, 5),
      note = sprintf(
        "Same n=%d subset. Negative delta means adding %s to the SVM-RBF kernel inputs degrades AUC — interpretable as the benchmark contributing noise relative to the LMM gate.",
        sum(valid_mask), benchmark_label
      )
    ),
    note = sprintf(
      "Logistic LOOCV: %s alone vs %s + %d gate marker(s). IDI/cNRI estimated with patient-level bootstrap (B=%d).",
      benchmark_label, benchmark_label, ncol(X_sub), n_boot
    ),
    plot_data = list(
      y_bin                      = y_bin,
      positive_label             = pos,
      p_bench_logistic           = p_bench,
      p_combined_logistic        = p_comb,
      p_svm_features_only        = res_only$predicted_probs,
      p_svm_features_with_bench  = res_with$predicted_probs,
      auc_svm_features_only      = round(m_only$auc, 4),
      auc_svm_features_with_bench= round(m_with$auc, 4)
    )
  )
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


# ==============================================================================
# NESTED LOO VALIDATION — fully-nested LMM feature selection inside outer fold
# ==============================================================================

#' @title LMM Gate Selection for a Single LOO Fold
#' @description
#' Fits per-marker LMMs on a longitudinal training subset (n-1 patients), applies
#' BH correction, and runs LOO sensitivity on FDR-passing markers to produce the
#' same two-stage gate used by Step 04, but trained on fold-local data only.
#'
#' @param df_lng_train Data.frame. Longitudinal data (meta + markers) for n-1 patients.
#' @param markers Character vector. All candidate marker names.
#' @param fdr_thr Numeric. FDR threshold (default 0.05).
#' @param loo_thr Numeric. Max LOO p-value threshold (default 0.05).
#' @return Character vector of markers that pass both gates. May be empty.
select_gate_for_fold <- function(df_lng_train, markers, fdr_thr = 0.05, loo_thr = 0.05) {
  results <- lapply(markers, function(mk) {
    vals   <- as.numeric(df_lng_train[[mk]])
    val_sd <- sd(vals, na.rm = TRUE)
    if (is.na(val_sd) || val_sd == 0) val_sd <- 1
    df_m <- data.frame(
      Value = vals / val_sd,
      Time  = as.factor(df_lng_train$Timepoint),
      Group = as.factor(df_lng_train$Group),
      ID    = as.factor(df_lng_train$Patient_ID)
    )
    df_m <- df_m[complete.cases(df_m), ]
    if (nrow(df_m) < 10 || length(unique(df_m$Group)) < 2)
      return(data.frame(Marker = mk, P = NA_real_))
    tryCatch({
      mod <- suppressMessages(suppressWarnings(
        lmerTest::lmer(Value ~ Time * Group + (1 | ID), data = df_m,
                       REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE))
      ))
      ct  <- summary(mod)$coefficients
      idx <- grep("Time.*:Group", rownames(ct))
      if (length(idx) != 1) return(data.frame(Marker = mk, P = NA_real_))
      p_col <- grep("Pr\\(>\\|t\\|\\)", colnames(ct))
      data.frame(Marker = mk, P = ct[idx, p_col])
    }, error = function(e) data.frame(Marker = mk, P = NA_real_))
  })

  res_df     <- do.call(rbind, results)
  res_df$FDR <- p.adjust(res_df$P, method = "BH")
  fdr_pass   <- res_df$Marker[!is.na(res_df$FDR) & res_df$FDR < fdr_thr]
  if (length(fdr_pass) == 0) return(character(0))

  loo_pass <- character(0)
  for (mk in fdr_pass) {
    unique_ids <- unique(df_lng_train$Patient_ID)
    max_p <- 0
    for (drop_id in unique_ids) {
      sub2  <- df_lng_train[df_lng_train$Patient_ID != drop_id, ]
      vals2 <- as.numeric(sub2[[mk]])
      sd2   <- sd(vals2, na.rm = TRUE)
      if (is.na(sd2) || sd2 == 0) sd2 <- 1
      df_m2 <- data.frame(
        Value = vals2 / sd2,
        Time  = as.factor(sub2$Timepoint),
        Group = as.factor(sub2$Group),
        ID    = as.factor(sub2$Patient_ID)
      )
      df_m2 <- df_m2[complete.cases(df_m2), ]
      if (nrow(df_m2) < 8 || length(unique(df_m2$Group)) < 2) next
      p_loo <- tryCatch({
        mod2 <- suppressMessages(suppressWarnings(
          lmerTest::lmer(Value ~ Time * Group + (1 | ID), data = df_m2,
                         REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE))
        ))
        ct2  <- summary(mod2)$coefficients
        idx2 <- grep("Time.*:Group", rownames(ct2))
        if (length(idx2) != 1) NA_real_
        else ct2[idx2, grep("Pr\\(>\\|t\\|\\)", colnames(ct2))]
      }, error = function(e) NA_real_)
      if (!is.na(p_loo) && p_loo > max_p) max_p <- p_loo
    }
    if (max_p < loo_thr) loo_pass <- c(loo_pass, mk)
  }
  loo_pass
}


#' @title Fully-Nested LOO Validation for SVM-RBF Classifier
#' @description
#' For each outer LOO fold, LMM-based feature selection (FDR + LOO sensitivity +
#' collinearity filter) is re-run on the n-1 training patients using the
#' longitudinal dataset, then an SVM-RBF is trained on T0 cross-sectional data
#' and tested on the held-out patient. This provides an unbiased AUC estimate
#' and quantifies the stability of the feature gate across folds.
#'
#' @param DATA_T0         Named list. Step 01 standard (T0) processed data object.
#' @param DATA_LONG       Named list. Step 01 longitudinal processed data object.
#' @param fdr_thr         Numeric. FDR threshold for LMM gate (default 0.05).
#' @param loo_thr         Numeric. LOO p-value threshold for LMM gate (default 0.05).
#' @param cor_threshold   Numeric. Collinearity filter threshold (default 0.85).
#' @param C_grid          Numeric vector. SVM cost values for inner CV.
#' @param gamma_grid      Numeric vector. SVM gamma values for inner CV.
#' @param inner_folds     Integer. Number of inner CV folds (default 5).
#' @param seed            Integer. Random seed for reproducibility.
#' @return Named list: metrics, gate_stability (data.frame), n_empty_folds, predicted_probs.
#' @export
run_nested_loocv_svm_validated <- function(DATA_T0, DATA_LONG,
                                           fdr_thr        = 0.05,
                                           loo_thr        = 0.05,
                                           cor_threshold  = 0.85,
                                           C_grid         = c(0.01, 0.1, 1, 10, 100),
                                           gamma_grid     = c(0.01, 0.1, 1, 10),
                                           inner_folds    = 5L,
                                           seed           = 2026L,
                                           positive_label = NULL) {

  META_COLS   <- c("Patient_ID", "Sample_ID", "Timepoint", "Group")
  std_markers <- setdiff(colnames(DATA_T0$hybrid_data_z),   META_COLS)
  lng_markers <- setdiff(colnames(DATA_LONG$hybrid_data_z), META_COLS)
  shared_mkrs <- intersect(std_markers, lng_markers)

  data_z_t0 <- as.matrix(DATA_T0$hybrid_data_z[, shared_mkrs, drop = FALSE])
  meta_t0   <- DATA_T0$metadata

  # Preserve existing factor level ordering (non-responder first, responder second)
  grp_fac   <- if (is.factor(meta_t0$Group)) meta_t0$Group
               else factor(meta_t0$Group)
  grp_lvls  <- levels(grp_fac)
  pos_label <- if (!is.null(positive_label)) positive_label else grp_lvls[2]
  neg_label <- setdiff(grp_lvls, pos_label)[1]

  pid_t0 <- meta_t0$Patient_ID
  n_t0   <- nrow(meta_t0)

  data_z_lng <- as.data.frame(DATA_LONG$hybrid_data_z[, shared_mkrs, drop = FALSE])
  meta_lng   <- DATA_LONG$metadata

  probs         <- numeric(n_t0)
  gate_per_fold <- vector("list", n_t0)

  for (i in seq_len(n_t0)) {
    pid_test  <- pid_t0[i]
    idx_train <- setdiff(seq_len(n_t0), i)
    X_train   <- data_z_t0[idx_train, , drop = FALSE]
    # Use group labels directly so SVM probability columns match pos_label / neg_label
    y_train   <- as.character(grp_fac[idx_train])
    X_test    <- data_z_t0[i, , drop = FALSE]

    lng_mask     <- meta_lng$Patient_ID != pid_test
    df_lng_train <- cbind(as.data.frame(meta_lng[lng_mask, ]),
                          as.data.frame(data_z_lng[lng_mask, , drop = FALSE]))

    gate <- select_gate_for_fold(df_lng_train, shared_mkrs, fdr_thr, loo_thr)

    if (length(gate) > 1) {
      filt <- filter_collinear_features(X_train[, gate, drop = FALSE], cor_threshold)
      gate <- filt$kept
    }
    gate_per_fold[[i]] <- gate

    if (length(gate) == 0) { probs[i] <- 0.5; next }

    X_tr   <- X_train[, gate, drop = FALSE]
    X_te   <- X_test[,  gate, drop = FALSE]
    mu     <- colMeans(X_tr, na.rm = TRUE)
    sds    <- apply(X_tr, 2, sd, na.rm = TRUE); sds[sds == 0] <- 1
    X_tr_s <- scale(X_tr, center = mu, scale = sds)
    X_te_s <- scale(X_te, center = mu, scale = sds)

    y_tr_fac   <- factor(y_train, levels = grp_lvls)
    safe_folds <- max(2L, min(inner_folds, min(table(y_train))))
    set.seed(seed + i)
    tune_res <- tryCatch(
      e1071::tune(e1071::svm, train.x = X_tr_s, train.y = y_tr_fac,
                  kernel = "radial", probability = TRUE,
                  ranges      = list(cost = C_grid, gamma = gamma_grid),
                  tunecontrol = e1071::tune.control(sampling = "cross", cross = safe_folds)),
      error = function(e) NULL
    )
    if (is.null(tune_res)) { probs[i] <- 0.5; next }

    fit <- tryCatch(
      e1071::svm(X_tr_s, y_tr_fac, kernel = "radial",
                 cost  = tune_res$best.parameters$cost,
                 gamma = tune_res$best.parameters$gamma,
                 probability = TRUE),
      error = function(e) NULL
    )
    if (is.null(fit)) { probs[i] <- 0.5; next }

    pred_obj <- predict(fit, X_te_s, probability = TRUE)
    pm       <- attr(pred_obj, "probabilities")
    probs[i] <- if (pos_label %in% colnames(pm)) pm[, pos_label] else 0.5
  }

  metrics <- compute_classification_metrics(grp_fac, probs, pos_label)

  all_gates <- unlist(gate_per_fold)
  freq_tab  <- sort(table(all_gates), decreasing = TRUE)
  gate_stab <- data.frame(
    Marker         = names(freq_tab),
    Folds_Selected = as.integer(freq_tab),
    Total_Folds    = n_t0,
    Pct            = round(100 * as.integer(freq_tab) / n_t0, 1),
    stringsAsFactors = FALSE
  )

  list(metrics        = metrics,
       gate_stability  = gate_stab,
       n_empty_folds   = sum(sapply(gate_per_fold, length) == 0L),
       predicted_probs = probs)
}


#' @title Plot Benchmark-Stratified Subgroup Analysis (2-panel)
#' @description
#' Renders a publication-ready 2-panel figure for run_pdl1_stratified() output:
#'   - Panel A: bar chart of response rate by benchmark bin (with N per bin
#'     labelled on bars) plus the chi-squared trend test p-value
#'   - Panel B: ROC curve of the model restricted to the "benchmark-low"
#'     subgroup (typically the clinically actionable population) with AUC and
#'     95% CI annotated. A single point at the benchmark binary-cut performance
#'     in the same subgroup is overlaid as a reference, when defined.
#'
#' Uses patchwork to lay out the two panels side by side. Pass the same
#' colors_viz used elsewhere to keep palette consistency.
#'
#' @param stratified_result List returned by run_pdl1_stratified().
#' @param colors_viz Named vector of clinical colours (optional; if missing,
#'   responder/non-responder default to standard journal palette).
#' @param positive_label Positive class label (e.g. "RP"). Defaults to inferred.
#' @param title Optional overall title for the figure.
#' @return A patchwork object combining the two panels, or NULL on failure.
#' @export
plot_benchmark_stratified <- function(stratified_result,
                                      colors_viz     = NULL,
                                      positive_label = NULL,
                                      title          = NULL) {
  if (is.null(stratified_result)) return(NULL)
  if (!requireNamespace("ggplot2",   quietly = TRUE)) return(NULL)
  if (!requireNamespace("patchwork", quietly = TRUE)) return(NULL)
  if (!requireNamespace("pROC",      quietly = TRUE)) return(NULL)

  pos_lbl <- if (!is.null(positive_label)) positive_label else "Responder"
  resp_color  <- if (!is.null(colors_viz) && pos_lbl %in% names(colors_viz)) colors_viz[[pos_lbl]] else "#2E8B57"
  nonresp_color <- if (!is.null(colors_viz)) {
    cn <- setdiff(names(colors_viz), pos_lbl)
    if (length(cn) > 0) colors_viz[[cn[1]]] else "#B2182B"
  } else "#B2182B"

  # Short label for titles — strip trailing "expression (%)" / "(%)" noise
  bench_short <- sub("\\s*(expression\\s*\\(%\\)|\\(%\\)|\\(.*\\))\\s*$", "",
                     stratified_result$label)
  bench_short <- trimws(bench_short)
  if (bench_short == "") bench_short <- stratified_result$label

  # ---- Panel A: Response rate bar chart ----
  bin_df <- stratified_result$bin_crosstab
  bin_df$Bin <- factor(bin_df$Bin, levels = stratified_result$bin_labels)
  bin_df$Label <- sprintf("%.1f%%\n(N=%d, R=%d)",
                          bin_df$Response_Rate, bin_df$N, bin_df$N_Responder)

  pA <- ggplot2::ggplot(bin_df, ggplot2::aes(x = Bin, y = Response_Rate)) +
    ggplot2::geom_col(fill = resp_color, alpha = 0.85, width = 0.65) +
    ggplot2::geom_text(ggplot2::aes(label = Label),
                       vjust = -0.4, size = 4.2, lineheight = 0.85) +
    ggplot2::geom_hline(yintercept = 50, linetype = "dashed", color = "grey60") +
    ggplot2::scale_y_continuous(limits = c(0, max(bin_df$Response_Rate) + 18),
                                 expand = c(0, 0),
                                 breaks = seq(0, 100, 25),
                                 labels = function(x) paste0(x, "%")) +
    ggplot2::labs(
      title = sprintf("Response rate by %s strata", bench_short),
      subtitle = sprintf("3-bin Fisher p = %.3f | Cochran-Armitage trend p = %.3f | N total = %d",
                         stratified_result$fisher_3bins_p,
                         stratified_result$ca_trend_p,
                         stratified_result$n_valid),
      x = stratified_result$label, y = "Responder rate"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "gray35", size = 11),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x   = ggplot2::element_text(face = "bold")
    )

  # ---- Panel B: ROC in benchmark-low subgroup ----
  sl <- stratified_result$subgroup_low
  pB <- NULL
  if (!is.null(sl$predicted_probs) && !is.na(sl$auc)) {
    y_bin   <- as.integer(sl$y_true == pos_lbl)
    roc_obj <- tryCatch(pROC::roc(y_bin, sl$predicted_probs, direction = "<", quiet = TRUE),
                        error = function(e) NULL)
    if (!is.null(roc_obj)) {
      roc_df <- data.frame(
        FPR = 1 - roc_obj$specificities,
        TPR = roc_obj$sensitivities
      )

      # Binary-cut reference point (TPS>=threshold within the full cohort)
      bc <- stratified_result$binary_cut
      bc_x <- 1 - bc$specificity; bc_y <- bc$sensitivity

      auc_lab <- sprintf("KI67-gate model (subgroup)\nAUC = %.3f [%.3f-%.3f]",
                         sl$auc, sl$auc_ci[1], sl$auc_ci[3])

      bench_threshold <- bc$threshold
      bench_lab <- sprintf("%s >=%g cut\n(full cohort)", bench_short, bench_threshold)

      pB <- ggplot2::ggplot(roc_df, ggplot2::aes(x = FPR, y = TPR)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                              color = "grey55", linewidth = 0.6) +
        ggplot2::geom_line(color = resp_color, linewidth = 1.3) +
        ggplot2::annotate("point",  x = bc_x, y = bc_y, color = nonresp_color,
                          shape = 18, size = 4.5) +
        ggplot2::annotate("text",   x = bc_x + 0.04, y = bc_y - 0.02,
                          label = bench_lab, hjust = 0, size = 3.4,
                          color = nonresp_color, lineheight = 0.9) +
        ggplot2::annotate("label",  x = 0.98, y = 0.06, label = auc_lab,
                          hjust = 1, size = 3.8, color = resp_color, fill = "white",
                          label.size = 0.4, fontface = "bold", lineheight = 0.9) +
        ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0.005, 0)) +
        ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0.005, 0)) +
        ggplot2::labs(
          title = sprintf("Model in %s < %g subgroup",
                          bench_short, bc$threshold),
          subtitle = sprintf("n = %d (responders=%d) | clinical context where %s alone is non-discriminative",
                             sl$n, sl$n_responder, bench_short),
          x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)"
        ) +
        ggplot2::theme_bw(base_size = 13) +
        ggplot2::theme(
          plot.title    = ggplot2::element_text(face = "bold"),
          plot.subtitle = ggplot2::element_text(color = "gray35", size = 11),
          panel.grid.minor = ggplot2::element_blank()
        )
    }
  }

  if (is.null(pB)) {
    out <- pA
  } else {
    out <- patchwork::wrap_plots(pA, pB, ncol = 2, widths = c(1, 1))
  }
  if (!is.null(title)) {
    out <- out + patchwork::plot_annotation(
      title = title,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14))
    )
  }
  out
}


#' @title Plot Combined-Model Information Gain (2-panel)
#' @description
#' Publication figure for run_combined_benchmark_model(). Panel A overlays the
#' three ROC curves of interest on the same n_valid subset (benchmark alone,
#' benchmark + LMM gate combined via logistic LOOCV, and the LMM-gate SVM-RBF
#' nested-LOOCV). Panel B is a forest-style plot of the three information-gain
#' point estimates (IDI, cNRI, ΔAUC) with their bootstrap 95% CIs, mapping
#' directly onto the Pencina reclassification framework expected by reviewers.
#'
#' @param combined_result List returned by run_combined_benchmark_model().
#' @param title Optional overall title.
#' @return A patchwork object, or NULL on failure.
#' @export
plot_combined_information_gain <- function(combined_result, title = NULL) {
  if (is.null(combined_result) || is.null(combined_result$plot_data)) return(NULL)
  if (!requireNamespace("ggplot2",   quietly = TRUE)) return(NULL)
  if (!requireNamespace("patchwork", quietly = TRUE)) return(NULL)
  if (!requireNamespace("pROC",      quietly = TRUE)) return(NULL)

  pd <- combined_result$plot_data
  lg <- combined_result$logistic
  ig <- combined_result$information_gain
  sv <- combined_result$svm_comparison

  y_bin <- pd$y_bin
  bench_label <- combined_result$label
  bench_short <- sub("\\s*(expression\\s*\\(%\\)|\\(%\\)|\\(.*\\))\\s*$", "", bench_label)
  bench_short <- trimws(bench_short)
  if (bench_short == "") bench_short <- bench_label

  build_roc <- function(p_vec, label, auc_val, auc_ci_lo, auc_ci_hi) {
    roc_obj <- tryCatch(pROC::roc(y_bin, p_vec, direction = "<", quiet = TRUE),
                        error = function(e) NULL)
    if (is.null(roc_obj)) return(NULL)
    data.frame(
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      Model = sprintf("%s (AUC=%.3f [%.3f-%.3f])",
                      label, auc_val, auc_ci_lo, auc_ci_hi),
      stringsAsFactors = FALSE
    )
  }

  # SVM CIs not stored explicitly; reconstruct quickly for labelling
  roc_only <- pROC::roc(y_bin, pd$p_svm_features_only,       direction = "<", quiet = TRUE)
  ci_only  <- tryCatch(as.numeric(pROC::ci.auc(roc_only, method = "delong")),
                       error = function(e) c(NA, pd$auc_svm_features_only, NA))
  roc_with <- pROC::roc(y_bin, pd$p_svm_features_with_bench, direction = "<", quiet = TRUE)
  ci_with  <- tryCatch(as.numeric(pROC::ci.auc(roc_with, method = "delong")),
                       error = function(e) c(NA, pd$auc_svm_features_with_bench, NA))

  rocs <- rbind(
    build_roc(pd$p_bench_logistic,    sprintf("%s alone (logistic)", bench_short),
              lg$auc_benchmark,
              lg$auc_benchmark_ci[[1]], lg$auc_benchmark_ci[[3]]),
    build_roc(pd$p_combined_logistic, sprintf("%s + KI67 gate (logistic)", bench_short),
              lg$auc_combined,
              lg$auc_combined_ci[[1]], lg$auc_combined_ci[[3]]),
    build_roc(pd$p_svm_features_only, "KI67 gate alone (SVM-RBF)",
              pd$auc_svm_features_only, ci_only[1], ci_only[3])
  )

  model_colors <- c("#7A5CA6", "#2E8B57", "#B2182B")
  names(model_colors) <- unique(rocs$Model)

  pA <- ggplot2::ggplot(rocs, ggplot2::aes(x = FPR, y = TPR, color = Model)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                          color = "grey55", linewidth = 0.6) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = model_colors) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0.005, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0.005, 0)) +
    ggplot2::labs(
      title    = "ROC: benchmark, combined, and gate-only",
      subtitle = sprintf("Same cohort subset (n = %d, non-NA %s)",
                         combined_result$n_valid, bench_short),
      x        = "1 - Specificity (FPR)", y = "Sensitivity (TPR)",
      color    = NULL
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "gray35", size = 11),
      legend.position = "bottom",
      legend.direction = "vertical",
      legend.text   = ggplot2::element_text(size = 9.5),
      legend.key.height = ggplot2::unit(8, "pt"),
      panel.grid.minor = ggplot2::element_blank()
    )

  # ---- Panel B: Information-gain forest ----
  forest_df <- data.frame(
    Metric = c(
      sprintf("IDI\n(bootstrap p = %.3f)",      ig$IDI_bootstrap_p),
      sprintf("Continuous NRI\n(bootstrap p = %.3f)", ig$cNRI_bootstrap_p),
      sprintf("Delta AUC (logistic)\n(DeLong p = %.3f)", lg$delong_p),
      sprintf("Delta AUC (SVM combined)\n(DeLong p = %.3f)", sv$delong_p)
    ),
    Estimate = c(ig$IDI, ig$cNRI, lg$delta_auc, sv$delta_auc),
    Lower    = c(ig$IDI_95CI[[1]], ig$cNRI_95CI[[1]], NA_real_, NA_real_),
    Upper    = c(ig$IDI_95CI[[2]], ig$cNRI_95CI[[2]], NA_real_, NA_real_),
    stringsAsFactors = FALSE
  )
  forest_df$Metric <- factor(forest_df$Metric, levels = rev(forest_df$Metric))
  forest_df$Significant <- ifelse(
    !is.na(forest_df$Lower) & forest_df$Lower > 0, "Positive (CI > 0)",
    ifelse(!is.na(forest_df$Lower) & forest_df$Upper < 0, "Negative (CI < 0)",
           ifelse(forest_df$Estimate > 0, "Positive (point)", "Negative (point)"))
  )
  sig_colors <- c(
    "Positive (CI > 0)"  = "#2E8B57",
    "Negative (CI < 0)"  = "#B2182B",
    "Positive (point)"   = "#7AC480",
    "Negative (point)"   = "#E08983"
  )

  pB <- ggplot2::ggplot(forest_df,
                        ggplot2::aes(x = Estimate, y = Metric, color = Significant)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = Lower, xmax = Upper),
                             height = 0.2, linewidth = 0.9, na.rm = TRUE) +
    ggplot2::geom_point(size = 4.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.3f", Estimate)),
                        hjust = -0.25, vjust = -0.7, size = 4, fontface = "bold",
                        color = "black") +
    ggplot2::scale_color_manual(values = sig_colors, na.translate = FALSE) +
    ggplot2::labs(
      title    = sprintf("Information gain (KI67 gate vs %s)", bench_short),
      subtitle = sprintf("Point estimates with 95%% bootstrap CIs (B=%d) where applicable",
                         ig$n_boot),
      x        = "Estimate (positive = KI67 gate adds value)",
      y        = NULL,
      color    = NULL
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "gray35", size = 11),
      legend.position = "bottom",
      axis.text.y   = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(8, 50, 8, 8))

  out <- patchwork::wrap_plots(pA, pB, ncol = 2, widths = c(1.05, 1))
  if (!is.null(title)) {
    out <- out + patchwork::plot_annotation(
      title = title,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14))
    )
  }
  out
}
