# R/modules_network.R
# ==============================================================================
# NETWORK INFERENCE MODULE
# Description: Partial Correlation (Shrinkage), Bootstrap, and Topology.
# Dependencies: corpcor, foreach, igraph, dplyr
# ==============================================================================

library(corpcor)
library(dplyr)
library(foreach)

#' @title Infer Robust Partial Correlation Network
#' @description 
#' Uses Schafer-Strimmer Shrinkage to directly estimate partial correlations.
#' 
#' @param mat Numeric matrix (Samples x Features).
#' @param fixed_lambda Optional numeric. If provided, enforces this shrinkage intensity.
#'        If NULL, lambda is estimated analytically (James-Stein).
#' @return Partial correlation matrix with diagonal = 1.
infer_network_pcor <- function(mat, fixed_lambda = NULL) {
  
  if (is.null(fixed_lambda)) {
    pcor_mat <- corpcor::pcor.shrink(mat, verbose = FALSE)
  } else {
    pcor_mat <- corpcor::pcor.shrink(mat, lambda = fixed_lambda, verbose = FALSE)
  }
  
  return(as.matrix(pcor_mat))
}

#' @title Single Bootstrap Worker
#' @description 
#' Worker function for parallel loops. Supports fixed lambda injection.
#' 
#' @param mat Numeric matrix (Samples x Features).
#' @param n Number of samples to resample.
#' @param lambda_val Optional fixed lambda value to pass to inference.
boot_worker_pcor <- function(mat, n, lambda_val = NULL) {
  idx <- sample(1:n, n, replace = TRUE)
  boot_sample <- mat[idx, ]

  # Discard resamples with constant columns (zero variance breaks pcor)
  col_vars <- apply(boot_sample, 2, var, na.rm = TRUE)
  if (any(col_vars == 0 | is.na(col_vars))) {
    return(NULL)
  }
  
  tryCatch({
    return(infer_network_pcor(boot_sample, fixed_lambda = lambda_val))
  }, error = function(e) return(NULL))
}

#' @title Check Bootstrap Yield
#' @description Warnings if too many bootstrap iterations were discarded due to zero variance.
#' @param res_obj The result list from aggregate_boot_results.
#' @param n_req The original number of requested bootstraps.
#' @param label String label for the group (e.g., "Control").
check_boot_yield <- function(res_obj, n_req, label) {
  if (is.null(res_obj) || is.null(res_obj$n_boot_valid)) return()
  
  if (res_obj$n_boot_valid < (n_req * 0.8)) {
    warning(sprintf(
      "[WARN] Low bootstrap yield for %s: %d/%d valid iterations. Results may be unstable.", 
      label, res_obj$n_boot_valid, n_req
    ))
  } else {
    message(sprintf("   -> %s Bootstrap Yield: %d/%d (OK)", label, res_obj$n_boot_valid, n_req))
  }
}

#' @title Aggregate Bootstrap Results
#' @description Calculates edge stability and mean weights based on CI.
aggregate_boot_results <- function(boot_list, alpha = 0.10) {
  valid_boots <- Filter(Negate(is.null), boot_list)
  n_boots <- length(valid_boots)
  if (n_boots == 0) stop("All bootstrap iterations failed.")
  
  p <- nrow(valid_boots[[1]])
  nodes <- rownames(valid_boots[[1]])
  
  adj_mat <- matrix(0, p, p, dimnames = list(nodes, nodes))
  mean_mat <- matrix(0, p, p, dimnames = list(nodes, nodes))
  stab_mat <- matrix(0, p, p, dimnames = list(nodes, nodes))
  
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      vals <- sapply(valid_boots, function(m) m[i, j])
      ci <- quantile(vals, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
      
      # Select edge if CI does not cross zero
      if ((ci[1] > 0) || (ci[2] < 0)) {
        adj_mat[i, j] <- 1; adj_mat[j, i] <- 1
        mu <- mean(vals, na.rm = TRUE)
        mean_mat[i, j] <- mu; mean_mat[j, i] <- mu
      }
      stab_mat[i, j] <- sum(sign(vals) == sign(mean(vals))) / n_boots
      stab_mat[j, i] <- stab_mat[i, j]
    }
  }
  return(list(adj = adj_mat, weights = mean_mat, stability = stab_mat, n_boot_valid = n_boots))
}

#' @title Select Network Features via sPLS-DA Loadings
#' @description
#' Reduces the marker set passed to network inference using sPLS-DA variable
#' importance scores produced by Step 03. Markers with absolute loading above
#' \code{min_abs_loading} in at least one sPLS-DA component are retained, plus
#' any markers listed in \code{always_include}.
#'
#' This two-stage design (discriminant pre-selection → partial correlation network)
#' keeps p/n < 1 for reliable James-Stein shrinkage estimation, and constrains the
#' network to co-regulatory patterns among markers already shown to be relevant to
#' the clinical stratification — a requirement for interpretability in a publication.
#'
#' @param available_markers Character vector of markers present in the data matrix.
#' @param json_path Path to the Machine_Metrics JSON produced by Step 03.
#' @param min_abs_loading Non-negative scalar. Minimum absolute sPLS-DA loading to
#'   retain a marker. Default 0 keeps any marker selected in at least one component.
#' @param always_include Character vector of markers to force-include regardless of
#'   their sPLS-DA importance (e.g. structural backbone markers).
#' @param min_features Minimum number of markers required after filtering. Falls back
#'   to the full panel with a warning if the selection is smaller.
#' @return Named list: \code{selected} (character), \code{n_splsda} (int),
#'   \code{n_forced} (int).
select_network_features <- function(available_markers,
                                    json_path,
                                    min_abs_loading = 0,
                                    always_include  = character(0),
                                    min_features    = 8) {
  fallback <- list(selected = available_markers, n_splsda = NA_integer_, n_forced = 0L)

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    warning("[FeatureSelection] jsonlite unavailable — using full panel.")
    return(fallback)
  }
  if (!file.exists(json_path)) {
    warning(sprintf("[FeatureSelection] JSON not found at '%s' — using full panel.", json_path))
    return(fallback)
  }

  metrics <- tryCatch(
    jsonlite::read_json(json_path, simplifyVector = TRUE),
    error = function(e) {
      warning(sprintf("[FeatureSelection] Could not parse JSON (%s) — using full panel.", e$message))
      NULL
    }
  )
  if (is.null(metrics) || is.null(metrics$top_drivers)) {
    warning("[FeatureSelection] 'top_drivers' absent from JSON — using full panel.")
    return(fallback)
  }

  drivers <- as.data.frame(metrics$top_drivers)
  if (!all(c("Marker", "Importance") %in% names(drivers))) {
    warning("[FeatureSelection] Required columns missing from top_drivers — using full panel.")
    return(fallback)
  }

  splsda_sel <- intersect(
    drivers$Marker[drivers$Importance > min_abs_loading],
    available_markers
  )
  forced   <- intersect(always_include, available_markers)
  selected <- union(splsda_sel, forced)

  if (length(selected) < min_features) {
    warning(sprintf(
      "[FeatureSelection] %d features pass filter (min=%d) — reverting to full panel.",
      length(selected), min_features
    ))
    return(list(selected = available_markers,
                n_splsda = length(splsda_sel),
                n_forced = length(forced)))
  }

  message(sprintf(
    "[FeatureSelection] %d features selected: %d sPLS-DA (|loading|>%.2f) + %d forced (%d unique additions).",
    length(selected), length(splsda_sel), min_abs_loading,
    length(forced), length(setdiff(forced, splsda_sel))
  ))
  message(sprintf("[FeatureSelection] Panel: %s", paste(sort(selected), collapse = ", ")))

  list(selected = selected, n_splsda = length(splsda_sel), n_forced = length(forced))
}

#' @title Compute Universal Baseline Network
#' @description
#' Calculates Pcor, Bootstrap Stability, and generates the finalized adjacency matrix
#' for a single macro-group. This operation is executed only once per group.
#'
#' @param mat Numeric matrix (Samples x Features).
#' @param label String label for the group.
#' @param n_boot Number of bootstrap iterations.
#' @param seed Random seed.
#' @param n_cores Number of cores.
#' @param threshold_type "percentile" or "absolute".
#' @param threshold_value Threshold magnitude for valid edges.
#' @return List containing stabilized topological networks and matrices.
compute_universal_baseline <- function(mat, label = "Group", n_boot = 100, seed = 123, n_cores = 1, threshold_type = "percentile", threshold_value = 0.85) {
  requireNamespace("parallel", quietly = TRUE)
  requireNamespace("doParallel", quietly = TRUE)
  requireNamespace("foreach", quietly = TRUE)
  requireNamespace("corpcor", quietly = TRUE)
  
  message(sprintf("   [%s Baseline] Start: Boot=%d, Cores=%d", label, n_boot, n_cores))
  
  cl <- parallel::makeCluster(n_cores)
  on.exit({
    if(!is.null(cl)) try(parallel::stopCluster(cl), silent=TRUE)
  }, add = TRUE)
  
  doParallel::registerDoParallel(cl)
  parallel::clusterEvalQ(cl, { library(corpcor) })
  parallel::clusterExport(cl, varlist = c("infer_network_pcor", "boot_worker_pcor"), envir = environment())
  parallel::clusterSetRNGStream(cl, seed)
  
  # 1. BOOTSTRAP STABILITY
  n_samp <- nrow(mat)
  res_list <- foreach::foreach(i = 1:n_boot, .packages = "corpcor") %dopar% {
    boot_worker_pcor(mat, n_samp, lambda_val = NULL)
  }
  
  boots <- Filter(Negate(is.null), res_list)
  agg <- aggregate_boot_results(boots, alpha = 0.10)
  check_boot_yield(agg, n_boot, label)
  
  # 2. PCOR & MAGNITUDE THRESHOLDING
  obs_pcor <- infer_network_pcor(mat, fixed_lambda = NULL)
  raw_cor <- suppressWarnings(cor(mat, method = "spearman", use = "pairwise.complete.obs"))
  
  if (threshold_value == 0) {
    actual_thresh <- 0
    message(sprintf("   [%s Baseline] Threshold set to 0.", label))
  } else if (threshold_type == "percentile") {
    actual_thresh <- quantile(abs(obs_pcor[upper.tri(obs_pcor)]), probs = threshold_value, na.rm = TRUE)
    message(sprintf("   [%s Baseline] Percentile threshold (%.2f) -> |rho| >= %.4f", label, threshold_value, actual_thresh))
  } else {
    actual_thresh <- threshold_value
    message(sprintf("   [%s Baseline] Absolute threshold -> |rho| >= %.4f", label, actual_thresh))
  }
  
  # 3. UNIVERSAL TOPOLOGICAL DEFINITION
  adj_final <- (agg$adj == 1) & (abs(obs_pcor) >= actual_thresh)
  n_edges <- sum(adj_final[upper.tri(adj_final)])
  message(sprintf("   [%s Baseline] Valid edges (Stable + Magnitude): %d", label, n_edges))
  
  return(list(
    label = label,
    mat = mat, 
    pcor = obs_pcor,
    raw_cor = raw_cor,
    stability = agg$adj,
    stability_freq = agg$stability,
    adj_final = adj_final,
    applied_threshold = actual_thresh
  ))
}

#' @title Compute Differential Overlay (Stability-Driven with Spearman Annotation)
#' @description 
#' Performs Permutation Testing leveraging pre-calculated Universal Baselines.
#' Bypasses redundant bootstrapping. Implements a Stability-Driven filtering 
#' approach, decoupling edge existence (Bootstrap stability + Delta Pcor) 
#' from biological interpretation (Spearman Rank Annotation).
#' 
#' @param base_ctrl Pre-calculated baseline object for the control group.
#' @param base_case Pre-calculated baseline object for the case group.
#' @param n_perm Number of permutations.
#' @param seed Random seed.
#' @param n_cores Number of cores.
#' @param pvalue_thresh Alpha threshold (Retained for continuous annotation).
#' @param min_marginal_cor Used as the threshold for Delta Pcor and Spearman magnitude.
#' @return List with structured edge table and network topologies.
compute_differential_overlay <- function(base_ctrl, base_case, n_perm = 1000, seed = 123, n_cores = 1, pvalue_thresh = 0.05, min_marginal_cor = 0.10) {
  requireNamespace("parallel", quietly = TRUE)
  requireNamespace("doParallel", quietly = TRUE)
  requireNamespace("foreach", quietly = TRUE)
  requireNamespace("corpcor", quietly = TRUE)
  
  label_ctrl <- base_ctrl$label
  label_case <- base_case$label
  mat_ctrl <- base_ctrl$mat
  mat_case <- base_case$mat
  
  message(sprintf("      [DiffNet] Overlaying %s vs %s (Perm=%d, Cores=%d)", label_ctrl, label_case, n_perm, n_cores))
  
  stable_mask <- base_ctrl$adj_final | base_case$adj_final
  n_stable_edges <- sum(stable_mask[upper.tri(stable_mask)])
  
  if (n_stable_edges == 0) {
    message("      [DiffNet] No valid edges found in either baseline. Returning empty results.")
    return(list(
      edges_table = data.frame(),
      networks = list(ctrl = base_ctrl$pcor, case = base_case$pcor),
      adj_final = list(ctrl = base_ctrl$adj_final, case = base_case$adj_final)
    ))
  }
  
  cl <- parallel::makeCluster(n_cores)
  on.exit({
    if(!is.null(cl)) try(parallel::stopCluster(cl), silent=TRUE)
  }, add = TRUE)
  doParallel::registerDoParallel(cl)
  parallel::clusterEvalQ(cl, { library(corpcor) })
  parallel::clusterExport(cl, varlist = c("infer_network_pcor"), envir = environment())
  parallel::clusterSetRNGStream(cl, seed)
  
  # 1. PERMUTATION TEST (Generates Empirical Null Distribution)
  pool_mat <- rbind(mat_ctrl, mat_case)
  n1 <- nrow(mat_ctrl)
  n_total <- nrow(pool_mat)
  
  null_diffs_raw <- foreach::foreach(i = 1:n_perm, .packages = "corpcor") %dopar% {
    shuffled_idx <- sample(1:n_total)
    p1 <- pool_mat[shuffled_idx[1:n1], , drop = FALSE]
    p2 <- pool_mat[shuffled_idx[(n1 + 1):n_total], , drop = FALSE]
    
    var1 <- apply(p1, 2, var, na.rm = TRUE)
    var2 <- apply(p2, 2, var, na.rm = TRUE)
    
    if (any(var1 == 0 | is.na(var1)) || any(var2 == 0 | is.na(var2))) {
      return(NULL)
    }
    
    tryCatch({
      r1 <- infer_network_pcor(p1, fixed_lambda = NULL)
      r2 <- infer_network_pcor(p2, fixed_lambda = NULL)
      return(abs(r1 - r2))
    }, error = function(e) return(NULL))
  }
  
  null_diffs <- Filter(Negate(is.null), null_diffs_raw)
  n_valid_perms <- length(null_diffs)
  
  if (n_valid_perms < (n_perm * 0.8)) {
    warning(sprintf("      [DiffNet] High invariant resample rate. %d/%d valid permutations used.", n_valid_perms, n_perm))
  }
  if (n_valid_perms == 0) {
    stop("      [DiffNet] All permutations failed due to zero variance.")
  }
  
  # 2. STATS & RAW DATA AGGREGATION
  obs_diff <- abs(base_ctrl$pcor - base_case$pcor)
  nodes <- colnames(mat_ctrl)
  edges_idx <- which(upper.tri(stable_mask) & stable_mask, arr.ind = TRUE)
  
  results_list <- list()
  denom <- n_valid_perms + 1
  
  for(k in 1:nrow(edges_idx)) {
    i <- edges_idx[k,1]
    j <- edges_idx[k,2]
    
    obs_val <- obs_diff[i,j]
    null_vals <- sapply(null_diffs, function(m) m[i, j])
    p_val <- (sum(null_vals >= obs_val) + 1) / denom
    
    is_valid_ctrl <- base_ctrl$adj_final[i,j]
    is_valid_case <- base_case$adj_final[i,j]
    category <- "Weak"
    
    if (is_valid_ctrl || is_valid_case) {
      if (is_valid_ctrl && is_valid_case) {
        same_sign <- (sign(base_ctrl$pcor[i,j]) == sign(base_case$pcor[i,j]))
        category <- ifelse(same_sign, "Conserved", "Inverted")
      } else {
        category <- "Specific"
      }
    }
    
    results_list[[k]] <- data.frame(
      Node1 = nodes[i],
      Node2 = nodes[j],
      Cor_Ctrl = base_ctrl$raw_cor[i,j], # Now strictly Spearman from baseline
      Weight_Ctrl = base_ctrl$pcor[i,j],
      Cor_Case = base_case$raw_cor[i,j], # Now strictly Spearman from baseline
      Weight_Case = base_case$pcor[i,j],
      Is_Stable_Ctrl = base_ctrl$stability[i,j] == 1,
      Is_Stable_Case = base_case$stability[i,j] == 1,
      StabFreq_Ctrl = base_ctrl$stability_freq[i,j], 
      StabFreq_Case = base_case$stability_freq[i,j],
      Is_Valid_Ctrl = is_valid_ctrl,
      Is_Valid_Case = is_valid_case,
      Diff_Score = obs_val,
      P_Value = p_val,
      Edge_Category = category,
      stringsAsFactors = FALSE
    )
  }
  
  res_df <- do.call(rbind, results_list)
  res_df$FDR <- stats::p.adjust(res_df$P_Value, method = "BH")
  
  # 3. STABILITY-DRIVEN FILTERING & SPEARMAN ANNOTATION PIPELINE
  res_df <- res_df %>%
    dplyr::mutate(
      # Mechanism Annotation: Compare Pcor with Spearman (Rank)
      Concordance_Ctrl = (sign(Weight_Ctrl) == sign(Cor_Ctrl)) & (abs(Cor_Ctrl) >= min_marginal_cor),
      Concordance_Case = (sign(Weight_Case) == sign(Cor_Case)) & (abs(Cor_Case) >= min_marginal_cor),
      
      Mechanism_Ctrl = dplyr::case_when(
        Is_Valid_Ctrl & Concordance_Ctrl ~ "Direct",
        Is_Valid_Ctrl & !Concordance_Ctrl ~ "Masked/Suppressor",
        TRUE ~ "None"
      ),
      
      Mechanism_Case = dplyr::case_when(
        Is_Valid_Case & Concordance_Case ~ "Direct",
        Is_Valid_Case & !Concordance_Case ~ "Masked/Suppressor",
        TRUE ~ "None"
      ),
      
      # Core Filtering Rule: Existence relies on Topological Stability AND Delta Magnitude.
      # Bypassing P-value/FDR constraints for structural existence to protect against small N power drop.
      Significant = (Is_Valid_Ctrl | Is_Valid_Case) & (Diff_Score >= min_marginal_cor),
      
      # Overwrite Edge_Category with precise dynamic labels for specific edges
      Edge_Category = dplyr::case_when(
        Edge_Category == "Specific" & Is_Valid_Ctrl ~ paste0("Specific_", label_ctrl),
        Edge_Category == "Specific" & Is_Valid_Case ~ paste0("Specific_", label_case),
        TRUE ~ Edge_Category
      )
    )
  
  # Order the output for readability
  dynamic_levels <- c(
    "Inverted", 
    paste0("Specific_", label_case), 
    paste0("Specific_", label_ctrl), 
    "Conserved", 
    "Weak"
  )
  
  res_df <- res_df %>%
    dplyr::arrange(factor(Edge_Category, levels = dynamic_levels), desc(Diff_Score))
  
  # Rename columns dynamically based on labels
  colnames(res_df)[colnames(res_df) == "Weight_Ctrl"] <- paste0("Pcor_", label_ctrl)
  colnames(res_df)[colnames(res_df) == "Weight_Case"] <- paste0("Pcor_", label_case)
  colnames(res_df)[colnames(res_df) == "Cor_Ctrl"] <- paste0("Spearman_", label_ctrl)
  colnames(res_df)[colnames(res_df) == "Cor_Case"] <- paste0("Spearman_", label_case)
  colnames(res_df)[colnames(res_df) == "Is_Stable_Ctrl"] <- paste0("Is_Stable_", label_ctrl)
  colnames(res_df)[colnames(res_df) == "Is_Stable_Case"] <- paste0("Is_Stable_", label_case)
  colnames(res_df)[colnames(res_df) == "Is_Valid_Ctrl"] <- paste0("Is_Valid_", label_ctrl)
  colnames(res_df)[colnames(res_df) == "Is_Valid_Case"] <- paste0("Is_Valid_", label_case)
  colnames(res_df)[colnames(res_df) == "Mechanism_Ctrl"] <- paste0("Mech_", label_ctrl)
  colnames(res_df)[colnames(res_df) == "Mechanism_Case"] <- paste0("Mech_", label_case)
  
  return(list(
    edges_table = res_df,
    networks = list(ctrl = base_ctrl$pcor, case = base_case$pcor),
    raw_cor = list(ctrl = base_ctrl$raw_cor, case = base_case$raw_cor),
    adj_final = list(ctrl = base_ctrl$adj_final, case = base_case$adj_final),
    stability = list(ctrl = base_ctrl$stability, case = base_case$stability),
    applied_threshold = base_ctrl$applied_threshold
  ))
}

#' @title Calculate Node Topology Metrics
#' @description 
#' Calculates centralities and topology metrics for a given adjacency matrix.
#' Returns a clean node-level dataframe.
#' 
#' @param adj_mat Adjacency matrix (0/1).
#' @return Dataframe with node metrics.
calculate_node_topology <- function(adj_mat) {
  
  requireNamespace("igraph", quietly = TRUE)
  
  if (sum(adj_mat) == 0) return(NULL)
  
  # Build Graph
  g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
  
  # Calculate Metrics
  metrics <- data.frame(
    Node = igraph::V(g)$name,
    Degree = igraph::degree(g),
    Betweenness = igraph::betweenness(g, normalized = TRUE),
    Closeness = igraph::closeness(g, normalized = TRUE),
    Eigen_Centrality = igraph::eigen_centrality(g)$vector,
    Cluster_Coeff = igraph::transitivity(g, type = "local"),
    stringsAsFactors = FALSE
  )
  
  # Order by Degree (Hubs first)
  metrics <- metrics[order(metrics$Degree, decreasing = TRUE), ]
  
  return(metrics)
}

#' @title Filter Baseline Network for Cytoscape Export (Backbone Strategy)
#' @description
#' Filters the raw baseline edges using stability and partial correlation thresholds.
#' Applies a maximum edge limit (Rank-Based) and returns both the main exported 
#' network and the remaining valid edges that exceeded the limit.
#'
#' @param edges_df Dataframe containing network edges (must include StabFreq and Pcor).
#' @param baseline_config List containing the baseline export configuration.
#' @return A list containing 'main' and 'remaining' dataframes.
filter_baseline_backbone <- function(edges_df, baseline_config) {
  
  if (is.null(edges_df) || nrow(edges_df) == 0) return(list(main = edges_df, remaining = NULL))
  
  min_stab <- if(!is.null(baseline_config$min_stability_freq)) baseline_config$min_stability_freq else 0.95
  min_pcor <- if(!is.null(baseline_config$min_abs_pcor)) baseline_config$min_abs_pcor else 0.10
  max_edges <- if(!is.null(baseline_config$max_export_edges)) baseline_config$max_export_edges else NULL
  
  df_filtered <- edges_df %>%
    dplyr::filter(StabFreq >= min_stab, abs(Pcor) >= min_pcor) %>%
    dplyr::arrange(dplyr::desc(abs(Pcor)))
  
  df_main <- df_filtered
  df_remaining <- NULL
  
  if (!is.null(max_edges) && nrow(df_filtered) > max_edges) {
    df_main <- utils::head(df_filtered, max_edges)
    df_remaining <- utils::tail(df_filtered, nrow(df_filtered) - max_edges)
  }
  
  return(list(main = df_main, remaining = df_remaining))
}

#' @title Filter Differential Network for Cytoscape Export
#' @description
#' Applies statistical and stability thresholds for differential network export.
#' Returns both the top N main edges and the remaining significant edges.
#'
#' @param edges_df Dataframe containing differential network edges.
#' @param diff_config List containing the differential export configuration.
#' @return A list containing 'main' and 'remaining' dataframes.
filter_differential_export <- function(edges_df, diff_config) {
  
  if (is.null(edges_df) || nrow(edges_df) == 0) return(list(main = edges_df, remaining = NULL))
  
  min_diff <- if(!is.null(diff_config$min_diff_score)) diff_config$min_diff_score else 0.15
  min_stab <- if(!is.null(diff_config$min_stability_freq_any)) diff_config$min_stability_freq_any else 0.95
  excl_cons <- if(!is.null(diff_config$exclude_conserved)) diff_config$exclude_conserved else TRUE
  max_edges <- if(!is.null(diff_config$max_export_edges)) diff_config$max_export_edges else NULL
  
  df_filtered <- edges_df %>%
    dplyr::filter(
      dplyr::if_any(dplyr::starts_with("Is_Valid_"), ~ .x == TRUE),
      Diff_Score >= min_diff,
      dplyr::if_any(dplyr::starts_with("StabFreq_"), ~ .x >= min_stab)
    )
  
  if (excl_cons) {
    df_filtered <- df_filtered %>% dplyr::filter(!Edge_Category %in% c("Conserved", "Weak"))
  } else {
    df_filtered <- df_filtered %>% dplyr::filter(Edge_Category != "Weak")
  }
  
  df_filtered <- df_filtered %>% dplyr::arrange(dplyr::desc(Diff_Score))
  
  df_main <- df_filtered
  df_remaining <- NULL
  
  if (!is.null(max_edges) && nrow(df_filtered) > max_edges) {
    df_main <- utils::head(df_filtered, max_edges)
    df_remaining <- utils::tail(df_filtered, nrow(df_filtered) - max_edges)
  }
  
  return(list(main = df_main, remaining = df_remaining))
}

#' @title Calculate Node Jaccard Rewiring
#' @description 
#' Quantifies how much a node's neighborhood changes between two networks.
#' Jaccard = Intersection(Neighbors) / Union(Neighbors).
#' 
#' @param adj_ctrl Adjacency matrix (0/1) for Control.
#' @param adj_case Adjacency matrix (0/1) for Case.
#' @param thresholds Vector of thresholds (lower bound, upper bound)
#' @return Dataframe with Jaccard index per node.
calculate_jaccard_rewiring <- function(adj_ctrl, adj_case, thresholds = c(0.3, 0.7)) {
  
  # Ensure node alignment
  nodes <- intersect(rownames(adj_ctrl), rownames(adj_case))
  if (length(nodes) == 0) return(NULL)
  
  res_df <- data.frame(
    Node = nodes,
    Degree_Ctrl = 0,
    Degree_Case = 0,
    Jaccard_Index = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(nodes)) {
    n <- nodes[i]
    
    neigh_ctrl <- names(which(adj_ctrl[n, ] != 0))
    neigh_case <- names(which(adj_case[n, ] != 0))
    
    # Sets
    intersect_set <- intersect(neigh_ctrl, neigh_case)
    union_set <- union(neigh_ctrl, neigh_case)
    
    # Store Degrees
    res_df$Degree_Ctrl[i] <- length(neigh_ctrl)
    res_df$Degree_Case[i] <- length(neigh_case)
    
    # Calculate Jaccard
    if (length(union_set) > 0) {
      res_df$Jaccard_Index[i] <- length(intersect_set) / length(union_set)
    } else {
      # Both isolated in both networks: no rewiring occurred
      res_df$Jaccard_Index[i] <- 1.0
    }
  }
  
  # Add interpretation
  res_df$Rewiring_Status <- cut(res_df$Jaccard_Index, 
                                breaks = c(-Inf, thresholds[1], thresholds[2], Inf), 
                                labels = c("High_Rewiring", "Moderate", "Stable"))
  
  return(res_df)
}

#' @title Integrate Drivers and Topology (Hub-Driver Analysis)
#' @description Merges PLS-DA drivers with Network Topology metrics to identify roles.
#' @param drivers_df Dataframe output from extract_plsda_loadings.
#' @param topo_df Dataframe output from calculate_node_topology.
#' @param topo_col String. Name of the topology metric to use for Y-axis classification (default "Degree").
#' @return Merged dataframe with "Role" column.
integrate_hub_drivers <- function(drivers_df, topo_df, topo_col = "Degree") {
  
  if (is.null(drivers_df) || nrow(drivers_df) == 0) return(NULL)
  if (is.null(topo_df) || nrow(topo_df) == 0) return(NULL)
  
  # Clean column names to ensure matching
  # Usually extract_plsda_loadings returns "Marker" and topo returns "Node"
  
  # Select only the highest importance per marker (if duplicates exist)
  df_drv <- drivers_df %>%
    dplyr::group_by(Marker) %>%
    dplyr::slice_max(order_by = Importance, n = 1) %>%
    dplyr::ungroup()
  
  # Merge
  merged <- merge(df_drv, topo_df, by.x = "Marker", by.y = "Node")
  
  if (nrow(merged) == 0) return(NULL)
  
  # Validation: check if topo_col exists in the merged dataframe
  if (!topo_col %in% names(merged)) {
    warning(sprintf("[Network] Topology column '%s' not found. Reverting to 'Degree'.", topo_col))
    topo_col <- "Degree"
  }
  
  # Calculate Quartiles/Medians for Categorization based on the SPECIFIC metric
  imp_med <- median(merged$Importance, na.rm = TRUE)
  topo_med <- median(merged[[topo_col]], na.rm = TRUE) # Dynamic Median calculation
  
  # Assign Roles dynamically
  merged$Role <- dplyr::case_when(
    merged$Importance >= imp_med & merged[[topo_col]] >= topo_med ~ "Master_Regulator",
    merged$Importance >= imp_med & merged[[topo_col]] < topo_med  ~ "Solo_Driver",
    merged$Importance < imp_med & merged[[topo_col]] >= topo_med  ~ "Structural_Connector",
    TRUE ~ "Background"
  )
  
  # Store the specific metric used in a generic column for easier plotting later
  merged$Topology_Metric_Value <- merged[[topo_col]] 
  
  return(merged)
}

#' @title Calculate Topology and Export Rich Cytoscape Table
#' @description 
#' Calculates node topology metrics and merges them into the edge list.
#' Creates a single "Rich Table" for Cytoscape (One-click import).
#' 
#' @param adj_mat Adjacency matrix (0/1).
#' @param weight_mat Partial correlation matrix.
#' @param group_label Label for the group (e.g., "Healthy", "LS").
#' @return A dataframe compatible with Cytoscape import.
get_rich_network_table <- function(adj_mat, weight_mat, group_label = "Unknown") {
  
  requireNamespace("igraph", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  
  # 1. Check if network is empty
  if (sum(adj_mat) == 0) return(NULL)
  
  # 2. Build Graph & Calculate Topology
  g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
  
  # Calculate Metrics
  node_metrics <- data.frame(
    Node = igraph::V(g)$name,
    Degree = igraph::degree(g),
    Betweenness = igraph::betweenness(g, normalized = TRUE),
    Closeness = igraph::closeness(g, normalized = TRUE),
    Eigen_Centrality = igraph::eigen_centrality(g)$vector,
    stringsAsFactors = FALSE
  )
  
  # 3. Build Edge List
  edge_list <- igraph::as_data_frame(g, what = "edges")
  
  # 4. Merge Logic (The "Rich" Table)
  # We join node metrics TWICE: once for Source, once for Target
  rich_table <- edge_list %>%
    rowwise() %>%
    mutate(
      Weight = weight_mat[from, to],
      Abs_Weight = abs(Weight),
      Sign = ifelse(Weight > 0, "Positive", "Negative"),
      Interaction = ifelse(Weight > 0, "Co-occurrence", "Mutual-Exclusion"),
      Group_ID = group_label
    ) %>%
    ungroup() %>%
    # Join Source Metrics
    left_join(node_metrics, by = c("from" = "Node")) %>%
    dplyr::rename(
      Source = from, 
      Source_Degree = Degree, 
      Source_Betweenness = Betweenness,
      Source_Closeness = Closeness
    ) %>%
    # Join Target Metrics
    left_join(node_metrics, by = c("to" = "Node")) %>%
    dplyr::rename(
      Target = to, 
      Target_Degree = Degree, 
      Target_Betweenness = Betweenness,
      Target_Closeness = Closeness
    ) %>%
    # Cleanup auxiliary columns from join if any
    select(-matches("Eigen_Centrality")) # Optional: keep or remove based on preference
  
  return(rich_table)
}

#' @title Integrate Global Hub and Driver Roles
#' @description Merges PLS-DA drivers with comprehensive Network Topology metrics.
#' @param drivers_df Dataframe output from extract_plsda_loadings.
#' @param topo_df Dataframe output from calculate_node_topology.
#' @return A unified dataframe with topological metrics and overall Role assignment.
integrate_global_hub_roles <- function(drivers_df, topo_df) {
  
  if (is.null(drivers_df) || nrow(drivers_df) == 0) return(NULL)
  if (is.null(topo_df) || nrow(topo_df) == 0) return(NULL)
  
  # Select highest importance per marker across components
  df_drv <- drivers_df %>%
    dplyr::group_by(Marker) %>%
    dplyr::slice_max(order_by = Importance, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  # Merge stats with topology
  merged <- merge(df_drv, topo_df, by.x = "Marker", by.y = "Node")
  if (nrow(merged) == 0) return(NULL)
  
  # Calculate Medians for generic classification
  imp_med <- median(merged$Importance, na.rm = TRUE)
  deg_med <- median(merged$Degree, na.rm = TRUE)
  bet_med <- median(merged$Betweenness, na.rm = TRUE)
  
  # Define a holistic biological role: 
  # Master_Regulator if it drives PLS-DA AND has high network centrality
  merged$Overall_Role <- dplyr::case_when(
    merged$Importance >= imp_med & (merged$Degree >= deg_med | merged$Betweenness >= bet_med) ~ "Master_Regulator",
    merged$Importance >= imp_med & merged$Degree < deg_med & merged$Betweenness < bet_med     ~ "Solo_Driver",
    merged$Importance < imp_med & (merged$Degree >= deg_med | merged$Betweenness >= bet_med)  ~ "Structural_Connector",
    TRUE ~ "Background"
  )
  
  # Order nicely for Excel
  merged <- merged %>% dplyr::arrange(dplyr::desc(Importance))

  return(merged)
}

#' @title Annotate Node Topology with Marker Domain
#' @description Labels nodes as FACS or Soluble based on config feature lists.
#'              Falls back to "Other" for unclassified nodes.
#' @param node_df Dataframe containing at least a 'Node' column.
#' @param config Global configuration list.
#' @return Dataframe with an appended 'Marker_Category' column.
annotate_marker_domain <- function(node_df, config) {
  if (is.null(node_df) || !("Node" %in% names(node_df))) return(node_df)
  if ("Marker_Category" %in% names(node_df)) node_df$Marker_Category <- NULL

  facs_feats <- if (!is.null(config$features$facs)) as.character(unlist(config$features$facs)) else character(0)
  sol_feats  <- if (!is.null(config$features$soluble)) as.character(unlist(config$features$soluble)) else character(0)

  node_df$Marker_Category <- dplyr::case_when(
    node_df$Node %in% facs_feats ~ "FACS",
    node_df$Node %in% sol_feats  ~ "Soluble",
    TRUE                          ~ "Other"
  )

  return(node_df)
}