# src/05_network_analysis.R
# ==============================================================================
# STEP 05: DIFFERENTIAL NETWORK ANALYSIS
# Description: Computes per-group partial correlation baselines (bootstrap) and
#              a differential overlay (permutation test) between the two clinical
#              groups. Exports Cytoscape-ready CSVs and PDF visualizations.
# Note: Runs on cross-sectional data only (standard pass). Longitudinal data is
#       not passed here — the T0/T1 mixed matrix would confound static networks.
# ==============================================================================

source("R/utils_io.R")
source("R/modules_network.R")
source("R/modules_viz.R")

message("\n=== PIPELINE STEP 5: DIFFERENTIAL NETWORK ANALYSIS ===")

# 1. Configuration & Data Loading
# ------------------------------------------------------------------------------
if (!exists("config")) stop("[FATAL] Global configuration object not detected in environment.")

input_rds <- file.path(config$output_root, "01_data_processing",
                       sprintf("data_processed_%s_standard.rds", config$project_name))
if (!file.exists(input_rds)) stop(sprintf("[FATAL] Step 01 standard output not found at: %s", input_rds))

DATA <- readRDS(input_rds)
message(sprintf("[Data] Loaded Step 01 output: %d samples x %d markers",
                nrow(DATA$metadata), length(DATA$hybrid_markers)))

out_dir <- file.path(config$output_root, "05_network_analysis")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 2. Group Setup
# ------------------------------------------------------------------------------
ctrl_label <- config$clinical$non_responder_label
case_label <- config$clinical$responder_label

message(sprintf("[Network] Clinical groups: Control = '%s' | Case = '%s'", ctrl_label, case_label))

meta <- DATA$metadata
mat_z_full <- as.matrix(DATA$hybrid_data_z[, DATA$hybrid_markers, drop = FALSE])
rownames(mat_z_full) <- meta$Sample_ID

ctrl_idx <- which(meta$Group == ctrl_label)
case_idx <- which(meta$Group == case_label)

if (length(ctrl_idx) == 0 || length(case_idx) == 0) {
  stop(sprintf("[FATAL] One or both clinical groups not found in metadata. Check config labels."))
}

mat_ctrl <- mat_z_full[ctrl_idx, , drop = FALSE]
mat_case <- mat_z_full[case_idx, , drop = FALSE]

message(sprintf("[Network] Group dimensions: %s = %d samples | %s = %d samples",
                ctrl_label, nrow(mat_ctrl), case_label, nrow(mat_case)))

# Safety: check minimum sample sizes
min_n <- if (!is.null(config$network$min_network_n)) config$network$min_network_n else 3
if (nrow(mat_ctrl) < min_n || nrow(mat_case) < min_n) {
  stop(sprintf("[FATAL] Insufficient samples for network analysis (min_network_n = %d).", min_n))
}

# Drop near-zero-variance features (per-group check — pooled filter)
var_thresh <- 1e-6
pool_vars <- apply(rbind(mat_ctrl, mat_case), 2, var, na.rm = TRUE)
valid_feats <- names(pool_vars)[!is.na(pool_vars) & pool_vars >= var_thresh]

if (length(valid_feats) < ncol(mat_z_full)) {
  message(sprintf("[Network] Dropping %d near-zero-variance features from network input.",
                  ncol(mat_z_full) - length(valid_feats)))
}
if (length(valid_feats) < 2) stop("[FATAL] Fewer than 2 features with sufficient variance for network analysis.")

mat_ctrl <- mat_ctrl[, valid_feats, drop = FALSE]
mat_case <- mat_case[, valid_feats, drop = FALSE]

# 3. Network Parameters
# ------------------------------------------------------------------------------
seed_val       <- if (!is.null(config$stats$seed)) config$stats$seed else 2026
n_boot         <- if (!is.null(config$network$n_boot)) config$network$n_boot else 500
n_perm         <- if (!is.null(config$network$n_perm)) config$network$n_perm else 1000
thresh_type    <- if (!is.null(config$network$threshold_type)) config$network$threshold_type else "absolute"
thresh_val     <- if (!is.null(config$network$edge_threshold)) config$network$edge_threshold else 0.10
min_marg_cor   <- if (!is.null(config$network$min_marginal_cor)) config$network$min_marginal_cor else 0.10

auto_cores <- max(1, parallel::detectCores(logical = TRUE) - 1, na.rm = TRUE)
n_cores <- if (!is.null(config$network$n_cores) && config$network$n_cores != "auto") {
  as.numeric(config$network$n_cores)
} else {
  auto_cores
}

message(sprintf("[Network] Parameters: n_boot=%d, n_perm=%d, threshold=%s(%.3f), cores=%d",
                n_boot, n_perm, thresh_type, thresh_val, n_cores))

# Helper: re-pad partial correlation matrices back to full feature space
# (needed if one group has fewer valid features than the other)
pad_matrix <- function(small_mat, full_names, default_val = 0) {
  p <- length(full_names)
  big_mat <- matrix(default_val, nrow = p, ncol = p, dimnames = list(full_names, full_names))
  valid_n <- intersect(rownames(small_mat), full_names)
  if (length(valid_n) > 0) big_mat[valid_n, valid_n] <- small_mat[valid_n, valid_n]
  return(big_mat)
}

# 4. Universal Baselines
# ------------------------------------------------------------------------------
message("\n--- PHASE 1: COMPUTING UNIVERSAL BASELINES ---")

compute_group_baseline <- function(mat, label) {
  base_raw <- compute_universal_baseline(
    mat = mat, label = label,
    n_boot = n_boot, seed = seed_val, n_cores = n_cores,
    threshold_type = thresh_type, threshold_value = thresh_val
  )

  all_feats <- valid_feats
  list(
    label          = base_raw$label,
    mat            = mat,
    pcor           = pad_matrix(base_raw$pcor, all_feats, 0),
    raw_cor        = pad_matrix(base_raw$raw_cor, all_feats, 0),
    stability      = pad_matrix(base_raw$stability, all_feats, 0),
    stability_freq = pad_matrix(base_raw$stability_freq, all_feats, 0),
    adj_final      = pad_matrix(base_raw$adj_final, all_feats, FALSE),
    applied_threshold = base_raw$applied_threshold
  )
}

base_ctrl <- compute_group_baseline(mat_ctrl, ctrl_label)
base_case <- compute_group_baseline(mat_case, case_label)

# 5. Cytoscape Baseline Export
# ------------------------------------------------------------------------------
cyto_base_dir <- file.path(out_dir, "cytoscape_baselines")
if (!dir.exists(cyto_base_dir)) dir.create(cyto_base_dir, recursive = TRUE)

b_conf <- config$network$cytoscape_export
if (is.null(b_conf)) b_conf <- list(min_stability_freq = 0.90, min_abs_pcor = 0.10, max_export_edges = 40)

export_baseline_cytoscape <- function(base_obj, conf, target_dir) {
  pcor_mat  <- base_obj$pcor
  stab_mat  <- base_obj$stability_freq
  nodes     <- colnames(pcor_mat)
  edges_idx <- which(upper.tri(pcor_mat), arr.ind = TRUE)
  if (nrow(edges_idx) == 0) return(invisible(NULL))

  base_edges_all <- data.frame(
    Source   = nodes[edges_idx[, 1]],
    Target   = nodes[edges_idx[, 2]],
    Pcor     = pcor_mat[edges_idx],
    StabFreq = stab_mat[edges_idx],
    stringsAsFactors = FALSE
  )

  res <- filter_baseline_backbone(base_edges_all, conf)

  write_edges <- function(edges_df, file_prefix) {
    if (is.null(edges_df) || nrow(edges_df) == 0) return(invisible(NULL))
    readr::write_csv(edges_df, file.path(target_dir, paste0(file_prefix, "_baseline_network.csv")))

    adj_exp <- matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
    for (k in seq_len(nrow(edges_df))) {
      adj_exp[edges_df$Source[k], edges_df$Target[k]] <- 1
      adj_exp[edges_df$Target[k], edges_df$Source[k]] <- 1
    }
    topo <- calculate_node_topology(adj_exp)
    if (!is.null(topo)) {
      topo <- annotate_marker_domain(topo, config)
      readr::write_csv(topo, file.path(target_dir, paste0(file_prefix, "_node_attributes.csv")))
    }
  }

  write_edges(res$main, base_obj$label)
  if (!is.null(res$remaining) && nrow(res$remaining) > 0) {
    rem_dir <- file.path(target_dir, "remaining_baselines")
    if (!dir.exists(rem_dir)) dir.create(rem_dir)
    write_edges(res$remaining, paste0(base_obj$label, "_remaining"))
  }
}

export_baseline_cytoscape(base_ctrl, b_conf, cyto_base_dir)
export_baseline_cytoscape(base_case, b_conf, cyto_base_dir)
message(sprintf("[Output] Baseline Cytoscape files saved to: %s", basename(cyto_base_dir)))

# 6. Differential Overlay
# ------------------------------------------------------------------------------
message("\n--- PHASE 2: DIFFERENTIAL NETWORK OVERLAY ---")

net_res <- compute_differential_overlay(
  base_ctrl        = base_ctrl,
  base_case        = base_case,
  n_perm           = n_perm,
  seed             = seed_val,
  n_cores          = n_cores,
  pvalue_thresh    = 0.05,
  min_marginal_cor = min_marg_cor
)

if (is.null(net_res) || nrow(net_res$edges_table) == 0) {
  message("[Network] No differential edges computed. Step 05 complete with empty output.")
} else {
  saveRDS(net_res, file.path(out_dir, sprintf("DifferentialNetwork_%s.rds", config$project_name)))

  # 7. Topology Metrics
  topo_ctrl <- if (sum(net_res$adj_final$ctrl) > 0) calculate_node_topology(net_res$adj_final$ctrl) else NULL
  topo_case <- if (sum(net_res$adj_final$case) > 0) calculate_node_topology(net_res$adj_final$case) else NULL

  # Write topology Excel
  wb_net <- openxlsx::createWorkbook()
  if (!is.null(topo_ctrl)) {
    topo_ctrl_ann <- annotate_marker_domain(topo_ctrl, config)
    openxlsx::addWorksheet(wb_net, safe_sheet_name(paste0("Topology_", ctrl_label)))
    openxlsx::writeData(wb_net, safe_sheet_name(paste0("Topology_", ctrl_label)), topo_ctrl_ann)
  }
  if (!is.null(topo_case)) {
    topo_case_ann <- annotate_marker_domain(topo_case, config)
    openxlsx::addWorksheet(wb_net, safe_sheet_name(paste0("Topology_", case_label)))
    openxlsx::writeData(wb_net, safe_sheet_name(paste0("Topology_", case_label)), topo_case_ann)
  }
  openxlsx::addWorksheet(wb_net, "Differential_Edges")
  openxlsx::writeData(wb_net, "Differential_Edges", net_res$edges_table)

  topo_path <- file.path(out_dir, sprintf("Network_Report_%s.xlsx", config$project_name))
  openxlsx::saveWorkbook(wb_net, topo_path, overwrite = TRUE)
  message(sprintf("[Output] Network topology and edges saved: %s", basename(topo_path)))

  # 8. Cytoscape Differential Export
  cyto_diff_dir <- file.path(out_dir, "cytoscape_differential")
  if (!dir.exists(cyto_diff_dir)) dir.create(cyto_diff_dir)

  d_conf <- config$network$cytoscape_export_diff
  if (is.null(d_conf)) d_conf <- list(min_diff_score = 0.10, min_stability_freq_any = 0.90,
                                      exclude_conserved = FALSE, max_export_edges = 40)

  diff_res <- filter_differential_export(net_res$edges_table, d_conf)

  export_diff_cyto <- function(edges_df, file_suffix, target_dir) {
    if (is.null(edges_df) || nrow(edges_df) == 0) return(invisible(NULL))
    if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)

    diff_net <- edges_df %>%
      dplyr::select(
        Source = Node1, Target = Node2, Weight = Diff_Score, Significant, P_Value, FDR,
        Edge_Category,
        dplyr::starts_with("Pcor_"), dplyr::starts_with("Spearman_"),
        dplyr::starts_with("Mech_"), dplyr::starts_with("Is_Stable_"), dplyr::starts_with("StabFreq_")
      ) %>%
      dplyr::rename(Interaction = Edge_Category)

    readr::write_csv(diff_net, file.path(target_dir,
                     sprintf("%s_%s_diff_network%s.csv", config$project_name, ctrl_label, file_suffix)))

    diff_nodes <- unique(c(edges_df$Node1, edges_df$Node2))
    adj_diff <- matrix(0, length(diff_nodes), length(diff_nodes), dimnames = list(diff_nodes, diff_nodes))
    for (k in seq_len(nrow(edges_df))) {
      adj_diff[edges_df$Node1[k], edges_df$Node2[k]] <- 1
      adj_diff[edges_df$Node2[k], edges_df$Node1[k]] <- 1
    }
    topo_diff <- calculate_node_topology(adj_diff)
    if (!is.null(topo_diff)) {
      topo_diff <- annotate_marker_domain(topo_diff, config)
      readr::write_csv(topo_diff, file.path(target_dir,
                       sprintf("%s_%s_diff_network%s_node_attributes.csv",
                               config$project_name, ctrl_label, file_suffix)))
    }
  }

  export_diff_cyto(diff_res$main, "", cyto_diff_dir)
  if (!is.null(diff_res$remaining) && nrow(diff_res$remaining) > 0) {
    export_diff_cyto(diff_res$remaining, "_remaining",
                     file.path(cyto_diff_dir, "remaining_differential"))
  }
  message(sprintf("[Output] Differential Cytoscape files saved to: %s", basename(cyto_diff_dir)))

  # 9. PDF Visualizations
  colors_viz <- get_clinical_colors(config)

  thresh_plot <- net_res$applied_threshold

  pdf(file.path(out_dir, sprintf("EdgeDensity_%s.pdf", config$project_name)), width = 8, height = 6)
  tryCatch({
    if (sum(abs(net_res$networks$ctrl), na.rm = TRUE) > 0)
      print(viz_plot_edge_density(net_res$networks$ctrl, adj_mat = net_res$adj_final$ctrl,
                                  threshold = thresh_plot, group_label = ctrl_label))
    if (sum(abs(net_res$networks$case), na.rm = TRUE) > 0)
      print(viz_plot_edge_density(net_res$networks$case, adj_mat = net_res$adj_final$case,
                                  threshold = thresh_plot, group_label = case_label))
  }, error = function(e) warning(paste("Edge density plot failed:", e$message)))
  dev.off()

  pdf(file.path(out_dir, sprintf("NetworkStructure_%s.pdf", config$project_name)), width = 12, height = 6)
  tryCatch({
    if (!is.null(topo_ctrl))
      print(plot_network_structure(net_res$adj_final$ctrl, net_res$networks$ctrl,
                                   title = paste("Control:", ctrl_label)))
    if (!is.null(topo_case))
      print(plot_network_structure(net_res$adj_final$case, net_res$networks$case,
                                   title = paste("Case:", case_label)))
  }, error = function(e) warning(paste("Network structure plot failed:", e$message)))
  dev.off()

  # 10. Hub-Driver Integration (requires sPLS-DA JSON from Step 03)
  json_path <- file.path(config$output_root, "03_statistical_analysis",
                         sprintf("Machine_Metrics_%s.json", config$project_name))

  if (file.exists(json_path) && requireNamespace("jsonlite", quietly = TRUE)) {
    tryCatch({
      metrics <- jsonlite::read_json(json_path, simplifyVector = TRUE)
      spls_drivers <- if (!is.null(metrics$top_drivers)) as.data.frame(metrics$top_drivers) else NULL

      if (!is.null(spls_drivers) && nrow(spls_drivers) > 0 &&
          "Marker" %in% names(spls_drivers) && "Importance" %in% names(spls_drivers)) {

        pdf(file.path(out_dir, sprintf("HubDriver_%s.pdf", config$project_name)), width = 11, height = 8)
        tryCatch({
          if (!is.null(topo_ctrl)) {
            hd_ctrl_deg <- integrate_hub_drivers(spls_drivers, topo_ctrl, "Degree")
            hd_ctrl_bet <- integrate_hub_drivers(spls_drivers, topo_ctrl, "Betweenness")
            if (!is.null(hd_ctrl_deg))
              print(plot_hub_driver_quadrant(hd_ctrl_deg, "Degree", paste("\nNetwork:", ctrl_label)))
            if (!is.null(hd_ctrl_bet))
              print(plot_hub_driver_quadrant(hd_ctrl_bet, "Betweenness", paste("\nNetwork:", ctrl_label)))
          }
          if (!is.null(topo_case)) {
            hd_case_deg <- integrate_hub_drivers(spls_drivers, topo_case, "Degree")
            hd_case_bet <- integrate_hub_drivers(spls_drivers, topo_case, "Betweenness")
            if (!is.null(hd_case_deg))
              print(plot_hub_driver_quadrant(hd_case_deg, "Degree", paste("\nNetwork:", case_label)))
            if (!is.null(hd_case_bet))
              print(plot_hub_driver_quadrant(hd_case_bet, "Betweenness", paste("\nNetwork:", case_label)))
          }
        }, error = function(e) warning(paste("Hub-driver plot failed:", e$message)))
        dev.off()
        message("[Output] Hub-Driver report exported.")
      }
    }, error = function(e) message(sprintf("[Network] Could not load sPLS-DA drivers: %s", e$message)))
  } else {
    message("[Network] Step 03 JSON not found — skipping Hub-Driver integration.")
  }
}

message("=== STEP 5 COMPLETE ===\n")
