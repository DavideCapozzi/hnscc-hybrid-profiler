# src/02_visualization.R
# ==============================================================================
# STEP 02: VISUALIZATION
# Description: Generates overview plots (Distributions, PCA, Heatmap)
# ==============================================================================

source("R/utils_io.R")          
source("R/modules_viz.R")

message("\n=== PIPELINE STEP 2: VISUALIZATION ===")

# 1. Load Config & Data
args <- commandArgs(trailingOnly = TRUE)
# 1. Load Config & Data
if (!exists("config")) {
  args <- commandArgs(trailingOnly = TRUE)
  config_path <- if (length(args) > 0) args[1] else "config/global_params.yml"
  config <- load_config(config_path)
}
input_file <- file.path(config$output_root, "01_data_processing", "data_processed.rds")

DATA <- readRDS(input_file)
meta_viz <- DATA$metadata
mat_z_global <- as.matrix(DATA$hybrid_data_z[, DATA$hybrid_markers])
raw_matrix <- DATA$raw_matrix

# 2. Colors Setup 
# ------------------------------------------------------------------------------
unique_groups <- levels(meta_viz$Group)

colors_viz <- c()
if (!is.null(config$colors$groups[[config$clinical$responder_label]])) {
  colors_viz[config$clinical$responder_label] <- config$colors$groups[[config$clinical$responder_label]]
} else { colors_viz[config$clinical$responder_label] <- "blue" }

if (!is.null(config$colors$groups[[config$clinical$non_responder_label]])) {
  colors_viz[config$clinical$non_responder_label] <- config$colors$groups[[config$clinical$non_responder_label]]
} else { colors_viz[config$clinical$non_responder_label] <- "red" }

out_dir <- file.path(config$output_root, "02_visualization")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 3. Raw Distributions
# ------------------------------------------------------------------------------
message("[Viz] Generating distribution plots...")
df_raw_viz <- cbind(meta_viz, as.data.frame(raw_matrix))

viz_save_distribution_report(
  data_df = df_raw_viz, 
  markers = DATA$markers, 
  file_path = file.path(out_dir, "Distributions_Raw_Hybrid.pdf"), 
  colors = colors_viz 
)

# 4. PCA Analysis
# ------------------------------------------------------------------------------
message("[PCA] Running Global PCA...")
res_pca <- FactoMineR::PCA(mat_z_global, scale.unit = FALSE, graph = FALSE)

pdf(file.path(out_dir, "PCA_Global_Individuals.pdf"), width = 11, height = 7) 

target_dims_list <- list(c(1,2), c(1,3), c(2,3))

for (dims in target_dims_list) {
  print(plot_pca_custom(res_pca, meta_viz, colors_viz, dims = dims, show_labels = TRUE))
}

message("   [PCA] Plotting Variance Scree Plot")
tryCatch({ print(plot_pca_variance_dashboard(res_pca, n_list = 10)) }, error = function(e) warning(e$message))

for (dims in target_dims_list) {
  p_var <- factoextra::fviz_pca_var(res_pca, axes = dims, col.var = "contrib", 
                                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) +
    ggplot2::labs(title = sprintf("Marker Contribution (PC%d vs PC%d)", dims[1], dims[2])) +
    theme_coda()
  print(p_var)
}
dev.off()

# 5. Heatmap Analysis
# ------------------------------------------------------------------------------
message("[Viz] Generating Heatmap...")
annotation_colors <- list(Group = colors_viz)

pdf(file.path(out_dir, "Heatmap_Stratification.pdf"), width = 10, height = 8)
tryCatch({
  ht_obj <- plot_stratification_heatmap(
    mat_z = mat_z_global,
    metadata = meta_viz,
    annotation_colors_list = annotation_colors,
    title = "Hybrid Clustering Analysis"
  )
  ComplexHeatmap::draw(ht_obj, merge_legend = TRUE)
}, error = function(e) warning(paste("Heatmap failed:", e$message)))
dev.off()

message("=== STEP 2 COMPLETE ===\n")
