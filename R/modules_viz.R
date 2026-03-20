# R/modules_viz.R
# ==============================================================================
# VISUALIZATION MODULE
# Description: Custom themes, PCA wrappers, and Plotting functions.
# Dependencies: ggplot2, FactoMineR, factoextra
# ==============================================================================

library(ggplot2)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(igraph)
library(tidygraph)
library(ggraph)
library(ComplexHeatmap)
library(circlize)

#' @title Standard Project Theme
#' @description A consistent ggplot2 theme for publication-quality figures.
#' @return A ggplot theme object.
theme_coda <- function() {
  theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "gray95")
    )
}

#' @title Generate Dynamic Group Palette (Root-Aware)
#' @description 
#' Centralized color logic with intelligent matching using a safe List-based approach.
#' Priority:
#' 1. Exact Match in Config (e.g., "NSCLC_LS" -> config "NSCLC_LS")
#' 2. Root Match (e.g., "NSCLC_EP" -> matches config "NSCLC")
#' 3. Control/Case Defaults
#' 
#' @param config The loaded configuration list.
#' @param match_groups Character vector of groups present in the data (optional). 
#' @return A named vector of hex codes.
get_palette <- function(config, match_groups = NULL) {
  
  # Initialize as LIST to allow safe NULL checks for missing keys
  # (Character vectors crash on [[missing_key]])
  final_palette <- list()
  
  # 1. Defaults setup (Safety unlist to handle YAML quirks)
  ctrl_col <- if(!is.null(config$colors$control)) unlist(config$colors$control) else "grey50"
  case_cols <- if(!is.null(config$colors$cases)) unlist(config$colors$cases) else c("#CD5C5C", "#4682B4")
  
  # Base lists from config
  defined_colors <- if(!is.null(config$colors$groups)) config$colors$groups else list()
  ctrl_groups_cfg <- unlist(config$control_group)
  case_groups_cfg <- unlist(config$case_groups)
  
  # Combine specific definitions into the palette first
  for (grp in names(defined_colors)) {
    final_palette[[grp]] <- defined_colors[[grp]]
  }
  
  # 2. Assign Generic Control/Case Colors for items in Config not yet colored
  if (!is.null(ctrl_groups_cfg)) {
    for(g in ctrl_groups_cfg) {
      if(is.null(final_palette[[g]])) final_palette[[g]] <- ctrl_col
    }
  }
  
  if (!is.null(case_groups_cfg)) {
    for (i in seq_along(case_groups_cfg)) {
      grp <- case_groups_cfg[i]
      if(is.null(final_palette[[grp]])) {
        # Cycle through available case colors
        col_idx <- (i - 1) %% length(case_cols) + 1
        final_palette[[grp]] <- case_cols[col_idx]
      }
    }
  }
  
  # 3. Data-Driven Extension (Root Matching)
  if (!is.null(match_groups)) {
    # Remove NAs and convert to character
    unique_data_groups <- unique(as.character(na.omit(match_groups)))
    
    for (g_data in unique_data_groups) {
      
      # Case A: Already has a color (Exact match)
      if (!is.null(final_palette[[g_data]])) next
      
      # Case B: Try Root Matching (Split by first underscore)
      # e.g., "NSCLC_EP" -> Root "NSCLC"
      root_name <- strsplit(g_data, "_")[[1]][1]
      
      # Safe lookup because final_palette is a list
      if (!is.null(final_palette[[root_name]])) {
        # Inherit color from Root
        final_palette[[g_data]] <- final_palette[[root_name]]
        
      } else {
        # Case C: Fallback
        if (g_data %in% ctrl_groups_cfg) {
          final_palette[[g_data]] <- ctrl_col
        } else {
          # Default fallback for unknown groups
          final_palette[[g_data]] <- case_cols[1] 
        }
      }
    }
  }
  
  # 4. Generic Fallback keys (required for legends sometimes)
  if (is.null(final_palette[["Case"]])) {
    final_palette[["Case"]] <- case_cols[1]
  }
  if (is.null(final_palette[["Control"]])) {
    final_palette[["Control"]] <- ctrl_col
  }
  
  # Convert back to named character vector for ggplot
  return(unlist(final_palette))
}

#' @title Plot Merged Raw Distribution with Highlights & Stats
#' @description 
#' Visualizes raw marker distribution.
#' Features: Violin (density), Jitter (points), Mean (Dashed Red), Median (Solid Black).
#' @param data_df Dataframe containing 'Group', metadata columns, and the marker.
#' @param marker_name Name of the marker column to plot.
#' @param colors Named vector of colors for the Groups.
#' @param highlight_pattern String pattern to search in 'Original_Source' (default: "_LS").
#' @param highlight_color Color for the border of highlighted points (default: "#FFD700").
#' @return A ggplot object.
plot_raw_distribution_merged <- function(data_df, marker_name, colors, 
                                         highlight_pattern = "_LS", 
                                         highlight_color = "#FFD700") {
  
  require(ggplot2)
  require(dplyr)
  
  # 1. Validation
  if (!marker_name %in% names(data_df)) return(NULL)
  
  # 2. Prepare Data & Highlighting Logic
  src_col <- names(data_df)[grepl("Original|Source", names(data_df), ignore.case = TRUE)][1]
  
  plot_data <- data_df
  plot_data$Value <- plot_data[[marker_name]]
  
  # Define Highlight Attributes defaults
  plot_data$Highlight_Type <- "Standard"
  plot_data$Stroke_Size <- 0.2
  
  if (!is.na(src_col)) {
    is_target <- grepl(highlight_pattern, plot_data[[src_col]])
    if (any(is_target)) {
      plot_data$Highlight_Type[is_target] <- "Target"
      plot_data$Stroke_Size[is_target] <- 1.5
    }
  }
  
  # 3. Calculate Stats for Subtitle
  n_na <- sum(is.na(plot_data$Value))
  pct_na <- round((n_na / nrow(plot_data)) * 100, 1)
  
  # 4. Explicitly remove NAs 
  plot_data_clean <- plot_data %>% filter(!is.na(Value))
  
  # 5. Plotting
  p <- ggplot(plot_data_clean, aes(x = Group, y = Value, fill = Group)) +
    
    # Layer 1: Violin
    geom_violin(alpha = 0.4, trim = FALSE, color = NA, scale = "width") +
    
    # Layer 2: Jittered Points
    geom_jitter(aes(color = Highlight_Type, stroke = Stroke_Size), 
                width = 0.2, size = 2.5, shape = 21, alpha = 0.8) +
    
    # Layer 3: Median (Solid Black Line)
    stat_summary(fun = median, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
                 width = 0.5, size = 0.8, color = "black") +
    
    # Layer 4: Mean (Dashed Red Line)
    stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)),
                 width = 0.5, size = 0.8, color = "darkred", linetype = "dashed") +
    
    # Colors
    scale_fill_manual(values = colors) +
    scale_color_manual(values = c("Standard" = "white", "Target" = highlight_color), 
                       guide = "none") +
    scale_continuous_identity(aesthetics = "stroke") +
    
    # Styling
    labs(
      title = paste("Raw Distribution:", marker_name),
      subtitle = sprintf("Stats: Black Line = Median | Red Dashed = Mean | Highlight: '%s'", 
                         highlight_pattern),
      y = "Raw Value",
      x = NULL
    ) +
    theme_coda() +
    theme(
      legend.position = "none", # Legend redundant as X axis labels exist
      axis.text.x = element_text(face = "bold", size = 11)
    )
  
  return(p)
}

#' @title Save Distribution Report (Multi-Page PDF)
#' @description 
#' Iterates through all markers and saves their raw distribution plots 
#' into a single PDF file. Encapsulates the looping logic.
#' 
#' @param data_df Merged dataframe with metadata and raw values.
#' @param markers Vector of marker names to plot.
#' @param file_path Output path for the PDF.
#' @param colors Named vector of colors for groups.
#' @param hl_pattern String pattern for highlighting (e.g. "_LS").
viz_save_distribution_report <- function(data_df, markers, file_path, colors, hl_pattern = "") {
  
  if (length(markers) == 0) {
    warning("[Viz] No markers provided for distribution report.")
    return(NULL)
  }
  
  message(sprintf("   [Viz] Saving distribution plots to: %s", basename(file_path)))
  
  pdf(file_path, width = 8, height = 6)
  # Ensure device is closed even if the loop crashes
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)
  
  for (marker in markers) {
    if (marker %in% names(data_df)) {
      # Try-catch ensures one bad plot doesn't crash the whole PDF generation
      tryCatch({
        p <- plot_raw_distribution_merged(
          data_df = data_df, 
          marker_name = marker,
          colors = colors,     
          highlight_pattern = hl_pattern,   
          highlight_color = "#FFD700"
        )
        if (!is.null(p)) print(p)
      }, error = function(e) {
        warning(sprintf("      [WARN] Failed to plot marker '%s': %s", marker, e$message))
      })
    }
  }
  
  dev.off()
}

#' @title Custom PCA Biplot with Multi-Pattern Highlighting
#' @description Plots Individuals with flexible highlighting based on metadata substrings.
#' @param pca_res Result from run_coda_pca().
#' @param metadata Dataframe containing 'Patient_ID', 'Group' and the subgroup column.
#' @param colors Named vector of colors for main Groups (fill).
#' @param dims Integer vector of length 2 indicating which PCs to plot.
#' @param show_labels Logical. If TRUE, adds Patient_ID labels.
#' @param highlight_patterns Named vector (Pattern -> Color). E.g., c("_LS" = "#FFD700").
#' @param find_col_keyword Keyword to identify the metadata column to search.
#' @return A ggplot object.
plot_pca_custom <- function(pca_res, metadata, colors, dims = c(1, 2), 
                            show_labels = FALSE, 
                            highlight_patterns = NULL,
                            find_col_keyword = "Original|Source") {
  
  # 1. Setup Coordinates & Metadata
  if (length(dims) != 2) stop("dims parameter must be a vector of length 2")
  ind_coords <- as.data.frame(pca_res$ind$coord[, dims])
  colnames(ind_coords) <- c("X_Coord", "Y_Coord")
  
  if (nrow(ind_coords) != nrow(metadata)) stop("Metadata/PCA dimension mismatch.")
  plot_data <- cbind(metadata, ind_coords)
  
  # 2. Logic for Highlights
  # Default state: No Highlight
  plot_data$Highlight_Type <- "Standard"
  plot_data$Stroke_Size <- 0.5 
  
  # Identify the column to search
  target_col <- names(metadata)[grepl(find_col_keyword, names(metadata), ignore.case = TRUE)][1]
  
  has_highlights <- !is.null(highlight_patterns) && length(highlight_patterns) > 0 && !is.na(target_col)
  
  if (has_highlights) {
    target_values <- as.character(plot_data[[target_col]])
    
    for (pattern in names(highlight_patterns)) {
      matches <- grepl(pattern, target_values)
      # Only overwrite if currently Standard (prevents overwriting higher priority matches)
      to_update <- matches & (plot_data$Highlight_Type == "Standard")
      
      if (any(to_update)) {
        plot_data$Highlight_Type[to_update] <- pattern 
        plot_data$Stroke_Size[to_update] <- 2.0  
      }
    }
  }
  
  # Prepare Color Scale for Borders
  border_colors <- c("Standard" = "white")
  if (has_highlights) {
    border_colors <- c(highlight_patterns, border_colors)
  }
  
  # 3. Variance Explained
  eig_val <- get_eigenvalue(pca_res)
  var_x <- round(eig_val[dims[1], 2], 1)
  var_y <- round(eig_val[dims[2], 2], 1)
  
  # 4. Plotting
  p <- ggplot(plot_data, aes(x = X_Coord, y = Y_Coord, fill = Group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray80") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray80")
  
  # Identify groups where ellipse calculation is mathematically possible
  valid_groups <- c()
  unique_grps <- unique(as.character(plot_data$Group))
  
  for (g in unique_grps) {
    sub_dat <- plot_data[plot_data$Group == g, c("X_Coord", "Y_Coord")]
    
    # Attempt covariance calculation as proxy for stat_ellipse feasibility
    is_computable <- tryCatch({
      if (nrow(sub_dat) < 3) stop("Insufficient points")
      # cov() will fail or return NA/Inf if variance is problematic
      cov_mat <- cov(sub_dat)
      if (any(is.na(cov_mat))) stop("NA in covariance")
      TRUE
    }, error = function(e) { FALSE })
    
    if (is_computable) valid_groups <- c(valid_groups, g)
  }
  
  # Only add the layer if we have valid groups
  if (length(valid_groups) > 0) {
    p <- p + stat_ellipse(
      data = plot_data[plot_data$Group %in% valid_groups, ],
      geom = "polygon", 
      alpha = 0.1, 
      show.legend = FALSE, 
      level = 0.95
    )
  }
  
  p <- p +
    # Points 
    geom_point(aes(color = Highlight_Type, stroke = Stroke_Size), 
               size = 3.5, shape = 21) +
    
    # Main Fill Colors (Clinical Group)
    scale_fill_manual(values = colors, guide = guide_legend(override.aes = list(shape = 21))) +
    
    # Border Colors (Highlights)
    scale_color_manual(name = "Subgroup", values = border_colors) +
    scale_continuous_identity(aesthetics = "stroke") +
    
    labs(
      x = sprintf("PC%d (%s%%)", dims[1], var_x),
      y = sprintf("PC%d (%s%%)", dims[2], var_y),
      title = sprintf("PCA (PC%d vs PC%d)", dims[1], dims[2]),
      subtitle = if(has_highlights) "Colored borders indicate specific subgroups" else NULL
    ) +
    theme_coda()
  
  # Hide legend if no highlights exist
  if (!has_highlights) {
    p <- p + guides(color = "none")
  }
  
  if (show_labels) {
    p <- p + geom_text_repel(aes(label = Patient_ID), size = 3, show.legend = FALSE, max.overlaps = 40)
  }
  
  return(p)
}

#' @title Variable Contribution Plot
#' @description Wraps fviz_pca_var for consistent styling.
plot_pca_variables <- function(pca_res) {
  fviz_pca_var(pca_res, 
               col.var = "contrib", 
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE) +
    labs(title = "Marker Contribution to Variance") +
    theme_coda()
}

#' @title Plot PCA Variance Dashboard
#' @description 
#' Visualizes explained variance per component with bars, curve, labels, 
#' and a side-list of cumulative variance values.
#' All labels are strictly in English.
#' 
#' @param pca_res Result object from FactoMineR::PCA.
#' @param n_list Integer. Number of top components to list on the side (default 10).
#' @return A ggplot object.
plot_pca_variance_dashboard <- function(pca_res, n_list = 10) {
  
  # 1. Extract Eigenvalues/Variance
  eig_df <- as.data.frame(pca_res$eig)
  # Standardize column names from FactoMineR
  colnames(eig_df) <- c("eigenvalue", "variance_percent", "cumulative_variance_percent")
  
  # Create PC identifiers (PC1, PC2...) and ensure factor order
  eig_df$PC <- factor(rownames(eig_df), levels = rownames(eig_df))
  eig_df$PC_Num <- 1:nrow(eig_df)
  
  # Limit to dimensions present
  n_pcs <- nrow(eig_df)
  # Dynamic limit for the plot (show max 15 bars for readability, or all if low dim)
  n_plot <- min(n_pcs, 15) 
  plot_data <- eig_df[1:n_plot, ]
  
  # 2. Prepare Side List Text
  n_list_actual <- min(n_pcs, n_list)
  list_data <- eig_df[1:n_list_actual, ]
  
  # Formatting the text table
  table_text <- paste0(
    "PC", list_data$PC_Num, ": ", 
    sprintf("%.1f%%", list_data$variance_percent), 
    " (Cum: ", sprintf("%.1f%%", list_data$cumulative_variance_percent), ")"
  )
  final_label <- paste(table_text, collapse = "\n")
  header_label <- paste0("Top ", n_list_actual, " Components\n(Var / Cumulative)")
  
  # 3. Plotting
  # We use a secondary axis scaling factor if needed, but for simplicity
  # since both variance and cumulative are %, we plot variance on Y.
  # The "Curve" usually represents the scree (variance), not cumulative, 
  # to match the histogram profile.
  
  p <- ggplot(plot_data, aes(x = PC, y = variance_percent)) +
    
    # A. Histograms (Bars)
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7, width = 0.7) +
    
    # B. The Curve (Scree line connecting bars)
    geom_line(aes(group = 1), color = "darkred", size = 1, linetype = "dashed") +
    geom_point(color = "darkred", size = 2) +
    
    # C. Labels above histograms
    geom_text(aes(label = sprintf("%.1f%%", variance_percent)), 
              vjust = -0.5, size = 3.5, fontface = "bold") +
    
    # D. The Side List (Annotation)
    # We create a text annotation on the right. 
    # Logic: Place it at x = n_plot + 0.5, y = max_variance.
    annotate("text", x = n_plot + 0.6, y = max(plot_data$variance_percent) * 0.9, 
             label = header_label, hjust = 0, vjust = 1, fontface = "bold", size = 4, color = "gray20") +
    annotate("text", x = n_plot + 0.6, y = max(plot_data$variance_percent) * 0.8, 
             label = final_label, hjust = 0, vjust = 1, size = 3.5, color = "gray30", lineheight = 1.2) +
    
    # Scales & Expansion
    scale_y_continuous(limits = c(0, max(plot_data$variance_percent) * 1.15), 
                       expand = c(0, 0)) +
    # Expand X axis to make room for the list on the right
    scale_x_discrete(expand = expansion(add = c(0.6, 4))) + 
    
    # Theme & Labels
    labs(
      title = "PCA Scree Plot & Variance Explained",
      subtitle = "Explained Variance per Principal Component",
      x = "Principal Component",
      y = "Explained Variance (%)"
    ) +
    theme_coda() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

#' @title Plot Stratification Heatmap (ComplexHeatmap)
#' @description 
#' Generates a clustered heatmap using Z-scored data.
#' @param mat_z A Z-scored numeric matrix (Samples x Markers).
#' @param metadata Dataframe containing patient metadata (aligned with mat_z rows).
#' @param annotation_colors_list A named list of color vectors for annotations.
#' @param title Plot title.
#' @return A ComplexHeatmap object (drawn).
plot_stratification_heatmap <- function(mat_z, metadata, annotation_colors_list, title = "Stratification") {
  
  # 1. Prepare Data
  mat_plot <- t(mat_z)
  
  if (!all(rownames(mat_z) == metadata$Patient_ID)) {
    stop("Mismatch between matrix rownames and metadata Patient_ID")
  }
  
  # Ensure metadata is a pure data.frame (not tibble) for HeatmapAnnotation
  metadata_clean <- as.data.frame(metadata)
  
  # 2. Setup Annotations
  ha <- HeatmapAnnotation(
    df = metadata_clean,
    col = annotation_colors_list,
    simple_anno_size = unit(0.3, "cm"),
    show_annotation_name = TRUE
  )
  
  # 3. Define Color Map for Z-Scores
  col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B"))
  
  # 4. Draw Heatmap
  hm <- Heatmap(
    mat_plot,
    name = "Z-Score",
    col = col_fun,
    
    # Clustering
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    
    # Annotations
    top_annotation = ha,
    
    # Aesthetics
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 6),
    show_column_names = FALSE, 
    row_dend_width = unit(2, "cm"),
    column_dend_height = unit(2, "cm"),
    
    column_title = paste0(title, " (n=", ncol(mat_plot), ")"),
    row_title = "Markers"
  )
  
  return(hm)
}

#' @title Generate Complex Heatmap Annotation Colors
#' @description 
#' Handles the logic for assigning colors to metadata columns in heatmaps.
#' Implements inheritance (Subgroup inherits from Group) and highlighting rules.
#' 
#' @param metadata Dataframe of metadata.
#' @param group_col Name of the main grouping column.
#' @param base_colors Named vector of base group colors.
#' @param hl_pattern Pattern to search for special highlighting.
#' @param hl_color Color to use for highlighted items.
#' @return A list of named vectors suitable for ComplexHeatmap.
viz_generate_complex_heatmap_colors <- function(metadata, group_col, base_colors, 
                                                hl_pattern = NULL, hl_color = "#FFD700") {
  
  annotation_colors <- list()
  
  # 1. Base Group Colors
  present_groups <- unique(metadata[[group_col]])
  # Intersect to avoid errors if config has more colors than data
  group_pal <- base_colors[intersect(names(base_colors), present_groups)]
  annotation_colors[[group_col]] <- group_pal
  
  # 2. Extra Columns Logic
  extra_cols <- setdiff(colnames(metadata), group_col)
  
  for (col_name in extra_cols) {
    # Get unique non-NA values
    raw_vals <- as.character(metadata[[col_name]])
    vals <- sort(unique(raw_vals[!is.na(raw_vals)]))
    
    col_pal <- character() # Force named vector
    
    for (v in vals) {
      # Priority 1: Exact Match to a known Group Name -> Inherit Color
      if (v %in% names(group_pal)) {
        col_pal[[v]] <- group_pal[[v]]
        
        # Priority 2: Matches Highlight Pattern -> Use Highlight Color
      } else if (!is.null(hl_pattern) && hl_pattern != "" && grepl(hl_pattern, v)) {
        col_pal[[v]] <- hl_color
        
        # Priority 3: Parent Group Containment (e.g. "NSCLC_LS" contains "NSCLC")
      } else {
        parent_match <- NA
        for (grp in names(group_pal)) {
          if (grepl(grp, v)) {
            parent_match <- grp
            break
          }
        }
        # Priority 4: Fallback Gray
        col_pal[[v]] <- if (!is.na(parent_match)) group_pal[[parent_match]] else "#D3D3D3"
      }
    }
    annotation_colors[[col_name]] <- col_pal
  }
  
  return(annotation_colors)
}

#' @title Plot Signature Boxplots (Top Drivers)
#' @description 
#' Helper function to visualize the distribution of the top discriminating markers.
#' Draws boxplots of Z-scores for the features with highest loadings.
#' 
#' @param data_matrix Numeric matrix of data used in PLS (usually Z-scored).
#' @param group_factor Factor vector defining groups.
#' @param loadings_df Dataframe of loadings from extract_plsda_loadings.
#' @param comp Integer. Which component to visualize?
#' @param n_top Integer. How many top markers to plot?
#' @param colors Named vector of colors.
#' @return A ggplot object.
viz_plot_signature_boxplots <- function(data_matrix, group_factor, loadings_df, 
                                        comp = 1, n_top = 6, colors) {
  
  col_name <- paste0("Comp", comp, "_Weight")
  if (!col_name %in% names(loadings_df)) return(NULL)
  
  # Filter top absolute weights
  top_mks <- loadings_df %>%
    arrange(desc(abs(!!sym(col_name)))) %>%
    head(n_top) %>%
    pull(Marker)
  
  if (length(top_mks) == 0) return(NULL)
  
  # Prepare long dataframe for plotting
  plot_df <- as.data.frame(data_matrix[, top_mks, drop=FALSE])
  plot_df$Group <- group_factor
  
  long_df <- plot_df %>%
    pivot_longer(-Group, names_to = "Marker", values_to = "Z_Score") %>%
    mutate(Marker = factor(Marker, levels = top_mks)) # Preserve order
  
  p <- ggplot(long_df, aes(x = Group, y = Z_Score, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~Marker, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = colors) +
    labs(
      title = sprintf("Top %d Drivers Distribution (Comp %d)", n_top, comp),
      subtitle = "Visual validation of PLS-DA weights (Z-Scores)",
      y = "Standardized Abundance (Z)", x = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", strip.background = element_rect(fill="gray95"))
  
  return(p)
}

#' @title Plot sPLS-DA Signature Heatmap
#' @description 
#' Custom heatmap for sPLS-DA markers with patients on rows and markers on columns.
#' @param mat_z Numeric matrix of Z-scored data (Patients x Markers).
#' @param metadata Dataframe with patient metadata.
#' @param group_col Column name defining the group.
#' @param colors Named vector of group colors.
#' @param title Plot title.
#' @return A ComplexHeatmap object.
viz_plot_splsda_heatmap <- function(mat_z, metadata, group_col, colors, title = "Global Signature Heatmap") {
  requireNamespace("ComplexHeatmap", quietly = TRUE)
  requireNamespace("circlize", quietly = TRUE)
  requireNamespace("grid", quietly = TRUE)
  
  if (!all(rownames(mat_z) == metadata$Patient_ID)) {
    stop("Mismatch between matrix rownames and metadata Patient_ID")
  }
  
  # Setup Row Annotations (Patients)
  meta_df <- as.data.frame(metadata[, group_col, drop = FALSE])
  anno_colors <- list()
  anno_colors[[group_col]] <- colors
  
  ra <- ComplexHeatmap::rowAnnotation(
    df = meta_df,
    col = anno_colors,
    show_annotation_name = FALSE,
    simple_anno_size = grid::unit(0.4, "cm")
  )
  
  # Standardized Z-Score color mapping
  col_fun <- circlize::colorRamp2(c(-3, 0, 3), c("#2166AC", "#F7F7F7", "#B2182B"))
  
  hm <- ComplexHeatmap::Heatmap(
    mat_z, 
    name = "Z-Score",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D2",
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "ward.D2",
    left_annotation = ra,
    row_names_gp = grid::gpar(fontsize = 7),
    column_names_gp = grid::gpar(fontsize = 8),
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_title = title,
    row_title = "Patients",
    column_title_side = "top"
  )
  
  return(hm)
}

#' @title Report sPLS-DA Visualization 
#' @description 
#' Generates a comprehensive PDF report for sPLS-DA results.
#' Includes robust error handling for PDF device closing and ggplot fallback.
#' 
#' @param pls_res The result object from run_splsda_model.
#' @param drivers_df Dataframe of extracted loadings.
#' @param metadata_viz Dataframe for plotting.
#' @param colors_viz Named vector of colors for groups.
#' @param out_path Path to save the PDF.
#' @param group_col The name of the metadata column to use for grouping/coloring (default: "Group").
#' @param n_top_boxplots Number of boxplots to plot. 
viz_report_plsda <- function(pls_res, drivers_df, metadata_viz, colors_viz, out_path, group_col = "Group", n_top_boxplots = 9) {
  
  if (is.null(pls_res$model)) return(NULL)
  
  requireNamespace("mixOmics", quietly = TRUE)
  requireNamespace("ComplexHeatmap", quietly = TRUE)
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(ggrepel)
  
  message(sprintf("   [Viz] Generating sPLS-DA graphical report: %s", basename(out_path)))
  
  if (!group_col %in% names(metadata_viz)) {
    stop(sprintf("Column '%s' not found in metadata for sPLS-DA visualization.", group_col))
  }
  
  pdf(out_path, width = 11, height = 8)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE) 
  
  plot_group_factor <- as.factor(metadata_viz[[group_col]])
  levels_present <- levels(plot_group_factor)
  plot_colors <- colors_viz[levels_present]
  
  if (any(is.na(plot_colors))) {
    missing_grps <- levels_present[is.na(plot_colors)]
    warning(paste("Missing colors for:", paste(missing_grps, collapse=", ")))
    plot_colors[is.na(plot_colors)] <- "gray50"
    names(plot_colors) <- levels_present
  }
  
  n_comps <- pls_res$model$ncomp
  
  # Helper for clean ggplot error messages
  plot_error_msg <- function(msg) {
    ggplot() + 
      annotate("text", x = 0, y = 0, label = msg, color = "darkred", size = 5, fontface = "bold") + 
      theme_void()
  }
  
  # 1. Native ggplot2 Biplot
  tryCatch({
    variates_df <- as.data.frame(pls_res$model$variates$X[, 1:2, drop = FALSE])
    colnames(variates_df) <- c("PC1", "PC2")
    variates_df$Group <- plot_group_factor
    variates_df$Patient_ID <- rownames(pls_res$model$variates$X)
    
    var_pct <- round(pls_res$model$prop_expl_var$X[1:2] * 100, 1)
    
    p_biplot <- ggplot(variates_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
      stat_ellipse(geom = "polygon", alpha = 0.1, level = 0.95, show.legend = FALSE) +
      geom_point(size = 3.5, shape = 21, color = "white", stroke = 0.5) +
      geom_text_repel(aes(label = Patient_ID), size = 3, show.legend = FALSE, max.overlaps = 20) +
      scale_fill_manual(values = plot_colors) +
      scale_color_manual(values = plot_colors) +
      labs(
        title = "sPLS-DA: Model Separation",
        subtitle = paste("Grouping by:", group_col),
        x = sprintf("X-variate 1: %s%% expl. var", var_pct[1]),
        y = sprintf("X-variate 2: %s%% expl. var", var_pct[2])
      ) +
      theme_coda()
    
    print(p_biplot)
  }, error = function(e) {
    print(plot_error_msg(paste("Biplot Error:", e$message)))
  })
  
  # 2. Loadings & Boxplots
  for (i in 1:n_comps) {
    comp_col <- paste0("Comp", i, "_Weight")
    if (!comp_col %in% names(drivers_df)) next
    
    variates <- pls_res$model$variates$X[, i]
    stat_groups <- plot_group_factor
    group_means <- tapply(variates, stat_groups, mean)
    
    pos_group_name <- names(group_means)[which.max(group_means)]
    neg_group_name <- names(group_means)[which.min(group_means)]
    
    label_sub <- sprintf("Direction: Positive -> %s | Negative -> %s", pos_group_name, neg_group_name)
    
    current_fill_colors <- c()
    if(pos_group_name %in% names(plot_colors)) current_fill_colors[[pos_group_name]] <- plot_colors[[pos_group_name]]
    if(neg_group_name %in% names(plot_colors)) current_fill_colors[[neg_group_name]] <- plot_colors[[neg_group_name]]
    if(is.null(current_fill_colors[[pos_group_name]])) current_fill_colors[[pos_group_name]] <- "gray"
    if(is.null(current_fill_colors[[neg_group_name]])) current_fill_colors[[neg_group_name]] <- "gray"
    
    df_comp <- drivers_df[abs(drivers_df[[comp_col]]) > 0, ]
    
    if (nrow(df_comp) > 0) {
      df_comp$Association <- ifelse(df_comp[[comp_col]] > 0, pos_group_name, neg_group_name)
      
      p_load <- ggplot(df_comp, aes(x = reorder(Marker, abs(.data[[comp_col]])), 
                                    y = .data[[comp_col]], fill = Association)) +
        geom_bar(stat = "identity", width = 0.7) +
        coord_flip() +
        scale_fill_manual(values = current_fill_colors) +
        labs(title = sprintf("sPLS-DA Loadings (Component %d)", i), 
             subtitle = label_sub, x = "Marker", y = "Weight Contribution", fill = "Associated Group") +
        theme_coda() + theme(legend.position = "bottom")
      print(p_load)
      
      if(exists("viz_plot_signature_boxplots")) {
        p_box <- viz_plot_signature_boxplots(
          data_matrix = pls_res$model$X, 
          group_factor = plot_group_factor,
          loadings_df = df_comp, 
          comp = i, 
          n_top = n_top_boxplots, 
          colors = plot_colors
        )
        if (!is.null(p_box)) print(p_box)
      }
    }
  }
  
  # 3. Transposed Signature Heatmap
  if (nrow(drivers_df) > 1) {
    tryCatch({
      mks <- unique(drivers_df$Marker)
      mat_sub <- pls_res$model$X[, mks, drop = FALSE]
      
      meta_hm <- metadata_viz[match(rownames(mat_sub), metadata_viz$Patient_ID), ]
      
      hm <- viz_plot_splsda_heatmap(
        mat_z = mat_sub,
        metadata = meta_hm,
        group_col = group_col,
        colors = plot_colors,
        title = "sPLS-DA Signature Heatmap"
      )
      
      ComplexHeatmap::draw(hm, merge_legend = TRUE)
    }, error = function(e) {
      print(plot_error_msg(paste("Heatmap Error:", e$message)))
    })
  }
}

#' @title Plot Network Graph (ggraph)
#' @description 
#' Visualizes the inferred network.
plot_network_structure <- function(adj_mat, weight_mat, title = "Network", 
                                   layout_type = "nicely", min_cor = 0) { 
  
  # 1. Build Graph Object
  g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
  
  # Add attributes
  E(g)$weight_raw <- NA
  E(g)$sign <- NA
  E(g)$weight <- NA 
  
  # Get edges
  el <- igraph::as_data_frame(g, what = "edges")
  
  if (nrow(el) > 0) {
    weights <- numeric(nrow(el))
    signs   <- character(nrow(el))
    keep_edge <- logical(nrow(el)) 
    
    for(k in 1:nrow(el)) {
      w <- weight_mat[el[k,1], el[k,2]]
      
      if (abs(w) >= min_cor) {
        keep_edge[k] <- TRUE
        weights[k] <- abs(w)
        signs[k]   <- ifelse(w > 0, "Positive", "Negative")
      } else {
        keep_edge[k] <- FALSE
      }
    }
    
    edges_to_remove <- E(g)[!keep_edge]
    g <- delete_edges(g, edges_to_remove)
    
    if (ecount(g) > 0) {
      E(g)$weight <- weights[keep_edge]
      E(g)$sign   <- signs[keep_edge]
    }
  }
  
  V(g)$degree <- igraph::degree(g)
  
  # 2. Plotting (Standard ggraph logic)
  tg <- tidygraph::as_tbl_graph(g)
  
  p <- ggraph(tg, layout = layout_type) + 
    geom_edge_link(aes(width = weight, color = sign), alpha = 0.6) +
    scale_edge_width(range = c(0.2, 1.5), guide = "none") +
    scale_edge_color_manual(values = c("Positive" = "#00BFC4", "Negative" = "#F8766D")) +
    
    geom_node_point(aes(size = degree), color = "gray20", fill = "white", shape=21) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3, max.overlaps = 20) +
    
    scale_size_continuous(range = c(2, 8)) +
    labs(title = title, 
         subtitle = sprintf("Nodes: %d | Edges: %d (Filter: |rho| > %.2f)", 
                            vcount(g), ecount(g), min_cor),
         edge_color = "Association") +
    theme_void() +
    theme(
      plot.title = element_text(face="bold", hjust=0.5),
      plot.subtitle = element_text(hjust=0.5, color="gray50"),
      legend.position = "bottom"
    )
  
  return(p)
}

#' @title Plot Differential Edge Overlap (Venn / UpSet)
#' @description 
#' Visualizes the intersection of differential edges across scenarios.
#' Uses ggVennDiagram for <= 3 sets and UpSet plots for > 3 sets.
#' 
#' @param edge_list A named list of character vectors (Edge IDs).
#' @param fill_colors Named vector of colors corresponding to names(edge_list).
#' @param title Plot title.
#' @return Returns the plot object.
viz_plot_differential_overlap <- function(edge_list, fill_colors = NULL, title = "Differential Overlap") {
  
  # Clean dependency loading without suppressPackageStartupMessages wrapper
  requireNamespace("ComplexHeatmap", quietly = TRUE)
  requireNamespace("grid", quietly = TRUE)
  requireNamespace("ggVennDiagram", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  
  if (length(edge_list) < 2) {
    warning("[Viz] Need at least 2 sets for overlap analysis.")
    return(NULL)
  }
  
  # --- Prepare Colors (Preserved Logic) ---
  final_colors <- NULL
  if (!is.null(fill_colors)) {
    common_names <- intersect(names(edge_list), names(fill_colors))
    if (length(common_names) > 0) {
      final_colors <- fill_colors[names(edge_list)]
      final_colors[is.na(final_colors)] <- "grey80" 
    }
  }
  
  # Default colors if mapping fails or not provided
  if (is.null(final_colors)) {
    defaults <- c("#4682B4", "#CD5C5C", "#E7B800", "#2E8B57", "#9932CC")
    n_needed <- length(edge_list)
    final_colors <- rep_len(defaults, n_needed)
    names(final_colors) <- names(edge_list)
  }
  
  # --- LOGIC: VENN (<= 3 sets) vs UPSET (> 3 sets) ---
  if (length(edge_list) <= 3) {
    
    # Render Venn
    p <- ggVennDiagram::ggVennDiagram(edge_list, label_alpha = 0, set_color = "black") +
      
      # Visualization Overrides
      ggplot2::scale_fill_gradient(low = "white", high = "white") +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(title = title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.2))
    
    print(p)
    return(p)
    
  } else {
    
    message("   [Viz] Drawing UpSet Plot (> 3 sets)...")
    
    # Create Combination Matrix
    m_comb <- ComplexHeatmap::make_comb_mat(edge_list)
    
    # Draw UpSet using ComplexHeatmap
    p <- ComplexHeatmap::draw(
      ComplexHeatmap::UpSet(
        m_comb,
        set_order = names(edge_list),
        comb_order = order(ComplexHeatmap::comb_size(m_comb), decreasing = TRUE),
        pt_size = grid::unit(3, "mm"), 
        lwd = 2,
        top_annotation = ComplexHeatmap::HeatmapAnnotation(
          "Intersection" = ComplexHeatmap::anno_barplot(
            ComplexHeatmap::comb_size(m_comb), 
            border = FALSE, 
            gp = grid::gpar(fill = "black"), 
            height = grid::unit(4, "cm")
          ), 
          annotation_name_side = "left"
        ),
        column_title = title
      )
    )
    return(p)
  }
}

#' @title Plot Hub-Driver Quadrant
#' @description Scatter plot of PLS-DA Importance vs Network Topology.
#' @param hub_driver_df Output from integrate_hub_drivers.
#' @param y_label String. Label for Y axis (e.g., "Degree" or "Betweenness").
#' @param title_suffix String to append to title.
#' @return ggplot object.
plot_hub_driver_quadrant <- function(hub_driver_df, y_label = "Degree", title_suffix = "") {
  
  require(ggplot2)
  require(ggrepel)
  
  if (is.null(hub_driver_df) || nrow(hub_driver_df) == 0) return(NULL)
  
  # Use the generic column created in integrate_hub_drivers
  y_col_name <- "Topology_Metric_Value" 
  
  # Fallback if the generic column is missing (e.g. using old integration)
  if (!y_col_name %in% names(hub_driver_df)) {
    y_col_name <- "Degree"
  }
  
  # Define Quadrant Lines
  x_mid <- median(hub_driver_df$Importance, na.rm = TRUE)
  y_mid <- median(hub_driver_df[[y_col_name]], na.rm = TRUE)
  
  p <- ggplot(hub_driver_df, aes(x = Importance, y = .data[[y_col_name]], fill = Role)) +
    # Quadrant Lines
    geom_vline(xintercept = x_mid, linetype = "dashed", color = "gray60") +
    geom_hline(yintercept = y_mid, linetype = "dashed", color = "gray60") +
    
    # Points
    geom_point(size = 4, shape = 21, alpha = 0.8) +
    
    # Labels
    geom_text_repel(aes(label = Marker), size = 3.5, max.overlaps = 20) +
    
    # Colors
    scale_fill_manual(values = c(
      "Master_Regulator" = "#B2182B",      # Red (High/High)
      "Solo_Driver" = "#D6604D",           # Light Red (High/Low)
      "Structural_Connector" = "#2166AC",  # Blue (Low/High)
      "Background" = "gray80"              # Gray (Low/Low)
    )) +
    
    # Scales
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    
    labs(
      title = paste("Hub-Driver Analysis", title_suffix),
      subtitle = paste("sPLS-DA Importance vs", y_label),
      x = "Statistical Importance (sPLS-DA)",
      y = paste("Topological Centrality (", y_label, ")", sep=""),
      caption = "Quadrants defined by median values"
    ) +
    theme_coda() +
    theme(legend.position = "bottom")
  
  return(p)
}

# R/modules_viz.R

#' @title Plot Partial Correlation Density
#' @description 
#' Visualizes the distribution of partial correlation values (Shrinkage).
#' Shows the cutoff threshold relative to the distribution of signal vs noise.
#' 
#' @param pcor_mat Numeric matrix of partial correlations.
#' @param adj_mat Optional. Adjacency matrix (0/1) of stable edges to count actual retained edges.
#' @param threshold Numeric. The magnitude threshold used for filtering.
#' @param group_label String. Label for the group (e.g., "Healthy").
#' @return A ggplot object.
viz_plot_edge_density <- function(pcor_mat, adj_mat = NULL, threshold = 0.15, group_label = "") {
  
  require(ggplot2)
  
  # Extract upper triangle values (exclude diagonal and duplicates)
  vals <- pcor_mat[upper.tri(pcor_mat)]
  n_total <- length(vals)
  
  df_plot <- data.frame(Value = vals)
  
  # Dynamic subtitle calculation based on provided adjacency matrix
  if (!is.null(adj_mat)) {
    n_stable <- sum(adj_mat[upper.tri(adj_mat)])
    sub_text <- sprintf("Threshold: |rho| > %.4f | Final Stable Edges: %d", threshold, n_stable)
  } else {
    n_kept <- sum(abs(vals) >= threshold)
    pct_kept <- round((n_kept / n_total) * 100, 1)
    sub_text <- sprintf("Threshold: |rho| > %.4f (Keeps %.1f%% of raw edges)", threshold, pct_kept)
  }
  
  p <- ggplot(df_plot, aes(x = Value)) +
    # Density Curve
    geom_density(fill = "steelblue", alpha = 0.3)
  
  # Conditionally add Threshold Lines only if threshold is greater than 0
  if (threshold > 0) {
    p <- p + geom_vline(xintercept = c(-threshold, threshold), 
                        linetype = "dashed", color = "red", size = 0.8)
  }
  
  p <- p +
    # Annotation
    labs(
      title = paste("Edge Weight Distribution:", group_label),
      subtitle = sub_text,
      x = "Partial Correlation (Shrinkage)",
      y = "Density"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40")
    ) +
    xlim(-1, 1) 
  
  return(p)
}

#' @title Plot Partial vs Raw Correlation Density Overlay
#' @description 
#' Visualizes the distribution of partial correlation values (Shrinkage) overlaid
#' on top of the raw Pearson correlation values to demonstrate the shrinkage effect.
#' 
#' @param pcor_mat Numeric matrix of partial correlations.
#' @param cor_mat Numeric matrix of raw Pearson correlations.
#' @param adj_mat Optional. Adjacency matrix (0/1) of stable edges.
#' @param threshold Numeric. The magnitude threshold used for filtering.
#' @param group_label String. Label for the group (e.g., "Healthy").
#' @return A ggplot object.
viz_plot_edge_density_overlay <- function(pcor_mat, cor_mat, adj_mat = NULL, threshold = 0.15, group_label = "") {
  
  require(ggplot2)
  
  vals_pcor <- pcor_mat[upper.tri(pcor_mat)]
  vals_cor <- cor_mat[upper.tri(cor_mat)]
  
  df_plot <- rbind(
    data.frame(Value = vals_cor, Metric = "Pearson (Raw)"),
    data.frame(Value = vals_pcor, Metric = "Partial (Shrinkage)")
  )
  
  # Set factor levels so Pearson renders behind Partial
  df_plot$Metric <- factor(df_plot$Metric, levels = c("Pearson (Raw)", "Partial (Shrinkage)"))
  
  if (!is.null(adj_mat)) {
    n_stable <- sum(adj_mat[upper.tri(adj_mat)])
    sub_text <- sprintf("Threshold: |rho| > %.4f | Final Stable Edges: %d", threshold, n_stable)
  } else {
    pct_kept <- round((sum(abs(vals_pcor) >= threshold) / length(vals_pcor)) * 100, 1)
    sub_text <- sprintf("Threshold: |rho| > %.4f (Keeps %.1f%% of Partial edges)", threshold, pct_kept)
  }
  
  p <- ggplot(df_plot, aes(x = Value, fill = Metric, color = Metric)) +
    geom_density(alpha = 0.4, size = 0.8)
  
  # Conditionally add Threshold Lines
  if (threshold > 0) {
    p <- p + geom_vline(xintercept = c(-threshold, threshold), 
                        linetype = "dashed", color = "red", size = 0.8)
  }
  
  p <- p + 
    scale_fill_manual(values = c("Pearson (Raw)" = "gray60", "Partial (Shrinkage)" = "steelblue")) +
    scale_color_manual(values = c("Pearson (Raw)" = "gray40", "Partial (Shrinkage)" = "#1F4E79")) +
    labs(
      title = paste("Correlation Distribution Overlay:", group_label),
      subtitle = sub_text,
      x = "Correlation Value",
      y = "Density",
      fill = "Metric",
      color = "Metric"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
      legend.position = "bottom"
    ) +
    xlim(-1, 1) 
  
  return(p)
}