# ============================================================================
# MSI-Specific Visualization and Analysis Functions
# Author: Xiaoyu Xu
# Purpose: Spatial analysis and visualization for msi_dataset objects
# ============================================================================

library(ggplot2)
library(viridis)
library(reshape2)
library(gridExtra)

# ============================================================================
# 1. Spatial Distribution Visualization
# ============================================================================

#' Plot spatial distribution of a single ion
#'
#' @param object An msi_dataset object
#' @param variable_id Variable ID to plot (alternative to mz)
#' @param mz m/z value to plot (alternative to variable_id)
#' @param mz_tolerance Tolerance for m/z matching (default: 0.01)
#' @param title Plot title (auto-generated if NULL)
#' @param color_scale Color scale: "viridis", "magma", "plasma", "inferno", "cividis"
#' @param transform Data transformation: "none", "log10", "sqrt"
#' @param add_grid Logical, add grid lines
#' @param point_size Size of points
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' plot_ion_spatial(msi_obj, variable_id = "var_1")
#' plot_ion_spatial(msi_obj, mz = 885.55, mz_tolerance = 0.05)
#' }
plot_ion_spatial <- function(object,
                             variable_id = NULL,
                             mz = NULL,
                             mz_tolerance = 0.01,
                             title = NULL,
                             color_scale = c("viridis", "magma", "plasma", "inferno", "cividis"),
                             transform = c("none", "log10", "sqrt"),
                             add_grid = TRUE,
                             point_size = 3) {
  
  # Validate inputs
  if (is.null(variable_id) && is.null(mz)) {
    stop("Either variable_id or mz must be provided")
  }
  
  color_scale <- match.arg(color_scale)
  transform <- match.arg(transform)
  
  # Get variable_id from mz if needed
  if (!is.null(mz) && is.null(variable_id)) {
    var_info <- variable_info(object)
    matched_vars <- which(abs(var_info$mz - mz) <= mz_tolerance)
    
    if (length(matched_vars) == 0) {
      stop(sprintf("No variable found with m/z = %.4f (tolerance = %.4f)", 
                   mz, mz_tolerance))
    }
    
    if (length(matched_vars) > 1) {
      warning(sprintf("Multiple variables found with m/z = %.4f, using first match", mz))
    }
    
    variable_id <- var_info$variable_id[matched_vars[1]]
    actual_mz <- var_info$mz[matched_vars[1]]
  } else {
    var_info <- variable_info(object)
    actual_mz <- var_info$mz[var_info$variable_id == variable_id]
  }
  
  # Extract data
  expr_data <- expression_data(object)
  pixel_inf <- pixel_info(object)
  
  # Check variable exists
  if (!variable_id %in% rownames(expr_data)) {
    stop(sprintf("Variable '%s' not found in expression data", variable_id))
  }
  
  # Get intensity values
  intensity <- as.numeric(expr_data[variable_id, ])
  
  # Apply transformation
  intensity_transformed <- switch(
    transform,
    "none" = intensity,
    "log10" = log10(intensity + 1),
    "sqrt" = sqrt(intensity)
  )
  
  # Create plotting data frame
  plot_data <- data.frame(
    x = pixel_inf$x,
    y = pixel_inf$y,
    intensity = intensity_transformed
  )
  
  # Generate title
  if (is.null(title)) {
    title <- sprintf("Spatial Distribution: m/z %.4f\n%s", 
                     actual_mz, variable_id)
    if (transform != "none") {
      title <- paste0(title, sprintf(" (%s transformed)", transform))
    }
  }
  
  # Create plot
  p <- ggplot(plot_data, aes(x = x, y = y, fill = intensity)) +
    geom_tile() +
    coord_equal() +
    labs(
      title = title,
      x = "X coordinate",
      y = "Y coordinate",
      fill = "Intensity"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = if(add_grid) element_line(color = "gray90") else element_blank(),
      legend.position = "right"
    )
  
  # Add color scale
  p <- p + switch(
    color_scale,
    "viridis" = scale_fill_viridis_c(option = "viridis"),
    "magma" = scale_fill_viridis_c(option = "magma"),
    "plasma" = scale_fill_viridis_c(option = "plasma"),
    "inferno" = scale_fill_viridis_c(option = "inferno"),
    "cividis" = scale_fill_viridis_c(option = "cividis")
  )
  
  return(p)
}


#' Plot multiple ions in a grid layout
#'
#' @param object An msi_dataset object
#' @param variable_ids Vector of variable IDs to plot
#' @param ncol Number of columns in grid
#' @param color_scale Color scale for all plots
#' @param transform Data transformation
#'
#' @return A grid of ggplot objects
#' @export
plot_ion_grid <- function(object,
                          variable_ids,
                          ncol = 3,
                          color_scale = "viridis",
                          transform = "none") {
  
  plots <- lapply(variable_ids, function(vid) {
    plot_ion_spatial(
      object,
      variable_id = vid,
      color_scale = color_scale,
      transform = transform,
      add_grid = FALSE,
      point_size = 2
    ) +
      theme(
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 8)
      )
  })
  
  do.call(grid.arrange, c(plots, ncol = ncol))
}


#' RGB overlay of three ions
#'
#' @param object An msi_dataset object
#' @param red_var Variable ID for red channel
#' @param green_var Variable ID for green channel
#' @param blue_var Variable ID for blue channel
#' @param title Plot title
#'
#' @return A ggplot object
#' @export
plot_rgb_overlay <- function(object,
                             red_var,
                             green_var,
                             blue_var,
                             title = "RGB Overlay") {
  
  expr_data <- expression_data(object)
  pixel_inf <- pixel_info(object)
  
  # Extract and normalize intensities
  normalize_channel <- function(var_id) {
    intensity <- as.numeric(expr_data[var_id, ])
    (intensity - min(intensity)) / (max(intensity) - min(intensity))
  }
  
  red <- normalize_channel(red_var)
  green <- normalize_channel(green_var)
  blue <- normalize_channel(blue_var)
  
  # Create RGB colors
  plot_data <- data.frame(
    x = pixel_inf$x,
    y = pixel_inf$y,
    color = rgb(red, green, blue)
  )
  
  # Plot
  ggplot(plot_data, aes(x = x, y = y, fill = color)) +
    geom_tile() +
    scale_fill_identity() +
    coord_equal() +
    labs(
      title = title,
      subtitle = sprintf("R: %s | G: %s | B: %s", red_var, green_var, blue_var),
      x = "X coordinate",
      y = "Y coordinate"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9)
    )
}


# ============================================================================
# 2. Quality Control Visualizations
# ============================================================================

#' Plot Total Ion Current (TIC) spatial distribution
#'
#' @param object An msi_dataset object
#' @param transform Transformation: "none", "log10", "sqrt"
#' @param color_scale Color scale
#'
#' @return A ggplot object
#' @export
plot_tic_spatial <- function(object,
                             transform = c("none", "log10", "sqrt"),
                             color_scale = "viridis") {
  
  transform <- match.arg(transform)
  
  expr_data <- expression_data(object)
  pixel_inf <- pixel_info(object)
  
  # Calculate TIC for each pixel
  tic <- colSums(expr_data)
  
  # Apply transformation
  tic_transformed <- switch(
    transform,
    "none" = tic,
    "log10" = log10(tic + 1),
    "sqrt" = sqrt(tic)
  )
  
  plot_data <- data.frame(
    x = pixel_inf$x,
    y = pixel_inf$y,
    tic = tic_transformed
  )
  
  title <- "Total Ion Current (TIC)"
  if (transform != "none") {
    title <- paste0(title, sprintf(" - %s transformed", transform))
  }
  
  ggplot(plot_data, aes(x = x, y = y, fill = tic)) +
    geom_tile() +
    coord_equal() +
    scale_fill_viridis_c(option = color_scale) +
    labs(
      title = title,
      x = "X coordinate",
      y = "Y coordinate",
      fill = "TIC"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}


#' Plot number of detected peaks per pixel
#'
#' @param object An msi_dataset object
#' @param threshold Intensity threshold for peak detection
#'
#' @return A ggplot object
#' @export
plot_peak_count_spatial <- function(object, threshold = 0) {
  
  expr_data <- expression_data(object)
  pixel_inf <- pixel_info(object)
  
  # Count peaks above threshold per pixel
  peak_counts <- colSums(expr_data > threshold)
  
  plot_data <- data.frame(
    x = pixel_inf$x,
    y = pixel_inf$y,
    n_peaks = peak_counts
  )
  
  ggplot(plot_data, aes(x = x, y = y, fill = n_peaks)) +
    geom_tile() +
    coord_equal() +
    scale_fill_viridis_c() +
    labs(
      title = sprintf("Number of Detected Peaks (threshold > %.1f)", threshold),
      x = "X coordinate",
      y = "Y coordinate",
      fill = "# Peaks"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}


# ============================================================================
# 3. Spatial Statistical Analysis
# ============================================================================

#' Calculate spatial autocorrelation (Moran's I)
#'
#' @param object An msi_dataset object
#' @param variable_id Variable ID to analyze
#' @param distance_threshold Maximum distance for neighbors
#'
#' @return A list with Moran's I statistic and p-value
#' @export
calculate_morans_i <- function(object, 
                               variable_id,
                               distance_threshold = NULL) {
  
  pixel_inf <- pixel_info(object)
  expr_data <- expression_data(object)
  
  # Get coordinates and intensity
  coords <- as.matrix(pixel_inf[, c("x", "y")])
  intensity <- as.numeric(expr_data[variable_id, ])
  
  n <- length(intensity)
  
  # Calculate distance matrix
  dist_matrix <- as.matrix(dist(coords))
  
  # Define spatial weights
  if (is.null(distance_threshold)) {
    # Use mean of minimum distances
    min_dists <- apply(dist_matrix + diag(Inf, n), 1, min)
    distance_threshold <- mean(min_dists) * 1.5
  }
  
  weights <- (dist_matrix <= distance_threshold) * 1
  diag(weights) <- 0
  
  # Standardize weights
  row_sums <- rowSums(weights)
  weights[row_sums > 0, ] <- weights[row_sums > 0, ] / row_sums[row_sums > 0]
  
  # Calculate Moran's I
  intensity_mean <- mean(intensity)
  intensity_centered <- intensity - intensity_mean
  
  numerator <- sum(weights * outer(intensity_centered, intensity_centered))
  denominator <- sum(intensity_centered^2)
  
  morans_i <- (n / sum(weights)) * (numerator / denominator)
  
  # Calculate expected value and variance (simplified)
  expected_i <- -1 / (n - 1)
  
  # Simple z-score
  variance_i <- 1 / (n - 1)  # Simplified
  z_score <- (morans_i - expected_i) / sqrt(variance_i)
  p_value <- 2 * pnorm(-abs(z_score))
  
  result <- list(
    morans_i = morans_i,
    expected_i = expected_i,
    z_score = z_score,
    p_value = p_value,
    distance_threshold = distance_threshold,
    interpretation = ifelse(
      morans_i > expected_i,
      "Positive spatial autocorrelation (clustering)",
      "Negative spatial autocorrelation (dispersion)"
    )
  )
  
  class(result) <- "morans_i_result"
  return(result)
}

#' Print method for Moran's I result
#' @export
print.morans_i_result <- function(x, ...) {
  cat("Moran's I Spatial Autocorrelation Test\n")
  cat("======================================\n\n")
  cat(sprintf("Moran's I:       %.4f\n", x$morans_i))
  cat(sprintf("Expected I:      %.4f\n", x$expected_i))
  cat(sprintf("Z-score:         %.4f\n", x$z_score))
  cat(sprintf("P-value:         %.4e\n", x$p_value))
  cat(sprintf("Distance:        %.2f\n", x$distance_threshold))
  cat("\nInterpretation:", x$interpretation, "\n")
  
  if (x$p_value < 0.001) {
    cat("Significance: *** (p < 0.001)\n")
  } else if (x$p_value < 0.01) {
    cat("Significance: ** (p < 0.01)\n")
  } else if (x$p_value < 0.05) {
    cat("Significance: * (p < 0.05)\n")
  } else {
    cat("Significance: Not significant (p >= 0.05)\n")
  }
}


#' Spatial clustering using k-means
#'
#' @param object An msi_dataset object
#' @param n_clusters Number of clusters
#' @param method Clustering method: "kmeans", "hierarchical"
#' @param plot Logical, create visualization
#'
#' @return An msi_dataset object with cluster assignments added to pixel_info
#' @export
spatial_clustering <- function(object,
                               n_clusters = 5,
                               method = c("kmeans", "hierarchical"),
                               plot = TRUE) {
  
  method <- match.arg(method)
  
  expr_data <- expression_data(object)
  pixel_inf <- pixel_info(object)
  
  # Transpose: pixels as rows, variables as columns
  data_matrix <- t(as.matrix(expr_data))
  
  # Perform clustering
  if (method == "kmeans") {
    clustering <- kmeans(data_matrix, centers = n_clusters, nstart = 25)
    clusters <- clustering$cluster
  } else if (method == "hierarchical") {
    hc <- hclust(dist(data_matrix))
    clusters <- cutree(hc, k = n_clusters)
  }
  
  # Add cluster assignments to pixel_info
  object@pixel_info$cluster <- as.factor(clusters)
  
  # Update pixel_info_note
  if (!"cluster" %in% object@pixel_info_note$name) {
    object@pixel_info_note <- rbind(
      object@pixel_info_note,
      data.frame(
        name = "cluster",
        meaning = sprintf("Spatial cluster assignment (k=%d, method=%s)", 
                          n_clusters, method),
        stringsAsFactors = FALSE
      )
    )
  }
  
  # Create visualization
  if (plot) {
    plot_data <- data.frame(
      x = pixel_inf$x,
      y = pixel_inf$y,
      cluster = clusters
    )
    
    p <- ggplot(plot_data, aes(x = x, y = y, fill = factor(cluster))) +
      geom_tile() +
      coord_equal() +
      scale_fill_viridis_d() +
      labs(
        title = sprintf("Spatial Clustering (%s, k=%d)", 
                        method, n_clusters),
        x = "X coordinate",
        y = "Y coordinate",
        fill = "Cluster"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
    print(p)
  }
  
  return(object)
}


# ============================================================================
# 4. ROI (Region of Interest) Analysis
# ============================================================================

#' Compare two regions of interest
#'
#' @param object An msi_dataset object
#' @param roi1_condition Condition for ROI 1 (e.g., "x <= 10 & y <= 10")
#' @param roi2_condition Condition for ROI 2
#' @param roi1_name Name for ROI 1
#' @param roi2_name Name for ROI 2
#' @param test_method Statistical test: "t.test", "wilcox.test"
#' @param top_n Number of top differential variables to display
#'
#' @return A data frame with differential analysis results
#' @export
compare_roi <- function(object,
                        roi1_condition,
                        roi2_condition,
                        roi1_name = "ROI1",
                        roi2_name = "ROI2",
                        test_method = c("t.test", "wilcox.test"),
                        top_n = 20) {
  
  test_method <- match.arg(test_method)
  
  pixel_inf <- pixel_info(object)
  expr_data <- expression_data(object)
  var_info <- variable_info(object)
  
  # Identify pixels in each ROI
  roi1_pixels <- with(pixel_inf, eval(parse(text = roi1_condition)))
  roi2_pixels <- with(pixel_inf, eval(parse(text = roi2_condition)))
  
  cat(sprintf("ROI 1 (%s): %d pixels\n", roi1_name, sum(roi1_pixels)))
  cat(sprintf("ROI 2 (%s): %d pixels\n", roi2_name, sum(roi2_pixels)))
  
  if (sum(roi1_pixels) == 0 || sum(roi2_pixels) == 0) {
    stop("One or both ROIs contain no pixels")
  }
  
  # Perform statistical test for each variable
  results <- data.frame(
    variable_id = var_info$variable_id,
    mz = var_info$mz,
    mean_roi1 = numeric(nrow(var_info)),
    mean_roi2 = numeric(nrow(var_info)),
    fold_change = numeric(nrow(var_info)),
    p_value = numeric(nrow(var_info)),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(results)) {
    var_id <- results$variable_id[i]
    
    roi1_values <- as.numeric(expr_data[var_id, roi1_pixels])
    roi2_values <- as.numeric(expr_data[var_id, roi2_pixels])
    
    results$mean_roi1[i] <- mean(roi1_values)
    results$mean_roi2[i] <- mean(roi2_values)
    results$fold_change[i] <- results$mean_roi2[i] / (results$mean_roi1[i] + 1e-10)
    
    # Perform test
    test_result <- switch(
      test_method,
      "t.test" = t.test(roi1_values, roi2_values),
      "wilcox.test" = wilcox.test(roi1_values, roi2_values)
    )
    
    results$p_value[i] <- test_result$p.value
  }
  
  # Adjust p-values
  results$p_adjusted <- p.adjust(results$p_value, method = "BH")
  results$log2_fc <- log2(results$fold_change)
  results$significant <- results$p_adjusted < 0.05
  
  # Sort by p-value
  results <- results[order(results$p_value), ]
  
  cat(sprintf("\nSignificant variables (p.adj < 0.05): %d\n", 
              sum(results$significant)))
  
  # Create volcano plot
  plot_data <- results
  plot_data$color <- ifelse(
    plot_data$significant & abs(plot_data$log2_fc) > 1,
    "Significant",
    "Not significant"
  )
  
  p <- ggplot(plot_data, aes(x = log2_fc, y = -log10(p_value), color = color)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray")) +
    labs(
      title = sprintf("ROI Comparison: %s vs %s", roi1_name, roi2_name),
      x = "log2(Fold Change)",
      y = "-log10(P-value)",
      color = ""
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "top"
    )
  
  print(p)
  
  # Display top results
  cat("\nTop", top_n, "most significant variables:\n")
  print(head(results[, c("variable_id", "mz", "log2_fc", "p_value", "p_adjusted")], 
             top_n))
  
  return(results)
}


# ============================================================================
# 5. Dimensionality Reduction with Spatial Context
# ============================================================================

#' Spatial PCA visualization
#'
#' @param object An msi_dataset object
#' @param n_components Number of PCA components
#' @param color_by How to color points: "spatial", "pc1", "pc2"
#'
#' @return A list containing PCA results and plots
#' @export
plot_spatial_pca <- function(object,
                             n_components = 2,
                             color_by = c("spatial", "pc1", "pc2", "cluster")) {
  
  color_by <- match.arg(color_by)
  
  expr_data <- expression_data(object)
  pixel_inf <- pixel_info(object)
  
  # Perform PCA (pixels as samples)
  data_matrix <- t(as.matrix(expr_data))
  pca_result <- prcomp(data_matrix, center = TRUE, scale. = TRUE)
  
  # Extract PC scores
  pc_scores <- pca_result$x[, 1:n_components]
  
  # Calculate variance explained
  var_explained <- summary(pca_result)$importance[2, 1:n_components] * 100
  
  # Create plotting data
  plot_data <- data.frame(
    x = pixel_inf$x,
    y = pixel_inf$y,
    PC1 = pc_scores[, 1],
    PC2 = if (n_components >= 2) pc_scores[, 2] else 0
  )
  
  # Choose coloring
  color_var <- switch(
    color_by,
    "spatial" = pixel_inf$x + pixel_inf$y,
    "pc1" = plot_data$PC1,
    "pc2" = plot_data$PC2,
    "cluster" = if ("cluster" %in% colnames(pixel_inf)) {
      as.numeric(pixel_inf$cluster)
    } else {
      pixel_inf$x  # fallback
    }
  )
  
  plot_data$color_value <- color_var
  
  # Create dual plot
  p1 <- ggplot(plot_data, aes(x = x, y = y, fill = color_value)) +
    geom_tile() +
    coord_equal() +
    scale_fill_viridis_c() +
    labs(
      title = "Spatial Map",
      subtitle = sprintf("Colored by: %s", color_by),
      x = "X coordinate",
      y = "Y coordinate",
      fill = color_by
    ) +
    theme_minimal()
  
  p2 <- ggplot(plot_data, aes(x = PC1, y = PC2, color = color_value)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_viridis_c() +
    labs(
      title = "PCA Space",
      x = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
      y = sprintf("PC2 (%.1f%% variance)", 
                  if (n_components >= 2) var_explained[2] else 0),
      color = color_by
    ) +
    theme_minimal()
  
  grid.arrange(p1, p2, ncol = 2)
  
  invisible(list(
    pca = pca_result,
    scores = pc_scores,
    variance_explained = var_explained
  ))
}