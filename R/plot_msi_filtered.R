# ============================================================================
# Visualization Functions for Filtered msi_dataset
# Author: Xiaoyu Xu
# Purpose: Plot spatial distribution of m/z ions in filtered MSI data
# ============================================================================

library(ggplot2)
library(viridis)

# ============================================================================
# 1. Plot Single Ion Spatial Distribution
# ============================================================================

#' Plot spatial distribution of a single ion in filtered msi_dataset
#'
#' @description
#' Visualizes the spatial distribution of a single m/z variable (ion) 
#' across the sampled pixels in a filtered msi_dataset object.
#'
#' @param msi_filtered A filtered msi_dataset object
#' @param variable_id Variable ID to plot (alternative to mz)
#' @param mz m/z value to plot (alternative to variable_id)
#' @param mz_tolerance Tolerance for m/z matching (default: 0.01)
#' @param title Plot title (auto-generated if NULL)
#' @param color_scale Color scale: "viridis" (default), "magma", "plasma", 
#'                     "inferno", "cividis"
#' @param transform Data transformation: "none" (default), "log10", "sqrt"
#' @param add_grid Logical, add grid lines (default: TRUE)
#'
#' @return A ggplot object
#'
#' @details
#' The function extracts intensity values for the specified variable and 
#' creates a heatmap-style visualization showing spatial distribution 
#' across the x,y coordinates.
#'
#' @examples
#' \dontrun{
#' # Plot by variable ID
#' plot_msi_filtered(msi_obj, variable_id = "var_1")
#' 
#' # Plot by m/z value with log transformation
#' plot_msi_filtered(msi_obj, mz = 885.55, transform = "log10")
#' }
#'
#' @export
plot_msi_filtered <- function(msi_filtered,
                              variable_id = NULL,
                              mz = NULL,
                              mz_tolerance = 0.01,
                              title = NULL,
                              color_scale = c("viridis", "magma", "plasma", "inferno", "cividis"),
                              transform = c("none", "log10", "sqrt"),
                              add_grid = TRUE) {
  
  # ========== 参数验证 ==========
  if (is.null(variable_id) && is.null(mz)) {
    stop("Either 'variable_id' or 'mz' must be provided")
  }
  
  color_scale <- match.arg(color_scale)
  transform <- match.arg(transform)
  
  # ========== 获取数据 ==========
  expr_data <- expression_data(msi_filtered)
  var_info <- variable_info(msi_filtered)
  sample_info <- sample_info(msi_filtered)
  
  # ========== 处理variable_id和mz的映射 ==========
  if (!is.null(mz) && is.null(variable_id)) {
    # 根据m/z查找variable_id
    matched_idx <- which(abs(var_info$mz - mz) <= mz_tolerance)
    
    if (length(matched_idx) == 0) {
      stop(sprintf("No variable found with m/z = %.4f (tolerance = %.4f)", 
                   mz, mz_tolerance))
    }
    
    if (length(matched_idx) > 1) {
      warning(sprintf("Multiple variables found with m/z = %.4f, using first match", mz))
      matched_idx <- matched_idx[1]
    }
    
    variable_id <- var_info$variable_id[matched_idx]
    actual_mz <- var_info$mz[matched_idx]
  } else {
    # 验证variable_id存在
    if (!variable_id %in% var_info$variable_id) {
      stop(sprintf("Variable '%s' not found in msi_dataset", variable_id))
    }
    actual_mz <- var_info$mz[var_info$variable_id == variable_id]
  }
  
  # ========== 检查变量是否在expression_data中 ==========
  if (!variable_id %in% rownames(expr_data)) {
    stop(sprintf("Variable '%s' not found in expression_data", variable_id))
  }
  
  # ========== 提取强度值 ==========
  intensity <- as.numeric(expr_data[variable_id, ])
  
  # ========== 数据转换 ==========
  intensity_transformed <- switch(
    transform,
    "none" = intensity,
    "log10" = log10(intensity + 1),
    "sqrt" = sqrt(intensity)
  )
  
  # ========== 创建绘图数据 ==========
  plot_data <- data.frame(
    x = sample_info$x,
    y = sample_info$y,
    intensity = intensity_transformed,
    stringsAsFactors = FALSE
  )
  
  # ========== 生成标题 ==========
  if (is.null(title)) {
    title <- sprintf("m/z %.4f | %s", actual_mz, variable_id)
    if (transform != "none") {
      title <- paste0(title, sprintf("\n(%s transformed)", transform))
    }
  }
  
  # ========== 绘制图形 ==========
  p <- ggplot(plot_data, aes(x = x, y = y, fill = intensity)) +
    geom_tile(color = NA) +
    coord_equal() +
    labs(
      title = title,
      x = "X coordinate",
      y = "Y coordinate",
      fill = "Intensity"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      panel.grid = if(add_grid) element_line(color = "gray95", size = 0.2) 
      else element_blank(),
      legend.position = "right",
      axis.title = element_text(size = 10)
    )
  
  # ========== 添加颜色映射 ==========
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


# ============================================================================
# 2. Plot Multiple Ions Grid
# ============================================================================

#' Plot multiple ions in a grid layout for filtered msi_dataset
#'
#' @description
#' Creates a grid of spatial distribution plots for multiple m/z variables,
#' useful for visualizing the distribution of several ions simultaneously.
#'
#' @param msi_filtered A filtered msi_dataset object
#' @param variable_ids Vector of variable IDs to plot
#' @param ncol Number of columns in the grid (default: 3)
#' @param color_scale Color scale for all plots (default: "viridis")
#' @param transform Data transformation for all plots (default: "none")
#'
#' @return A grid of ggplot objects
#'
#' @examples
#' \dontrun{
#' # Plot top 6 variables in 2x3 grid
#' top_vars <- head(variable_info(msi_obj)$variable_id, 6)
#' plot_msi_grid(msi_obj, top_vars, ncol = 3)
#' }
#'
#' @export
plot_msi_grid <- function(msi_filtered,
                          variable_ids,
                          ncol = 3,
                          color_scale = "viridis",
                          transform = "none") {
  
  # ========== 参数验证 ==========
  if (!is.character(variable_ids) || length(variable_ids) == 0) {
    stop("'variable_ids' must be a non-empty character vector")
  }
  
  if (ncol < 1) {
    stop("'ncol' must be at least 1")
  }
  
  # ========== 创建每个离子的单独图 ==========
  plots <- lapply(variable_ids, function(vid) {
    plot_msi_filtered(
      msi_filtered,
      variable_id = vid,
      color_scale = color_scale,
      transform = transform,
      add_grid = FALSE
    ) +
      theme(
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8)
      )
  })
  
  # ========== 组合成网格 ==========
  gridExtra::grid.arrange(grobs = plots, ncol = ncol)
}


# ============================================================================
# 3. Plot Spatial Correlation with Spatial Coordinates
# ============================================================================

#' Analyze and visualize spatial correlation of m/z variables
#'
#' @description
#' Calculates correlation between m/z variable intensities and spatial 
#' coordinates (x, y) to identify ions with spatial patterns. Returns 
#' visualization and correlation statistics.
#'
#' @param msi_filtered A filtered msi_dataset object
#' @param n_vars Number of top variables to analyze (default: 20)
#' @param method Correlation method: "pearson" (default), "spearman", "kendall"
#' @param show_labels Logical, show m/z labels on points (default: FALSE)
#' @param label_top_n Number of top correlated variables to label (default: 5)
#'
#' @return A list containing:
#'   - plot_x: ggplot object for X-coordinate correlation
#'   - plot_y: ggplot object for Y-coordinate correlation
#'   - data: Data frame with correlation results
#'   - summary: Summary statistics
#'
#' @details
#' This function identifies variables whose intensity patterns correlate 
#' with spatial location, which can indicate biologically meaningful 
#' spatial distribution patterns.
#'
#' @examples
#' \dontrun{
#' result <- plot_spatial_correlation(msi_obj, n_vars = 30)
#' # View plots
#' result$plot_x
#' result$plot_y
#' # Check correlation data
#' head(result$data)
#' }
#'
#' @export
plot_spatial_correlation <- function(msi_filtered,
                                     n_vars = 20,
                                     method = c("pearson", "spearman", "kendall"),
                                     show_labels = FALSE,
                                     label_top_n = 5) {
  
  # ========== 参数验证 ==========
  method <- match.arg(method)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required")
  }
  
  # ========== 提取数据 ==========
  expr_data <- expression_data(msi_filtered)
  var_info <- variable_info(msi_filtered)
  sample_info <- sample_info(msi_filtered)
  
  # 检查空间坐标
  if (!all(c("x", "y") %in% colnames(sample_info))) {
    stop("sample_info must contain 'x' and 'y' columns for spatial coordinates")
  }
  
  # ========== 计算平均强度并筛选 ==========
  mean_intensity <- rowMeans(expr_data, na.rm = TRUE)
  
  # 按平均强度排序，选择前n_vars个
  top_idx <- order(mean_intensity, decreasing = TRUE)[1:min(n_vars, nrow(expr_data))]
  selected_vars <- top_idx
  n_selected <- length(selected_vars)
  
  # ========== 计算与空间坐标的相关性 ==========
  x_coords <- sample_info$x
  y_coords <- sample_info$y
  
  cors_x <- numeric(n_selected)
  cors_y <- numeric(n_selected)
  pvals_x <- numeric(n_selected)
  pvals_y <- numeric(n_selected)
  
  for (i in seq_along(selected_vars)) {
    var_idx <- selected_vars[i]
    intensity_values <- as.numeric(expr_data[var_idx, ])
    
    # X坐标相关性
    cor_x_test <- cor.test(intensity_values, x_coords, 
                           method = method, use = "complete.obs")
    cors_x[i] <- cor_x_test$estimate
    pvals_x[i] <- cor_x_test$p.value
    
    # Y坐标相关性
    cor_y_test <- cor.test(intensity_values, y_coords, 
                           method = method, use = "complete.obs")
    cors_y[i] <- cor_y_test$estimate
    pvals_y[i] <- cor_y_test$p.value
  }
  
  # ========== 准备绘图数据 ==========
  plot_data <- data.frame(
    variable_id = var_info$variable_id[selected_vars],
    mz = var_info$mz[selected_vars],
    mean_intensity = mean_intensity[selected_vars],
    cor_x = cors_x,
    cor_y = cors_y,
    pval_x = pvals_x,
    pval_y = pvals_y,
    abs_cor_x = abs(cors_x),
    abs_cor_y = abs(cors_y),
    stringsAsFactors = FALSE
  )
  
  # 添加显著性标记
  plot_data$sig_x <- ifelse(plot_data$pval_x < 0.05, "Significant", "Not Sig.")
  plot_data$sig_y <- ifelse(plot_data$pval_y < 0.05, "Significant", "Not Sig.")
  
  # ========== 准备标签数据 ==========
  label_data_x <- NULL
  label_data_y <- NULL
  
  if (show_labels) {
    top_x_idx <- order(plot_data$abs_cor_x, decreasing = TRUE)[1:min(label_top_n, nrow(plot_data))]
    label_data_x <- plot_data[top_x_idx, ]
    
    top_y_idx <- order(plot_data$abs_cor_y, decreasing = TRUE)[1:min(label_top_n, nrow(plot_data))]
    label_data_y <- plot_data[top_y_idx, ]
  }
  
  # ========== 创建X相关性图 ==========
  p1 <- ggplot(plot_data, aes(x = mz, y = cor_x)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
    geom_hline(yintercept = c(-0.5, 0.5), linetype = "dotted", 
               color = "gray70", size = 0.3, alpha = 0.5) +
    geom_point(aes(color = abs_cor_x, 
                   size = mean_intensity,
                   shape = sig_x),
               alpha = 0.7) +
    scale_color_gradient(
      low = "lightblue", 
      high = "darkred",
      name = "|Correlation|",
      limits = c(0, 1)
    ) +
    scale_size_continuous(
      name = "Mean Intensity",
      range = c(2, 6)
    ) +
    scale_shape_manual(
      values = c("Significant" = 16, "Not Sig." = 1),
      name = "p < 0.05"
    ) +
    labs(
      title = "m/z vs X-coordinate Correlation",
      subtitle = paste0("Method: ", method, " | n=", n_selected),
      x = "m/z",
      y = paste0(tools::toTitleCase(method), " Correlation with X")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  # 添加标签
  if (!is.null(label_data_x) && requireNamespace("ggrepel", quietly = TRUE)) {
    p1 <- p1 + 
      ggrepel::geom_text_repel(
        data = label_data_x,
        aes(label = sprintf("%.1f", mz)),
        size = 3,
        max.overlaps = 20
      )
  }
  
  # ========== 创建Y相关性图 ==========
  p2 <- ggplot(plot_data, aes(x = mz, y = cor_y)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
    geom_hline(yintercept = c(-0.5, 0.5), linetype = "dotted", 
               color = "gray70", size = 0.3, alpha = 0.5) +
    geom_point(aes(color = abs_cor_y, 
                   size = mean_intensity,
                   shape = sig_y),
               alpha = 0.7) +
    scale_color_gradient(
      low = "lightgreen", 
      high = "darkblue",
      name = "|Correlation|",
      limits = c(0, 1)
    ) +
    scale_size_continuous(
      name = "Mean Intensity",
      range = c(2, 6)
    ) +
    scale_shape_manual(
      values = c("Significant" = 16, "Not Sig." = 1),
      name = "p < 0.05"
    ) +
    labs(
      title = "m/z vs Y-coordinate Correlation",
      subtitle = paste0("Method: ", method, " | n=", n_selected),
      x = "m/z",
      y = paste0(tools::toTitleCase(method), " Correlation with Y")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  # 添加标签
  if (!is.null(label_data_y) && requireNamespace("ggrepel", quietly = TRUE)) {
    p2 <- p2 + 
      ggrepel::geom_text_repel(
        data = label_data_y,
        aes(label = sprintf("%.1f", mz)),
        size = 3,
        max.overlaps = 20
      )
  }
  
  # ========== 编译结果 ==========
  result <- list(
    plot_x = p1,
    plot_y = p2,
    data = plot_data,
    summary = list(
      n_vars = n_selected,
      n_pixels = nrow(sample_info),
      method = method,
      sig_x = sum(plot_data$sig_x == "Significant"),
      sig_y = sum(plot_data$sig_y == "Significant"),
      cor_x_range = range(cors_x),
      cor_y_range = range(cors_y)
    )
  )
  
  class(result) <- c("spatial_correlation", "list")
  
  return(result)
}


# ============================================================================
# S3 Methods for spatial_correlation class
# ============================================================================

#' Print spatial_correlation object
#' @export
print.spatial_correlation <- function(x, ...) {
  cat("\n========== Spatial Correlation Analysis ==========\n")
  cat(sprintf("Variables analyzed: %d\n", x$summary$n_vars))
  cat(sprintf("Total pixels: %d\n", x$summary$n_pixels))
  cat(sprintf("Correlation method: %s\n\n", x$summary$method))
  
  cat("Significant correlations (p < 0.05):\n")
  cat(sprintf("  X-coordinate: %d variables\n", x$summary$sig_x))
  cat(sprintf("  Y-coordinate: %d variables\n", x$summary$sig_y))
  
  cat("\nCorrelation ranges:\n")
  cat(sprintf("  X: [%.3f, %.3f]\n", x$summary$cor_x_range[1], x$summary$cor_x_range[2]))
  cat(sprintf("  Y: [%.3f, %.3f]\n", x$summary$cor_y_range[1], x$summary$cor_y_range[2]))
  cat("==================================================\n\n")
  
  invisible(x)
}

#' Summary for spatial_correlation object
#' @export
summary.spatial_correlation <- function(object, top_n = 5, ...) {
  cat("\n========== Detailed Spatial Correlation Summary ==========\n\n")
  
  cat(sprintf("Analysis Parameters:\n"))
  cat(sprintf("  - Variables: %d\n", object$summary$n_vars))
  cat(sprintf("  - Pixels: %d\n", object$summary$n_pixels))
  cat(sprintf("  - Method: %s\n\n", object$summary$method))
  
  # Top X correlations
  cat(sprintf("Top %d variables by X-correlation:\n", top_n))
  top_x <- object$data[order(abs(object$data$cor_x), decreasing = TRUE), ]
  top_x <- head(top_x[, c("mz", "cor_x", "pval_x")], top_n)
  print(top_x)
  
  cat(sprintf("\nTop %d variables by Y-correlation:\n", top_n))
  top_y <- object$data[order(abs(object$data$cor_y), decreasing = TRUE), ]
  top_y <- head(top_y[, c("mz", "cor_y", "pval_y")], top_n)
  print(top_y)
  
  cat("\n========================================================\n\n")
  
  invisible(object)
}