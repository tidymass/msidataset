#' Plot m/z vs Spatial Correlation for MSI Data
#'
#' @description
#' 分析 m/z 变量与空间坐标(x, y)的相关性,识别具有空间模式的离子。
#' 这对于发现与组织结构相关的代谢物/脂质非常有用。
#'
#' @param mass_obj A mass_dataset object (converted from msi_dataset)
#' @param n_vars Integer, number of top variables to analyze (default: 20)
#' @param method Correlation method: "pearson", "spearman", "kendall" (default: "pearson")
#' @param intensity_threshold Minimum mean intensity threshold (default: 0)
#' @param show_labels Logical, whether to show m/z labels on points (default: FALSE)
#' @param label_top_n Integer, number of top correlated variables to label (default: 5)
#'
#' @return A patchwork object containing two correlation plots
#' @export
#'
#' @examples
#' plot_mz_spatial_correlation(my_mass, n_vars = 50)
#' plot_mz_spatial_correlation(my_mass, n_vars = 100, method = "spearman", show_labels = TRUE)
#'
plot_mz_spatial_correlation <- function(mass_obj, 
                                        n_vars = 20,
                                        method = c("pearson", "spearman", "kendall"),
                                        intensity_threshold = 0,
                                        show_labels = FALSE,
                                        label_top_n = 5) {
  
  # 参数验证
  method <- match.arg(method)
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.")
  }
  
  # ========== 步骤 1: 提取数据 ==========
  
  # 获取表达数据矩阵
  expr_data <- mass_obj@expression_data
  
  # 获取变量信息(包含 m/z)
  variable_info <- mass_obj@variable_info
  
  # 获取样本信息(包含空间坐标 x, y)
  sample_info <- mass_obj@sample_info
  
  # 检查是否有空间坐标
  if (!all(c("x", "y") %in% colnames(sample_info))) {
    stop("sample_info must contain 'x' and 'y' columns. ",
         "Did you use add_spatial_to_sample_info = TRUE when converting?")
  }
  
  cat("Analyzing", nrow(expr_data), "variables in", 
      ncol(expr_data), "pixels...\n")
  
  
  # ========== 步骤 2: 计算平均强度并筛选变量 ==========
  
  # 计算每个变量的平均强度
  mean_intensity <- rowMeans(expr_data, na.rm = TRUE)
  
  # 过滤低强度变量
  high_intensity_idx <- which(mean_intensity > intensity_threshold)
  
  if (length(high_intensity_idx) == 0) {
    stop("No variables pass the intensity threshold. Lower the threshold.")
  }
  
  # 按平均强度排序,选择前 n_vars 个
  top_intensity_idx <- order(mean_intensity, decreasing = TRUE)
  top_intensity_idx <- intersect(top_intensity_idx, high_intensity_idx)
  
  # 限制变量数量
  n_vars <- min(n_vars, length(top_intensity_idx))
  selected_vars <- top_intensity_idx[1:n_vars]
  
  cat("Selected top", n_vars, "variables based on mean intensity\n")
  
  
  # ========== 步骤 3: 计算与空间坐标的相关性 ==========
  
  cat("Calculating correlations with spatial coordinates...\n")
  
  # 提取空间坐标
  x_coords <- sample_info$x
  y_coords <- sample_info$y
  
  # 初始化结果向量
  cors_x <- numeric(n_vars)
  cors_y <- numeric(n_vars)
  pvals_x <- numeric(n_vars)
  pvals_y <- numeric(n_vars)
  
  # 对每个选定的变量计算相关性
  for (i in seq_along(selected_vars)) {
    var_idx <- selected_vars[i]
    
    # 提取该变量在所有像素的强度值
    intensity_values <- as.numeric(expr_data[var_idx, ])
    
    # 计算与 x 坐标的相关性
    cor_test_x <- cor.test(intensity_values, x_coords, 
                           method = method, 
                           use = "complete.obs")
    cors_x[i] <- cor_test_x$estimate
    pvals_x[i] <- cor_test_x$p.value
    
    # 计算与 y 坐标的相关性
    cor_test_y <- cor.test(intensity_values, y_coords, 
                           method = method,
                           use = "complete.obs")
    cors_y[i] <- cor_test_y$estimate
    pvals_y[i] <- cor_test_y$p.value
  }
  
  
  # ========== 步骤 4: 准备绘图数据 ==========
  
  # 创建数据框
  plot_data <- data.frame(
    variable_id = variable_info$variable_id[selected_vars],
    mz = variable_info$mz[selected_vars],
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
  plot_data$sig_x <- ifelse(plot_data$pval_x < 0.05, "Significant", "Not Significant")
  plot_data$sig_y <- ifelse(plot_data$pval_y < 0.05, "Significant", "Not Significant")
  
  # 计算总体相关性(用于着色)
  plot_data$total_cor <- sqrt(cors_x^2 + cors_y^2)
  
  cat("Found", sum(plot_data$sig_x == "Significant"), 
      "variables significantly correlated with X\n")
  cat("Found", sum(plot_data$sig_y == "Significant"), 
      "variables significantly correlated with Y\n")
  
  
  # ========== 步骤 5: 准备标签数据 ==========
  
  if (show_labels) {
    # 选择 X 方向相关性最强的变量
    top_x_idx <- order(plot_data$abs_cor_x, decreasing = TRUE)[1:label_top_n]
    label_data_x <- plot_data[top_x_idx, ]
    
    # 选择 Y 方向相关性最强的变量
    top_y_idx <- order(plot_data$abs_cor_y, decreasing = TRUE)[1:label_top_n]
    label_data_y <- plot_data[top_y_idx, ]
  }
  
  
  # ========== 步骤 6: 创建 X 相关性图 ==========
  
  p1 <- ggplot(plot_data, aes(x = mz, y = cor_x)) +
    # 添加显著性阈值线
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
    geom_hline(yintercept = c(-0.5, 0.5), linetype = "dotted", 
               color = "gray70", size = 0.5) +
    
    # 添加散点,按显著性和相关性强度着色和调整大小
    geom_point(aes(color = abs_cor_x, 
                   size = mean_intensity,
                   shape = sig_x),
               alpha = 0.7) +
    
    # 颜色映射:从灰色(低相关)到红色(高相关)
    scale_color_gradient(
      low = "lightblue", 
      high = "darkred",
      name = "|Correlation|",
      limits = c(0, 1)
    ) +
    
    # 大小映射:强度
    scale_size_continuous(
      name = "Mean Intensity",
      range = c(2, 8)
    ) +
    
    # 形状映射:显著性
    scale_shape_manual(
      values = c("Significant" = 16, "Not Significant" = 1),
      name = paste0("Significance\n(p < 0.05)")
    ) +
    
    # 标题和坐标轴标签
    labs(
      title = "m/z vs Correlation with X-coordinate",
      subtitle = paste0("Method: ", method, " | Top ", n_vars, " variables"),
      x = "m/z",
      y = paste0(tools::toTitleCase(method), " Correlation with X")
    ) +
    
    # 主题设置
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  # 添加标签(如果需要)
  if (show_labels) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p1 <- p1 + 
        ggrepel::geom_text_repel(
          data = label_data_x,
          aes(label = sprintf("%.2f", mz)),
          size = 3,
          max.overlaps = 20
        )
    } else {
      warning("Package 'ggrepel' not available. Labels will not be shown.")
    }
  }
  
  
  # ========== 步骤 7: 创建 Y 相关性图 ==========
  
  p2 <- ggplot(plot_data, aes(x = mz, y = cor_y)) +
    # 添加显著性阈值线
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
    geom_hline(yintercept = c(-0.5, 0.5), linetype = "dotted", 
               color = "gray70", size = 0.5) +
    
    # 添加散点
    geom_point(aes(color = abs_cor_y, 
                   size = mean_intensity,
                   shape = sig_y),
               alpha = 0.7) +
    
    # 颜色映射
    scale_color_gradient(
      low = "lightgreen", 
      high = "darkblue",
      name = "|Correlation|",
      limits = c(0, 1)
    ) +
    
    # 大小映射
    scale_size_continuous(
      name = "Mean Intensity",
      range = c(2, 8)
    ) +
    
    # 形状映射
    scale_shape_manual(
      values = c("Significant" = 16, "Not Significant" = 1),
      name = paste0("Significance\n(p < 0.05)")
    ) +
    
    # 标题和标签
    labs(
      title = "m/z vs Correlation with Y-coordinate",
      subtitle = paste0("Method: ", method, " | Top ", n_vars, " variables"),
      x = "m/z",
      y = paste0(tools::toTitleCase(method), " Correlation with Y")
    ) +
    
    # 主题
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  # 添加标签(如果需要)
  if (show_labels) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p2 <- p2 + 
        ggrepel::geom_text_repel(
          data = label_data_y,
          aes(label = sprintf("%.2f", mz)),
          size = 3,
          max.overlaps = 20
        )
    }
  }
  
  
  # ========== 步骤 8: 组合图形并返回 ==========
  
  # 使用 patchwork 组合两个图
  combined_plot <- p1 / p2
  
  # 添加总标题
  combined_plot <- combined_plot + 
    plot_annotation(
      title = "Spatial Correlation Analysis of m/z Variables",
      subtitle = paste0("Analyzed ", n_vars, " variables across ", 
                        ncol(expr_data), " pixels"),
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
  
  # 打印统计摘要
  cat("\n========== Correlation Summary ==========\n")
  cat("X-correlation range:", 
      sprintf("[%.3f, %.3f]", min(cors_x), max(cors_x)), "\n")
  cat("Y-correlation range:", 
      sprintf("[%.3f, %.3f]", min(cors_y), max(cors_y)), "\n")
  cat("Strong X-correlations (|r| > 0.5):", 
      sum(abs(cors_x) > 0.5), "\n")
  cat("Strong Y-correlations (|r| > 0.5):", 
      sum(abs(cors_y) > 0.5), "\n")
  cat("=========================================\n\n")
  
  # 返回图形和数据
  result <- list(
    plot = combined_plot,
    data = plot_data,
    summary = list(
      n_vars = n_vars,
      n_pixels = ncol(expr_data),
      method = method,
      sig_x = sum(plot_data$sig_x == "Significant"),
      sig_y = sum(plot_data$sig_y == "Significant")
    )
  )
  
  class(result) <- c("mz_spatial_cor", "list")
  
  return(result)
}


#' Print method for mz_spatial_cor object
#' @export
print.mz_spatial_cor <- function(x, ...) {
  print(x$plot)
}


#' Summary method for mz_spatial_cor object
#' @export
summary.mz_spatial_cor <- function(object, ...) {
  cat("\n========== Spatial Correlation Analysis Summary ==========\n\n")
  cat("Variables analyzed:", object$summary$n_vars, "\n")
  cat("Total pixels:", object$summary$n_pixels, "\n")
  cat("Correlation method:", object$summary$method, "\n\n")
  
  cat("Significant correlations (p < 0.05):\n")
  cat("  X-coordinate:", object$summary$sig_x, "\n")
  cat("  Y-coordinate:", object$summary$sig_y, "\n\n")
  
  cat("Top 5 variables by X-correlation:\n")
  top_x <- object$data[order(abs(object$data$cor_x), decreasing = TRUE)[1:5], ]
  print(top_x[, c("mz", "cor_x", "pval_x")])
  
  cat("\nTop 5 variables by Y-correlation:\n")
  top_y <- object$data[order(abs(object$data$cor_y), decreasing = TRUE)[1:5], ]
  print(top_y[, c("mz", "cor_y", "pval_y")])
  
  cat("\n==========================================================\n")
}