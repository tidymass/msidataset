# ============================================================================
# Supplementary filter_msi_dataset Functions - UPDATED
# 移除：频率筛选（freq已从variable_info中删除）
# 保留：高频筛选（基于零值计数），零值过多筛选
# ============================================================================

# ════════════════════════════════════════════════════════════════════════════════
# 补充筛选函数 1: 排除零值过多的 M/Z（> 60% 零值）
# ════════════════════════════════════════════════════════════════════════════════

#' Filter Out M/Z Features with Excessive Zero Values (> 60% zeros)
#'
#' @description
#' Remove m/z features where more than 60% of the values are zero.
#' This is useful for removing low-quality, sparse features.
#'
#' @param object A msi_dataset object
#' @param zero_threshold Threshold for percentage of zeros (default: 0.6)
#'
#' @return Filtered msi_dataset with low-zero features only
#' @export
#'
#' @examples
#' \dontrun{
#' msi_no_zeros <- filter_out_zero_heavy_mz(msi)
#' }

filter_out_zero_heavy_mz <- function(object, zero_threshold = 0.6) {
  
  if (!methods::is(object, "msi_dataset"))
    stop("object must be a msi_dataset.")
  
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("Filtering Out Zero-Heavy M/Z Features (zeros > ", zero_threshold*100, "%)\n")
  cat("───────────────────────────────────────────────────────────────\n\n")
  
  # Calculate percentage of zeros for each feature
  zero_percentages <- apply(object@expression_data, 1, function(x) {
    sum(x == 0, na.rm = TRUE) / length(x)
  })
  
  # Keep features with zero percentage <= threshold
  keep_idx <- zero_percentages <= zero_threshold
  keep_ids <- object@variable_info$variable_id[keep_idx]
  
  num_original <- nrow(object@variable_info)
  num_remaining <- sum(keep_idx)
  pct_retained <- 100 * num_remaining / num_original
  
  cat(sprintf("Zero value threshold: <= %.1f%%\n", zero_threshold * 100))
  cat(sprintf("\nFeatures with <= %.1f%% zeros: %d\n", zero_threshold*100, num_remaining))
  cat(sprintf("Features removed (> %.1f%% zeros): %d\n", zero_threshold*100, 
              num_original - num_remaining))
  cat(sprintf("Retention rate: %.1f%%\n\n", pct_retained))
  
  # Show distribution of zero percentages
  zero_stats <- zero_percentages[keep_idx] * 100
  if (length(zero_stats) > 0) {
    cat(sprintf("Zero percentage statistics (remaining features):\n"))
    cat(sprintf("  Min:    %.2f%%\n", min(zero_stats)))
    cat(sprintf("  Q1:     %.2f%%\n", quantile(zero_stats, 0.25)))
    cat(sprintf("  Median: %.2f%%\n", median(zero_stats)))
    cat(sprintf("  Q3:     %.2f%%\n", quantile(zero_stats, 0.75)))
    cat(sprintf("  Max:    %.2f%%\n\n", max(zero_stats)))
  }
  
  # Update object
  new_expr <- object@expression_data[keep_ids, , drop = FALSE]
  new_vinfo <- object@variable_info[keep_idx, ]
  
  object@expression_data <- new_expr
  object@variable_info <- new_vinfo
  
  # Record in process_info
  process_key <- paste0("filter_zero_heavy_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "filter_out_zero_heavy_mz()",
    parameter = list(zero_threshold = zero_threshold,
                     num_retained = num_remaining,
                     num_removed = num_original - num_remaining),
    time = Sys.time()
  )
  object@process_info[[process_key]] <- parameter
  
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  return(object)
}


# ════════════════════════════════════════════════════════════════════════════════
# 补充筛选函数 2: 高频 M/Z 筛选（基于非零像素数量）
# ════════════════════════════════════════════════════════════════════════════════

#' Filter High-Frequency M/Z Features (Based on Non-Zero Pixel Count)
#'
#' @description
#' Retain only m/z features that have non-zero values in at least N pixels.
#' This is useful for selecting robust, abundant features.
#' 
#' Note: This calculates frequency dynamically from the expression data,
#' rather than reading from a stored 'freq' column.
#'
#' @param object A msi_dataset object
#' @param min_pixels Minimum number of pixels where feature must be non-zero
#'
#' @return Filtered msi_dataset with high-frequency features only
#' @export
#'
#' @examples
#' \dontrun{
#' msi_high_freq <- filter_high_frequency_mz_by_count(msi, min_pixels = 10)
#' }

filter_high_frequency_mz_by_count <- function(object, min_pixels = NULL) {
  
  if (!methods::is(object, "msi_dataset"))
    stop("object must be a msi_dataset.")
  
  # If min_pixels not specified, use default of 50% of pixels
  if (is.null(min_pixels)) {
    min_pixels <- ceiling(ncol(object@expression_data) * 0.5)
  }
  
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("Filtering High-Frequency M/Z Features (min_pixels >= ", min_pixels, ")\n")
  cat("───────────────────────────────────────────────────────────────\n\n")
  
  # Calculate number of non-zero pixels for each feature
  non_zero_counts <- apply(object@expression_data, 1, function(x) {
    sum(x != 0, na.rm = TRUE)
  })
  
  # Keep features with non-zero count >= threshold
  keep_idx <- non_zero_counts >= min_pixels
  keep_ids <- object@variable_info$variable_id[keep_idx]
  
  num_original <- nrow(object@variable_info)
  num_remaining <- sum(keep_idx)
  pct_retained <- 100 * num_remaining / num_original
  
  cat(sprintf("Minimum non-zero pixels: %d out of %d (%.1f%%)\n", 
              min_pixels, ncol(object@expression_data),
              100 * min_pixels / ncol(object@expression_data)))
  cat(sprintf("\nFeatures with >= %d non-zero pixels: %d\n", min_pixels, num_remaining))
  cat(sprintf("Features removed (< %d non-zero pixels): %d\n", min_pixels, 
              num_original - num_remaining))
  cat(sprintf("Retention rate: %.1f%%\n\n", pct_retained))
  
  # Show distribution of non-zero counts
  count_stats <- non_zero_counts[keep_idx]
  if (length(count_stats) > 0) {
    cat(sprintf("Non-zero pixel count statistics (remaining features):\n"))
    cat(sprintf("  Min:    %d pixels\n", min(count_stats)))
    cat(sprintf("  Q1:     %d pixels\n", quantile(count_stats, 0.25)))
    cat(sprintf("  Median: %d pixels\n", median(count_stats)))
    cat(sprintf("  Q3:     %d pixels\n", quantile(count_stats, 0.75)))
    cat(sprintf("  Max:    %d pixels\n\n", max(count_stats)))
  }
  
  # Update object
  new_expr <- object@expression_data[keep_ids, , drop = FALSE]
  new_vinfo <- object@variable_info[keep_idx, ]
  
  object@expression_data <- new_expr
  object@variable_info <- new_vinfo
  
  # Record in process_info
  process_key <- paste0("filter_high_freq_count_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "filter_high_frequency_mz_by_count()",
    parameter = list(min_pixels = min_pixels,
                     num_retained = num_remaining,
                     num_removed = num_original - num_remaining),
    time = Sys.time()
  )
  object@process_info[[process_key]] <- parameter
  
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  return(object)
}


# ════════════════════════════════════════════════════════════════════════════════
# 组合函数：同时应用两种补充筛选
# ════════════════════════════════════════════════════════════════════════════════

#' Apply Both Supplementary Filters
#'
#' @description
#' Apply both zero-heavy and high-frequency filters in sequence.
#' This is a convenience function combining the two supplementary filters.
#'
#' @param object A msi_dataset object
#' @param min_pixels Minimum number of pixels for high-frequency filter (default: NULL, uses 50%)
#' @param zero_threshold Zero percentage threshold (default: 0.6)
#'
#' @return Filtered msi_dataset
#' @export

filter_supplementary <- function(object, min_pixels = NULL, zero_threshold = 0.6) {
  
  if (!methods::is(object, "msi_dataset"))
    stop("object must be a msi_dataset.")
  
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════════╗\n")
  cat("║              Applying Supplementary Filters                   ║\n")
  cat("╚═══════════════════════════════════════════════════════════════╝\n\n")
  
  # Step 1: Filter zero-heavy m/z
  cat("Step 1: Zero-Heavy M/Z Filter\n")
  cat("─────────────────────────────────────────────────────────────\n")
  object <- filter_out_zero_heavy_mz(object, zero_threshold = zero_threshold)
  
  # Step 2: Filter high-frequency m/z
  cat("Step 2: High-Frequency M/Z Filter (by non-zero count)\n")
  cat("─────────────────────────────────────────────────────────────\n")
  object <- filter_high_frequency_mz_by_count(object, min_pixels = min_pixels)
  
  cat("╔═══════════════════════════════════════════════════════════════╗\n")
  cat("║              Supplementary Filtering Complete                 ║\n")
  cat("╚═══════════════════════════════════════════════════════════════╝\n\n")
  
  return(object)
}


# ════════════════════════════════════════════════════════════════════════════════
# 增强的 filter_msi_dataset 函数：支持新的筛选类型
# ════════════════════════════════════════════════════════════════════════════════

#' Extended filter_msi_dataset with Supplementary Filter Types
#'
#' @description
#' Extended version of filter_msi_dataset that supports two supplementary filter types:
#' - "high_frequency_count": Filter features by minimum non-zero pixel count
#' - "zero_heavy": Filter out features with > 60% zero values
#'
#' @param object A msi_dataset object
#' @param target Which component to filter ("variable_info")
#' @param filter_type The type of filter to apply
#' @param ... Additional arguments
#'
#' @return Filtered msi_dataset object
#' @export

filter_msi_dataset_supplementary <- function(object,
                                             target = "variable_info",
                                             filter_type = NULL,
                                             ...) {
  
  if (!methods::is(object, "msi_dataset"))
    stop("object must be a msi_dataset.")
  
  if (is.null(filter_type))
    stop("filter_type must be specified.")
  
  args <- list(...)
  
  switch(filter_type,
         
         # ────────────────────────────────────────────────────────────────────
         # 补充类型 1: 高频 M/Z 筛选（基于非零像素计数）
         # ────────────────────────────────────────────────────────────────────
         "high_frequency_count" = {
           min_pixels <- args$min_pixels %||% NULL
           object <- filter_high_frequency_mz_by_count(object, min_pixels = min_pixels)
         },
         
         # ────────────────────────────────────────────────────────────────────
         # 补充类型 2: 零值过多的 M/Z 筛选
         # ────────────────────────────────────────────────────────────────────
         "zero_heavy" = {
           zero_threshold <- args$zero_threshold %||% 0.6
           object <- filter_out_zero_heavy_mz(object, zero_threshold = zero_threshold)
         },
         
         # ────────────────────────────────────────────────────────────────────
         # 补充类型 3: 组合补充筛选
         # ────────────────────────────────────────────────────────────────────
         "supplementary" = {
           min_pixels <- args$min_pixels %||% NULL
           zero_threshold <- args$zero_threshold %||% 0.6
           object <- filter_supplementary(object, min_pixels, zero_threshold)
         },
         
         {
           stop(sprintf("Unknown supplementary filter_type: '%s'.
                   Available types are: 'high_frequency_count', 'zero_heavy', 'supplementary'",
                        filter_type))
         }
  )
  
  return(object)
}


# ════════════════════════════════════════════════════════════════════════════════
# 辅助函数：显示 M/Z 统计信息（不依赖 freq 列）
# ════════════════════════════════════════════════════════════════════════════════

#' Show M/Z Statistics (Calculated from Expression Data)
#'
#' @description
#' Display detailed statistics about m/z features by calculating from expression data.
#' No longer requires 'freq' column in variable_info.
#'
#' @param object A msi_dataset object
#'
#' @return Data frame with statistics (invisibly)
#' @export

show_mz_statistics <- function(object) {
  
  if (!methods::is(object, "msi_dataset"))
    stop("object must be a msi_dataset.")
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("                      M/Z Statistics\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  # Calculate statistics
  num_features <- nrow(object@variable_info)
  num_pixels <- ncol(object@expression_data)
  
  # Calculate non-zero counts for each feature
  non_zero_counts <- apply(object@expression_data, 1, function(x) {
    sum(x != 0, na.rm = TRUE)
  })
  
  non_zero_freq <- non_zero_counts / num_pixels
  
  cat("Non-Zero Frequency Statistics (calculated from expression data):\n")
  cat(sprintf("  Count:  %d features\n", num_features))
  cat(sprintf("  Min:    %.4f (%.1f%% of %d pixels)\n", 
              min(non_zero_freq), min(non_zero_freq)*100, num_pixels))
  cat(sprintf("  Q1:     %.4f (%.1f%% of %d pixels)\n", 
              quantile(non_zero_freq, 0.25), 
              quantile(non_zero_freq, 0.25)*100, num_pixels))
  cat(sprintf("  Median: %.4f (%.1f%% of %d pixels)\n", 
              median(non_zero_freq), median(non_zero_freq)*100, num_pixels))
  cat(sprintf("  Q3:     %.4f (%.1f%% of %d pixels)\n", 
              quantile(non_zero_freq, 0.75), 
              quantile(non_zero_freq, 0.75)*100, num_pixels))
  cat(sprintf("  Max:    %.4f (%.1f%% of %d pixels)\n\n", 
              max(non_zero_freq), max(non_zero_freq)*100, num_pixels))
  
  # Count features in different frequency ranges
  cat(sprintf("  Features with non-zero freq > 0.8: %d (%.1f%%)\n",
              sum(non_zero_freq > 0.8), 100*sum(non_zero_freq > 0.8)/num_features))
  cat(sprintf("  Features with non-zero freq > 0.5: %d (%.1f%%)\n",
              sum(non_zero_freq > 0.5), 100*sum(non_zero_freq > 0.5)/num_features))
  cat(sprintf("  Features with non-zero freq > 0.2: %d (%.1f%%)\n",
              sum(non_zero_freq > 0.2), 100*sum(non_zero_freq > 0.2)/num_features))
  cat(sprintf("  Features with non-zero freq > 0.1: %d (%.1f%%)\n\n",
              sum(non_zero_freq > 0.1), 100*sum(non_zero_freq > 0.1)/num_features))
  
  # Zero percentage statistics
  zero_percentages <- apply(object@expression_data, 1, function(x) {
    sum(x == 0, na.rm = TRUE) / length(x)
  }) * 100
  
  cat("Zero Percentage Statistics (% of zero values in each feature):\n")
  cat(sprintf("  Min:    %.2f%%\n", min(zero_percentages)))
  cat(sprintf("  Q1:     %.2f%%\n", quantile(zero_percentages, 0.25)))
  cat(sprintf("  Median: %.2f%%\n", median(zero_percentages)))
  cat(sprintf("  Q3:     %.2f%%\n", quantile(zero_percentages, 0.75)))
  cat(sprintf("  Max:    %.2f%%\n\n", max(zero_percentages)))
  
  # Count features in different zero percentage ranges
  cat(sprintf("  Features with < 20%% zeros: %d (%.1f%%)\n",
              sum(zero_percentages < 20), 100*sum(zero_percentages < 20)/num_features))
  cat(sprintf("  Features with < 40%% zeros: %d (%.1f%%)\n",
              sum(zero_percentages < 40), 100*sum(zero_percentages < 40)/num_features))
  cat(sprintf("  Features with < 60%% zeros: %d (%.1f%%)\n",
              sum(zero_percentages < 60), 100*sum(zero_percentages < 60)/num_features))
  cat(sprintf("  Features with < 80%% zeros: %d (%.1f%%)\n",
              sum(zero_percentages < 80), 100*sum(zero_percentages < 80)/num_features))
  cat(sprintf("  Features with > 60%% zeros: %d (%.1f%%)\n\n",
              sum(zero_percentages > 60), 100*sum(zero_percentages > 60)/num_features))
  
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  # Return statistics invisibly
  stats_df <- data.frame(
    metric = c("num_features", "num_pixels", 
               "non_zero_freq_min", "non_zero_freq_median", "non_zero_freq_max",
               "zero_pct_min", "zero_pct_median", "zero_pct_max"),
    value = c(num_features, num_pixels,
              min(non_zero_freq), median(non_zero_freq), max(non_zero_freq),
              min(zero_percentages), median(zero_percentages), max(zero_percentages))
  )
  
  invisible(stats_df)
}


# ════════════════════════════════════════════════════════════════════════════════
# 辅助函数：建议最佳筛选参数
# ════════════════════════════════════════════════════════════════════════════════

#' Recommend Optimal Filter Parameters
#'
#' @description
#' Analyze the data distribution and recommend suitable filter parameters.
#' Now based on non-zero counts instead of stored freq column.
#'
#' @param object A msi_dataset object
#'
#' @return List with recommended parameters
#' @export

recommend_filter_params <- function(object) {
  
  if (!methods::is(object, "msi_dataset"))
    stop("object must be a msi_dataset.")
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("              Recommended Filter Parameters\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  # Analyze non-zero count distribution
  non_zero_counts <- apply(object@expression_data, 1, function(x) {
    sum(x != 0, na.rm = TRUE)
  })
  
  num_pixels <- ncol(object@expression_data)
  non_zero_median <- median(non_zero_counts)
  non_zero_q3 <- quantile(non_zero_counts, 0.75)
  
  cat("Non-Zero Count Analysis:\n")
  cat(sprintf("  Median non-zero count: %d out of %d pixels (%.1f%%)\n", 
              non_zero_median, num_pixels, 100*non_zero_median/num_pixels))
  cat(sprintf("  Q3 non-zero count:     %d out of %d pixels (%.1f%%)\n", 
              non_zero_q3, num_pixels, 100*non_zero_q3/num_pixels))
  
  if (non_zero_median < 0.3 * num_pixels) {
    cat("  Recommendation: Data is sparse. Use min_pixels = ", 
        ceiling(0.2 * num_pixels), "\n\n")
    min_pixels_rec <- ceiling(0.2 * num_pixels)
  } else if (non_zero_median < 0.6 * num_pixels) {
    cat("  Recommendation: Data is moderately sparse. Use min_pixels = ", 
        ceiling(0.3 * num_pixels), "\n\n")
    min_pixels_rec <- ceiling(0.3 * num_pixels)
  } else {
    cat("  Recommendation: Data is dense. Use min_pixels = ", 
        ceiling(0.5 * num_pixels), "\n\n")
    min_pixels_rec <- ceiling(0.5 * num_pixels)
  }
  
  # Analyze zero distribution
  zero_percentages <- apply(object@expression_data, 1, function(x) {
    sum(x == 0, na.rm = TRUE) / length(x)
  })
  
  zero_median <- median(zero_percentages)
  zero_q3 <- quantile(zero_percentages, 0.75)
  
  cat("Zero Value Analysis:\n")
  cat(sprintf("  Median zero %%: %.1f%%\n", zero_median*100))
  cat(sprintf("  Q3 zero %%:     %.1f%%\n", zero_q3*100))
  
  if (zero_median > 0.7) {
    cat("  Recommendation: Many sparse features. Use zero_threshold = 0.7-0.8\n\n")
    zero_rec <- 0.7
  } else if (zero_median > 0.5) {
    cat("  Recommendation: Moderate sparsity. Use zero_threshold = 0.5-0.6\n\n")
    zero_rec <- 0.5
  } else {
    cat("  Recommendation: Low sparsity. Use zero_threshold = 0.4-0.5\n\n")
    zero_rec <- 0.4
  }
  
  cat("Suggested Usage:\n")
  cat(sprintf("  # Zero-heavy filter\n"))
  cat(sprintf("  msi <- filter_msi_dataset_supplementary(\n"))
  cat(sprintf("    msi, filter_type = 'zero_heavy', zero_threshold = %.1f)\n\n", zero_rec))
  
  cat(sprintf("  # High-frequency filter (by non-zero count)\n"))
  cat(sprintf("  msi <- filter_msi_dataset_supplementary(\n"))
  cat(sprintf("    msi, filter_type = 'high_frequency_count', min_pixels = %d)\n\n", 
              min_pixels_rec))
  
  cat(sprintf("  # Or both together\n"))
  cat(sprintf("  msi <- filter_msi_dataset_supplementary(\n"))
  cat(sprintf("    msi, filter_type = 'supplementary',\n"))
  cat(sprintf("    min_pixels = %d, zero_threshold = %.1f)\n\n", 
              min_pixels_rec, zero_rec))
  
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  return(invisible(list(
    min_pixels = min_pixels_rec,
    zero_threshold = zero_rec
  )))
}


# ════════════════════════════════════════════════════════════════════════════════
# Helper: %||% operator (null coalescing)
# ════════════════════════════════════════════════════════════════════════════════

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}