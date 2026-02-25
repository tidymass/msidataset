# ============================================================================
# filter_msi_dataset 补充筛选功能
# 补充 1: 筛选在超过80%像素都存在的m/z（高频特征）
# 补充 2: 排除掉超过60%的特征都是0的m/z（低质量特征）
# ============================================================================

# ════════════════════════════════════════════════════════════════════════════════
# 补充筛选函数 1: 高频 M/Z 筛选（freq > 0.8）
# ════════════════════════════════════════════════════════════════════════════════

#' Filter High-Frequency M/Z Features (> 80% pixels)
#'
#' @description
#' Retain only m/z features that are present in more than 80% of pixels.
#' This is useful for selecting robust, abundant features.
#'
#' @param object A msi_dataset object
#' @param freq_threshold Frequency threshold (default: 0.8)
#'
#' @return Filtered msi_dataset with high-frequency features only
#' @export
#'
#' @examples
#' \dontrun{
#' msi_high_freq <- filter_high_frequency_mz(msi)
#' }

filter_high_frequency_mz <- function(object, freq_threshold = 0.8) {
  
  if (!methods::is(object, "msi_dataset"))
    stop("object must be a msi_dataset.")
  
  if (!("freq" %in% colnames(object@variable_info))) {
    stop("freq column not found in variable_info.
         Please ensure variable_info has a 'freq' column.")
  }
  
  cat("═══════════════════════════════════════════════════════════════════\n")
  cat("Filtering High-Frequency M/Z Features (freq > ", freq_threshold, ")\n")
  cat("───────────────────────────────────────────────────────────────────\n\n")
  
  # Find features with frequency > threshold
  keep_idx <- object@variable_info$freq > freq_threshold
  keep_ids <- object@variable_info$variable_id[keep_idx]
  
  num_original <- nrow(object@variable_info)
  num_remaining <- sum(keep_idx)
  pct_retained <- 100 * num_remaining / num_original
  
  cat(sprintf("Frequency threshold: > %.1f%%\n", freq_threshold * 100))
  cat(sprintf("\nFeatures with frequency > %.1f%%: %d\n", freq_threshold * 100, num_remaining))
  cat(sprintf("Features removed: %d\n", num_original - num_remaining))
  cat(sprintf("Retention rate: %.1f%%\n\n", pct_retained))
  
  # Show distribution of frequencies
  freq_stats <- object@variable_info$freq[keep_idx]
  if (length(freq_stats) > 0) {
    cat(sprintf("Frequency statistics (remaining features):\n"))
    cat(sprintf("  Min:    %.4f (%.1f%%)\n", min(freq_stats), min(freq_stats)*100))
    cat(sprintf("  Q1:     %.4f (%.1f%%)\n", quantile(freq_stats, 0.25), 
                quantile(freq_stats, 0.25)*100))
    cat(sprintf("  Median: %.4f (%.1f%%)\n", median(freq_stats), median(freq_stats)*100))
    cat(sprintf("  Q3:     %.4f (%.1f%%)\n", quantile(freq_stats, 0.75), 
                quantile(freq_stats, 0.75)*100))
    cat(sprintf("  Max:    %.4f (%.1f%%)\n\n", max(freq_stats), max(freq_stats)*100))
  }
  
  # Update object
  new_expr <- object@expression_data[keep_ids, , drop = FALSE]
  new_vinfo <- object@variable_info[keep_idx, ]
  
  object@expression_data <- new_expr
  object@variable_info <- new_vinfo
  
  # Record in process_info
  process_key <- paste0("filter_high_freq_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "filter_high_frequency_mz()",
    parameter = list(freq_threshold = freq_threshold,
                     num_retained = num_remaining,
                     num_removed = num_original - num_remaining),
    time = Sys.time()
  )
  object@process_info[[process_key]] <- parameter
  
  cat("═══════════════════════════════════════════════════════════════════\n\n")
  
  return(object)
}


# ════════════════════════════════════════════════════════════════════════════════
# 补充筛选函数 2: 排除零值过多的 M/Z（> 60% 零值）
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
  
  cat("═══════════════════════════════════════════════════════════════════\n")
  cat("Filtering Out Zero-Heavy M/Z Features (zeros > ", zero_threshold*100, "%)\n")
  cat("───────────────────────────────────────────────────────────────────\n\n")
  
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
  
  cat("═══════════════════════════════════════════════════════════════════\n\n")
  
  return(object)
}


# ════════════════════════════════════════════════════════════════════════════════
# 组合函数：同时应用两种补充筛选
# ════════════════════════════════════════════════════════════════════════════════

#' Apply Both Supplementary Filters
#'
#' @description
#' Apply both high-frequency and zero-heavy filters in sequence.
#' This is a convenience function combining the two supplementary filters.
#'
#' @param object A msi_dataset object
#' @param freq_threshold Frequency threshold for high-frequency filter (default: 0.8)
#' @param zero_threshold Zero percentage threshold (default: 0.6)
#'
#' @return Filtered msi_dataset
#' @export

filter_supplementary <- function(object, freq_threshold = 0.8, zero_threshold = 0.6) {
  
  if (!methods::is(object, "msi_dataset"))
    stop("object must be a msi_dataset.")
  
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════════════╗\n")
  cat("║              Applying Supplementary Filters                       ║\n")
  cat("╚═══════════════════════════════════════════════════════════════════╝\n\n")
  
  # Step 1: Filter high-frequency m/z
  cat("Step 1: High-Frequency M/Z Filter\n")
  cat("─────────────────────────────────────────────────────────────────\n")
  object <- filter_high_frequency_mz(object, freq_threshold = freq_threshold)
  
  # Step 2: Filter zero-heavy m/z
  cat("Step 2: Zero-Heavy M/Z Filter\n")
  cat("─────────────────────────────────────────────────────────────────\n")
  object <- filter_out_zero_heavy_mz(object, zero_threshold = zero_threshold)
  
  cat("╔═══════════════════════════════════════════════════════════════════╗\n")
  cat("║              Supplementary Filtering Complete                     ║\n")
  cat("╚═══════════════════════════════════════════════════════════════════╝\n\n")
  
  return(object)
}


# ════════════════════════════════════════════════════════════════════════════════
# 增强的 filter_msi_dataset 函数：支持新的筛选类型
# ════════════════════════════════════════════════════════════════════════════════

#' Enhanced filter_msi_dataset with Supplementary Filter Types
#'
#' @description
#' Extended version of filter_msi_dataset that supports two new filter types:
#' - "high_frequency": Filter features present in > 80% of pixels
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
         # 补充类型 1: 高频 M/Z 筛选
         # ────────────────────────────────────────────────────────────────────
         "high_frequency" = {
           freq_threshold <- args$freq_threshold %||% 0.8
           object <- filter_high_frequency_mz(object, freq_threshold = freq_threshold)
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
           freq_threshold <- args$freq_threshold %||% 0.8
           zero_threshold <- args$zero_threshold %||% 0.6
           object <- filter_supplementary(object, freq_threshold, zero_threshold)
         },
         
         {
           stop(sprintf("Unknown supplementary filter_type: '%s'.
                   Available types are: 'high_frequency', 'zero_heavy', 'supplementary'",
                        filter_type))
         }
  )
  
  return(object)
}


# ════════════════════════════════════════════════════════════════════════════════
# 辅助函数：显示 M/Z 统计信息
# ════════════════════════════════════════════════════════════════════════════════

#' Show M/Z Statistics Before and After Filtering
#'
#' @description
#' Display detailed statistics about m/z features, including frequency and zero distributions.
#'
#' @param object A msi_dataset object
#' @param show_distribution Logical, whether to show distribution plots
#'
#' @return Data frame with statistics (invisibly)
#' @export

show_mz_statistics <- function(object, show_distribution = TRUE) {
  
  if (!methods::is(object, "msi_dataset"))
    stop("object must be a msi_dataset.")
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════════\n")
  cat("                      M/Z Statistics\n")
  cat("═══════════════════════════════════════════════════════════════════\n\n")
  
  # Calculate statistics
  num_features <- nrow(object@variable_info)
  num_pixels <- ncol(object@expression_data)
  
  # Frequency statistics
  if ("freq" %in% colnames(object@variable_info)) {
    freqs <- object@variable_info$freq
    cat("Frequency Statistics (proportion of pixels where m/z is present):\n")
    cat(sprintf("  Count:  %d\n", num_features))
    cat(sprintf("  Min:    %.4f (%.1f%%)\n", min(freqs), min(freqs)*100))
    cat(sprintf("  Q1:     %.4f (%.1f%%)\n", quantile(freqs, 0.25), 
                quantile(freqs, 0.25)*100))
    cat(sprintf("  Median: %.4f (%.1f%%)\n", median(freqs), median(freqs)*100))
    cat(sprintf("  Q3:     %.4f (%.1f%%)\n", quantile(freqs, 0.75), 
                quantile(freqs, 0.75)*100))
    cat(sprintf("  Max:    %.4f (%.1f%%)\n\n", max(freqs), max(freqs)*100))
    
    # Count features in different frequency ranges
    cat(sprintf("  Features with freq > 0.8: %d (%.1f%%)\n",
                sum(freqs > 0.8), 100*sum(freqs > 0.8)/num_features))
    cat(sprintf("  Features with freq > 0.5: %d (%.1f%%)\n",
                sum(freqs > 0.5), 100*sum(freqs > 0.5)/num_features))
    cat(sprintf("  Features with freq > 0.2: %d (%.1f%%)\n",
                sum(freqs > 0.2), 100*sum(freqs > 0.2)/num_features))
    cat(sprintf("  Features with freq > 0.1: %d (%.1f%%)\n\n",
                sum(freqs > 0.1), 100*sum(freqs > 0.1)/num_features))
  }
  
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
  
  cat("═══════════════════════════════════════════════════════════════════\n\n")
  
  # Return statistics invisibly
  stats_df <- data.frame(
    metric = c("num_features", "num_pixels", "freq_min", "freq_median", "freq_max",
               "zero_min", "zero_median", "zero_max"),
    value = c(num_features, num_pixels,
              min(freqs), median(freqs), max(freqs),
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
#'
#' @param object A msi_dataset object
#'
#' @return List with recommended parameters
#' @export

recommend_filter_params <- function(object) {
  
  if (!methods::is(object, "msi_dataset"))
    stop("object must be a msi_dataset.")
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════════\n")
  cat("              Recommended Filter Parameters\n")
  cat("═══════════════════════════════════════════════════════════════════\n\n")
  
  # Analyze frequency distribution
  if ("freq" %in% colnames(object@variable_info)) {
    freqs <- object@variable_info$freq
    freq_median <- median(freqs)
    freq_q3 <- quantile(freqs, 0.75)
    
    cat("Frequency Analysis:\n")
    cat(sprintf("  Median frequency: %.3f (%.1f%%)\n", freq_median, freq_median*100))
    cat(sprintf("  Q3 frequency:     %.3f (%.1f%%)\n", freq_q3, freq_q3*100))
    
    if (freq_median < 0.3) {
      cat("  Recommendation: Data is sparse. Use freq_threshold = 0.2-0.3\n\n")
      freq_rec <- 0.2
    } else if (freq_median < 0.6) {
      cat("  Recommendation: Data is moderately sparse. Use freq_threshold = 0.3-0.5\n\n")
      freq_rec <- 0.3
    } else {
      cat("  Recommendation: Data is dense. Use freq_threshold = 0.5-0.8\n\n")
      freq_rec <- 0.5
    }
  } else {
    freq_rec <- 0.2
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
  cat(sprintf("  # High-frequency filter\n"))
  cat(sprintf("  msi <- filter_msi_dataset_supplementary(\n"))
  cat(sprintf("    msi, filter_type = 'high_frequency', freq_threshold = %.1f)\n\n", freq_rec))
  
  cat(sprintf("  # Zero-heavy filter\n"))
  cat(sprintf("  msi <- filter_msi_dataset_supplementary(\n"))
  cat(sprintf("    msi, filter_type = 'zero_heavy', zero_threshold = %.1f)\n\n", zero_rec))
  
  cat(sprintf("  # Or both together\n"))
  cat(sprintf("  msi <- filter_msi_dataset_supplementary(\n"))
  cat(sprintf("    msi, filter_type = 'supplementary',\n"))
  cat(sprintf("    freq_threshold = %.1f, zero_threshold = %.1f)\n\n", freq_rec, zero_rec))
  
  cat("═══════════════════════════════════════════════════════════════════\n\n")
  
  return(invisible(list(
    freq_threshold = freq_rec,
    zero_threshold = zero_rec
  )))
}