#' Filter a msi_dataset Object with Multiple Filter Options
#'
#' @description
#' A comprehensive filtering function for msi_dataset objects that provides
#' multiple filtering strategies for both samples (pixels) and variables (m/z features).
#'
#' @param object A msi_dataset object to be filtered
#' @param target Which component to filter. One of "sample_info" (default) or "variable_info"
#' @param filter_type The type of filter to apply
#' @param ... Additional arguments passed to specific filter functions
#'
#' @return A filtered msi_dataset object with updated process_info
#' @export
#'
#' @examples
#' \dontrun{
#' # Filter samples by spatial coordinates
#' obj <- filter_msi_dataset(obj, target = "sample_info", 
#'                           filter_type = "spatial",
#'                           x_min = 2, x_max = 8, y_min = 1, y_max = 5)
#' 
#' # Filter variables by m/z range
#' obj <- filter_msi_dataset(obj, target = "variable_info",
#'                           filter_type = "mz_range",
#'                           mz_min = 100, mz_max = 1000)
#' }

filter_msi_dataset <- function(object,
                               target = "sample_info",
                               filter_type = NULL,
                               ...) {
  
  if (!methods::is(object, "msi_dataset"))
    stop("object must be a msi_dataset.")
  
  if (!(target %in% c("sample_info", "variable_info")))
    stop("target must be either 'sample_info' or 'variable_info'.")
  
  if (is.null(filter_type))
    stop("filter_type must be specified.")
  
  # Route to appropriate filter function based on target
  if (target == "sample_info") {
    object <- filter_samples_comprehensive(object, filter_type = filter_type, ...)
  } else if (target == "variable_info") {
    object <- filter_variables_comprehensive(object, filter_type = filter_type, ...)
  }
  
  return(object)
}


# ============================================================================
# SAMPLE (PIXEL) FILTERING FUNCTIONS
# ============================================================================

#' Comprehensive Sample Filtering Function
#'
#' @description
#' Provides multiple filtering strategies for samples/pixels:
#' - "spatial": Filter by spatial coordinates (x, y, z)
#' - "section": Filter by section_id
#' - "run": Filter by experimental run/batch
#' - "intensity": Filter by signal intensity (median, mean)
#' - "custom": Apply custom filter function
#'
#' @keywords internal

filter_samples_comprehensive <- function(object,
                                         filter_type = "spatial",
                                         ...) {
  
  args <- list(...)
  
  switch(filter_type,
         
         # ────────────────────────────────────────────────────────────────────
         # Type 1: Spatial filtering
         # ────────────────────────────────────────────────────────────────────
         "spatial" = {
           cat("Applying spatial filtering...\n")
           object <- filter_spatial(object, args = args)
         },
         
         # ────────────────────────────────────────────────────────────────────
         # Type 2: Section filtering
         # ────────────────────────────────────────────────────────────────────
         "section" = {
           cat("Applying section filtering...\n")
           object <- filter_by_section(object, args = args)
         },
         
         # ────────────────────────────────────────────────────────────────────
         # Type 3: Run/Batch filtering
         # ────────────────────────────────────────────────────────────────────
         "run" = {
           cat("Applying run/batch filtering...\n")
           object <- filter_by_run(object, args = args)
         },
         
         # ────────────────────────────────────────────────────────────────────
         # Type 4: Intensity-based filtering
         # ────────────────────────────────────────────────────────────────────
         "intensity" = {
           cat("Applying intensity-based filtering...\n")
           object <- filter_by_intensity(object, args = args)
         },
         
         # ────────────────────────────────────────────────────────────────────
         # Type 5: Combined filtering
         # ────────────────────────────────────────────────────────────────────
         "combined" = {
           cat("Applying combined filtering...\n")
           object <- filter_samples_combined(object, args = args)
         },
         
         # ────────────────────────────────────────────────────────────────────
         # Type 6: Custom filter function
         # ────────────────────────────────────────────────────────────────────
         "custom" = {
           cat("Applying custom filter function...\n")
           object <- filter_by_custom_function(object, args = args)
         },
         
         {
           stop(sprintf("Unknown filter_type: '%s'. Available types are:
                   'spatial', 'section', 'run', 'intensity', 'combined', 'custom'", 
                        filter_type))
         }
  )
  
  return(object)
}


# ────────────────────────────────────────────────────────────────────────────
# Sample Filter 1: Spatial Filtering
# ────────────────────────────────────────────────────────────────────────────

#' @keywords internal
filter_spatial <- function(object, args) {
  
  # Extract parameters with defaults
  x_min <- args$x_min %||% min(object@sample_info$x, na.rm = TRUE)
  x_max <- args$x_max %||% max(object@sample_info$x, na.rm = TRUE)
  y_min <- args$y_min %||% min(object@sample_info$y, na.rm = TRUE)
  y_max <- args$y_max %||% max(object@sample_info$y, na.rm = TRUE)
  
  # Optional z filtering if available
  z_min <- args$z_min %||% NULL
  z_max <- args$z_max %||% NULL
  
  # Apply filtering
  keep_idx <- object@sample_info$x >= x_min & 
    object@sample_info$x <= x_max &
    object@sample_info$y >= y_min & 
    object@sample_info$y <= y_max
  
  # Filter z if specified
  if (!is.null(z_min) && !is.null(z_max)) {
    if ("z" %in% colnames(object@sample_info)) {
      keep_idx <- keep_idx & 
        object@sample_info$z >= z_min & 
        object@sample_info$z <= z_max
    }
  }
  
  keep_ids <- object@sample_info$sample_id[keep_idx]
  
  cat(sprintf("  Spatial bounds: x[%.1f-%.1f], y[%.1f-%.1f]\n",
              x_min, x_max, y_min, y_max))
  if (!is.null(z_min) && !is.null(z_max)) {
    cat(sprintf("                   z[%.1f-%.1f]\n", z_min, z_max))
  }
  cat(sprintf("  Pixels remaining: %d / %d (%.1f%%)\n",
              length(keep_ids), ncol(object@expression_data),
              100 * length(keep_ids) / ncol(object@expression_data)))
  
  # Update object
  new_expr <- object@expression_data[, keep_ids, drop = FALSE]
  new_sinfo <- object@sample_info[keep_idx, ]
  
  object@expression_data <- new_expr
  object@sample_info <- new_sinfo
  
  # Record in process_info
  process_key <- paste0("filter_spatial_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "filter_msi_dataset(spatial)",
    parameter = list(x_min = x_min, x_max = x_max, 
                     y_min = y_min, y_max = y_max,
                     z_min = z_min, z_max = z_max),
    time = Sys.time()
  )
  object@process_info[[process_key]] <- parameter
  
  return(object)
}


# ────────────────────────────────────────────────────────────────────────────
# Sample Filter 2: Section Filtering
# ────────────────────────────────────────────────────────────────────────────

#' @keywords internal
filter_by_section <- function(object, args) {
  
  if (!("section_id" %in% colnames(object@sample_info))) {
    stop("section_id column not found in sample_info")
  }
  
  section_ids <- args$section_ids
  
  if (is.null(section_ids)) {
    stop("section_ids must be provided. Example: section_ids = c('slice_01', 'slice_02')")
  }
  
  if (is.character(section_ids)) {
    section_ids <- list(section_ids)
  }
  
  keep_idx <- object@sample_info$section_id %in% unlist(section_ids)
  keep_ids <- object@sample_info$sample_id[keep_idx]
  
  cat(sprintf("  Sections: %s\n", paste(unlist(section_ids), collapse = ", ")))
  cat(sprintf("  Pixels remaining: %d / %d (%.1f%%)\n",
              length(keep_ids), ncol(object@expression_data),
              100 * length(keep_ids) / ncol(object@expression_data)))
  
  # Update object
  new_expr <- object@expression_data[, keep_ids, drop = FALSE]
  new_sinfo <- object@sample_info[keep_idx, ]
  
  object@expression_data <- new_expr
  object@sample_info <- new_sinfo
  
  # Record in process_info
  process_key <- paste0("filter_section_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "filter_msi_dataset(section)",
    parameter = list(section_ids = unlist(section_ids)),
    time = Sys.time()
  )
  object@process_info[[process_key]] <- parameter
  
  return(object)
}


# ────────────────────────────────────────────────────────────────────────────
# Sample Filter 3: Run/Batch Filtering
# ────────────────────────────────────────────────────────────────────────────

#' @keywords internal
filter_by_run <- function(object, args) {
  
  if (!("run" %in% colnames(object@sample_info))) {
    stop("run column not found in sample_info")
  }
  
  run_ids <- args$run_ids
  
  if (is.null(run_ids)) {
    stop("run_ids must be provided. Example: run_ids = c('run1', 'run2')")
  }
  
  keep_idx <- object@sample_info$run %in% run_ids
  keep_ids <- object@sample_info$sample_id[keep_idx]
  
  cat(sprintf("  Runs: %s\n", paste(run_ids, collapse = ", ")))
  cat(sprintf("  Pixels remaining: %d / %d (%.1f%%)\n",
              length(keep_ids), ncol(object@expression_data),
              100 * length(keep_ids) / ncol(object@expression_data)))
  
  # Update object
  new_expr <- object@expression_data[, keep_ids, drop = FALSE]
  new_sinfo <- object@sample_info[keep_idx, ]
  
  object@expression_data <- new_expr
  object@sample_info <- new_sinfo
  
  # Record in process_info
  process_key <- paste0("filter_run_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "filter_msi_dataset(run)",
    parameter = list(run_ids = run_ids),
    time = Sys.time()
  )
  object@process_info[[process_key]] <- parameter
  
  return(object)
}


# ────────────────────────────────────────────────────────────────────────────
# Sample Filter 4: Intensity-based Filtering
# ────────────────────────────────────────────────────────────────────────────

#' @keywords internal
filter_by_intensity <- function(object, args) {
  
  intensity_type <- args$intensity_type %||% "median"  # median or mean
  intensity_threshold <- args$intensity_threshold %||% 0
  
  # Calculate intensity for each sample
  if (intensity_type == "median") {
    intensities <- apply(object@expression_data, 2, median, na.rm = TRUE)
  } else if (intensity_type == "mean") {
    intensities <- apply(object@expression_data, 2, mean, na.rm = TRUE)
  } else {
    stop("intensity_type must be 'median' or 'mean'")
  }
  
  keep_idx <- intensities >= intensity_threshold
  keep_ids <- object@sample_info$sample_id[keep_idx]
  
  cat(sprintf("  Intensity type: %s\n", intensity_type))
  cat(sprintf("  Threshold: >= %.2f\n", intensity_threshold))
  cat(sprintf("  Pixels remaining: %d / %d (%.1f%%)\n",
              length(keep_ids), ncol(object@expression_data),
              100 * length(keep_ids) / ncol(object@expression_data)))
  
  # Update object
  new_expr <- object@expression_data[, keep_idx, drop = FALSE]
  new_sinfo <- object@sample_info[keep_idx, ]
  
  object@expression_data <- new_expr
  object@sample_info <- new_sinfo
  
  # Record in process_info
  process_key <- paste0("filter_intensity_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "filter_msi_dataset(intensity)",
    parameter = list(intensity_type = intensity_type, 
                     intensity_threshold = intensity_threshold),
    time = Sys.time()
  )
  object@process_info[[process_key]] <- parameter
  
  return(object)
}


# ────────────────────────────────────────────────────────────────────────────
# Sample Filter 5: Combined Filtering
# ────────────────────────────────────────────────────────────────────────────

#' @keywords internal
filter_samples_combined <- function(object, args) {
  
  cat("  Applying combined filters in sequence...\n")
  
  # Step 1: Spatial filtering
  if (!is.null(args$spatial)) {
    cat("  Step 1: Spatial filtering\n")
    spatial_args <- c(list(), args$spatial)
    object <- filter_spatial(object, args = spatial_args)
  }
  
  # Step 2: Section filtering
  if (!is.null(args$section)) {
    cat("  Step 2: Section filtering\n")
    section_args <- c(list(), args$section)
    object <- filter_by_section(object, args = section_args)
  }
  
  # Step 3: Run filtering
  if (!is.null(args$run)) {
    cat("  Step 3: Run filtering\n")
    run_args <- c(list(), args$run)
    object <- filter_by_run(object, args = run_args)
  }
  
  # Step 4: Intensity filtering
  if (!is.null(args$intensity)) {
    cat("  Step 4: Intensity filtering\n")
    intensity_args <- c(list(), args$intensity)
    object <- filter_by_intensity(object, args = intensity_args)
  }
  
  return(object)
}


# ────────────────────────────────────────────────────────────────────────────
# Sample Filter 6: Custom Function Filtering
# ────────────────────────────────────────────────────────────────────────────

#' @keywords internal
filter_by_custom_function <- function(object, args) {
  
  filter_func <- args$filter_function
  
  if (is.null(filter_func)) {
    stop("filter_function must be provided")
  }
  
  if (!is.function(filter_func)) {
    stop("filter_function must be a function")
  }
  
  # Apply function to each sample
  keep_idx <- sapply(1:ncol(object@expression_data), function(i) {
    filter_func(object@expression_data[, i])
  })
  
  keep_ids <- object@sample_info$sample_id[keep_idx]
  
  cat(sprintf("  Custom function applied\n"))
  cat(sprintf("  Pixels remaining: %d / %d (%.1f%%)\n",
              length(keep_ids), ncol(object@expression_data),
              100 * length(keep_ids) / ncol(object@expression_data)))
  
  # Update object
  new_expr <- object@expression_data[, keep_idx, drop = FALSE]
  new_sinfo <- object@sample_info[keep_idx, ]
  
  object@expression_data <- new_expr
  object@sample_info <- new_sinfo
  
  # Record in process_info
  process_key <- paste0("filter_custom_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "filter_msi_dataset(custom)",
    parameter = list(filter_function = "user_defined"),
    time = Sys.time()
  )
  object@process_info[[process_key]] <- parameter
  
  return(object)
}


# ============================================================================
# VARIABLE (M/Z FEATURE) FILTERING FUNCTIONS
# ============================================================================

#' Comprehensive Variable Filtering Function
#'
#' @description
#' Provides multiple filtering strategies for variables/m/z features:
#' - "mz_range": Filter by m/z range
#' - "count": Filter by count (number of pixels)
#' - "intensity": Filter by feature intensity (mean, median)
#' - "custom": Apply custom filter function
#'
#' Note: "frequency" filtering has been removed as freq column is no longer in variable_info
#'
#' @keywords internal

filter_variables_comprehensive <- function(object,
                                           filter_type = "mz_range",
                                           ...) {
  
  args <- list(...)
  
  switch(filter_type,
         
         # ────────────────────────────────────────────────────────────────────
         # Type 1: M/Z Range filtering
         # ────────────────────────────────────────────────────────────────────
         "mz_range" = {
           cat("Applying m/z range filtering...\n")
           object <- filter_mz_range(object, args = args)
         },
         
         # ────────────────────────────────────────────────────────────────────
         # Type 2: Count filtering
         # ────────────────────────────────────────────────────────────────────
         "count" = {
           cat("Applying count filtering...\n")
           object <- filter_by_count(object, args = args)
         },
         
         # ────────────────────────────────────────────────────────────────────
         # Type 3: Intensity-based filtering
         # ────────────────────────────────────────────────────────────────────
         "intensity" = {
           cat("Applying feature intensity filtering...\n")
           object <- filter_variable_intensity(object, args = args)
         },
         
         # ────────────────────────────────────────────────────────────────────
         # Type 4: Combined filtering
         # ────────────────────────────────────────────────────────────────────
         "combined" = {
           cat("Applying combined variable filtering...\n")
           object <- filter_variables_combined(object, args = args)
         },
         
         # ────────────────────────────────────────────────────────────────────
         # Type 5: Custom filter function
         # ────────────────────────────────────────────────────────────────────
         "custom" = {
           cat("Applying custom filter function...\n")
           object <- filter_variable_by_custom(object, args = args)
         },
         
         {
           stop(sprintf("Unknown filter_type: '%s'. Available types are:
                   'mz_range', 'count', 'intensity', 'combined', 'custom'
                   Note: 'frequency' filtering has been removed as freq column is no longer in variable_info", 
                        filter_type))
         }
  )
  
  return(object)
}


# ────────────────────────────────────────────────────────────────────────────
# Variable Filter 1: M/Z Range Filtering
# ────────────────────────────────────────────────────────────────────────────

#' @keywords internal
filter_mz_range <- function(object, args) {
  
  mz_min <- args$mz_min %||% min(object@variable_info$mz, na.rm = TRUE)
  mz_max <- args$mz_max %||% max(object@variable_info$mz, na.rm = TRUE)
  
  keep_idx <- object@variable_info$mz >= mz_min & 
    object@variable_info$mz <= mz_max
  keep_ids <- object@variable_info$variable_id[keep_idx]
  
  cat(sprintf("  M/Z range: [%.4f - %.4f]\n", mz_min, mz_max))
  cat(sprintf("  Variables remaining: %d / %d (%.1f%%)\n",
              length(keep_ids), nrow(object@variable_info),
              100 * length(keep_ids) / nrow(object@variable_info)))
  
  # Update object
  new_expr <- object@expression_data[keep_ids, , drop = FALSE]
  new_vinfo <- object@variable_info[keep_idx, ]
  
  object@expression_data <- new_expr
  object@variable_info <- new_vinfo
  
  # Record in process_info
  process_key <- paste0("filter_mz_range_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "filter_msi_dataset(mz_range)",
    parameter = list(mz_min = mz_min, mz_max = mz_max),
    time = Sys.time()
  )
  object@process_info[[process_key]] <- parameter
  
  return(object)
}


# ────────────────────────────────────────────────────────────────────────────
# Variable Filter 2: Count Filtering
# ────────────────────────────────────────────────────────────────────────────

#' @keywords internal
filter_by_count <- function(object, args) {
  
  if (!("count" %in% colnames(object@variable_info))) {
    stop("count column not found in variable_info.
         If you need to filter by frequency/count, add these columns first.")
  }
  
  count_min <- args$count_min %||% 0
  count_max <- args$count_max %||% ncol(object@expression_data)
  
  keep_idx <- object@variable_info$count >= count_min & 
    object@variable_info$count <= count_max
  keep_ids <- object@variable_info$variable_id[keep_idx]
  
  cat(sprintf("  Count range: [%d - %d]\n", count_min, count_max))
  cat(sprintf("  Variables remaining: %d / %d (%.1f%%)\n",
              length(keep_ids), nrow(object@variable_info),
              100 * length(keep_ids) / nrow(object@variable_info)))
  
  # Update object
  new_expr <- object@expression_data[keep_ids, , drop = FALSE]
  new_vinfo <- object@variable_info[keep_idx, ]
  
  object@expression_data <- new_expr
  object@variable_info <- new_vinfo
  
  # Record in process_info
  process_key <- paste0("filter_count_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "filter_msi_dataset(count)",
    parameter = list(count_min = count_min, count_max = count_max),
    time = Sys.time()
  )
  object@process_info[[process_key]] <- parameter
  
  return(object)
}


# ────────────────────────────────────────────────────────────────────────────
# Variable Filter 3: Feature Intensity Filtering
# ────────────────────────────────────────────────────────────────────────────

#' @keywords internal
filter_variable_intensity <- function(object, args) {
  
  intensity_type <- args$intensity_type %||% "median"
  intensity_threshold <- args$intensity_threshold %||% 0
  
  # Calculate intensity for each feature
  if (intensity_type == "median") {
    intensities <- apply(object@expression_data, 1, median, na.rm = TRUE)
  } else if (intensity_type == "mean") {
    intensities <- apply(object@expression_data, 1, mean, na.rm = TRUE)
  } else if (intensity_type == "max") {
    intensities <- apply(object@expression_data, 1, max, na.rm = TRUE)
  } else {
    stop("intensity_type must be 'median', 'mean', or 'max'")
  }
  
  keep_idx <- intensities >= intensity_threshold
  keep_ids <- object@variable_info$variable_id[keep_idx]
  
  cat(sprintf("  Intensity type: %s\n", intensity_type))
  cat(sprintf("  Threshold: >= %.2f\n", intensity_threshold))
  cat(sprintf("  Variables remaining: %d / %d (%.1f%%)\n",
              length(keep_ids), nrow(object@variable_info),
              100 * length(keep_ids) / nrow(object@variable_info)))
  
  # Update object
  new_expr <- object@expression_data[keep_idx, , drop = FALSE]
  new_vinfo <- object@variable_info[keep_idx, ]
  
  object@expression_data <- new_expr
  object@variable_info <- new_vinfo
  
  # Record in process_info
  process_key <- paste0("filter_var_intensity_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "filter_msi_dataset(intensity)",
    parameter = list(intensity_type = intensity_type,
                     intensity_threshold = intensity_threshold),
    time = Sys.time()
  )
  object@process_info[[process_key]] <- parameter
  
  return(object)
}


# ────────────────────────────────────────────────────────────────────────────
# Variable Filter 4: Combined Filtering
# ────────────────────────────────────────────────────────────────────────────

#' @keywords internal
filter_variables_combined <- function(object, args) {
  
  cat("  Applying combined filters in sequence...\n")
  
  # Step 1: M/Z range filtering
  if (!is.null(args$mz_range)) {
    cat("  Step 1: M/Z range filtering\n")
    mz_args <- c(list(), args$mz_range)
    object <- filter_mz_range(object, args = mz_args)
  }
  
  # Step 2: Count filtering
  if (!is.null(args$count)) {
    cat("  Step 2: Count filtering\n")
    count_args <- c(list(), args$count)
    object <- filter_by_count(object, args = count_args)
  }
  
  # Step 3: Intensity filtering
  if (!is.null(args$intensity)) {
    cat("  Step 3: Intensity filtering\n")
    intensity_args <- c(list(), args$intensity)
    object <- filter_variable_intensity(object, args = intensity_args)
  }
  
  return(object)
}


# ────────────────────────────────────────────────────────────────────────────
# Variable Filter 5: Custom Function Filtering
# ────────────────────────────────────────────────────────────────────────────

#' @keywords internal
filter_variable_by_custom <- function(object, args) {
  
  filter_func <- args$filter_function
  
  if (is.null(filter_func)) {
    stop("filter_function must be provided")
  }
  
  if (!is.function(filter_func)) {
    stop("filter_function must be a function")
  }
  
  # Apply function to each variable
  keep_idx <- sapply(1:nrow(object@variable_info), function(i) {
    filter_func(object@expression_data[i, ])
  })
  
  keep_ids <- object@variable_info$variable_id[keep_idx]
  
  cat(sprintf("  Custom function applied\n"))
  cat(sprintf("  Variables remaining: %d / %d (%.1f%%)\n",
              length(keep_ids), nrow(object@variable_info),
              100 * length(keep_ids) / nrow(object@variable_info)))
  
  # Update object
  new_expr <- object@expression_data[keep_idx, , drop = FALSE]
  new_vinfo <- object@variable_info[keep_idx, ]
  
  object@expression_data <- new_expr
  object@variable_info <- new_vinfo
  
  # Record in process_info
  process_key <- paste0("filter_var_custom_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "filter_msi_dataset(custom)",
    parameter = list(filter_function = "user_defined"),
    time = Sys.time()
  )
  object@process_info[[process_key]] <- parameter
  
  return(object)
}


# ============================================================================
# Helper: %||% operator (null coalescing)
# ============================================================================

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


# ============================================================================
# Additional Convenience Functions
# ============================================================================

#' Print Available Filter Types
#' @export
show_filter_types <- function() {
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("              Available Filter Types\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  cat("FOR SAMPLES (target = 'sample_info'):\n")
  cat("  1. 'spatial'    - Filter by spatial coordinates (x, y, z)\n")
  cat("  2. 'section'    - Filter by section_id\n")
  cat("  3. 'run'        - Filter by experimental run/batch\n")
  cat("  4. 'intensity'  - Filter by signal intensity (median, mean)\n")
  cat("  5. 'combined'   - Apply multiple filters in sequence\n")
  cat("  6. 'custom'     - Apply custom filter function\n\n")
  
  cat("FOR VARIABLES (target = 'variable_info'):\n")
  cat("  1. 'mz_range'   - Filter by m/z range\n")
  cat("  2. 'count'      - Filter by count (number of pixels with value)\n")
  cat("  3. 'intensity'  - Filter by feature intensity\n")
  cat("  4. 'combined'   - Apply multiple filters in sequence\n")
  cat("  5. 'custom'     - Apply custom filter function\n\n")
  
  cat("NOTE: 'frequency' filtering has been removed as the freq column\n")
  cat("is no longer in the standard variable_info structure.\n\n")
  
  cat("═══════════════════════════════════════════════════════════════\n")
}