# ============================================================================
# Core S4 Methods for msi_dataset (show, [, dim, dimnames, update)
# Author: Xiaoyu Xu <XIAOYU021@e.ntu.edu.sg>
# ============================================================================

# ---------------------------------------------------------------------------
# show()
# ---------------------------------------------------------------------------

#' Display a msi_dataset object
#'
#' @description
#' Prints a formatted summary of the object including dimensions, spatial
#' extent, section information, and processing history.
#'
#' @param object A \code{msi_dataset} object.
#' @return Invisibly returns the object.
#'
#' @docType methods
#' @name show,msi_dataset-method
#' @rdname show-methods
#' @aliases show,msi_dataset-method
#' @export

setMethod(
  f          = "show",
  signature  = "msi_dataset",
  definition = function(object) {
    
    ## ---- re-validate ----
    check_result <- check_msi_dataset_class(object)
    if (!isTRUE(check_result)) {
      cat(crayon::red("Validation errors:\n"))
      cat(crayon::red(paste(check_result, collapse = "\n"), "\n"))
      cat(crayon::red(
        "Slots may have been manually modified. ",
        "Try update_msi_dataset().\n"
      ))
    }
    
    sep <- paste(rep("-", 20), collapse = "")
    cat(crayon::yellow(sep, "\n"))
    cat(crayon::green("msidataset version:", object@version, "\n"))
    cat(crayon::yellow("Structure: massdataset-compatible\n"))
    cat(crayon::yellow("Activated slot:", object@activated, "\n"))
    cat(crayon::yellow(sep, "\n"))
    
    ## 1. expression_data
    cat(crayon::green("1. expression_data: "))
    data_type <- if (methods::is(object@expression_data, "DelayedArray")) {
      "DelayedArray"
    } else if (is.matrix(object@expression_data)) {
      "matrix"
    } else {
      "data.frame"
    }
    cat("[", nrow(object@expression_data), "x",
        ncol(object@expression_data), data_type, "]\n")
    
    ## 2. variable_info
    cat(crayon::green("2. variable_info: "))
    cat("[", nrow(object@variable_info), "x",
        ncol(object@variable_info), "data.frame]\n")
    nv <- nrow(object@variable_info)
    vid <- if (nv <= 5) {
      object@variable_info$variable_id
    } else {
      c(head(object@variable_info$variable_id, 3), "...",
        tail(object@variable_info$variable_id, 2))
    }
    cat(nv, "variables:", paste(vid, collapse = " "), "\n")
    
    ## m/z range
    if ("mz" %in% colnames(object@variable_info)) {
      mz_vals <- object@variable_info$mz
      cat("  m/z range: [",
          round(min(mz_vals), 4), " - ", round(max(mz_vals), 4), "]\n",
          sep = "")
    }
    
    ## 3. sample_info
    cat(crayon::green("3. sample_info: "))
    cat("[", nrow(object@sample_info), "x",
        ncol(object@sample_info), "data.frame]\n")
    np <- nrow(object@sample_info)
    
    pid <- if (np <= 5) {
      object@sample_info$sample_id
    } else {
      c(head(object@sample_info$sample_id, 3), "...",
        tail(object@sample_info$sample_id, 2))
    }
    cat(np, "pixels:", paste(pid, collapse = " "), "\n")
    
    ## spatial extent
    if (all(c("x", "y") %in% colnames(object@sample_info))) {
      cat("  Spatial extent: x[",
          min(object@sample_info$x), "-", max(object@sample_info$x), "], ",
          "y[", min(object@sample_info$y), "-", max(object@sample_info$y), "]",
          sep = "")
      if ("z" %in% colnames(object@sample_info)) {
        cat(", ",
            "z[", min(object@sample_info$z), "-", max(object@sample_info$z), "]",
            sep = "")
      }
      cat("\n")
    }
    
    ## 4. notes
    has_sinfo_note <- nrow(object@sample_info_note) > 0
    has_vinfo_note <- nrow(object@variable_info_note) > 0
    if (has_sinfo_note | has_vinfo_note) {
      cat(crayon::green("4. notes: "))
      note_parts <- c()
      if (has_sinfo_note) note_parts <- c(note_parts, "sample_info_note")
      if (has_vinfo_note) note_parts <- c(note_parts, "variable_info_note")
      cat(paste(note_parts, collapse = ", "), "\n")
    }

    process_label <- "5"
    roi_summary <- try(list_roi(object), silent = TRUE)
    if (!inherits(roi_summary, "try-error") && nrow(roi_summary) > 0) {
      cat(crayon::green("5. roi: "))
      cat(nrow(roi_summary), "ROI(s), ",
          sum(roi_summary$n_pixels), "pixel-memberships\n", sep = "")
      process_label <- "6"
    }
    
    ## processing history
    if (length(object@process_info) > 0) {
      cat(crayon::green(paste0(process_label, ". process_info: ")))
      cat(length(object@process_info), "steps\n")
      
      # ✅ IMPROVED: Show ALL steps (not just last 3)
      all_steps <- names(object@process_info)
      
      for (step_index in seq_along(all_steps)) {
        step_name <- all_steps[step_index]
        param <- object@process_info[[step_name]]
        
        # Get function name
        func_name <- if (methods::is(param, "tidymass_parameter")) {
          param@function_name
        } else {
          "unknown"
        }
        
        # ✅ IMPROVED: Get timestamp (YYYY-MM-DD HH:MM:SS format)
        timestamp <- if (methods::is(param, "tidymass_parameter")) {
          format(param@time, "%Y-%m-%d %H:%M:%S")
        } else {
          "N/A"
        }
        
        # ✅ IMPROVED: Extract key parameters
        params_summary <- ""
        if (methods::is(param, "tidymass_parameter") && 
            length(param@parameter) > 0) {
          
          key_params <- c()
          params <- param@parameter
          
          # Retention rate
          if (!is.null(params$retention_rate)) {
            key_params <- c(key_params, params$retention_rate)
          } else if (!is.null(params$num_retained) && !is.null(params$num_removed)) {
            total <- params$num_retained + params$num_removed
            if (total > 0) {
              rate <- sprintf("%.1f%%", 100 * params$num_retained / total)
              key_params <- c(key_params, rate)
            }
          }
          
          # Threshold
          if (!is.null(params$zero_threshold_percent)) {
            key_params <- c(key_params, sprintf("threshold=%s", params$zero_threshold_percent))
          } else if (!is.null(params$zero_threshold)) {
            key_params <- c(key_params, sprintf("threshold=%.1f%%", params$zero_threshold * 100))
          }
          
          # Min pixels
          if (!is.null(params$min_pixels)) {
            key_params <- c(key_params, sprintf("min_pixels=%d", params$min_pixels))
          }
          
          # Method
          if (!is.null(params$method)) {
            key_params <- c(key_params, sprintf("method=%s", params$method))
          }
          
          # Format parameter summary
          if (length(key_params) > 0) {
            params_summary <- paste0(" [", paste(key_params, collapse=", "), "]")
          }
        }
        
        # ✅ IMPROVED: Output with step number, function name, timestamp, and parameters
        cat(sprintf("  Step %d: %s\n", step_index, func_name))
        cat(sprintf("    └─ Time: %s%s\n", timestamp, params_summary))
      }
    }
    
    invisible(object)
  }
)


# ---------------------------------------------------------------------------
# [ (subset)
# ---------------------------------------------------------------------------

#' Subset a msi_dataset object
#'
#' @description
#' Subset by variables (rows) and/or pixels (columns).
#'
#' @param x A \code{msi_dataset} object.
#' @param i Row indices (variables): numeric, logical, or character.
#' @param j Column indices (pixels): numeric, logical, or character.
#' @param drop Logical; default \code{FALSE}.
#' @return A subsetted \code{msi_dataset} object.
#'
#' @docType methods
#' @rdname extract-methods
#' @aliases [,msi_dataset,ANY,ANY,ANY-method
#' @export

setMethod(
  "[",
  signature = "msi_dataset",
  function(x, i, j, drop = FALSE) {
    
    if (missing(i)) i <- seq_len(nrow(x@expression_data))
    if (missing(j)) j <- seq_len(ncol(x@expression_data))
    
    new_expr  <- x@expression_data[i, j, drop = FALSE]
    new_vinfo <- x@variable_info[i, , drop = FALSE]
    new_sinfo <- x@sample_info[j, , drop = FALSE]
    
    ## Rebuild notes if needed (keep existing columns only)
    new_snote <- x@sample_info_note
    new_vnote <- x@variable_info_note
    
    new_obj <- new(
      Class              = "msi_dataset",
      expression_data    = new_expr,
      ms2_data           = x@ms2_data,
      annotation_table   = x@annotation_table,
      sample_info        = new_sinfo,
      variable_info      = new_vinfo,
      sample_info_note   = new_snote,
      variable_info_note = new_vnote,
      process_info       = x@process_info,
      other_files        = x@other_files,
      version            = x@version,
      activated          = x@activated
    )
    new_obj
  }
)


# ---------------------------------------------------------------------------
# dim() and dimnames()
# ---------------------------------------------------------------------------

#' Dimensions of a msi_dataset
#' @param x A \code{msi_dataset} object.
#' @return Integer vector \code{c(n_variables, n_pixels)}.
#' @docType methods
#' @rdname dim-methods
#' @aliases dim,msi_dataset-method
#' @export
setMethod("dim", "msi_dataset", function(x) {
  c(nrow(x@expression_data), ncol(x@expression_data))
})

#' Dimension names of a msi_dataset
#' @param x A \code{msi_dataset} object.
#' @return A list: variable IDs (rows) and pixel IDs (columns).
#' @docType methods
#' @rdname dimnames-methods
#' @aliases dimnames,msi_dataset-method
#' @export
setMethod("dimnames", "msi_dataset", function(x) {
  list(rownames(x@expression_data), colnames(x@expression_data))
})


# ---------------------------------------------------------------------------
# update_msi_dataset()
# ---------------------------------------------------------------------------

#' Update a msi_dataset Object
#'
#' @description
#' Re-synchronises \code{expression_data} dimnames with the id columns in
#' \code{variable_info} / \code{sample_info} and regenerates note tables
#' when they are missing or have the wrong number of rows.
#'
#' @param object A \code{msi_dataset} object.
#' @return The updated \code{msi_dataset} object (invisibly).
#' @export

update_msi_dataset <- function(object) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }

  old_sample_info <- as.data.frame(object@sample_info, stringsAsFactors = FALSE,
                                   check.names = FALSE)
  old_sample_id <- if ("sample_id" %in% colnames(old_sample_info)) {
    as.character(old_sample_info$sample_id)
  } else {
    character()
  }

  sample_info <- .normalize_sample_info(old_sample_info)
  new_sample_id <- as.character(sample_info$sample_id)
  variable_info <- .normalize_variable_info(object@variable_info)
  expression_data <- object@expression_data

  if (length(old_sample_id) == length(new_sample_id) && length(new_sample_id) > 0) {
    rename_map <- setNames(new_sample_id, old_sample_id)
    old_colnames <- colnames(expression_data)
    if (!is.null(old_colnames)) {
      hit <- old_colnames %in% names(rename_map)
      old_colnames[hit] <- unname(rename_map[old_colnames[hit]])
      colnames(expression_data) <- old_colnames
    }

    roi_set <- object@other_files$roi
    if (!is.null(roi_set) && methods::is(roi_set, "msi_roi_set")) {
      for (i in seq_along(roi_set@rois)) {
        roi_ids <- roi_set@rois[[i]]@sample_id
        hit <- roi_ids %in% names(rename_map)
        roi_ids[hit] <- unname(rename_map[roi_ids[hit]])
        roi_set@rois[[i]]@sample_id <- unique(as.character(roi_ids))
      }
      object@other_files$roi <- roi_set
    }
  }

  intersect_sample_id <- intersect(sample_info$sample_id, colnames(expression_data))
  intersect_variable_id <- intersect(variable_info$variable_id, rownames(expression_data))

  # If IDs are fully mismatched but dimensions still align, repair dimnames by position.
  if (length(intersect_sample_id) == 0 &&
      ncol(expression_data) == nrow(sample_info) &&
      nrow(sample_info) > 0) {
    colnames(expression_data) <- as.character(sample_info$sample_id)
    intersect_sample_id <- as.character(sample_info$sample_id)
  }
  if (length(intersect_variable_id) == 0 &&
      nrow(expression_data) == nrow(variable_info) &&
      nrow(variable_info) > 0) {
    rownames(expression_data) <- as.character(variable_info$variable_id)
    intersect_variable_id <- as.character(variable_info$variable_id)
  }

  expression_data <- expression_data[intersect_variable_id, intersect_sample_id, drop = FALSE]
  sample_info <- sample_info[match(intersect_sample_id, sample_info$sample_id), , drop = FALSE]
  variable_info <- variable_info[match(intersect_variable_id, variable_info$variable_id), , drop = FALSE]

  object@sample_info <- sample_info
  object@variable_info <- variable_info
  object@expression_data <- expression_data
  
  object@sample_info_note <- .sync_info_note(object@sample_info_note,
                                             object@sample_info)
  object@variable_info_note <- .sync_info_note(object@variable_info_note,
                                               object@variable_info)
  
  object <- add_process_log(
    object = object,
    function_name = "update_msi_dataset()",
    parameter = list(
      sample_id_standardized = !identical(old_sample_id, new_sample_id)
    )
  )
  
  ## validate
  check_result <- check_msi_dataset_class(object)
  if (!isTRUE(check_result))
    warning("Validation errors after update:\n",
            paste(check_result, collapse = "\n"))
  
  invisible(object)
}
