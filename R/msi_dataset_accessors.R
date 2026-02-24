# ============================================================================
# Accessor Methods for msi_dataset (getters + setters)
# Author: Xiaoyu Xu <XIAOYU021@e.ntu.edu.sg>
# ============================================================================

# ---------------------------------------------------------------------------
# Slot GETTERS
# ---------------------------------------------------------------------------

#' Get expression data from a msi_dataset object
#' @param object A \code{msi_dataset} object.
#' @return A \code{data.frame} or \code{matrix}.
#' @export
setMethod("expression_data", "msi_dataset", function(object) {
  object@expression_data
})

#' Get sample information (pixel metadata)
#' @param object A \code{msi_dataset} object.
#' @return A \code{data.frame}.
#' @export
setMethod("sample_info", "msi_dataset", function(object) {
  object@sample_info
})

#' Get variable information (m/z feature metadata)
#' @param object A \code{msi_dataset} object.
#' @return A \code{data.frame}.
#' @export
setMethod("variable_info", "msi_dataset", function(object) {
  object@variable_info
})

#' Get sample_info_note
#' @param object A \code{msi_dataset} object.
#' @return A \code{data.frame}.
#' @export
setMethod("sample_info_note", "msi_dataset", function(object) {
  object@sample_info_note
})

#' Get variable_info_note
#' @param object A \code{msi_dataset} object.
#' @return A \code{data.frame}.
#' @export
setMethod("variable_info_note", "msi_dataset", function(object) {
  object@variable_info_note
})

#' Get MS2 data
#' @param object A \code{msi_dataset} object.
#' @return A \code{list}.
#' @export
setMethod("ms2_data", "msi_dataset", function(object) {
  object@ms2_data
})

#' Get process information
#' @param object A \code{msi_dataset} object.
#' @return A \code{list}.
#' @export
setMethod("process_info", "msi_dataset", function(object) {
  object@process_info
})

#' Get annotation table
#' @param object A \code{msi_dataset} object.
#' @return A \code{data.frame}.
#' @export
setMethod("annotation_table", "msi_dataset", function(object) {
  object@annotation_table
})

# ---------------------------------------------------------------------------
# Slot SETTERS (replacement methods)
# ---------------------------------------------------------------------------

#' Set expression data
#'
#' @description
#' Replaces the \code{expression_data} slot, ensuring row/column names stay
#' consistent with \code{variable_info} and \code{sample_info}.
#'
#' @param object A \code{msi_dataset} object.
#' @param value A \code{data.frame} or \code{matrix} with matching dimensions.
#' @return Updated \code{msi_dataset} object.
#' @export
setReplaceMethod("expression_data", "msi_dataset", function(object, value) {
  value <- as.data.frame(value, check.names = FALSE)
  rownames(value) <- object@variable_info$variable_id
  colnames(value) <- object@sample_info$sample_id
  object@expression_data <- value
  methods::validObject(object)
  object
})

#' Set sample_info
#' @param object A \code{msi_dataset} object.
#' @param value A \code{data.frame}; must contain \code{sample_id}, \code{x},
#'   \code{y}.
#' @return Updated \code{msi_dataset} object.
#' @export
setReplaceMethod("sample_info", "msi_dataset", function(object, value) {
  object@sample_info <- as.data.frame(value, stringsAsFactors = FALSE,
                                       check.names = FALSE)
  # rebuild note if dimension changed
  if (nrow(object@sample_info_note) != ncol(object@sample_info)) {
    object@sample_info_note <- data.frame(
      name    = colnames(object@sample_info),
      meaning = colnames(object@sample_info),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }
  colnames(object@expression_data) <- object@sample_info$sample_id
  methods::validObject(object)
  object
})

#' Set variable_info
#' @param object A \code{msi_dataset} object.
#' @param value A \code{data.frame}; must contain \code{variable_id}.
#' @return Updated \code{msi_dataset} object.
#' @export
setReplaceMethod("variable_info", "msi_dataset", function(object, value) {
  object@variable_info <- as.data.frame(value, stringsAsFactors = FALSE,
                                         check.names = FALSE)
  if (nrow(object@variable_info_note) != ncol(object@variable_info)) {
    object@variable_info_note <- data.frame(
      name    = colnames(object@variable_info),
      meaning = colnames(object@variable_info),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }
  rownames(object@expression_data) <- object@variable_info$variable_id
  methods::validObject(object)
  object
})

#' Set ms2_data
#' @param object A \code{msi_dataset} object.
#' @param value A \code{list}.
#' @return Updated \code{msi_dataset} object.
#' @export
setReplaceMethod("ms2_data", "msi_dataset", function(object, value) {
  object@ms2_data <- value
  object
})

#' Set annotation_table
#' @param object A \code{msi_dataset} object.
#' @param value A \code{data.frame}.
#' @return Updated \code{msi_dataset} object.
#' @export
setReplaceMethod("annotation_table", "msi_dataset", function(object, value) {
  object@annotation_table <- as.data.frame(value, stringsAsFactors = FALSE,
                                            check.names = FALSE)
  object
})

#' Set sample_info_note
#' @param object A \code{msi_dataset} object.
#' @param value A \code{data.frame} with \code{name} and \code{meaning}.
#' @return Updated \code{msi_dataset} object.
#' @export
setReplaceMethod("sample_info_note", "msi_dataset", function(object, value) {
  object@sample_info_note <- as.data.frame(value, stringsAsFactors = FALSE,
                                            check.names = FALSE)
  methods::validObject(object)
  object
})

#' Set variable_info_note
#' @param object A \code{msi_dataset} object.
#' @param value A \code{data.frame} with \code{name} and \code{meaning}.
#' @return Updated \code{msi_dataset} object.
#' @export
setReplaceMethod("variable_info_note", "msi_dataset", function(object, value) {
  object@variable_info_note <- as.data.frame(value, stringsAsFactors = FALSE,
                                              check.names = FALSE)
  methods::validObject(object)
  object
})

# ---------------------------------------------------------------------------
# Convenience COUNT accessors
# ---------------------------------------------------------------------------

#' Number of variables in a msi_dataset
#' @param object A \code{msi_dataset} object.
#' @return Integer.
#' @export
setMethod("n_variables", "msi_dataset", function(object) {
  nrow(object@variable_info)
})

#' Number of samples/pixels in a msi_dataset
#' @param object A \code{msi_dataset} object.
#' @return Integer.
#' @export
setMethod("n_samples", "msi_dataset", function(object) {
  nrow(object@sample_info)
})

#' Number of pixels in a msi_dataset (alias for \code{n_samples})
#' @param object A \code{msi_dataset} object.
#' @return Integer.
#' @export
setMethod("n_pixels", "msi_dataset", function(object) {
  nrow(object@sample_info)
})

#' Extract m/z values from a msi_dataset
#' @param object A \code{msi_dataset} object.
#' @return Numeric vector of m/z values, or \code{NULL} with a warning.
#' @export
setMethod("mz", "msi_dataset", function(object) {
  if ("mz" %in% colnames(object@variable_info)) {
    return(object@variable_info$mz)
  }
  warning("No 'mz' column found in variable_info.")
  return(NULL)
})

# ---------------------------------------------------------------------------
# Backward-compatibility: pixel_info == sample_info
# ---------------------------------------------------------------------------

#' Get pixel_info (backward compatibility — returns sample_info)
#' @param object A \code{msi_dataset} object.
#' @return A \code{data.frame} (same as \code{sample_info(object)}).
#' @export
setMethod("pixel_info", "msi_dataset", function(object) {
  message("Note: 'pixel_info' is now stored as 'sample_info'. ",
          "Use sample_info() instead.")
  object@sample_info
})

#' Get pixel_info_note (backward compatibility — returns sample_info_note)
#' @param object A \code{msi_dataset} object.
#' @return A \code{data.frame}.
#' @export
setMethod("pixel_info_note", "msi_dataset", function(object) {
  message("Note: 'pixel_info_note' is now stored as 'sample_info_note'. ",
          "Use sample_info_note() instead.")
  object@sample_info_note
})

# ---------------------------------------------------------------------------
# Spatial convenience accessors
# ---------------------------------------------------------------------------

#' Get spatial coordinates (x, y, z) from a msi_dataset
#'
#' @param object A \code{msi_dataset} object.
#' @return A \code{data.frame} with columns \code{sample_id}, \code{x},
#'   \code{y}, \code{z}.
#' @export
setMethod("get_spatial_coords", "msi_dataset", function(object) {
  required_cols <- c("sample_id", "x", "y")
  if (!all(required_cols %in% colnames(object@sample_info)))
    stop("sample_info must contain: sample_id, x, y.")

  out <- object@sample_info[, required_cols, drop = FALSE]
  out$z <- if ("z" %in% colnames(object@sample_info)) object@sample_info$z else 1L
  out
})

#' Get unique section identifiers
#' @param object A \code{msi_dataset} object.
#' @return Character vector of section IDs, or \code{NULL}.
#' @export
setMethod("get_section_ids", "msi_dataset", function(object) {
  if ("section_id" %in% colnames(object@sample_info)) {
    return(unique(object@sample_info$section_id))
  }
  warning("No 'section_id' column found in sample_info.")
  return(NULL)
})

#' Get pixels for a specific tissue section
#'
#' @param object A \code{msi_dataset} object.
#' @param section_id Character scalar — the desired section identifier.
#' @return A \code{data.frame} subset of \code{sample_info}.
#' @export
setMethod(
  "get_section_pixels",
  signature = c(object = "msi_dataset", section_id = "character"),
  function(object, section_id) {
    if (!"section_id" %in% colnames(object@sample_info))
      stop("No 'section_id' column found in sample_info.")
    object@sample_info[object@sample_info$section_id == section_id, ,
                       drop = FALSE]
  }
)

#' Get sample (pixel) IDs
#' @param object A \code{msi_dataset} object.
#' @return Character vector.
#' @export
setMethod("get_sample_id", "msi_dataset", function(object) {
  object@sample_info$sample_id
})

#' Get variable IDs
#' @param object A \code{msi_dataset} object.
#' @return Character vector.
#' @export
setMethod("get_variable_id", "msi_dataset", function(object) {
  object@variable_info$variable_id
})
