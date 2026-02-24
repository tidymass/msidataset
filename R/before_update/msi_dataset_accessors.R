# ============================================================================
# Accessor Methods for msi_dataset Class (Restructured for massdataset compatibility)
# Author: Xiaoyu Xu
# ============================================================================

# ---------------------------------------------------------------------------
# Basic statistics accessors
# ---------------------------------------------------------------------------

#' Get number of variables from msi_dataset
#'
#' @param object A msi_dataset object
#' @return Integer number of variables
#' @export
setMethod(
  "n_variables",
  "msi_dataset",
  function(object) {
    nrow(object@variable_info)
  }
)

#' Get number of samples/pixels from msi_dataset
#'
#' @param object A msi_dataset object
#' @return Integer number of samples (pixels)
#' @export
setMethod(
  "n_samples",
  "msi_dataset",
  function(object) {
    nrow(object@sample_info)
  }
)

#' Get number of pixels from msi_dataset (alias for n_samples)
#'
#' @param object A msi_dataset object
#' @return Integer number of pixels
#' @export
setMethod(
  "n_pixels",
  "msi_dataset",
  function(object) {
    nrow(object@sample_info)
  }
)

#' Get m/z values from msi_dataset
#'
#' @param object A msi_dataset object
#' @return Numeric vector of m/z values
#' @export
setMethod(
  "mz",
  "msi_dataset",
  function(object) {
    if ("mz" %in% colnames(object@variable_info)) {
      return(object@variable_info$mz)
    } else {
      warning("No 'mz' column found in variable_info")
      return(NULL)
    }
  }
)

# ---------------------------------------------------------------------------
# Slot accessors (matching massdataset structure)
# ---------------------------------------------------------------------------

#' Get expression data
#'
#' @param object A msi_dataset object
#' @return matrix, data.frame, or DelayedArray
#' @export
setMethod(
  "expression_data",
  "msi_dataset",
  function(object) {
    object@expression_data
  }
)

#' Get sample information (now contains pixel spatial data)
#'
#' @param object A msi_dataset object
#' @return data.frame
#' @export
setMethod(
  "sample_info",
  "msi_dataset",
  function(object) {
    object@sample_info
  }
)

#' Get variable information
#'
#' @param object A msi_dataset object
#' @return data.frame
#' @export
setMethod(
  "variable_info",
  "msi_dataset",
  function(object) {
    object@variable_info
  }
)

#' Get sample info note (now contains pixel info notes)
#'
#' @param object A msi_dataset object
#' @return data.frame
#' @export
setMethod(
  "sample_info_note",
  "msi_dataset",
  function(object) {
    object@sample_info_note
  }
)

#' Get variable info note
#'
#' @param object A msi_dataset object
#' @return data.frame
#' @export
setMethod(
  "variable_info_note",
  "msi_dataset",
  function(object) {
    object@variable_info_note
  }
)

#' Get MS2 data
#'
#' @param object A msi_dataset object
#' @return list
#' @export
setMethod(
  "ms2_data",
  "msi_dataset",
  function(object) {
    object@ms2_data
  }
)

#' Get process information
#'
#' @param object A msi_dataset object
#' @return list
#' @export
setMethod(
  "process_info",
  "msi_dataset",
  function(object) {
    object@process_info
  }
)

#' Get annotation table
#'
#' @param object A msi_dataset object
#' @return data.frame
#' @export
setMethod(
  "annotation_table",
  "msi_dataset",
  function(object) {
    object@annotation_table
  }
)

# ---------------------------------------------------------------------------
# Backward compatibility accessors (for old pixel_info)
# ---------------------------------------------------------------------------

#' Get pixel information (backward compatibility)
#'
#' @description
#' Returns sample_info for backward compatibility.
#' In the new structure, pixel information is stored in sample_info.
#'
#' @param object A msi_dataset object
#' @return data.frame
#' @export
setMethod(
  "pixel_info",
  "msi_dataset",
  function(object) {
    message("Note: pixel_info is now part of sample_info in the restructured msi_dataset.")
    object@sample_info
  }
)

#' Get pixel info note (backward compatibility)
#'
#' @description
#' Returns sample_info_note for backward compatibility.
#' In the new structure, pixel information notes are in sample_info_note.
#'
#' @param object A msi_dataset object
#' @return data.frame
#' @export
setMethod(
  "pixel_info_note",
  "msi_dataset",
  function(object) {
    message("Note: pixel_info_note is now part of sample_info_note in the restructured msi_dataset.")
    object@sample_info_note
  }
)

# ---------------------------------------------------------------------------
# Spatial data accessors (convenience functions)
# ---------------------------------------------------------------------------

#' Get spatial coordinates from msi_dataset
#'
#' @param object A msi_dataset object
#' @return data.frame with sample_id, x, y, z columns
#' @export
setMethod(
  "get_spatial_coords",
  "msi_dataset",
  function(object) {
    required_cols <- c("sample_id", "x", "y")
    if (!all(required_cols %in% colnames(object@sample_info))) {
      stop("sample_info must contain: sample_id, x, y columns")
    }
    
    out <- object@sample_info[, required_cols, drop = FALSE]
    if ("z" %in% colnames(object@sample_info)) {
      out$z <- object@sample_info$z
    } else {
      out$z <- 1
    }
    out
  }
)

#' Get unique section IDs
#'
#' @param object A msi_dataset object
#' @return Character vector of section IDs
#' @export
setMethod(
  "get_section_ids",
  "msi_dataset",
  function(object) {
    if ("section_id" %in% colnames(object@sample_info)) {
      return(unique(object@sample_info$section_id))
    } else {
      warning("No 'section_id' column found in sample_info")
      return(NULL)
    }
  }
)

#' Get pixels for a specific section
#'
#' @param object A msi_dataset object
#' @param section_id Section identifier
#' @return data.frame of sample_info for specified section
#' @export
setMethod(
  "get_section_pixels",
  signature = c(object = "msi_dataset", section_id = "character"),
  function(object, section_id) {
    if (!"section_id" %in% colnames(object@sample_info)) {
      stop("No 'section_id' column found in sample_info")
    }
    object@sample_info[object@sample_info$section_id == section_id, , drop = FALSE]
  }
)

#' Get sample IDs (pixel IDs)
#'
#' @param object A msi_dataset object
#' @return Character vector of sample IDs
#' @export
setMethod(
  "get_sample_id",
  "msi_dataset",
  function(object) {
    object@sample_info$sample_id
  }
)

#' Get variable IDs
#'
#' @param object A msi_dataset object
#' @return Character vector of variable IDs
#' @export
setMethod(
  "get_variable_id",
  "msi_dataset",
  function(object) {
    object@variable_info$variable_id
  }
)