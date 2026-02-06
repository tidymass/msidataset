# ============================================================================
# Accessor Methods for msi_dataset Class
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

#' Get number of pixels from msi_dataset
#'
#' @param object A msi_dataset object
#' @return Integer number of pixels
#' @export
setMethod(
  "n_pixels",
  "msi_dataset",
  function(object) {
    nrow(object@pixel_info)
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
    object@variable_info$mz
  }
)

# ---------------------------------------------------------------------------
# Slot accessors
# ---------------------------------------------------------------------------

#' Get expression data
#'
#' @param object A msi_dataset object
#' @return data.frame
#' @export
setMethod(
  "expression_data",
  "msi_dataset",
  function(object) {
    object@expression_data
  }
)

#' Get pixel information
#'
#' @param object A msi_dataset object
#' @return data.frame
#' @export
setMethod(
  "pixel_info",
  "msi_dataset",
  function(object) {
    object@pixel_info
  }
)

#' Get sample information
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

#' Get pixel info note
#'
#' @param object A msi_dataset object
#' @return data.frame
#' @export
setMethod(
  "pixel_info_note",
  "msi_dataset",
  function(object) {
    object@pixel_info_note
  }
)

#' Get sample info note
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
