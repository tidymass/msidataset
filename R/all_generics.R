# ============================================================================
# S4 generic definitions for msi_dataset
# ============================================================================

#' Get sample information
#' @export
setGeneric(
  name = "sample_info",
  def = function(object, ...)
    standardGeneric("sample_info")
)

#' Get pixel information
#' @export
setGeneric(
  name = "pixel_info",
  def = function(object, ...)
    standardGeneric("pixel_info")
)

#' Get expression data
#' @export
setGeneric(
  name = "expression_data",
  def = function(object, ...)
    standardGeneric("expression_data")
)

#' Get variable information
#' @export
setGeneric(
  name = "variable_info",
  def = function(object, ...)
    standardGeneric("variable_info")
)

#' Get process information
#' @export
setGeneric(
  name = "process_info",
  def = function(object, ...)
    standardGeneric("process_info")
)

#' Get MS2 data
#' @export
setGeneric(
  name = "ms2_data",
  def = function(object, ...)
    standardGeneric("ms2_data")
)

#' Get sample info note
#' @export
setGeneric(
  name = "sample_info_note",
  def = function(object, ...)
    standardGeneric("sample_info_note")
)

#' Get variable info note
#' @export
setGeneric(
  name = "variable_info_note",
  def = function(object, ...)
    standardGeneric("variable_info_note")
)

# ---------------------------------------------------------------------------
# Additional generics used by accessor methods
# ---------------------------------------------------------------------------

#' Get number of variables
#' @export
setGeneric(
  name = "n_variables",
  def = function(object, ...)
    standardGeneric("n_variables")
)

#' Get number of pixels
#' @export
setGeneric(
  name = "n_pixels",
  def = function(object, ...)
    standardGeneric("n_pixels")
)

#' Get m/z values
#' @export
setGeneric(
  name = "mz",
  def = function(object, ...)
    standardGeneric("mz")
)

#' Get pixel info note
#' @export
setGeneric(
  name = "pixel_info_note",
  def = function(object, ...)
    standardGeneric("pixel_info_note")
)

#' Get annotation table
#' @export
setGeneric(
  name = "annotation_table",
  def = function(object, ...)
    standardGeneric("annotation_table")
)
