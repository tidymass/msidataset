# ============================================================================
# Generic Functions for msi_dataset (Restructured)
# Author: Xiaoyu Xu
# ============================================================================

# ---------------------------------------------------------------------------
# Basic generics (compatible with massdataset)
# ---------------------------------------------------------------------------

#' @title n_variables
#' @description Get number of variables
#' @param object An msi_dataset object
#' @return Integer
#' @export
setGeneric(
  name = "n_variables",
  def = function(object) {
    standardGeneric("n_variables")
  }
)

#' @title n_samples
#' @description Get number of samples (pixels in MSI context)
#' @param object An msi_dataset object
#' @return Integer
#' @export
setGeneric(
  name = "n_samples",
  def = function(object) {
    standardGeneric("n_samples")
  }
)

#' @title n_pixels
#' @description Get number of pixels (alias for n_samples)
#' @param object An msi_dataset object
#' @return Integer
#' @export
setGeneric(
  name = "n_pixels",
  def = function(object) {
    standardGeneric("n_pixels")
  }
)

#' @title mz
#' @description Get m/z values
#' @param object An msi_dataset object
#' @return Numeric vector
#' @export
setGeneric(
  name = "mz",
  def = function(object) {
    standardGeneric("mz")
  }
)

# ---------------------------------------------------------------------------
# Slot accessor generics (massdataset-compatible)
# ---------------------------------------------------------------------------

#' @title expression_data
#' @description Get expression data
#' @param object An msi_dataset object
#' @return Matrix, data.frame, or DelayedArray
#' @export
setGeneric(
  name = "expression_data",
  def = function(object) {
    standardGeneric("expression_data")
  }
)

#' @title sample_info
#' @description Get sample information (contains pixel spatial data)
#' @param object An msi_dataset object
#' @return data.frame
#' @export
setGeneric(
  name = "sample_info",
  def = function(object) {
    standardGeneric("sample_info")
  }
)

#' @title variable_info
#' @description Get variable information
#' @param object An msi_dataset object
#' @return data.frame
#' @export
setGeneric(
  name = "variable_info",
  def = function(object) {
    standardGeneric("variable_info")
  }
)

#' @title sample_info_note
#' @description Get sample information notes (contains pixel info notes)
#' @param object An msi_dataset object
#' @return data.frame
#' @export
setGeneric(
  name = "sample_info_note",
  def = function(object) {
    standardGeneric("sample_info_note")
  }
)

#' @title variable_info_note
#' @description Get variable information notes
#' @param object An msi_dataset object
#' @return data.frame
#' @export
setGeneric(
  name = "variable_info_note",
  def = function(object) {
    standardGeneric("variable_info_note")
  }
)

#' @title ms2_data
#' @description Get MS2 data
#' @param object An msi_dataset object
#' @return list
#' @export
setGeneric(
  name = "ms2_data",
  def = function(object) {
    standardGeneric("ms2_data")
  }
)

#' @title process_info
#' @description Get processing information
#' @param object An msi_dataset object
#' @return list
#' @export
setGeneric(
  name = "process_info",
  def = function(object) {
    standardGeneric("process_info")
  }
)

#' @title annotation_table
#' @description Get annotation table
#' @param object An msi_dataset object
#' @return data.frame
#' @export
setGeneric(
  name = "annotation_table",
  def = function(object) {
    standardGeneric("annotation_table")
  }
)

# ---------------------------------------------------------------------------
# Backward compatibility generics
# ---------------------------------------------------------------------------

#' @title pixel_info
#' @description Get pixel information (backward compatibility - returns sample_info)
#' @param object An msi_dataset object
#' @return data.frame
#' @export
setGeneric(
  name = "pixel_info",
  def = function(object) {
    standardGeneric("pixel_info")
  }
)

#' @title pixel_info_note
#' @description Get pixel information notes (backward compatibility - returns sample_info_note)
#' @param object An msi_dataset object
#' @return data.frame
#' @export
setGeneric(
  name = "pixel_info_note",
  def = function(object) {
    standardGeneric("pixel_info_note")
  }
)

# ---------------------------------------------------------------------------
# Spatial data generics
# ---------------------------------------------------------------------------

#' @title get_spatial_coords
#' @description Get spatial coordinates (x, y, z)
#' @param object An msi_dataset object
#' @return data.frame
#' @export
setGeneric(
  name = "get_spatial_coords",
  def = function(object) {
    standardGeneric("get_spatial_coords")
  }
)

#' @title get_section_ids
#' @description Get unique section IDs
#' @param object An msi_dataset object
#' @return Character vector
#' @export
setGeneric(
  name = "get_section_ids",
  def = function(object) {
    standardGeneric("get_section_ids")
  }
)

#' @title get_section_pixels
#' @description Get pixels for a specific section
#' @param object An msi_dataset object
#' @param section_id Section identifier
#' @return data.frame
#' @export
setGeneric(
  name = "get_section_pixels",
  def = function(object, section_id) {
    standardGeneric("get_section_pixels")
  }
)

#' @title get_sample_id
#' @description Get sample IDs (pixel IDs)
#' @param object An msi_dataset object
#' @return Character vector
#' @export
setGeneric(
  name = "get_sample_id",
  def = function(object) {
    standardGeneric("get_sample_id")
  }
)

#' @title get_variable_id
#' @description Get variable IDs
#' @param object An msi_dataset object
#' @return Character vector
#' @export
setGeneric(
  name = "get_variable_id",
  def = function(object) {
    standardGeneric("get_variable_id")
  }
)