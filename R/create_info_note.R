# ============================================================================
# Enhanced sample_info_note and variable_info_note Documentation - UPDATED
# Updated to match new structure:
# - sample_info: sample_id, x, y, z, subject_id, section_id, resolution, matrix (8 columns)
# - variable_info: variable_id, mz, ms2 (3 columns)
# ============================================================================

#' Create Enhanced sample_info_note
#' @description Generates comprehensive documentation for sample_info columns
#' @param msi_object An msi_dataset object
#' @return A data.frame with columns: name, meaning
#' @export

create_enhanced_sample_info_note <- function(msi_object) {
  
  # Define meanings for standard MSI spatial and metadata columns
  standard_meanings <- list(
    x = "The X-axis coordinate of the pixel in the imaging space. Usually corresponds to the horizontal position.",
    
    y = "The Y-axis coordinate of the pixel in the imaging space. Usually corresponds to the vertical position.",
    
    z = "The Z-axis coordinate of the pixel. Represents the depth of a tissue section or different slice layers in a 3D dataset.",
    
    sample_id = "Unique pixel identifier in the image (one id per pixel).",
    
    subject_id = "Unique identifier for the biological sample/subject.",
    
    section_id = "Unique section identifier. If a sample is cut into multiple slices for imaging, this ID distinguishes different slices.",
    
    resolution = "Imaging resolution. Usually expressed in micrometers per pixel (µm/pixel), indicating the actual tissue size represented by one pixel.",
    
    matrix = "The matrix used in MALDI (Matrix-Assisted Laser Desorption/Ionization). Records which matrix was applied to assist ionization for this sample."
    
  )
  
  # Build note table
  sample_info_cols <- colnames(msi_object@sample_info)
  
  meanings <- sapply(sample_info_cols, function(col) {
    if (col %in% names(standard_meanings)) {
      standard_meanings[[col]]
    } else {
      col  # Default: use column name as meaning if not predefined
    }
  }, USE.NAMES = FALSE)
  
  sample_info_note <- data.frame(
    name = sample_info_cols,
    meaning = meanings,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  return(sample_info_note)
}


#' Create Enhanced variable_info_note
#' @description Generates comprehensive documentation for variable_info columns
#' @param msi_object An msi_dataset object
#' @return A data.frame with columns: name, meaning
#' @export

create_enhanced_variable_info_note <- function(msi_object) {
  
  # Define meanings for standard m/z and feature-related columns
  standard_meanings <- list(
    variable_id = "Unique variable identifier. Used to uniquely index each m/z peak within the internal data structure, facilitating data matching and manipulation.",
    
    mz = "Mass-to-charge ratio. The core variable in mass spectrometry data, representing the ratio of an ion's mass to its charge number.",
    
    ms2 = "MS2 availability indicator (yes or no). Indicates whether tandem mass spectrometry (MS2) fragmentation data is available for this m/z feature."
  )
  
  # Build note table
  variable_info_cols <- colnames(msi_object@variable_info)
  
  meanings <- sapply(variable_info_cols, function(col) {
    if (col %in% names(standard_meanings)) {
      standard_meanings[[col]]
    } else {
      col  # Default: use column name as meaning if not predefined
    }
  }, USE.NAMES = FALSE)
  
  variable_info_note <- data.frame(
    name = variable_info_cols,
    meaning = meanings,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  return(variable_info_note)
}


#' Update both notes in msi_dataset object
#' @description Convenience function to update both sample_info_note and variable_info_note
#' @param msi_object An msi_dataset object
#' @return Updated msi_dataset object with enhanced notes
#' @export

update_msi_notes <- function(msi_object) {
  msi_object@sample_info_note <- create_enhanced_sample_info_note(msi_object)
  msi_object@variable_info_note <- create_enhanced_variable_info_note(msi_object)
  return(msi_object)
}


# ============================================================================
# USAGE EXAMPLE 1: Direct application after conversion
# ============================================================================

# After converting from Cardinal object:
# my_msi_full <- convert_cardinal2msi_dataset(cardinal_object, ...)

# Update with enhanced notes:
# my_msi_full <- update_msi_notes(my_msi_full)

# Then proceed with filtering:
# my_msi_clean <- my_msi_full %>%
#   filter_pixels(x > 2) %>%
#   filter_variables(freq > 0.2)


# ============================================================================
# FIELD DEFINITIONS - Updated Structure
# ============================================================================

# SAMPLE_INFO COLUMNS (8 required/common columns):
# ─────────────────────────────────────────────────
# 
# pixel_id (REQUIRED)
#   Purpose: Unique identifier for each pixel in the imaging space
#   Type: character
#   Usage: Column names of expression_data; primary index for pixels
#   Example: "pixel_1", "pixel_2", "pixel_54"
#
# x (REQUIRED)
#   Purpose: X-axis spatial coordinate in the imaging space
#   Type: numeric
#   Usage: Spatial mapping, visualization, filtering
#   Range: Depends on imaging resolution and sample size (e.g., 1-10)
#
# y (REQUIRED)
#   Purpose: Y-axis spatial coordinate in the imaging space
#   Type: numeric
#   Usage: Spatial mapping, visualization, filtering
#   Range: Depends on imaging resolution and sample size (e.g., 1-6)
#
# z (OPTIONAL)
#   Purpose: Z-axis spatial coordinate; depth of tissue section or slice layer
#   Type: numeric
#   Usage: 3D reconstruction, multi-slice datasets
#   Example: 1 for single slice, 1-5 for multi-slice
#
# sample_id (REQUIRED)
#   Purpose: Unique identifier for the biological sample
#   Type: character
#   Usage: Distinguish different individuals or biological replicates
#   Example: "renal_cancer_1", "control_2"
#
# section_id (OPTIONAL)
#   Purpose: Unique section/slice identifier within a sample
#   Type: character
#   Usage: Distinguish multiple sections cut from the same tissue
#   Example: "slice_01", "slice_02"
#
# resolution (OPTIONAL)
#   Purpose: Spatial resolution of the imaging
#   Type: numeric
#   Unit: micrometers per pixel (µm/pixel)
#   Usage: Recording imaging resolution for downstream analysis
#   Example: 25 (meaning 25 µm/pixel)
#
# matrix (OPTIONAL)
#   Purpose: MALDI matrix type used for ionization
#   Type: character
#   Common values: "DHB", "CHCA", "NEDC", "noMatrix"
#   Usage: Method documentation and cross-sample comparisons


# VARIABLE_INFO COLUMNS (3 columns):
# ─────────────────────────────────
#
# variable_id (REQUIRED)
#   Purpose: Unique identifier for each m/z feature
#   Type: character
#   Usage: Row names of expression_data; primary index for features
#   Example: "m50.0038_1", "m200.1234_15", "m500.5678_3"
#   Generated by: generate_variable_id() function
#
# mz (REQUIRED)
#   Purpose: Mass-to-charge ratio
#   Type: numeric
#   Usage: Chemical composition identification, feature filtering
#   Range: Depends on mass spectrometer (e.g., 50-1000)
#   Precision: Usually 4-5 decimal places
#   Example: 200.1234, 500.5678
#
# ms2 (OPTIONAL)
#   Purpose: Indicates availability of MS2 fragmentation data
#   Type: character (yes/no) or logical (TRUE/FALSE)
#   Usage: Identifying features with MS2 fragmentation spectra
#   Example: "yes" or "no"; or TRUE or FALSE
#
# (Optional additional columns that may be added later)
# count (OPTIONAL)
#   Purpose: Number of non-zero pixels for this m/z feature
#   Type: numeric (integer)
#   Usage: Feature quality filtering
#
# freq (OPTIONAL)
#   Purpose: Proportion of pixels where this m/z is present
#   Type: numeric (0-1)
#   Calculation: count / total_pixels
#   Usage: Feature filtering by prevalence
#   Example: 0.2 means feature present in 20% of pixels


# ============================================================================
# Summary of Changes from Previous Version
# ============================================================================

# 1. sample_info structure updated:
#    - Added pixel_id as REQUIRED (used for column names of expression_data)
#    - sample_id is now a separate column (not used for column names)
#    - Total columns: 8 (pixel_id, x, y, z, sample_id, section_id, resolution, matrix)
#
# 2. variable_info structure updated:
#    - Simplified to 3 core columns: variable_id, mz, ms2
#    - ms2 is now a yes/no indicator (character) not logical
#    - Removed optional count/freq from standard structure
#
# 3. expression_data naming:
#    - Column names now use pixel_id (not sample_id)
#    - Row names use variable_id (unchanged)
#
# 4. Updated notes to reflect new structure:
#    - All 8 sample_info columns documented
#    - All 3 variable_info columns documented
