# ============================================================================
# Enhanced sample_info_note and variable_info_note Documentation
# ============================================================================

#' Create Enhanced sample_info_note
#' @description Generates comprehensive documentation for sample_info columns
#' @param msi_object An msi_dataset object
#' @return A data.frame with columns: name, meaning
#' @export

create_enhanced_sample_info_note <- function(msi_object) {
  
  # Define meanings for standard MSI spatial and metadata columns
  standard_meanings <- list(
    sample_id = "Unique sample identifier used to distinguish different biological samples (e.g., tissue from different individuals).",
    
    x = "The X-axis coordinate of the sample/pixel in the imaging space. Usually corresponds to the horizontal position.",
    
    y = "The Y-axis coordinate of the sample/pixel in the imaging space. Usually corresponds to the vertical position.",
    
    z = "The Z-axis coordinate of the sample/pixel. May represent the depth of a tissue section or different slice layers in a 3D dataset.",
    
    run = "Experimental run or batch number. Used to distinguish data acquired at different times or under different instrument conditions; important for batch effect correction.",
    
    X3DPositionX = "The X coordinate in 3D space (if the data includes 3D reconstruction information).",
    
    X3DPositionY = "The Y coordinate in 3D space.",
    
    X3DPositionZ = "The Z coordinate in 3D space.",
    
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
    
    count = "Count. Refers to the number of non-zero pixels for this m/z value in the dataset, or the total frequency of this feature's appearance (depending on the data processing pipeline).",
    
    freq = "Frequency. Usually refers to the proportion of pixels where this m/z value is present (i.e., count divided by the total number of pixels).",
    
    ms2 = "Logical indicator. Indicates whether MS2 (tandem mass spectrometry) data is available for this m/z feature."
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
# USAGE EXAMPLE 2: Manual enhancement with custom meanings
# ============================================================================

# create_enhanced_sample_info_note <- function(msi_object) {
#   sample_info_note <- create_enhanced_sample_info_note(msi_object)
#   
#   # Add custom meanings for additional columns
#   # sample_info_note$meaning[sample_info_note$name == "custom_col"] <- 
#   #   "My custom explanation..."
#   
#   return(sample_info_note)
# }


# ============================================================================
# Detailed Field Definitions Reference
# ============================================================================

# SAMPLE_INFO COLUMNS (Spatial & Metadata):
# ──────────────────────────────────────────
# 
# sample_id
#   Purpose: Unique identifier for each pixel/sample
#   Type: character
#   Usage: Internal indexing, traceability
#   Example: "renal_cancer_pixel_1", "renal_cancer_pixel_2"
#
# x, y
#   Purpose: 2D spatial coordinates in the imaging space
#   Type: numeric
#   Usage: Spatial mapping, visualization
#   Range: Depends on imaging resolution and sample size
#
# z (optional)
#   Purpose: 3D depth coordinate or slice indicator
#   Type: numeric
#   Usage: 3D reconstruction, multi-slice datasets
#   Common in: Tissue sectioning, volumetric imaging
#
# run (optional)
#   Purpose: Batch effect tracking
#   Type: character/numeric
#   Usage: Quality control, batch correction
#   Example: "run1", "run2", or 1, 2, 3
#
# section_id (optional)
#   Purpose: Slice identifier for multi-section samples
#   Type: character
#   Usage: Distinguishing multiple sections of same tissue
#   Example: "slice_01", "slice_02"
#
# resolution (optional)
#   Purpose: Recording imaging resolution
#   Type: numeric
#   Unit: micrometers per pixel (µm/pixel)
#   Usage: Spatial scale information
#   Example: 25 (25 µm/pixel)
#
# matrix (optional)
#   Purpose: MALDI matrix recording
#   Type: character
#   Common values: "DHB", "CHCA", "MALDI", "noMatrix"
#   Usage: Method documentation, comparison studies


# VARIABLE_INFO COLUMNS (m/z & Feature Data):
# ──────────────────────────────────────────
#
# variable_id
#   Purpose: Unique identifier for each m/z feature
#   Type: character
#   Usage: Internal indexing, feature traceability
#   Example: "m50_1", "m200_15", "m500_3"
#   Generated by: generate_variable_id() function
#
# mz
#   Purpose: Mass-to-charge ratio
#   Type: numeric
#   Usage: Chemical composition identification, filtering
#   Range: Depends on mass spectrometer (e.g., 50-1000 m/z)
#   Precision: Usually 4-5 decimal places
#
# count (optional)
#   Purpose: Feature prevalence measure
#   Type: numeric (integer)
#   Interpretation: Number of pixels where this m/z is detected (non-zero)
#   Usage: Quality filtering, spurious signal removal
#   Example: If count=50 and total pixels=100, freq would be 0.5
#
# freq (optional)
#   Purpose: Feature prevalence as proportion
#   Type: numeric (0-1)
#   Calculation: count / total_pixels
#   Usage: Filtering common/rare features
#   Example: freq=0.2 means feature present in 20% of pixels
#
# ms2 (optional)
#   Purpose: MS2 spectra availability
#   Type: logical
#   Usage: Identifying features with fragmentation data
#   Example: TRUE if MS2 available, FALSE otherwise
