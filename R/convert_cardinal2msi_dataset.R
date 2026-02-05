#' Convert Cardinal Object to msi_dataset
#' 
#' @author Xiaoyu Xu
#' @description
#' Converts a Cardinal MSImagingExperiment or MSImagingArrays object to an msi_dataset object.
#' For MSImagingArrays (Processed mode), automatic peak alignment is performed at 10ppm 
#' to establish a common reference m/z axis.
#' 
#' @param cardinal_object A Cardinal object (MSImagingExperiment or MSImagingArrays)
#' @param tolerance Alignment tolerance, default is 10
#' @param units Tolerance unit, default is "ppm"
#' @param sample_name Sample name. If NULL, extracted from Cardinal object or set to "sample1"
#' @param section_name Section/slice name, default is "section1"
#' @param resolution Spatial resolution of the acquisition (e.g., "50um"). Default is "unknown"
#' @param matrix_name Chemical matrix used for ionization. Default is "unknown"
#' 
#' @return An msi_dataset object
#' 
#' @details
#' This function performs the following steps:
#' 1. Automatic peak alignment for MSImagingArrays objects (10ppm tolerance by default)
#' 2. Extraction of expression matrix (variables × pixels)
#' 3. Extraction of common m/z reference axis
#' 4. Extraction of spatial coordinates (supports x, y, z for 3D data)
#' 5. Construction of pixel_info and sample_info data frames
#' 6. Creation of the final msi_dataset object
#' 
#' The function ensures compatibility with the msi_dataset class structure:
#' - expression_data: matrix with variables as rows and pixels as columns
#' - pixel_info: contains pixel_id, x, y, z coordinates
#' - sample_info: contains sample_id, section_id, z, resolution, matrix
#' - variable_info: contains variable_id, mz, and ms2 flag
#' 
#' @export
#' @examples
#' \dontrun{
#' # Load Cardinal data
#' library(Cardinal)
#' data(MSImagingExperiment)
#' 
#' # Convert to msi_dataset
#' msi_obj <- convert_cardinal2msi_dataset(
#'   cardinal_object = MSImagingExperiment,
#'   sample_name = "tissue_sample1",
#'   section_name = "slice_001",
#'   resolution = "50um",
#'   matrix_name = "NEDC"
#' )
#' 
#' # View the object
#' msi_obj
#' }
#' 
convert_cardinal2msi_dataset <- function(cardinal_object,
                                         tolerance = 10,
                                         units = "ppm",
                                         sample_name = NULL,
                                         section_name = "section1",
                                         resolution = "unknown",
                                         matrix_name = "unknown") {
  
  cat("\n========================================================\n")
  cat("  Cardinal → msi_dataset Conversion (auto-aligned)\n")
  cat("========================================================\n\n")
  
  # ---- Step 1: Peak Alignment ----
  # MSImagingArrays are unaligned and must establish a common m/z axis
  if (is(cardinal_object, "MSImagingArrays")) {
    cat("[1/7] MSImagingArrays detected, performing automatic alignment...\n")
    cat("      Tolerance:", tolerance, units, "\n")
    
    # This step converts Arrays to aligned Experiment object
    cardinal_object <- Cardinal::peakAlign(cardinal_object, 
                                           tolerance = tolerance, 
                                           units = units)
    cat("  ✓ Alignment completed, converted to MSImagingExperiment\n\n")
  } else {
    cat("[1/7] Input object is already aligned (MSImagingExperiment)\n")
    cat("      Skipping alignment step\n\n")
  }
  
  # ---- Step 2: Extract Expression Matrix ----
  cat("[2/7] Extracting expression matrix...\n")
  # After alignment, spectra() returns a regular matrix
  expr_matrix <- as.matrix(Cardinal::spectra(cardinal_object))
  
  cat("      Dimensions:", nrow(expr_matrix), "variables ×", 
      ncol(expr_matrix), "pixels\n")
  cat("  ✓ Expression matrix extracted\n\n")
  
  # ---- Step 3: Extract Common m/z Reference Axis ----
  cat("[3/7] Extracting common m/z reference axis...\n")
  # After alignment, mz(object) is a unique feature vector
  mz_values <- as.numeric(Cardinal::mz(cardinal_object))
  
  if (length(mz_values) != nrow(expr_matrix)) {
    stop("Length of m/z axis does not match number of matrix rows. Check alignment.")
  }
  
  # Create variable_info data frame
  variable_info <- data.frame(
    variable_id = paste0("mz_", sprintf("%05d", 1:length(mz_values))),
    mz = mz_values,
    ms2 = "no",  # Default to "no", can be updated later
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat("      m/z range: [", round(min(mz_values), 4), ",", 
      round(max(mz_values), 4), "]\n")
  cat("      Number of variables:", nrow(variable_info), "\n")
  cat("  ✓ variable_info created\n\n")
  
  # ---- Step 4: Extract Spatial Coordinates (3D supported) ----
  cat("[4/7] Extracting spatial coordinates...\n")
  pdata <- Cardinal::pixelData(cardinal_object)
  
  # Automatically detect x, y, z coordinates
  coord_names <- intersect(c("x", "y", "z"), names(pdata))
  coord_data <- as.data.frame(pdata[, coord_names, drop = FALSE])
  
  # Force conversion of S4 vectors to standard numeric
  for (col in names(coord_data)) {
    coord_data[[col]] <- as.numeric(as.vector(coord_data[[col]]))
  }
  
  # Ensure z column exists (required by msi_dataset structure)
  if (!"z" %in% colnames(coord_data)) {
    coord_data$z <- 1  # Default z slice
  }
  
  cat("      Spatial dimensions:", length(coord_names), "D\n")
  cat("      x range: [", min(coord_data$x), ",", max(coord_data$x), "]\n")
  cat("      y range: [", min(coord_data$y), ",", max(coord_data$y), "]\n")
  if ("z" %in% coord_names) {
    cat("      z range: [", min(coord_data$z), ",", max(coord_data$z), "]\n")
  }
  cat("  ✓ Spatial coordinates extracted\n\n")
  
  # ---- Step 5: Extract Sample Information ----
  cat("[5/7] Extracting sample information...\n")
  
  # Determine sample_id
  if (!is.null(sample_name)) {
    sample_ids <- rep(sample_name, ncol(expr_matrix))
  } else if ("run" %in% names(pdata)) {
    sample_ids <- as.character(pdata$run)
  } else {
    sample_ids <- rep("sample1", ncol(expr_matrix))
  }
  
  # Create pixel_info data frame
  pixel_info <- data.frame(
    pixel_id = paste0("pixel_", sprintf("%05d", 1:ncol(expr_matrix))),
    x = coord_data$x,
    y = coord_data$y,
    z = coord_data$z,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat("      Number of pixels:", nrow(pixel_info), "\n")
  cat("      Number of unique samples:", length(unique(sample_ids)), "\n")
  cat("  ✓ pixel_info created\n\n")
  
  # ---- Step 6: Create sample_info ----
  cat("[6/7] Creating sample_info...\n")
  
  # Get unique z slices to create sample_info
  unique_z <- unique(coord_data$z)
  
  sample_info <- data.frame(
    sample_id = unique(sample_ids),
    section_id = section_name,
    z = unique_z,
    resolution = resolution,
    matrix = matrix_name,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # If there are multiple z slices, create multiple sample entries
  if (length(unique_z) > 1 && length(unique(sample_ids)) == 1) {
    sample_info <- data.frame(
      sample_id = paste0(sample_info$sample_id[1], "_z", unique_z),
      section_id = paste0(section_name, "_z", unique_z),
      z = unique_z,
      resolution = resolution,
      matrix = matrix_name,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  
  cat("      Number of samples:", nrow(sample_info), "\n")
  cat("  ✓ sample_info created\n\n")
  
  # ---- Step 7: Construct msi_dataset Object ----
  cat("[7/7] Constructing msi_dataset object...\n")
  
  # Set row and column names for expression_data
  rownames(expr_matrix) <- variable_info$variable_id
  colnames(expr_matrix) <- pixel_info$pixel_id
  
  # Convert to data.frame
  expr_data_df <- as.data.frame(expr_matrix)
  
  # Create info_note data frames
  pixel_info_note <- data.frame(
    name = colnames(pixel_info),
    meaning = c(
      "Pixel ID",
      "Horizontal spatial coordinate in the matrix",
      "Vertical spatial coordinate in the image",
      "Tissue slice depth / relative to a reference slice"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  sample_info_note <- data.frame(
    name = colnames(sample_info),
    meaning = c(
      "Unique identifier for the biological sample",
      "Identifier for the specific tissue section/slice",
      "Tissue slice depth / relative to a reference slice",
      "Spatial resolution of the acquisition (e.g., 50µm)",
      "Chemical matrix used for ionization"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  variable_info_note <- data.frame(
    name = colnames(variable_info),
    meaning = c(
      "Variable ID",
      "Mass to charge ratio",
      "Yes or no"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Create the msi_dataset object
  msi_obj <- create_msi_dataset(
    expression_data = expr_data_df,
    pixel_info = pixel_info,
    sample_info = sample_info,
    variable_info = variable_info,
    pixel_info_note = pixel_info_note,
    sample_info_note = sample_info_note,
    variable_info_note = variable_info_note
  )
  
  cat("========================================================\n")
  cat(" ✓ Conversion successful!\n")
  cat("========================================================\n")
  cat("\nSummary:\n")
  cat("  - Variables:", nrow(variable_info), "\n")
  cat("  - Pixels:", nrow(pixel_info), "\n")
  cat("  - Samples:", nrow(sample_info), "\n")
  cat("  - m/z range: [", round(min(mz_values), 4), ",", 
      round(max(mz_values), 4), "]\n\n")
  
  return(msi_obj)
}