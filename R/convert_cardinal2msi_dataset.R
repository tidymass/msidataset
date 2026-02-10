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
#' @param variable_id_mode Character string specifying variable ID naming mode:
#'   - "mz_based" (default): IDs based on m/z integer values (e.g., m50_1, m50_2, m51_1)
#'   - "sequential": Sequential numbering (e.g., mz_00001, mz_00002)
#' @param variable_id_prefix Prefix for variable IDs when using mz_based mode. Default is "m"
#' @param variable_id_separator Separator between mass and count. Default is "_"
#' @param variable_id_zero_pad Logical, whether to zero-pad count numbers. Default is FALSE
#' 
#' @return An msi_dataset object
#' 
#' @details
#' This function performs the following steps:
#' 1. Automatic peak alignment for MSImagingArrays objects (10ppm tolerance by default)
#' 2. Extraction of expression matrix (variables × pixels)
#' 3. Extraction of common m/z reference axis
#' 4. Generation of variable IDs based on m/z values
#' 5. Extraction of spatial coordinates (supports x, y, z for 3D data)
#' 6. Construction of pixel_info and sample_info data frames
#' 7. Creation of the final msi_dataset object
#' 
#' Variable ID naming modes:
#' - mz_based: IDs reflect the integer m/z value with sequential numbering within each mass
#'   Example: m/z values 50.001, 50.002, 51.001 become m50_1, m50_2, m51_1
#' - sequential: Simple sequential numbering (legacy mode)
#'   Example: mz_00001, mz_00002, mz_00003
#' 
#' @export
#' @examples
#' \dontrun{
#' # Load Cardinal data
#' library(Cardinal)
#' data(MSImagingExperiment)
#' 
#' # Convert with m/z-based variable IDs (recommended)
#' msi_obj <- convert_cardinal2msi_dataset_v1(
#'   cardinal_object = MSImagingExperiment,
#'   sample_name = "tissue_sample1",
#'   section_name = "slice_001",
#'   resolution = "50um",
#'   matrix_name = "DHB",
#'   variable_id_mode = "mz_based"
#' )
#' 
#' # Convert with sequential variable IDs (legacy)
#' msi_obj_seq <- convert_cardinal2msi_dataset_v1(
#'   cardinal_object = MSImagingExperiment,
#'   variable_id_mode = "sequential"
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
                                         matrix_name = "unknown",
                                         variable_id_mode = c("mz_based", "sequential"),
                                         variable_id_prefix = "m",
                                         variable_id_separator = "_",
                                         variable_id_zero_pad = FALSE) {
  
  # Match argument
  variable_id_mode <- match.arg(variable_id_mode)
  
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
  
  # ---- Generate variable_id based on selected mode ----
  cat("      Variable ID mode:", variable_id_mode, "\n")
  
  if (variable_id_mode == "mz_based") {
    # Generate IDs based on m/z values
    variable_ids <- generate_variable_id_extended(
      mz_values = mz_values,
      prefix = variable_id_prefix,
      separator = variable_id_separator,
      zero_pad = variable_id_zero_pad,
      pad_width = 3
    )
    
    # Show some examples
    n_examples <- min(5, length(variable_ids))
    cat("      Variable ID examples:\n")
    for (i in 1:n_examples) {
      cat("        m/z", round(mz_values[i], 5), "->", variable_ids[i], "\n")
    }
    
  } else {
    # Sequential mode (legacy)
    variable_ids <- paste0("mz_", sprintf("%05d", 1:length(mz_values)))
    cat("      Using sequential numbering: mz_00001, mz_00002, ...\n")
  }
  
  # Validate uniqueness
  validate_variable_id(variable_ids, mz_values)
  
  # Create variable_info data frame
  variable_info <- data.frame(
    variable_id = variable_ids,
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
      round(max(mz_values), 4), "]\n")
  cat("  - Variable ID mode:", variable_id_mode, "\n\n")
  
  return(msi_obj)
}