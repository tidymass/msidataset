#' Convert Cardinal MSImageSet to msi_dataset
#'
#' @description
#' Converts a Cardinal MSImageSet object to the new msi_dataset structure
#' that is compatible with massdataset. Pixels are stored in sample_info.
#'
#' @param cardinal_object A Cardinal MSImageSet or MSImagingExperiment object
#' @param sample_id Optional sample identifier for the dataset
#' @param section_id Optional section identifier
#' @param resolution Optional spatial resolution in micrometers
#' @param matrix Optional matrix used for ionization (e.g., "DHB", "CHCA")
#' @return An msi_dataset object
#' @export
#' @examples
#' \dontrun{
#' library(Cardinal)
#' # Load Cardinal data
#' data <- readMSIData("path/to/imzML")
#' 
#' # Convert to msi_dataset
#' msi_obj <- convert_cardinal2msi_dataset(
#'   cardinal_object = data,
#'   sample_id = "mouse_brain_1",
#'   section_id = "section_1",
#'   resolution = 50,
#'   matrix = "DHB"
#' )
#' }

convert_cardinal2msi_dataset <- function(cardinal_object,
                                         sample_id = "sample_1",
                                         section_id = "section_1", 
                                         resolution = NA,
                                         matrix = NA) {
  
  # Check Cardinal package
  if (!requireNamespace("Cardinal", quietly = TRUE)) {
    stop("Package 'Cardinal' is required but not installed.")
  }
  
  # ========== Extract expression data ==========
  # Get intensity matrix
  if (methods::is(cardinal_object, "MSImageSet")) {
    expression_data <- Cardinal::iData(cardinal_object)
  } else if (methods::is(cardinal_object, "MSImagingExperiment")) {
    expression_data <- Cardinal::spectra(cardinal_object)
  } else {
    stop("cardinal_object must be MSImageSet or MSImagingExperiment")
  }
  
  # Convert to matrix if needed
  if (methods::is(expression_data, "matter_mat")) {
    expression_data <- as.matrix(expression_data)
  }
  
  # ========== Extract spatial coordinates (now goes into sample_info) ==========
  coord_data <- Cardinal::coord(cardinal_object)
  
  # Create pixel IDs
  pixel_ids <- paste0("pixel_", 
                      sprintf("%04d", coord_data$x), "_",
                      sprintf("%04d", coord_data$y), "_",
                      sprintf("%02d", coord_data$z))
  
  # ========== Create sample_info (contains pixel spatial data) ==========
  sample_info <- data.frame(
    sample_id = pixel_ids,
    x = coord_data$x,
    y = coord_data$y,
    z = if ("z" %in% colnames(coord_data)) coord_data$z else 1,
    stringsAsFactors = FALSE
  )
  
  # Add optional metadata
  sample_info$section_id <- section_id
  if (!is.na(resolution)) {
    sample_info$resolution <- resolution
  }
  if (!is.na(matrix)) {
    sample_info$matrix <- matrix
  }
  
  # Add run/acquisition info if available
  if (methods::is(cardinal_object, "MSImageSet")) {
    if ("run" %in% colnames(Cardinal::pixelData(cardinal_object))) {
      sample_info$run <- Cardinal::pixelData(cardinal_object)$run
    }
  } else if (methods::is(cardinal_object, "MSImagingExperiment")) {
    pdata <- Cardinal::pixelData(cardinal_object)
    if ("run" %in% colnames(pdata)) {
      sample_info$run <- pdata$run
    }
  }
  
  # ========== Create variable_info ==========
  # Get m/z values
  if (methods::is(cardinal_object, "MSImageSet")) {
    mz_values <- Cardinal::mz(cardinal_object)
  } else {
    mz_values <- Cardinal::mz(cardinal_object)
  }
  
  # Create variable IDs
  variable_ids <- paste0("mz_", sprintf("%.4f", mz_values))
  
  variable_info <- data.frame(
    variable_id = variable_ids,
    mz = mz_values,
    ms2 = FALSE,  # Assume MS1 data unless specified
    stringsAsFactors = FALSE
  )
  
  # ========== Set dimnames ==========
  colnames(expression_data) <- sample_info$sample_id
  rownames(expression_data) <- variable_info$variable_id
  
  # ========== Create sample_info_note ==========
  sample_info_note <- data.frame(
    name = c("sample_id", "x", "y", "z", "section_id"),
    meaning = c(
      "Pixel ID",
      "Horizontal spatial coordinate in the matrix",
      "Vertical spatial coordinate in the image", 
      "Tissue slice depth / relative to a reference slice",
      "Identifier for the specific tissue section/slice"
    ),
    stringsAsFactors = FALSE
  )
  
  # Add optional metadata notes
  if (!is.na(resolution)) {
    sample_info_note <- rbind(
      sample_info_note,
      data.frame(
        name = "resolution",
        meaning = paste0("Spatial resolution of the acquisition (e.g., ", resolution, "μm)"),
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (!is.na(matrix)) {
    sample_info_note <- rbind(
      sample_info_note,
      data.frame(
        name = "matrix",
        meaning = "Chemical matrix used for ionization",
        stringsAsFactors = FALSE
      )
    )
  }
  
  if ("run" %in% colnames(sample_info)) {
    sample_info_note <- rbind(
      sample_info_note,
      data.frame(
        name = "run",
        meaning = "Acquisition run identifier",
        stringsAsFactors = FALSE
      )
    )
  }
  
  # ========== Create variable_info_note ==========
  variable_info_note <- data.frame(
    name = c("variable_id", "mz", "ms2"),
    meaning = c(
      "Variable ID",
      "Mass to charge ratio",
      "Whether MS2 data is available (TRUE/FALSE)"
    ),
    stringsAsFactors = FALSE
  )
  
  # ========== Create msi_dataset object ==========
  msi_obj <- create_msi_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info,
    sample_info_note = sample_info_note,
    variable_info_note = variable_info_note
  )
  
  # ========== Add Cardinal metadata to process_info ==========
  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "convert_cardinal2msi_dataset()",
    parameter = list(
      sample_id = sample_id,
      section_id = section_id,
      resolution = resolution,
      matrix = matrix,
      cardinal_class = class(cardinal_object)[1]
    ),
    time = Sys.time()
  )
  
  msi_obj@process_info$convert_from_cardinal <- parameter
  
  return(msi_obj)
}


#' Convert msi_dataset to Cardinal MSImagingExperiment
#'
#' @description
#' Converts an msi_dataset object to a Cardinal MSImagingExperiment object.
#' Useful for using Cardinal-specific visualization and analysis functions.
#'
#' @param msi_object An msi_dataset object
#' @return A Cardinal MSImagingExperiment object
#' @export
#' @examples
#' \dontrun{
#' # Convert msi_dataset to Cardinal
#' cardinal_obj <- convert_msi_dataset2cardinal(msi_obj)
#' 
#' # Now use Cardinal functions
#' Cardinal::image(cardinal_obj, mz = 500)
#' }

convert_msi_dataset2cardinal <- function(msi_object) {
  
  # Check Cardinal package
  if (!requireNamespace("Cardinal", quietly = TRUE)) {
    stop("Package 'Cardinal' is required but not installed.")
  }
  
  # Extract data
  expr_data <- expression_data(msi_object)
  samp_info <- sample_info(msi_object)
  var_info <- variable_info(msi_object)
  
  # Create coordinate data frame
  coord_df <- data.frame(
    x = samp_info$x,
    y = samp_info$y,
    z = samp_info$z
  )
  
  # Get m/z values
  if ("mz" %in% colnames(var_info)) {
    mz_values <- var_info$mz
  } else {
    stop("variable_info must contain 'mz' column")
  }
  
  # Create MSImagingExperiment
  cardinal_obj <- Cardinal::MSImagingExperiment(
    imageData = expr_data,
    featureData = Cardinal::MassDataFrame(mz = mz_values),
    pixelData = Cardinal::PositionDataFrame(coord = coord_df)
  )
  
  return(cardinal_obj)
}


#' Batch convert multiple Cardinal objects to msi_dataset
#'
#' @description
#' Converts multiple Cardinal objects (e.g., multiple tissue sections) 
#' to a single msi_dataset object.
#'
#' @param cardinal_list List of Cardinal MSImageSet or MSImagingExperiment objects
#' @param sample_ids Character vector of sample identifiers (one per object)
#' @param section_ids Character vector of section identifiers (one per object)
#' @param resolution Numeric vector of spatial resolutions
#' @param matrix Character vector of matrix names
#' @return An msi_dataset object containing all sections
#' @export

batch_convert_cardinal2msi_dataset <- function(cardinal_list,
                                               sample_ids = NULL,
                                               section_ids = NULL,
                                               resolution = NA,
                                               matrix = NA) {
  
  n_objects <- length(cardinal_list)
  
  # Set default IDs if not provided
  if (is.null(sample_ids)) {
    sample_ids <- paste0("sample_", seq_len(n_objects))
  }
  if (is.null(section_ids)) {
    section_ids <- paste0("section_", seq_len(n_objects))
  }
  
  # Replicate resolution and matrix if single values
  if (length(resolution) == 1) {
    resolution <- rep(resolution, n_objects)
  }
  if (length(matrix) == 1) {
    matrix <- rep(matrix, n_objects)
  }
  
  # Convert each object
  msi_list <- lapply(seq_len(n_objects), function(i) {
    convert_cardinal2msi_dataset(
      cardinal_object = cardinal_list[[i]],
      sample_id = sample_ids[i],
      section_id = section_ids[i],
      resolution = resolution[i],
      matrix = matrix[i]
    )
  })
  
  # Combine all objects
  # This would require a combine_msi_dataset function
  # For now, return the list
  message("Multiple msi_dataset objects created. Consider using combine_msi_dataset() to merge them.")
  return(msi_list)
}