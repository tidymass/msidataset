#' Convert Cardinal MSImagingArrays/MSImageSet to msi_dataset
#'
#' @description
#' Converts a Cardinal MSImageSet, MSImagingExperiment, or MSImagingArrays object 
#' to the new msi_dataset structure that is compatible with massdataset. 
#' 
#' @param cardinal_object A Cardinal object
#' @param sample_id Optional sample identifier
#' @param section_id Optional section identifier
#' @param resolution Optional spatial resolution in micrometers
#' @param matrix Optional matrix used for ionization
#' @return An msi_dataset object
#' @export

convert_cardinal2msi_dataset <- function(cardinal_object,
                                         sample_id = "sample_1",
                                         section_id = "section_1", 
                                         resolution = NA,
                                         matrix = NA) {
  
  if (!requireNamespace("Cardinal", quietly = TRUE)) {
    stop("Package 'Cardinal' is required but not installed.")
  }
  
  # ========== Check object type ==========
  is_msimageset <- inherits(cardinal_object, "MSImageSet")
  is_msimagingexp <- inherits(cardinal_object, "MSImagingExperiment")
  is_msimarrays <- inherits(cardinal_object, "MSImagingArrays")
  
  is_valid_cardinal <- is_msimageset || is_msimagingexp || is_msimarrays
  
  if (!is_valid_cardinal) {
    stop(
      "cardinal_object must be a Cardinal imaging object.\n",
      "  Supported types: MSImageSet, MSImagingExperiment, MSImagingArrays\n",
      "  Provided object class: ", paste(class(cardinal_object), collapse = ", "), "\n",
      "  Please ensure you are using a valid Cardinal object from readMSIData()."
    )
  }
  
  if (is_msimarrays) {
    cat("Converting MSImagingArrays object (Cardinal 3.x)...\n")
  } else if (is_msimagingexp) {
    cat("Converting MSImagingExperiment object (Cardinal 2.x)...\n")
  } else {
    cat("Converting MSImageSet object (Cardinal < 2.0)...\n")
  }
  
  # ========== Extract expression data ==========
  if (is_msimarrays || is_msimagingexp) {
    expression_data <- Cardinal::spectra(cardinal_object)
  } else {
    expression_data <- Cardinal::iData(cardinal_object)
  }
  
  # Convert to matrix if needed
  if (methods::is(expression_data, "matter_mat") || 
      methods::is(expression_data, "matter_matc") ||
      methods::is(expression_data, "matter_arr")) {
    cat("Converting matter array to standard matrix...\n")
    expression_data <- as.matrix(expression_data)
  }
  
  if (methods::is(expression_data, "sparseMatrix") ||
      methods::is(expression_data, "dgCMatrix") ||
      methods::is(expression_data, "dgTMatrix")) {
    cat("Converting sparse matrix to dense matrix...\n")
    expression_data <- as.matrix(expression_data)
  }
  
  if (methods::is(expression_data, "SimpleList") || 
      (is.list(expression_data) && !is.data.frame(expression_data))) {
    cat("Converting list of matrices to single matrix...\n")
    expression_data <- expression_data[[1]]
    if (!is.matrix(expression_data)) {
      expression_data <- as.matrix(expression_data)
    }
  }
  
  cat("Expression data dimensions:", nrow(expression_data), "x", ncol(expression_data), "\n")
  
  # ========== Extract spatial coordinates ==========
  coord_data <- Cardinal::coord(cardinal_object)
  
  if (is.data.frame(coord_data)) {
    coords_df <- coord_data
  } else if (methods::is(coord_data, "DataFrame") || 
             methods::is(coord_data, "DFrame")) {
    coords_df <- as.data.frame(coord_data)
  } else if (is.matrix(coord_data)) {
    coords_df <- as.data.frame(coord_data)
  } else {
    stop("Unexpected coordinate data format: ", paste(class(coord_data), collapse = ", "))
  }
  
  if (!"x" %in% colnames(coords_df)) {
    stop("Coordinate data must contain 'x' column")
  }
  if (!"y" %in% colnames(coords_df)) {
    stop("Coordinate data must contain 'y' column")
  }
  
  z_vals <- if ("z" %in% colnames(coords_df)) coords_df$z else rep(1, nrow(coords_df))
  
  pixel_ids <- paste0("pixel_", 
                      sprintf("%04d", coords_df$x), "_",
                      sprintf("%04d", coords_df$y), "_",
                      sprintf("%02d", z_vals))
  
  # ========== Create sample_info ==========
  sample_info <- data.frame(
    sample_id = pixel_ids,
    x = coords_df$x,
    y = coords_df$y,
    z = z_vals,
    stringsAsFactors = FALSE
  )
  
  sample_info$section_id <- section_id
  if (!is.na(resolution)) {
    sample_info$resolution <- resolution
  }
  if (!is.na(matrix)) {
    sample_info$matrix <- matrix
  }
  
  tryCatch({
    pdata <- Cardinal::pixelData(cardinal_object)
    if (methods::is(pdata, "DataFrame") || methods::is(pdata, "DFrame")) {
      pdata <- as.data.frame(pdata)
    }
    if ("run" %in% colnames(pdata)) {
      sample_info$run <- pdata$run
    }
  }, error = function(e) {
    cat("Note: Could not extract pixelData (not critical)\n")
  })
  
  # ========== Extract m/z values (FIXED: use isS4 correctly) ==========
  cat("Extracting m/z values...\n")
  
  mz_values <- Cardinal::mz(cardinal_object)
  
  # Handle NumericList
  if (methods::is(mz_values, "NumericList")) {
    cat("  m/z format: NumericList\n")
    mz_values <- mz_values[[1]]
  }
  
  # Handle SimpleList
  if (methods::is(mz_values, "SimpleList")) {
    cat("  m/z format: SimpleList\n")
    mz_values <- mz_values[[1]]
  }
  
  # Handle Rle (run-length encoding)
  if (methods::is(mz_values, "Rle")) {
    cat("  m/z format: Rle (run-length encoded)\n")
    mz_values <- as.vector(mz_values)
  }
  
  # Check if it's an S4 object - FIXED: use isS4() without namespace
  if (isS4(mz_values)) {
    cat("  m/z format: S4 object (", paste(class(mz_values), collapse = ", "), ")\n", sep = "")
    
    # Try to access .Data slot (common in S4Vectors)
    has_data_slot <- tryCatch({
      methods::.hasSlot(mz_values, ".Data")
    }, error = function(e) {
      FALSE
    })
    
    if (has_data_slot) {
      cat("  Extracting from .Data slot\n")
      mz_values <- slot(mz_values, ".Data")
    } else {
      # Try generic conversion methods
      cat("  Attempting generic S4 conversion\n")
      converted <- FALSE
      
      # Try as.vector
      if (!converted) {
        tryCatch({
          mz_values <- as.vector(mz_values)
          converted <- TRUE
          cat("  Converted using as.vector()\n")
        }, error = function(e) {})
      }
      
      # Try as.matrix then as.vector
      if (!converted) {
        tryCatch({
          mz_values <- as.vector(as.matrix(mz_values))
          converted <- TRUE
          cat("  Converted using as.vector(as.matrix())\n")
        }, error = function(e) {})
      }
      
      # Try direct slot access
      if (!converted) {
        slot_names <- methods::slotNames(mz_values)
        if (length(slot_names) > 0) {
          cat("  Available slots:", paste(slot_names, collapse = ", "), "\n")
        }
        stop(
          "Cannot extract m/z values from S4 object.\n",
          "  Class: ", paste(class(mz_values), collapse = ", "), "\n",
          "  Please report this issue with your Cardinal version."
        )
      }
    }
  }
  
  # Final conversion to numeric vector
  if (!is.numeric(mz_values)) {
    tryCatch({
      mz_values <- as.numeric(mz_values)
    }, error = function(e) {
      stop(
        "Failed to convert m/z values to numeric.\n",
        "  Current class: ", paste(class(mz_values), collapse = ", "), "\n",
        "  Error: ", e$message
      )
    })
  }
  
  # Ensure it's a plain vector
  if (is.matrix(mz_values)) {
    mz_values <- as.vector(mz_values)
  }
  
  cat("Number of m/z values:", length(mz_values), "\n")
  cat("m/z range:", round(min(mz_values), 2), "-", round(max(mz_values), 2), "\n")
  
  # ========== Create variable_info ==========
  variable_ids <- paste0("mz_", sprintf("%.4f", mz_values))
  
  variable_info <- data.frame(
    variable_id = variable_ids,
    mz = mz_values,
    ms2 = FALSE,
    stringsAsFactors = FALSE
  )
  
  # ========== Verify dimensions ==========
  if (nrow(expression_data) != nrow(variable_info)) {
    warning(
      "Dimension mismatch!\n",
      "  Expression data rows: ", nrow(expression_data), "\n",
      "  Variable info rows: ", nrow(variable_info), "\n",
      "  Adjusting to match..."
    )
    min_vars <- min(nrow(expression_data), nrow(variable_info))
    expression_data <- expression_data[1:min_vars, , drop = FALSE]
    variable_info <- variable_info[1:min_vars, , drop = FALSE]
  }
  
  if (ncol(expression_data) != nrow(sample_info)) {
    warning(
      "Dimension mismatch!\n",
      "  Expression data cols: ", ncol(expression_data), "\n",
      "  Sample info rows: ", nrow(sample_info), "\n",
      "  Adjusting to match..."
    )
    min_samples <- min(ncol(expression_data), nrow(sample_info))
    expression_data <- expression_data[, 1:min_samples, drop = FALSE]
    sample_info <- sample_info[1:min_samples, , drop = FALSE]
  }
  
  # ========== Set dimnames ==========
  colnames(expression_data) <- sample_info$sample_id
  rownames(expression_data) <- variable_info$variable_id
  
  # ========== Create info notes ==========
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
  
  if (!is.na(resolution)) {
    sample_info_note <- rbind(
      sample_info_note,
      data.frame(
        name = "resolution",
        meaning = paste0("Spatial resolution of the acquisition (", resolution, "μm)"),
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
  cat("Creating msi_dataset object...\n")
  msi_obj <- create_msi_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info,
    sample_info_note = sample_info_note,
    variable_info_note = variable_info_note
  )
  
  # ========== Add metadata ==========
  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "convert_cardinal2msi_dataset()",
    parameter = list(
      sample_id = sample_id,
      section_id = section_id,
      resolution = resolution,
      matrix = matrix,
      cardinal_class = paste(class(cardinal_object), collapse = ", ")
    ),
    time = Sys.time()
  )
  
  msi_obj@process_info$convert_from_cardinal <- parameter
  
  cat("\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("✓ Conversion completed successfully!\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("  Cardinal object type:", class(cardinal_object)[1], "\n")
  cat("  Variables:", nrow(variable_info), "\n")
  cat("  Samples (pixels):", nrow(sample_info), "\n")
  cat("  m/z range:", round(min(mz_values), 2), "-", round(max(mz_values), 2), "\n")
  cat("  Spatial extent: x[", min(sample_info$x), "-", max(sample_info$x), "], ",
      "y[", min(sample_info$y), "-", max(sample_info$y), "]\n", sep = "")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  return(msi_obj)
}


#' Convert msi_dataset to Cardinal MSImagingExperiment
#' @param msi_object An msi_dataset object
#' @return A Cardinal MSImagingExperiment object
#' @export

convert_msi_dataset2cardinal <- function(msi_object) {
  
  if (!requireNamespace("Cardinal", quietly = TRUE)) {
    stop("Package 'Cardinal' is required but not installed.")
  }
  
  expr_data <- expression_data(msi_object)
  samp_info <- sample_info(msi_object)
  var_info <- variable_info(msi_object)
  
  if (!all(c("x", "y") %in% colnames(samp_info))) {
    stop("sample_info must contain 'x' and 'y' columns.")
  }
  
  z_vals <- if ("z" %in% colnames(samp_info)) samp_info$z else rep(1, nrow(samp_info))
  
  coord_df <- data.frame(
    x = samp_info$x,
    y = samp_info$y,
    z = z_vals
  )
  
  if ("mz" %in% colnames(var_info)) {
    mz_values <- var_info$mz
  } else {
    stop("variable_info must contain 'mz' column")
  }
  
  cardinal_obj <- Cardinal::MSImagingExperiment(
    imageData = expr_data,
    featureData = Cardinal::MassDataFrame(mz = mz_values),
    pixelData = Cardinal::PositionDataFrame(coord = coord_df)
  )
  
  return(cardinal_obj)
}


#' Batch convert multiple Cardinal objects
#' @param cardinal_list List of Cardinal objects
#' @param sample_ids Character vector of sample identifiers
#' @param section_ids Character vector of section identifiers
#' @param resolution Numeric vector of spatial resolutions
#' @param matrix Character vector of matrix names
#' @return List of msi_dataset objects
#' @export

batch_convert_cardinal2msi_dataset <- function(cardinal_list,
                                               sample_ids = NULL,
                                               section_ids = NULL,
                                               resolution = NA,
                                               matrix = NA) {
  
  n_objects <- length(cardinal_list)
  
  if (is.null(sample_ids)) {
    sample_ids <- paste0("sample_", seq_len(n_objects))
  }
  if (is.null(section_ids)) {
    section_ids <- paste0("section_", seq_len(n_objects))
  }
  
  if (length(resolution) == 1) {
    resolution <- rep(resolution, n_objects)
  }
  if (length(matrix) == 1) {
    matrix <- rep(matrix, n_objects)
  }
  
  msi_list <- lapply(seq_len(n_objects), function(i) {
    cat("\n--- Converting object", i, "of", n_objects, "---\n")
    convert_cardinal2msi_dataset(
      cardinal_object = cardinal_list[[i]],
      sample_id = sample_ids[i],
      section_id = section_ids[i],
      resolution = resolution[i],
      matrix = matrix[i]
    )
  })
  
  message("\nMultiple msi_dataset objects created.")
  message("Consider using combine_msi_dataset() to merge them.")
  return(msi_list)
}