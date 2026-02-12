#' Convert mass_dataset back to msi_dataset
#'
#' @description
#' Converts a mass_dataset object back to an msi_dataset object, recovering
#' MSI-specific information (pixel coordinates, sample metadata) that was
#' saved during the original conversion.
#'
#' @param mass_object A mass_dataset object (originally converted from msi_dataset)
#' @param metadata_source Source of MSI metadata: "auto" (automatic detection),
#'   "env" (from environment cache), or "file" (from RDS file)
#' @param metadata_path Path to metadata RDS file (required if metadata_source = "file")
#' @param object_id Conversion object ID (auto-detected from process_info if NULL)
#' @param create_annotation_table Logical, whether to create annotation_table from
#'   mass_dataset results (default: TRUE)
#'
#' @return An msi_dataset object
#'
#' @details
#' This function performs the reverse conversion:
#' - Each "sample" becomes a pixel in msi_dataset
#' - sample_id → pixel_id
#' - Spatial information (x, y, z) is recovered from saved metadata
#' - MSI-specific sample_info is restored
#' - Any annotations from mass_dataset analysis are incorporated
#'
#' The function will attempt to automatically find the metadata from:
#' 1. Environment cache (.msi_conversion_cache)
#' 2. process_info in the mass_dataset object
#' 3. Specified file path
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Automatic metadata recovery (same R session)
#' my_msi_restored <- convert_mass2msi_dataset(my_mass)
#'
#' # Load from file (different R session)
#' my_msi_restored <- convert_mass2msi_dataset(
#'   my_mass,
#'   metadata_source = "file",
#'   metadata_path = "msi_metadata.rds"
#' )
#'
#' # Specify object ID explicitly
#' my_msi_restored <- convert_mass2msi_dataset(
#'   my_mass,
#'   object_id = "msi_20250209_143022"
#' )
#' }
#'
convert_mass2msi_dataset <- function(mass_object,
                                     metadata_source = c("auto", "env", "file"),
                                     metadata_path = NULL,
                                     object_id = NULL,
                                     create_annotation_table = TRUE) {
  
  cat("\n========================================\n")
  cat("  Mass Dataset → MSI Dataset\n")
  cat("========================================\n\n")
  
  # ---- Validate input ----
  if (!inherits(mass_object, "mass_dataset")) {
    stop("Input must be a mass_dataset object")
  }
  
  metadata_source <- match.arg(metadata_source)
  
  cat("[1/7] Validating input...\n")
  cat("  Input class: mass_dataset\n")
  cat("  Variables:", nrow(mass_object@variable_info), "\n")
  cat("  Samples:", nrow(mass_object@sample_info), "\n\n")
  
  # ---- Retrieve metadata ----
  cat("[2/7] Retrieving MSI metadata...\n")
  
  msi_metadata <- NULL
  
  # Try to find object_id from process_info
  if (is.null(object_id)) {
    if ("msi_conversion" %in% names(mass_object@process_info)) {
      conv_info <- mass_object@process_info$msi_conversion
      if (!is.null(conv_info@parameter$conversion_id)) {
        object_id <- conv_info@parameter$conversion_id
        cat("  ✓ Found conversion ID in process_info:", object_id, "\n")
      }
    }
  }
  
  # Retrieve metadata based on source
  if (metadata_source == "auto") {
    # Try environment first
    if (exists(".msi_conversion_cache", envir = .GlobalEnv)) {
      cache_env <- get(".msi_conversion_cache", envir = .GlobalEnv)
      if (!is.null(object_id) && exists(object_id, envir = cache_env)) {
        msi_metadata <- cache_env[[object_id]]
        cat("  ✓ Metadata retrieved from environment cache\n")
      }
    }
    
    # If not found in environment, try file
    if (is.null(msi_metadata) && !is.null(metadata_path)) {
      if (file.exists(metadata_path)) {
        msi_metadata <- readRDS(metadata_path)
        cat("  ✓ Metadata loaded from file:", metadata_path, "\n")
      }
    }
    
  } else if (metadata_source == "env") {
    if (!exists(".msi_conversion_cache", envir = .GlobalEnv)) {
      stop("No conversion cache found in environment")
    }
    cache_env <- get(".msi_conversion_cache", envir = .GlobalEnv)
    
    if (is.null(object_id)) {
      # List available IDs
      available_ids <- ls(cache_env)
      if (length(available_ids) == 0) {
        stop("No metadata found in cache")
      }
      cat("  Available object IDs:", paste(available_ids, collapse = ", "), "\n")
      stop("Please specify object_id parameter")
    }
    
    if (!exists(object_id, envir = cache_env)) {
      stop("Object ID '", object_id, "' not found in cache")
    }
    
    msi_metadata <- cache_env[[object_id]]
    cat("  ✓ Metadata retrieved from environment\n")
    
  } else if (metadata_source == "file") {
    if (is.null(metadata_path)) {
      stop("metadata_path must be provided when metadata_source = 'file'")
    }
    if (!file.exists(metadata_path)) {
      stop("Metadata file not found: ", metadata_path)
    }
    msi_metadata <- readRDS(metadata_path)
    cat("  ✓ Metadata loaded from file:", metadata_path, "\n")
  }
  
  if (is.null(msi_metadata)) {
    stop("Failed to retrieve MSI metadata. Please provide metadata_path or ensure ",
         "metadata is in environment cache.")
  }
  
  cat("\n")
  
  # ---- Validate metadata ----
  cat("[3/7] Validating metadata...\n")
  
  # Check dimensions match
  if (nrow(mass_object@sample_info) != msi_metadata$n_pixels) {
    warning("Sample count mismatch. Expected ", msi_metadata$n_pixels, 
            " pixels, got ", nrow(mass_object@sample_info), " samples. ",
            "This may occur if samples were filtered in mass_dataset.")
  }
  
  if (nrow(mass_object@variable_info) != msi_metadata$n_variables) {
    warning("Variable count mismatch. Expected ", msi_metadata$n_variables,
            " variables, got ", nrow(mass_object@variable_info), " variables. ",
            "This may occur if variables were filtered in mass_dataset.")
  }
  
  cat("  ✓ Metadata validated\n\n")
  
  # ---- Extract data from mass_dataset ----
  cat("[4/7] Extracting data from mass_dataset...\n")
  
  expr_data <- mass_object@expression_data
  mass_sample_info <- mass_object@sample_info
  mass_variable_info <- mass_object@variable_info
  
  cat("  ✓ Expression data extracted\n")
  cat("  ✓ Sample info extracted\n")
  cat("  ✓ Variable info extracted\n\n")
  
  # ---- Reconstruct pixel_info ----
  cat("[5/7] Reconstructing pixel_info...\n")
  
  # Get original pixel_info
  pixel_inf <- msi_metadata$pixel_info
  
  # Map sample_id back to pixel_id
  # Handle case where samples may have been filtered
  pixel_ids_in_mass <- mass_sample_info$sample_id
  
  # Filter pixel_info to match current samples
  pixel_inf_filtered <- pixel_inf[pixel_inf$pixel_id %in% pixel_ids_in_mass, ]
  
  # Ensure order matches expression_data
  pixel_inf_filtered <- pixel_inf_filtered[
    match(pixel_ids_in_mass, pixel_inf_filtered$pixel_id), 
  ]
  
  cat("  Original pixels:", nrow(msi_metadata$pixel_info), "\n")
  cat("  Current pixels:", nrow(pixel_inf_filtered), "\n")
  cat("  ✓ pixel_info reconstructed\n\n")
  
  # ---- Reconstruct sample_info ----
  cat("[6/7] Reconstructing sample_info...\n")
  
  # Use original sample_info from metadata
  sample_inf <- msi_metadata$sample_info
  
  cat("  ✓ sample_info reconstructed\n\n")
  
  # ---- Reconstruct variable_info ----
  cat("[7/7] Reconstructing variable_info...\n")
  
  # Create variable_info for msi_dataset
  msi_variable_info <- data.frame(
    variable_id = mass_variable_info$variable_id,
    mz = mass_variable_info$mz,
    ms2 = "no",  # Default value
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat("  ✓ variable_info reconstructed\n\n")
  
  # ---- Create annotation_table ----
  if (create_annotation_table) {
    cat("Creating annotation_table from mass_dataset results...\n")
    
    # Check if annotation exists in variable_info
    annotation_cols <- setdiff(
      colnames(mass_variable_info),
      c("variable_id", "mz", "rt")
    )
    
    if (length(annotation_cols) > 0) {
      annotation_table <- mass_variable_info[, c("variable_id", annotation_cols), 
                                             drop = FALSE]
      cat("  ✓ Annotation columns found:", paste(annotation_cols, collapse = ", "), "\n")
    } else {
      annotation_table <- data.frame(
        variable_id = character(),
        stringsAsFactors = FALSE
      )
      cat("  No annotation columns found, created empty table\n")
    }
  } else {
    annotation_table <- data.frame(
      variable_id = character(),
      stringsAsFactors = FALSE
    )
  }
  
  cat("\n")
  
  # ---- Restore info_note data frames ----
  pixel_inf_note <- msi_metadata$pixel_info_note
  sample_inf_note <- msi_metadata$sample_info_note
  
  # Reconstruct variable_info_note
  variable_inf_note <- data.frame(
    name = c("variable_id", "mz", "ms2"),
    meaning = c(
      "Variable ID",
      "Mass to charge ratio",
      "Yes or no"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # ---- Create msi_dataset object ----
  cat("Creating msi_dataset object...\n")
  
  msi_obj <- create_msi_dataset(
    expression_data = expr_data,
    pixel_info = pixel_inf_filtered,
    sample_info = sample_inf,
    variable_info = msi_variable_info,
    pixel_info_note = pixel_inf_note,
    sample_info_note = sample_inf_note,
    variable_info_note = variable_inf_note
  )
  
  # Add annotation_table
  msi_obj@annotation_table <- annotation_table
  
  # Merge process_info
  if (length(mass_object@process_info) > 0) {
    msi_obj@process_info <- c(
      msi_obj@process_info,
      mass_object@process_info
    )
  }
  
  # Add conversion info
  conversion_parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "convert_mass2msi_dataset()",
    parameter = list(
      conversion_id = msi_metadata$object_id,
      metadata_source = metadata_source
    ),
    time = Sys.time()
  )
  
  msi_obj@process_info$mass_conversion <- conversion_parameter
  
  cat("\n========================================\n")
  cat(" ✓ Conversion complete!\n")
  cat("========================================\n\n")
  
  cat("Summary:\n")
  cat("  Input:  mass_dataset (", nrow(mass_variable_info), "variables ×",
      nrow(mass_sample_info), "samples)\n")
  cat("  Output: msi_dataset  (", n_variables(msi_obj), "variables ×",
      n_pixels(msi_obj), "pixels)\n")
  cat("  Conversion ID:", msi_metadata$object_id, "\n")
  if (nrow(annotation_table) > 0 && ncol(annotation_table) > 1) {
    cat("  Annotation columns:", ncol(annotation_table) - 1, "\n")
  }
  cat("\n")
  
  return(msi_obj)
}


#' Extract Annotation Table from mass_dataset
#'
#' @description
#' Extracts annotation information from a mass_dataset object and formats it
#' as an annotation_table suitable for msi_dataset.
#'
#' @param mass_object A mass_dataset object
#' @param annotation_cols Optional vector of column names to extract.
#'   If NULL, all non-standard columns are extracted.
#'
#' @return A data frame with variable_id and annotation columns
#' @export
extract_annotation_from_mass <- function(mass_object, 
                                         annotation_cols = NULL) {
  
  variable_info <- mass_object@variable_info
  
  # Standard columns to exclude
  standard_cols <- c("variable_id", "mz", "rt")
  
  if (is.null(annotation_cols)) {
    # Extract all non-standard columns
    annotation_cols <- setdiff(colnames(variable_info), standard_cols)
  }
  
  if (length(annotation_cols) == 0) {
    warning("No annotation columns found")
    return(data.frame(
      variable_id = variable_info$variable_id,
      stringsAsFactors = FALSE
    ))
  }
  
  # Create annotation table
  annotation_table <- variable_info[, c("variable_id", annotation_cols), 
                                    drop = FALSE]
  
  cat("Extracted annotation columns:\n")
  cat("  ", paste(annotation_cols, collapse = ", "), "\n")
  cat("  Total annotations:", nrow(annotation_table), "\n")
  
  return(annotation_table)
}


#' Get Conversion Summary
#'
#' @description
#' Provides a summary of the conversion between msi_dataset and mass_dataset.
#'
#' @param object Either an msi_dataset or mass_dataset object
#'
#' @return A list with conversion information
#' @export
get_conversion_summary <- function(object) {
  
  if (!inherits(object, c("msi_dataset", "mass_dataset"))) {
    stop("Object must be either msi_dataset or mass_dataset")
  }
  
  # Check for conversion info in process_info
  proc_info <- object@process_info
  
  conversions <- list()
  
  if ("msi_conversion" %in% names(proc_info)) {
    conversions$msi_to_mass <- list(
      conversion_id = proc_info$msi_conversion@parameter$conversion_id,
      time = proc_info$msi_conversion@time,
      metadata_saved = proc_info$msi_conversion@parameter$metadata_saved,
      spatial_added = proc_info$msi_conversion@parameter$spatial_added
    )
  }
  
  if ("mass_conversion" %in% names(proc_info)) {
    conversions$mass_to_msi <- list(
      conversion_id = proc_info$mass_conversion@parameter$conversion_id,
      time = proc_info$mass_conversion@time,
      metadata_source = proc_info$mass_conversion@parameter$metadata_source
    )
  }
  
  if (length(conversions) == 0) {
    cat("No conversion history found in this object\n")
    return(invisible(NULL))
  }
  
  cat("Conversion History:\n")
  cat("===================\n\n")
  
  if (!is.null(conversions$msi_to_mass)) {
    cat("MSI → Mass conversion:\n")
    cat("  Conversion ID:", conversions$msi_to_mass$conversion_id, "\n")
    cat("  Time:", as.character(conversions$msi_to_mass$time), "\n")
    cat("  Metadata saved:", conversions$msi_to_mass$metadata_saved, "\n")
    cat("  Spatial info added:", conversions$msi_to_mass$spatial_added, "\n\n")
  }
  
  if (!is.null(conversions$mass_to_msi)) {
    cat("Mass → MSI conversion:\n")
    cat("  Conversion ID:", conversions$mass_to_msi$conversion_id, "\n")
    cat("  Time:", as.character(conversions$mass_to_msi$time), "\n")
    cat("  Metadata source:", conversions$mass_to_msi$metadata_source, "\n\n")
  }
  
  invisible(conversions)
}