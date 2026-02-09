#' Convert msi_dataset to mass_dataset
#'
#' @description
#' Converts an msi_dataset object to a mass_dataset object for use with
#' massdataset package functions. MSI-specific information (pixel coordinates,
#' sample metadata) is preserved for later conversion back to msi_dataset.
#'
#' @param msi_object An msi_dataset object
#' @param save_metadata Logical, whether to save MSI metadata to environment (default: TRUE)
#' @param metadata_path Optional path to save metadata as RDS file
#' @param add_spatial_to_sample_info Logical, whether to add x, y, z coordinates 
#'   to sample_info as additional columns (default: FALSE)
#' @param object_id Optional unique identifier for this conversion (auto-generated if NULL)
#'
#' @return A mass_dataset object
#'
#' @details
#' This function performs the following conversions:
#' - Each pixel becomes a "sample" in mass_dataset
#' - pixel_id → sample_id
#' - pixel_info is stored for later recovery
#' - sample_info (MSI-specific) is stored for later recovery
#' - variable_info is mostly preserved (ms2 column dropped, rt column added as 0)
#'
#' The MSI-specific metadata is stored in:
#' 1. A package environment variable (.msi_conversion_cache) by default
#' 2. Optionally in an RDS file if metadata_path is provided
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic conversion
#' my_mass <- convert_msi2mass_dataset(my_msi)
#'
#' # Save metadata to file for later
#' my_mass <- convert_msi2mass_dataset(
#'   my_msi,
#'   metadata_path = "msi_metadata.rds"
#' )
#'
#' # Add spatial coordinates to sample_info
#' my_mass <- convert_msi2mass_dataset(
#'   my_msi,
#'   add_spatial_to_sample_info = TRUE
#' )
#' }
#'
convert_msi2mass_dataset <- function(msi_object,
                                     save_metadata = TRUE,
                                     metadata_path = NULL,
                                     add_spatial_to_sample_info = FALSE,
                                     object_id = NULL) {
  
  cat("\n========================================\n")
  cat("  MSI Dataset → Mass Dataset\n")
  cat("========================================\n\n")
  
  # ---- Validate input ----
  if (!inherits(msi_object, "msi_dataset")) {
    stop("Input must be an msi_dataset object")
  }
  
  # Generate unique ID for this conversion
  if (is.null(object_id)) {
    object_id <- paste0("msi_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  
  cat("[1/6] Validating input...\n")
  cat("  Input class: msi_dataset\n")
  cat("  Variables:", n_variables(msi_object), "\n")
  cat("  Pixels:", n_pixels(msi_object), "\n")
  cat("  Conversion ID:", object_id, "\n\n")
  
  # ---- Extract data from msi_dataset ----
  cat("[2/6] Extracting data from msi_dataset...\n")
  
  expr_data <- expression_data(msi_object)
  pixel_inf <- pixel_info(msi_object)
  sample_inf <- sample_info(msi_object)
  variable_inf <- variable_info(msi_object)
  
  pixel_inf_note <- pixel_info_note(msi_object)
  sample_inf_note <- sample_info_note(msi_object)
  variable_inf_note <- variable_info_note(msi_object)
  
  cat("  ✓ Expression data extracted\n")
  cat("  ✓ Pixel info extracted\n")
  cat("  ✓ Sample info extracted\n")
  cat("  ✓ Variable info extracted\n\n")
  
  # ---- Create mass_dataset sample_info ----
  cat("[3/6] Creating mass_dataset sample_info...\n")
  
  # Map pixel_id to sample_id
  mass_sample_info <- data.frame(
    sample_id = pixel_inf$pixel_id,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Optionally add spatial coordinates
  if (add_spatial_to_sample_info) {
    mass_sample_info$x <- pixel_inf$x
    mass_sample_info$y <- pixel_inf$y
    mass_sample_info$z <- pixel_inf$z
    cat("  ✓ Spatial coordinates added to sample_info\n")
  }
  
  # Add placeholder columns expected by mass_dataset
  mass_sample_info$subject_id <- NA
  mass_sample_info$injection_order <- seq_len(nrow(pixel_inf))
  mass_sample_info$class <- "pixel"
  mass_sample_info$group <- "default"
  
  cat("  Samples created:", nrow(mass_sample_info), "\n")
  cat("  ✓ sample_info created\n\n")
  
  # ---- Create mass_dataset variable_info ----
  cat("[4/6] Creating mass_dataset variable_info...\n")
  
  # Create variable_info for mass_dataset
  mass_variable_info <- data.frame(
    variable_id = variable_inf$variable_id,
    mz = variable_inf$mz,
    rt = 0,  # MSI data has no retention time
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat("  Variables created:", nrow(mass_variable_info), "\n")
  cat("  ✓ variable_info created\n\n")
  
  # ---- Create info_note data frames ----
  cat("[5/6] Creating info_note tables...\n")
  
  # sample_info_note
  mass_sample_info_note <- data.frame(
    name = colnames(mass_sample_info),
    meaning = c(
      "Sample ID (originally pixel_id from MSI data)",
      if (add_spatial_to_sample_info) c(
        "Horizontal spatial coordinate",
        "Vertical spatial coordinate",
        "Tissue slice depth"
      ) else NULL,
      "Subject ID",
      "Injection order",
      "Sample class",
      "Sample group"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # variable_info_note
  mass_variable_info_note <- data.frame(
    name = c("variable_id", "mz", "rt"),
    meaning = c(
      "Variable ID",
      "Mass to charge ratio",
      "Retention time (set to 0 for MSI data)"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat("  ✓ info_note tables created\n\n")
  
  # ---- Save MSI metadata ----
  if (save_metadata) {
    cat("[6/6] Saving MSI metadata...\n")
    
    # Create metadata object
    msi_metadata <- list(
      object_id = object_id,
      conversion_time = Sys.time(),
      pixel_info = pixel_inf,
      sample_info = sample_inf,
      pixel_info_note = pixel_inf_note,
      sample_info_note = sample_inf_note,
      original_variable_info_note = variable_inf_note,
      n_pixels = nrow(pixel_inf),
      n_variables = nrow(variable_inf),
      spatial_added = add_spatial_to_sample_info
    )
    
    # Save to environment
    if (!exists(".msi_conversion_cache", envir = .GlobalEnv)) {
      assign(".msi_conversion_cache", new.env(), envir = .GlobalEnv)
    }
    
    cache_env <- get(".msi_conversion_cache", envir = .GlobalEnv)
    cache_env[[object_id]] <- msi_metadata
    
    cat("  ✓ Metadata saved to environment (.msi_conversion_cache$", object_id, ")\n")
    
    # Optionally save to file
    if (!is.null(metadata_path)) {
      saveRDS(msi_metadata, metadata_path)
      cat("  ✓ Metadata saved to file:", metadata_path, "\n")
    }
    
    cat("\n")
  }
  
  # ---- Create mass_dataset object ----
  cat("Creating mass_dataset object...\n")
  
  # Use massdataset::create_mass_dataset if available
  if (requireNamespace("massdataset", quietly = TRUE)) {
    mass_obj <- massdataset::create_mass_dataset(
      expression_data = expr_data,
      sample_info = mass_sample_info,
      variable_info = mass_variable_info,
      sample_info_note = mass_sample_info_note,
      variable_info_note = mass_variable_info_note
    )
  } else {
    stop("massdataset package is required for this conversion.\n",
         "Please install it with: install.packages('massdataset')")
  }
  
  # Add conversion info to process_info
  conversion_parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "convert_msi2mass_dataset()",
    parameter = list(
      conversion_id = object_id,
      metadata_saved = save_metadata,
      spatial_added = add_spatial_to_sample_info
    ),
    time = Sys.time()
  )
  
  if (length(mass_obj@process_info) == 0) {
    mass_obj@process_info <- list(msi_conversion = conversion_parameter)
  } else {
    mass_obj@process_info$msi_conversion <- conversion_parameter
  }
  
  cat("\n========================================\n")
  cat(" ✓ Conversion complete!\n")
  cat("========================================\n\n")
  
  cat("Summary:\n")
  cat("  Input:  msi_dataset  (", n_variables(msi_object), "variables ×", 
      n_pixels(msi_object), "pixels)\n")
  cat("  Output: mass_dataset (", nrow(mass_variable_info), "variables ×", 
      nrow(mass_sample_info), "samples)\n")
  cat("  Conversion ID:", object_id, "\n")
  if (save_metadata) {
    cat("  Metadata saved: Yes\n")
    if (!is.null(metadata_path)) {
      cat("  Metadata file:", metadata_path, "\n")
    }
  }
  cat("\n")
  
  return(mass_obj)
}


#' Save MSI Metadata to File
#'
#' @description
#' Saves MSI-specific metadata from the conversion cache to an RDS file.
#'
#' @param object_id Conversion object ID
#' @param file_path Path to save the metadata RDS file
#'
#' @export
save_msi_metadata <- function(object_id, file_path = "msi_metadata.rds") {
  
  if (!exists(".msi_conversion_cache", envir = .GlobalEnv)) {
    stop("No conversion cache found in environment")
  }
  
  cache_env <- get(".msi_conversion_cache", envir = .GlobalEnv)
  
  if (!exists(object_id, envir = cache_env)) {
    stop("Object ID '", object_id, "' not found in conversion cache")
  }
  
  metadata <- cache_env[[object_id]]
  saveRDS(metadata, file_path)
  
  cat("✓ MSI metadata saved to:", file_path, "\n")
  cat("  Object ID:", object_id, "\n")
  cat("  Pixels:", metadata$n_pixels, "\n")
  cat("  Variables:", metadata$n_variables, "\n")
}


#' Load MSI Metadata from File
#'
#' @description
#' Loads MSI-specific metadata from an RDS file into the conversion cache.
#'
#' @param file_path Path to the metadata RDS file
#'
#' @return Invisibly returns the loaded metadata
#' @export
load_msi_metadata <- function(file_path) {
  
  if (!file.exists(file_path)) {
    stop("Metadata file not found: ", file_path)
  }
  
  metadata <- readRDS(file_path)
  
  # Validate metadata structure
  required_fields <- c("object_id", "pixel_info", "sample_info")
  missing_fields <- setdiff(required_fields, names(metadata))
  
  if (length(missing_fields) > 0) {
    stop("Invalid metadata file. Missing fields: ", 
         paste(missing_fields, collapse = ", "))
  }
  
  # Add to cache
  if (!exists(".msi_conversion_cache", envir = .GlobalEnv)) {
    assign(".msi_conversion_cache", new.env(), envir = .GlobalEnv)
  }
  
  cache_env <- get(".msi_conversion_cache", envir = .GlobalEnv)
  cache_env[[metadata$object_id]] <- metadata
  
  cat("✓ MSI metadata loaded from:", file_path, "\n")
  cat("  Object ID:", metadata$object_id, "\n")
  cat("  Pixels:", metadata$n_pixels, "\n")
  cat("  Variables:", metadata$n_variables, "\n")
  cat("  Added to conversion cache\n")
  
  invisible(metadata)
}


#' List Available MSI Metadata in Cache
#'
#' @description
#' Lists all available MSI metadata objects in the conversion cache.
#'
#' @return A data frame with conversion information
#' @export
list_msi_metadata <- function() {
  
  if (!exists(".msi_conversion_cache", envir = .GlobalEnv)) {
    cat("No conversion cache found\n")
    return(invisible(NULL))
  }
  
  cache_env <- get(".msi_conversion_cache", envir = .GlobalEnv)
  object_ids <- ls(cache_env)
  
  if (length(object_ids) == 0) {
    cat("Conversion cache is empty\n")
    return(invisible(NULL))
  }
  
  # Create summary table
  summary_list <- lapply(object_ids, function(id) {
    metadata <- cache_env[[id]]
    data.frame(
      object_id = id,
      conversion_time = as.character(metadata$conversion_time),
      n_pixels = metadata$n_pixels,
      n_variables = metadata$n_variables,
      spatial_added = metadata$spatial_added,
      stringsAsFactors = FALSE
    )
  })
  
  summary_df <- do.call(rbind, summary_list)
  
  cat("Available MSI metadata in cache:\n")
  print(summary_df)
  
  invisible(summary_df)
}