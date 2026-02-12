# ============================================================================
# Improved Conversion Functions for msidataset Package
# Author: Xiaoyu Xu
# Purpose: Enhanced conversion between msi_dataset and mass_dataset
# ============================================================================

#' Convert msi_dataset to mass_dataset (Enhanced Version)
#'
#' @description
#' Enhanced version with better metadata preservation and flexibility.
#'
#' @param msi_object An msi_dataset object
#' @param save_metadata Logical, whether to save MSI metadata (default: TRUE)
#' @param metadata_path Optional path to save metadata as RDS file
#' @param add_spatial_to_sample_info Logical, add x, y, z to sample_info (default: FALSE)
#' @param preserve_pixel_metadata Logical, preserve extra pixel_info columns (default: TRUE)
#' @param preserve_variable_metadata Logical, preserve extra variable_info columns (default: TRUE)
#' @param object_id Optional unique identifier (auto-generated if NULL)
#'
#' @return A mass_dataset object
#' @export
convert_msi2mass_dataset_v2 <- function(msi_object,
                                        save_metadata = TRUE,
                                        metadata_path = NULL,
                                        add_spatial_to_sample_info = FALSE,
                                        preserve_pixel_metadata = TRUE,
                                        preserve_variable_metadata = TRUE,
                                        object_id = NULL) {
  
  cat("\n========================================\n")
  cat("  MSI Dataset → Mass Dataset (v2)\n")
  cat("========================================\n\n")
  
  # ---- Validate input ----
  if (!inherits(msi_object, "msi_dataset")) {
    stop("Input must be an msi_dataset object")
  }
  
  # Generate unique ID
  if (is.null(object_id)) {
    object_id <- paste0("msi_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  
  cat("[1/7] Validating input...\n")
  cat("  Input class: msi_dataset\n")
  cat("  Variables:", n_variables(msi_object), "\n")
  cat("  Pixels:", n_pixels(msi_object), "\n")
  cat("  Conversion ID:", object_id, "\n\n")
  
  # ---- Extract data ----
  cat("[2/7] Extracting data from msi_dataset...\n")
  
  expr_data <- expression_data(msi_object)
  pixel_inf <- pixel_info(msi_object)
  sample_inf <- sample_info(msi_object)
  variable_inf <- variable_info(msi_object)
  
  pixel_inf_note <- pixel_info_note(msi_object)
  sample_inf_note <- sample_info_note(msi_object)
  variable_inf_note <- variable_info_note(msi_object)
  
  cat("  ✓ All data extracted\n\n")
  
  # ---- Create mass_dataset sample_info ----
  cat("[3/7] Creating mass_dataset sample_info...\n")
  
  # Base columns
  mass_sample_info <- data.frame(
    sample_id = pixel_inf$pixel_id,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Add spatial coordinates if requested
  if (add_spatial_to_sample_info) {
    mass_sample_info$x <- pixel_inf$x
    mass_sample_info$y <- pixel_inf$y
    mass_sample_info$z <- pixel_inf$z
    cat("  ✓ Spatial coordinates added\n")
  }
  
  # Preserve extra pixel_info columns
  if (preserve_pixel_metadata) {
    base_pixel_cols <- c("pixel_id", "x", "y", "z")
    extra_pixel_cols <- setdiff(colnames(pixel_inf), base_pixel_cols)
    
    if (length(extra_pixel_cols) > 0) {
      mass_sample_info <- cbind(
        mass_sample_info,
        pixel_inf[, extra_pixel_cols, drop = FALSE]
      )
      cat("  ✓ Preserved", length(extra_pixel_cols), 
          "extra pixel metadata columns:", 
          paste(extra_pixel_cols, collapse = ", "), "\n")
    }
  }
  
  # Add standard mass_dataset columns
  mass_sample_info$subject_id <- NA
  mass_sample_info$injection_order <- seq_len(nrow(pixel_inf))
  mass_sample_info$class <- "pixel"
  mass_sample_info$group <- "default"
  
  cat("  Samples created:", nrow(mass_sample_info), "\n")
  cat("  ✓ sample_info created\n\n")
  
  # ---- Create mass_dataset variable_info ----
  cat("[4/7] Creating mass_dataset variable_info...\n")
  
  # Base columns
  mass_variable_info <- data.frame(
    variable_id = variable_inf$variable_id,
    mz = variable_inf$mz,
    rt = 0,  # MSI data has no retention time
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Preserve extra variable_info columns (except ms2)
  if (preserve_variable_metadata) {
    base_var_cols <- c("variable_id", "mz", "ms2")
    extra_var_cols <- setdiff(colnames(variable_inf), base_var_cols)
    
    if (length(extra_var_cols) > 0) {
      mass_variable_info <- cbind(
        mass_variable_info,
        variable_inf[, extra_var_cols, drop = FALSE]
      )
      cat("  ✓ Preserved", length(extra_var_cols), 
          "extra variable metadata columns:", 
          paste(extra_var_cols, collapse = ", "), "\n")
    }
  }
  
  cat("  Variables created:", nrow(mass_variable_info), "\n")
  cat("  ✓ variable_info created\n\n")
  
  # ---- Create info_note tables ----
  cat("[5/7] Creating info_note tables...\n")
  
  # sample_info_note
  sample_note_meanings <- c(
    sample_id = "Sample ID (originally pixel_id from MSI data)"
  )
  
  if (add_spatial_to_sample_info) {
    sample_note_meanings <- c(
      sample_note_meanings,
      x = "Horizontal spatial coordinate",
      y = "Vertical spatial coordinate",
      z = "Tissue slice depth"
    )
  }
  
  if (preserve_pixel_metadata) {
    extra_pixel_cols <- setdiff(
      colnames(pixel_inf), 
      c("pixel_id", "x", "y", "z")
    )
    for (col in extra_pixel_cols) {
      # Try to get meaning from pixel_info_note
      meaning_idx <- which(pixel_inf_note$name == col)
      if (length(meaning_idx) > 0) {
        sample_note_meanings[col] <- pixel_inf_note$meaning[meaning_idx]
      } else {
        sample_note_meanings[col] <- paste("Pixel metadata:", col)
      }
    }
  }
  
  sample_note_meanings <- c(
    sample_note_meanings,
    subject_id = "Subject ID",
    injection_order = "Injection order",
    class = "Sample class",
    group = "Sample group"
  )
  
  mass_sample_info_note <- data.frame(
    name = names(sample_note_meanings),
    meaning = as.character(sample_note_meanings),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # variable_info_note
  var_note_meanings <- c(
    variable_id = "Variable ID",
    mz = "Mass to charge ratio",
    rt = "Retention time (set to 0 for MSI data)"
  )
  
  if (preserve_variable_metadata) {
    extra_var_cols <- setdiff(
      colnames(variable_inf),
      c("variable_id", "mz", "ms2")
    )
    for (col in extra_var_cols) {
      meaning_idx <- which(variable_inf_note$name == col)
      if (length(meaning_idx) > 0) {
        var_note_meanings[col] <- variable_inf_note$meaning[meaning_idx]
      } else {
        var_note_meanings[col] <- paste("Variable metadata:", col)
      }
    }
  }
  
  mass_variable_info_note <- data.frame(
    name = names(var_note_meanings),
    meaning = as.character(var_note_meanings),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat("  ✓ info_note tables created\n\n")
  
  # ---- Save MSI metadata ----
  cat("[6/7] Preparing MSI metadata...\n")
  
  msi_metadata <- list(
    object_id = object_id,
    conversion_time = Sys.time(),
    conversion_version = "v2",
    
    # Core MSI data
    pixel_info = pixel_inf,
    sample_info = sample_inf,
    
    # Notes
    pixel_info_note = pixel_inf_note,
    sample_info_note = sample_inf_note,
    original_variable_info_note = variable_inf_note,
    
    # Original variable_info (for ms2 recovery)
    original_variable_info = variable_inf,
    
    # Statistics
    n_pixels = nrow(pixel_inf),
    n_variables = nrow(variable_inf),
    
    # Conversion parameters
    spatial_added = add_spatial_to_sample_info,
    pixel_metadata_preserved = preserve_pixel_metadata,
    variable_metadata_preserved = preserve_variable_metadata,
    
    # Column tracking
    original_pixel_columns = colnames(pixel_inf),
    original_variable_columns = colnames(variable_inf),
    extra_pixel_columns = setdiff(
      colnames(pixel_inf),
      c("pixel_id", "x", "y", "z")
    ),
    extra_variable_columns = setdiff(
      colnames(variable_inf),
      c("variable_id", "mz", "ms2")
    )
  )
  
  if (save_metadata) {
    # Save to environment
    if (!exists(".msi_conversion_cache", envir = .GlobalEnv)) {
      assign(".msi_conversion_cache", new.env(), envir = .GlobalEnv)
    }
    
    cache_env <- get(".msi_conversion_cache", envir = .GlobalEnv)
    cache_env[[object_id]] <- msi_metadata
    
    cat("  ✓ Metadata saved to environment\n")
    
    # Optionally save to file
    if (!is.null(metadata_path)) {
      saveRDS(msi_metadata, metadata_path)
      cat("  ✓ Metadata saved to file:", metadata_path, "\n")
    }
  }
  
  cat("\n")
  
  # ---- Create mass_dataset object ----
  cat("[7/7] Creating mass_dataset object...\n")
  
  if (requireNamespace("massdataset", quietly = TRUE)) {
    mass_obj <- massdataset::create_mass_dataset(
      expression_data = expr_data,
      sample_info = mass_sample_info,
      variable_info = mass_variable_info,
      sample_info_note = mass_sample_info_note,
      variable_info_note = mass_variable_info_note
    )
  } else {
    stop("massdataset package is required.\n",
         "Install with: remotes::install_github('tidymass/massdataset')")
  }
  
  # Add conversion info to process_info
  conversion_parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "convert_msi2mass_dataset_v2()",
    parameter = list(
      conversion_id = object_id,
      version = "v2",
      metadata_saved = save_metadata,
      spatial_added = add_spatial_to_sample_info,
      pixel_metadata_preserved = preserve_pixel_metadata,
      variable_metadata_preserved = preserve_variable_metadata
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
  cat("  Version: v2 (enhanced)\n")
  
  if (preserve_pixel_metadata && length(msi_metadata$extra_pixel_columns) > 0) {
    cat("  Extra pixel metadata preserved:",
        length(msi_metadata$extra_pixel_columns), "columns\n")
  }
  
  if (preserve_variable_metadata && length(msi_metadata$extra_variable_columns) > 0) {
    cat("  Extra variable metadata preserved:",
        length(msi_metadata$extra_variable_columns), "columns\n")
  }
  
  cat("\n")
  
  return(mass_obj)
}


#' Convert mass_dataset back to msi_dataset (Enhanced Version)
#'
#' @description
#' Enhanced version with better annotation handling and metadata recovery.
#'
#' @param mass_object A mass_dataset object
#' @param metadata_source Source: "auto", "env", "file"
#' @param metadata_path Path to metadata RDS file
#' @param object_id Conversion object ID (auto-detected if NULL)
#' @param create_annotation_table Logical, create annotation_table (default: TRUE)
#' @param annotation_cols Specific columns to treat as annotations (NULL = auto-detect)
#' @param restore_original_columns Logical, restore original ms2 and other columns (default: TRUE)
#'
#' @return An msi_dataset object
#' @export
convert_mass2msi_dataset_v2 <- function(mass_object,
                                        metadata_source = c("auto", "env", "file"),
                                        metadata_path = NULL,
                                        object_id = NULL,
                                        create_annotation_table = TRUE,
                                        annotation_cols = NULL,
                                        restore_original_columns = TRUE) {
  
  cat("\n========================================\n")
  cat("  Mass Dataset → MSI Dataset (v2)\n")
  cat("========================================\n\n")
  
  # ---- Validate input ----
  if (!inherits(mass_object, "mass_dataset")) {
    stop("Input must be a mass_dataset object")
  }
  
  metadata_source <- match.arg(metadata_source)
  
  cat("[1/8] Validating input...\n")
  cat("  Input class: mass_dataset\n")
  cat("  Variables:", nrow(mass_object@variable_info), "\n")
  cat("  Samples:", nrow(mass_object@sample_info), "\n\n")
  
  # ---- Retrieve metadata ----
  cat("[2/8] Retrieving MSI metadata...\n")
  
  msi_metadata <- NULL
  
  # Try to find object_id from process_info
  if (is.null(object_id)) {
    if ("msi_conversion" %in% names(mass_object@process_info)) {
      conv_info <- mass_object@process_info$msi_conversion
      if (!is.null(conv_info@parameter$conversion_id)) {
        object_id <- conv_info@parameter$conversion_id
        cat("  ✓ Found conversion ID:", object_id, "\n")
      }
    }
  }
  
  # Retrieve metadata
  if (metadata_source == "auto") {
    # Try environment first
    if (exists(".msi_conversion_cache", envir = .GlobalEnv)) {
      cache_env <- get(".msi_conversion_cache", envir = .GlobalEnv)
      if (!is.null(object_id) && exists(object_id, envir = cache_env)) {
        msi_metadata <- cache_env[[object_id]]
        cat("  ✓ Metadata from environment cache\n")
      }
    }
    
    # Try file if not found
    if (is.null(msi_metadata) && !is.null(metadata_path)) {
      if (file.exists(metadata_path)) {
        msi_metadata <- readRDS(metadata_path)
        cat("  ✓ Metadata loaded from file\n")
      }
    }
    
  } else if (metadata_source == "env") {
    if (!exists(".msi_conversion_cache", envir = .GlobalEnv)) {
      stop("No conversion cache found in environment")
    }
    cache_env <- get(".msi_conversion_cache", envir = .GlobalEnv)
    
    if (is.null(object_id)) {
      available_ids <- ls(cache_env)
      stop("Please specify object_id. Available: ",
           paste(available_ids, collapse = ", "))
    }
    
    msi_metadata <- cache_env[[object_id]]
    cat("  ✓ Metadata from environment\n")
    
  } else if (metadata_source == "file") {
    if (is.null(metadata_path)) {
      stop("metadata_path required when metadata_source = 'file'")
    }
    msi_metadata <- readRDS(metadata_path)
    cat("  ✓ Metadata from file\n")
  }
  
  if (is.null(msi_metadata)) {
    stop("Failed to retrieve MSI metadata")
  }
  
  # Check metadata version
  metadata_version <- if (!is.null(msi_metadata$conversion_version)) {
    msi_metadata$conversion_version
  } else {
    "v1"
  }
  cat("  Metadata version:", metadata_version, "\n\n")
  
  # ---- Validate metadata ----
  cat("[3/8] Validating metadata...\n")
  
  if (nrow(mass_object@sample_info) != msi_metadata$n_pixels) {
    warning("Sample count mismatch (may be due to filtering)")
  }
  
  if (nrow(mass_object@variable_info) != msi_metadata$n_variables) {
    warning("Variable count mismatch (may be due to filtering)")
  }
  
  cat("  ✓ Metadata validated\n\n")
  
  # ---- Extract data ----
  cat("[4/8] Extracting data from mass_dataset...\n")
  
  expr_data <- mass_object@expression_data
  mass_sample_info <- mass_object@sample_info
  mass_variable_info <- mass_object@variable_info
  
  cat("  ✓ Data extracted\n\n")
  
  # ---- Reconstruct pixel_info ----
  cat("[5/8] Reconstructing pixel_info...\n")
  
  pixel_inf <- msi_metadata$pixel_info
  pixel_ids_in_mass <- mass_sample_info$sample_id
  
  # Filter to match current samples
  pixel_inf_filtered <- pixel_inf[pixel_inf$pixel_id %in% pixel_ids_in_mass, ]
  pixel_inf_filtered <- pixel_inf_filtered[
    match(pixel_ids_in_mass, pixel_inf_filtered$pixel_id),
  ]
  
  cat("  Original pixels:", nrow(msi_metadata$pixel_info), "\n")
  cat("  Current pixels:", nrow(pixel_inf_filtered), "\n")
  cat("  ✓ pixel_info reconstructed\n\n")
  
  # ---- Reconstruct sample_info ----
  cat("[6/8] Reconstructing sample_info...\n")
  
  sample_inf <- msi_metadata$sample_info
  cat("  ✓ sample_info reconstructed\n\n")
  
  # ---- Reconstruct variable_info ----
  cat("[7/8] Reconstructing variable_info...\n")
  
  # Determine annotation columns
  if (is.null(annotation_cols)) {
    # Auto-detect based on metadata
    base_cols <- c("variable_id", "mz", "rt")
    
    if (metadata_version == "v2" && 
        !is.null(msi_metadata$original_variable_columns)) {
      # Exclude original columns (except ms2)
      original_cols_no_ms2 <- setdiff(
        msi_metadata$original_variable_columns,
        "ms2"
      )
      annotation_cols <- setdiff(
        colnames(mass_variable_info),
        c(base_cols, original_cols_no_ms2)
      )
    } else {
      # Simple detection
      annotation_cols <- setdiff(colnames(mass_variable_info), base_cols)
    }
  }
  
  # Create base variable_info
  msi_variable_info <- data.frame(
    variable_id = mass_variable_info$variable_id,
    mz = mass_variable_info$mz,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Restore ms2 column
  if (restore_original_columns && 
      !is.null(msi_metadata$original_variable_info)) {
    
    original_var_info <- msi_metadata$original_variable_info
    
    if ("ms2" %in% colnames(original_var_info)) {
      match_idx <- match(
        msi_variable_info$variable_id,
        original_var_info$variable_id
      )
      msi_variable_info$ms2 <- original_var_info$ms2[match_idx]
      cat("  ✓ Restored ms2 column from original data\n")
    } else {
      msi_variable_info$ms2 <- "no"
    }
    
    # Restore other original columns if they exist
    if (metadata_version == "v2" && 
        !is.null(msi_metadata$extra_variable_columns)) {
      
      extra_cols <- intersect(
        msi_metadata$extra_variable_columns,
        colnames(mass_variable_info)
      )
      
      if (length(extra_cols) > 0) {
        msi_variable_info <- cbind(
          msi_variable_info,
          mass_variable_info[, extra_cols, drop = FALSE]
        )
        cat("  ✓ Restored", length(extra_cols), 
            "extra variable columns\n")
      }
    }
  } else {
    msi_variable_info$ms2 <- "no"
  }
  
  cat("  ✓ variable_info reconstructed\n\n")
  
  # ---- Create annotation_table ----
  cat("[8/8] Creating annotation_table...\n")
  
  if (create_annotation_table && length(annotation_cols) > 0) {
    annotation_table <- mass_variable_info[, 
                                           c("variable_id", annotation_cols),
                                           drop = FALSE]
    cat("  ✓ Annotation columns:", 
        paste(annotation_cols, collapse = ", "), "\n")
  } else {
    annotation_table <- data.frame(
      variable_id = character(),
      stringsAsFactors = FALSE
    )
    cat("  No annotation columns found\n")
  }
  
  cat("\n")
  
  # ---- Restore info_notes ----
  pixel_inf_note <- msi_metadata$pixel_info_note
  sample_inf_note <- msi_metadata$sample_info_note
  
  # Reconstruct variable_info_note
  var_note_names <- colnames(msi_variable_info)
  var_note_meanings <- rep("", length(var_note_names))
  
  # Try to get meanings from original note
  if (!is.null(msi_metadata$original_variable_info_note)) {
    orig_note <- msi_metadata$original_variable_info_note
    for (i in seq_along(var_note_names)) {
      name <- var_note_names[i]
      meaning_idx <- which(orig_note$name == name)
      if (length(meaning_idx) > 0) {
        var_note_meanings[i] <- orig_note$meaning[meaning_idx]
      } else {
        var_note_meanings[i] <- name
      }
    }
  } else {
    var_note_meanings <- var_note_names
  }
  
  variable_inf_note <- data.frame(
    name = var_note_names,
    meaning = var_note_meanings,
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
    function_name = "convert_mass2msi_dataset_v2()",
    parameter = list(
      conversion_id = msi_metadata$object_id,
      version = "v2",
      metadata_source = metadata_source,
      annotation_cols = if(create_annotation_table) annotation_cols else NULL
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
  cat("  Version: v2 (enhanced)\n")
  
  if (nrow(annotation_table) > 0 && ncol(annotation_table) > 1) {
    cat("  Annotation columns:", ncol(annotation_table) - 1, "\n")
  }
  
  cat("\n")
  
  return(msi_obj)
}