#' @title Convert Cardinal Object to msi_dataset
#' @description 
#' Convert Cardinal object to msi_dataset with new structure:
#' - sample_info: sample_id, x, y, z, subject_id, section_id, resolution, matrix
#' - variable_info: variable_id, mz, ms2
#'
#' @param cardinal_object Cardinal object
#' @param sample_name Sample name (character)
#' @param section_name Section/slice name (character)
#' @param resolution Spatial resolution in micrometers per pixel (numeric, optional)
#' @param matrix_name MALDI matrix name: CHCA, DHB, NEDC, noMatrix (character, optional)
#' @export

.extract_cardinal_mz_numeric <- function(cardinal_object) {
  raw_mz <- Cardinal::mz(cardinal_object)

  extract_numeric <- function(x) {
    if (is.null(x)) return(numeric(0))
    if (is.numeric(x) || is.integer(x)) return(as.numeric(x))
    if (is.character(x)) {
      out <- suppressWarnings(as.numeric(x))
      return(out[is.finite(out)])
    }
    if (is.matrix(x) || is.data.frame(x)) {
      out <- suppressWarnings(as.numeric(as.matrix(x)))
      return(out[is.finite(out)])
    }
    if (isS4(x)) {
      for (nm in c("mz", "mzs", "mass", "values", "data")) {
        if (nm %in% methods::slotNames(x)) {
          out <- extract_numeric(methods::slot(x, nm))
          if (length(out) > 0) return(out)
        }
      }
      as_list <- try(as.list(x), silent = TRUE)
      if (!inherits(as_list, "try-error")) {
        out <- extract_numeric(as_list)
        if (length(out) > 0) return(out)
      }
      return(numeric(0))
    }
    if (is.list(x)) {
      if (!is.null(x$mz)) {
        out <- extract_numeric(x$mz)
        if (length(out) > 0) return(out)
      }
      for (i in seq_along(x)) {
        out <- extract_numeric(x[[i]])
        if (length(out) > 0) return(out)
      }
      pooled <- suppressWarnings(as.numeric(unlist(x, recursive = TRUE, use.names = FALSE)))
      return(pooled[is.finite(pooled)])
    }
    out <- suppressWarnings(as.numeric(as.vector(x)))
    out[is.finite(out)]
  }

  extract_numeric(raw_mz)
}

.coerce_mz_column_numeric <- function(x) {
  if (is.null(x)) return(numeric(0))
  if (is.numeric(x) || is.integer(x)) return(as.numeric(x))
  if (is.character(x)) return(suppressWarnings(as.numeric(x)))
  if (is.matrix(x) || is.data.frame(x)) {
    return(suppressWarnings(as.numeric(as.matrix(x))))
  }
  if (isS4(x)) {
    as_list <- try(as.list(x), silent = TRUE)
    if (!inherits(as_list, "try-error")) {
      return(.coerce_mz_column_numeric(as_list))
    }
    return(numeric(0))
  }
  if (is.list(x)) {
    return(vapply(x, function(el) {
      v <- suppressWarnings(as.numeric(el))
      v <- v[is.finite(v)]
      if (length(v) == 0) NA_real_ else v[[1]]
    }, numeric(1)))
  }
  suppressWarnings(as.numeric(x))
}

convert_cardinal2msi_dataset <- function(cardinal_object,
                                         sample_name = "subject_1",
                                         section_name = "section_1",
                                         resolution = NA,
                                         matrix_name = NA) {
  
  # 1. Check Cardinal package
  if (!requireNamespace("Cardinal", quietly = TRUE)) {
    stop("please install Cardinal package: BiocManager::install('Cardinal')")
  }
  
  # 2. Extract Expression Data (intensity matrix)
  cat(crayon::blue("Extracting intensity matrix...\n"))
  expression_data <- as.matrix(Cardinal::spectra(cardinal_object))
  
  # 3. Extract and Prepare Sample Info
  cat(crayon::blue("Preparing sample information...\n"))
  sample_info <- as.data.frame(Cardinal::pData(cardinal_object))
  
  # Extract pixel IDs (x, y coordinates)
  pixel_ids <- rownames(Cardinal::pData(cardinal_object))
  
  # Create sample_id column (pixel-level unique id)
  if (is.null(pixel_ids) || all(pixel_ids == seq_len(nrow(sample_info)))) {
    sample_info$sample_id <- paste0("pixel_", seq_len(nrow(sample_info)))
  } else {
    sample_info$sample_id <- pixel_ids
  }
  
  # Ensure x and y columns exist
  if (!("x" %in% colnames(sample_info)) || !("y" %in% colnames(sample_info))) {
    # Extract coordinates from Cardinal object
    coords <- Cardinal::coord(cardinal_object)
    if (ncol(coords) >= 2) {
      sample_info$x <- coords[, 1]
      sample_info$y <- coords[, 2]
    } else {
      stop("Cannot extract x, y coordinates from Cardinal object")
    }
  }
  
  # Add z coordinate if 3D data
  if (!("z" %in% colnames(sample_info))) {
    coords <- Cardinal::coord(cardinal_object)
    if (ncol(coords) >= 3) {
      sample_info$z <- coords[, 3]
    } else {
      # For 2D data, set z to 1
      sample_info$z <- 1
    }
  }
  
  # Add metadata columns
  sample_info$subject_id <- sample_name
  sample_info$section_id <- section_name
  
  if (!is.na(resolution)) {
    sample_info$resolution <- resolution
  }
  
  if (!is.na(matrix_name)) {
    sample_info$matrix <- matrix_name
  }
  
  # Reorder columns to match expected structure
  expected_cols <- c("sample_id", "x", "y", "z", "subject_id", "section_id",
                     "resolution", "matrix")
  existing_cols <- intersect(expected_cols, colnames(sample_info))
  
  # Reorder to have core columns first
  sample_info <- sample_info[, existing_cols, drop = FALSE]
  
  # 4. Extract and Prepare Variable Info
  cat(crayon::blue("Preparing variable (m/z) information...\n"))
  variable_info <- tryCatch(
    {
      as.data.frame(Cardinal::fData(cardinal_object))
    },
    error = function(e) {
      mz_vals <- .extract_cardinal_mz_numeric(cardinal_object)
      if (length(mz_vals) == 0) {
        stop("Cannot extract feature metadata: both fData() and mz() are unavailable/empty.")
      }
      data.frame(mz = mz_vals, stringsAsFactors = FALSE)
    }
  )
  
  # Ensure mz column exists
  if (!"mz" %in% colnames(variable_info)) {
    mz_vals <- .extract_cardinal_mz_numeric(cardinal_object)
    if (length(mz_vals) == 0) {
      stop("Cannot extract mz values from Cardinal object.")
    }
    variable_info$mz <- mz_vals
  }

  # Ensure mz column is pure numeric (handles S4/list-like Cardinal metadata)
  mz_col <- .coerce_mz_column_numeric(variable_info$mz)
  if (length(mz_col) == nrow(variable_info) && sum(is.finite(mz_col)) > 0) {
    variable_info$mz <- mz_col
  } else {
    mz_vals <- .extract_cardinal_mz_numeric(cardinal_object)
    if (length(mz_vals) == nrow(variable_info)) {
      variable_info$mz <- mz_vals
    }
  }
  
  # Keep variable_info rows aligned with expression_data rows
  if (nrow(variable_info) != nrow(expression_data)) {
    mz_vals <- .extract_cardinal_mz_numeric(cardinal_object)
    if (length(mz_vals) == nrow(expression_data)) {
      variable_info <- data.frame(mz = mz_vals, stringsAsFactors = FALSE)
    } else if (nrow(variable_info) >= nrow(expression_data)) {
      variable_info <- variable_info[seq_len(nrow(expression_data)), , drop = FALSE]
    } else {
      stop(
        "Feature dimension mismatch: nrow(expression_data) = ", nrow(expression_data),
        ", nrow(variable_info) = ", nrow(variable_info), "."
      )
    }
  }

  # Final mz sanity check before ID generation
  if (sum(is.finite(variable_info$mz)) < nrow(variable_info)) {
    mz_vals <- .extract_cardinal_mz_numeric(cardinal_object)
    if (length(mz_vals) == nrow(variable_info)) {
      variable_info$mz <- mz_vals
    } else {
      stop("Failed to normalize mz column to numeric values for all variables.")
    }
  }
  
  # Check for generate_variable_id function
  if (!exists("generate_variable_id")) {
    stop("Function generate_variable_id not found. Please ensure source('generate_variable_id.R') has been run.")
  }
  
  # Generate variable_id
  variable_info$variable_id <- generate_variable_id(variable_info$mz)
  
  # Add ms2 column (yes/no indicator)
  if (!("ms2" %in% colnames(variable_info))) {
    # Default: set all to "no" unless you have specific information
    variable_info$ms2 <- "no"
  } else {
    # Convert to yes/no if it's boolean
    if (is.logical(variable_info$ms2)) {
      variable_info$ms2 <- ifelse(variable_info$ms2, "Yes", "no")
    } else if (is.numeric(variable_info$ms2)) {
      variable_info$ms2 <- ifelse(variable_info$ms2 > 0, "Yes", "no")
    }
    # Otherwise assume it's already yes/no character
  }
  
  # Keep only required columns: variable_id, mz, ms2
  variable_info <- variable_info[, c("variable_id", "mz", "ms2"), drop = FALSE]
  
  # 5. Set Expression Data Dimensions
  cat(crayon::blue("Setting expression data dimensions...\n"))
  rownames(expression_data) <- as.character(variable_info$variable_id)
  colnames(expression_data) <- as.character(sample_info$sample_id)
  
  # 6. Create msi_dataset Object
  cat(crayon::blue("Creating msi_dataset object...\n"))
  object <- create_msi_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )
  
  # 7. Log Processing Activity
  object <- add_process_log(
    object = object,
    function_name = "convert_cardinal2msi_dataset()",
    parameter = list(
      sample_name = sample_name,
      section_name = section_name,
      cardinal_class = class(cardinal_object)[1],
      num_pixels = ncol(expression_data),
      num_variables = nrow(expression_data),
      resolution = resolution,
      matrix_name = matrix_name
    )
  )
  
  # 8. Summary Output
  cat(crayon::green("Successfully converted!\n\n"))
  cat(crayon::cyan("Summary:\n"))
  cat(sprintf("  Pixels (samples): %d\n", ncol(expression_data)))
  cat(sprintf("  Variables (m/z features): %d\n", nrow(expression_data)))
  cat(sprintf("  Sample name: %s\n", sample_name))
  cat(sprintf("  Section name: %s\n", section_name))
  
  if (!is.na(resolution)) {
    cat(sprintf("  Resolution: %.2f µm/pixel\n", resolution))
  }
  if (!is.na(matrix_name)) {
    cat(sprintf("  MALDI matrix: %s\n", matrix_name))
  }
  cat("\n")
  
  return(object)
}
