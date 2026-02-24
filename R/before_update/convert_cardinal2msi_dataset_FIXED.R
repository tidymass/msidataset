# ---------------------------------------------------------------------------
# Internal: m/z bin-centre calculation
#   Input : sorted numeric vector of raw m/z values
#   Output: sorted numeric vector of bin centres (one per bin)
# ---------------------------------------------------------------------------
#' @keywords internal
.bin_mz_values <- function(mz_sorted, tolerance, units = "ppm") {
  if (length(mz_sorted) == 0L) return(numeric(0))
  
  centres <- numeric(length(mz_sorted))
  bin_sum  <- mz_sorted[1]
  bin_n    <- 1L
  ref      <- mz_sorted[1]
  out_idx  <- 1L
  
  for (i in seq_along(mz_sorted)[-1L]) {
    delta <- if (units == "ppm") {
      abs(mz_sorted[i] - ref) / ref * 1e6
    } else {
      abs(mz_sorted[i] - ref)
    }
    
    if (delta > tolerance) {
      # close current bin
      centres[out_idx] <- bin_sum / bin_n
      out_idx <- out_idx + 1L
      # start new bin
      bin_sum <- mz_sorted[i]
      bin_n   <- 1L
      ref     <- mz_sorted[i]
    } else {
      bin_sum <- bin_sum + mz_sorted[i]
      bin_n   <- bin_n   + 1L
    }
  }
  # close last bin
  centres[out_idx] <- bin_sum / bin_n
  
  sort(centres[seq_len(out_idx)])
}


# ---------------------------------------------------------------------------
# Internal: assign each raw m/z to the nearest bin centre
#   Returns an integer index into `bin_centres` for each element of `mz_vec`,
#   or NA when the nearest bin is farther than `tolerance`.
# ---------------------------------------------------------------------------
#' @keywords internal
.match_to_bins <- function(mz_vec, bin_centres, tolerance, units = "ppm") {
  if (length(mz_vec) == 0L) return(integer(0))
  
  # findInterval gives the left boundary index; check both neighbours
  lo  <- findInterval(mz_vec, bin_centres, all.inside = FALSE)
  lo  <- pmax(lo, 1L)
  hi  <- pmin(lo + 1L, length(bin_centres))
  
  dist_lo <- if (units == "ppm") {
    abs(bin_centres[lo] - mz_vec) / mz_vec * 1e6
  } else {
    abs(bin_centres[lo] - mz_vec)
  }
  dist_hi <- if (units == "ppm") {
    abs(bin_centres[hi] - mz_vec) / mz_vec * 1e6
  } else {
    abs(bin_centres[hi] - mz_vec)
  }
  
  idx      <- ifelse(dist_lo <= dist_hi, lo, hi)
  best_d   <- pmin(dist_lo, dist_hi)
  idx[best_d > tolerance] <- NA_integer_
  idx
}


# ---------------------------------------------------------------------------
# Internal: coerce exotic spectra objects to a plain base matrix
# ---------------------------------------------------------------------------
#' @keywords internal
.coerce_to_matrix <- function(x) {
  if (methods::is(x, "matter_mat")   ||
      methods::is(x, "matter_matc")  ||
      methods::is(x, "matter_arr")   ||
      methods::is(x, "sparseMatrix") ||
      methods::is(x, "dgCMatrix")    ||
      methods::is(x, "dgTMatrix"))   return(as.matrix(x))
  
  if (methods::is(x, "SimpleList") ||
      (is.list(x) && !is.data.frame(x)))  return(as.matrix(x[[1]]))
  
  if (!is.matrix(x) && !is.data.frame(x)) return(as.matrix(x))
  x
}


# ---------------------------------------------------------------------------
# Internal: extract a plain numeric mz vector from exotic S4 / list types
# 修复版本 v2：添加 matter_list 支持和改进的列表处理
# ---------------------------------------------------------------------------
#' @keywords internal
.extract_mz_vector <- function(mz_raw) {
  # 处理 matter_list（Cardinal/matter 包中的特殊列表类型）
  if (methods::is(mz_raw, "matter_list")) {
    tryCatch({
      # matter_list 可以转换为标准 list
      mz_std_list <- as.list(mz_raw)
      return(sort(unique(as.numeric(unlist(
        lapply(mz_std_list, as.numeric)
      )))))
    }, error = function(e) {
      # 备用方案：直接尝试转换
      tryCatch({
        return(sort(unique(as.numeric(unlist(mz_raw)))))
      }, error = function(e2) {
        stop("Failed to convert matter_list to numeric. ",
             "Original error: ", conditionMessage(e))
      })
    })
  }
  
  # 处理特殊的 S4 列表类型 (NumericList, SimpleList)
  if (methods::is(mz_raw, "NumericList"))  return(as.numeric(mz_raw[[1]]))
  
  if (methods::is(mz_raw, "SimpleList")) {
    # 如果是 SimpleList，尝试提取第一个元素
    tryCatch({
      return(sort(unique(as.numeric(mz_raw[[1]]))))
    }, error = function(e) {
      # 如果失败，尝试展平整个列表
      return(sort(unique(as.numeric(unlist(mz_raw)))))
    })
  }
  
  # 处理标准 R list（非 DataFrame）
  if (is.list(mz_raw) && !is.data.frame(mz_raw)) {
    # 使用 lapply + unlist 的组合处理嵌套列表
    return(sort(unique(as.numeric(unlist(
      lapply(mz_raw, as.numeric)
    )))))
  }
  
  # 处理 Rle 对象
  if (methods::is(mz_raw, "Rle")) return(as.numeric(as.vector(mz_raw)))
  
  # 处理一般 S4 对象
  if (isS4(mz_raw)) {
    if (tryCatch(methods::.hasSlot(mz_raw, ".Data"), error = function(e) FALSE)) {
      return(as.numeric(slot(mz_raw, ".Data")))
    }
    return(as.numeric(as.vector(mz_raw)))
  }
  
  # 最后的后备方案
  as.numeric(as.vector(mz_raw))
}


# ===========================================================================
# convert_cardinal2msi_dataset()
# ===========================================================================

#' Convert a Cardinal Imaging Object to msi_dataset
#'
#' @description
#' Converts \code{MSImagingArrays} (Cardinal 3.x), \code{MSImagingExperiment}
#' (Cardinal 2.x), or \code{MSImageSet} (Cardinal < 2) into an
#' \code{msi_dataset}.
#'
#' \strong{MSImagingArrays}: each pixel stores its own variable-length m/z
#' vector.  This function collects all unique m/z values across pixels, bins
#' them using \code{tolerance} + \code{units}, and assembles a dense
#' intensity matrix using fully vectorised matching (no nested R loops).
#'
#' @param cardinal_object A Cardinal imaging object
#'   (\code{MSImagingArrays}, \code{MSImagingExperiment}, or
#'   \code{MSImageSet}).
#' @param sample_name  Character.  Biological sample identifier.
#'   Default \code{"sample_1"}.
#' @param section_name Character.  Tissue section identifier.
#'   Default \code{"section_1"}.
#' @param resolution   Spatial resolution string (e.g. \code{"50um"}) or
#'   \code{NA}.
#' @param matrix_name  MALDI matrix name (e.g. \code{"CHCA"}) or \code{NA}.
#' @param tolerance    Numeric.  m/z binning tolerance applied when assembling
#'   the common m/z axis from \code{MSImagingArrays} per-pixel spectra.
#'   Typical values: \code{10} (tight), \code{50} (loose), \code{100} (very
#'   loose).  Set to \code{NA} to disable binning and use exact matching
#'   (only recommended when spectra are already on a common grid).
#'   Default \code{10}.
#' @param units        Tolerance units: \code{"ppm"} (default) or
#'   \code{"mz"} (absolute Da).
#'
#' @return An \code{msi_dataset} object.
#' @export
#'
#' @examples
#' \dontrun{
#' # tight binning (10 ppm)
#' msi_10 <- convert_cardinal2msi_dataset(obj, tolerance = 10,  units = "ppm")
#'
#' # loose binning (50 ppm)
#' msi_50 <- convert_cardinal2msi_dataset(obj, tolerance = 50,  units = "ppm")
#'
#' # absolute Da window (0.02 Da)
#' msi_da <- convert_cardinal2msi_dataset(obj, tolerance = 0.02, units = "mz")
#'
#' # no binning (exact m/z match)
#' msi_ex <- convert_cardinal2msi_dataset(obj, tolerance = NA)
#' }

convert_cardinal2msi_dataset <- function(cardinal_object,
                                         sample_name  = "sample_1",
                                         section_name = "section_1",
                                         resolution   = NA,
                                         matrix_name  = NA,
                                         tolerance    = 10,
                                         units        = "ppm") {
  
  if (!requireNamespace("Cardinal", quietly = TRUE))
    stop("Package 'Cardinal' is required but not installed.")
  
  # ---- validate tolerance --------------------------------------------------
  use_binning <- !is.na(tolerance) && is.numeric(tolerance) && tolerance > 0
  units       <- match.arg(units, c("ppm", "mz"))
  
  if (use_binning) {
    cat(sprintf("Tolerance: %g %s\n", tolerance, units))
  } else {
    cat("Tolerance: NA (exact match)\n")
  }
  
  # ---- type check ----------------------------------------------------------
  is_msimageset   <- inherits(cardinal_object, "MSImageSet")
  is_msimagingexp <- inherits(cardinal_object, "MSImagingExperiment")
  is_msimarrays   <- inherits(cardinal_object, "MSImagingArrays")
  
  if (!any(is_msimageset, is_msimagingexp, is_msimarrays))
    stop("cardinal_object must be MSImageSet, MSImagingExperiment, or ",
         "MSImagingArrays. Got: ", paste(class(cardinal_object), collapse = ", "))
  
  cat("Converting", class(cardinal_object)[1], "...\n")
  
  # ---- extract expression data & metadata ---------------------------------
  
  if (is_msimarrays) {
    # Cardinal 3.x: MSImagingArrays
    cat("Expression data:", nrow(cardinal_object), "x", ncol(cardinal_object), "\n")
    
    expr_raw    <- cardinal_object@.Data
    pixel_names <- colnames(cardinal_object)
    
    # for MSImagingArrays, coord slots store the spatial positions
    pixelData   <- S4Vectors::as.data.frame(Cardinal::pixelData(cardinal_object))
    
  } else if (is_msimagingexp) {
    # Cardinal 2.x: MSImagingExperiment
    cat("Expression data:", nrow(cardinal_object), "x", ncol(cardinal_object), "\n")
    
    expr_raw    <- SummarizedExperiment::assay(cardinal_object)
    pixel_names <- colnames(cardinal_object)
    pixelData   <- S4Vectors::as.data.frame(SummarizedExperiment::colData(cardinal_object))
    
  } else {
    # Cardinal < 2: MSImageSet (older Cardinal API)
    cat("Expression data:", nrow(cardinal_object), "x", ncol(cardinal_object), "\n")
    
    expr_raw    <- Cardinal::spectra(cardinal_object)
    pixel_names <- colnames(cardinal_object)
    pixelData   <- S4Vectors::as.data.frame(Cardinal::pixelData(cardinal_object))
  }
  
  # ---- extract m/z values -------------------------------------------------
  
  cat("Extracting m/z values...\n")
  
  mz_raw <- Cardinal::mz(cardinal_object)
  mz_extracted <- tryCatch({
    .extract_mz_vector(mz_raw)
  }, error = function(e) {
    cat("Error during .extract_mz_vector():\n")
    cat("  Class of mz_raw:", class(mz_raw), "\n")
    cat("  Error message:", conditionMessage(e), "\n")
    stop("Failed to extract m/z vector. Check that your Cardinal object has valid m/z data.")
  })
  
  if (length(mz_extracted) == 0) {
    stop("No m/z values extracted from cardinal_object. ",
         "Check that the object contains valid spectra.")
  }
  
  mz_values <- sort(unique(mz_extracted))
  
  cat("  Found", length(mz_values), "unique m/z values\n")
  
  # ---- bin m/z values (if requested) --------------------------------------
  
  if (use_binning) {
    mz_values <- .bin_mz_values(mz_values, tolerance = tolerance, units = units)
    cat("  Binned to", length(mz_values), "m/z values\n")
  }
  
  # ---- extract spatial coordinates ----------------------------------------
  
  cat("Extracting spatial coordinates...\n")
  
  if ("x" %in% colnames(pixelData) && "y" %in% colnames(pixelData)) {
    coord_x <- pixelData[["x"]]
    coord_y <- pixelData[["y"]]
  } else if ("coord" %in% colnames(pixelData)) {
    # PositionDataFrame may pack coordinates into 'coord' column
    coord_packed <- pixelData[["coord"]]
    if (is.matrix(coord_packed) && ncol(coord_packed) >= 2) {
      coord_x <- coord_packed[, 1]
      coord_y <- coord_packed[, 2]
    } else {
      stop("Cannot extract x/y from 'coord' column in pixelData.")
    }
  } else {
    stop("pixelData must contain 'x' and 'y' columns (or 'coord' matrix).")
  }
  
  # ---- coerce expression to matrix ----------------------------------------
  
  cat("Coercing expression data to matrix...\n")
  expr_matrix <- .coerce_to_matrix(expr_raw)
  
  if (nrow(expr_matrix) == 0 || ncol(expr_matrix) == 0) {
    stop("Expression matrix is empty after coercion.")
  }
  
  # ---- handle MSImagingArrays per-pixel m/z matching ----------------------
  
  if (is_msimarrays) {
    cat("Matching per-pixel spectra to common m/z axis...\n")
    
    # expr_raw is a list of per-pixel matrices
    # Assemble into a shared dense matrix using vectorised matching
    npixels   <- length(expr_raw)
    nmz_out   <- length(mz_values)
    
    # Validate dimensions
    if (npixels == 0) {
      stop("expr_raw is empty - no pixels found")
    }
    if (nmz_out == 0) {
      stop("mz_values is empty - no m/z values found")
    }
    
    # Pre-allocate output matrix
    dense_mat <- matrix(0, nrow = nmz_out, ncol = npixels,
                        dimnames = list(NULL, pixel_names))
    
    # Track statistics
    pixels_processed <- 0
    pixels_empty <- 0
    pixels_failed <- 0
    
    for (pix in seq_len(npixels)) {
      spec  <- expr_raw[[pix]]
      
      # Check 1: spec is not NULL
      if (is.null(spec)) {
        pixels_empty <- pixels_empty + 1
        next
      }
      
      # Check 2: spec is not empty
      if (length(spec) == 0) {
        pixels_empty <- pixels_empty + 1
        next
      }
      
      # Extract m/z and intensity from this pixel's spectrum
      spec_mz <- tryCatch({
        .extract_mz_vector(Cardinal::mz(spec))
      }, error = function(e) {
        # If Cardinal::mz() fails on the per-pixel spectrum,
        # assume its row indices are m/z bin indices
        seq_len(nrow(spec))
      })
      
      # Check 3: spec_mz is not NULL/empty
      if (is.null(spec_mz) || length(spec_mz) == 0) {
        pixels_failed <- pixels_failed + 1
        next
      }
      
      # Check 4: Safely convert intensity data
      spec_int <- tryCatch({
        result <- as.numeric(spec)
        
        # If result is NULL/empty, try flattening
        if (is.null(result) || length(result) == 0) {
          result <- as.numeric(unlist(spec))
        }
        
        result
      }, error = function(e) {
        NULL
      })
      
      # Check 5: spec_int is not NULL/empty
      if (is.null(spec_int) || length(spec_int) == 0) {
        pixels_failed <- pixels_failed + 1
        next
      }
      
      # Check 6: Length match verification
      if (length(spec_int) != length(spec_mz)) {
        min_len <- min(length(spec_int), length(spec_mz))
        spec_int <- spec_int[seq_len(min_len)]
        spec_mz <- spec_mz[seq_len(min_len)]
        
        if (length(spec_int) == 0 || length(spec_mz) == 0) {
          pixels_failed <- pixels_failed + 1
          next
        }
      }
      
      # Vectorised matching
      idx <- tryCatch({
        .match_to_bins(spec_mz, mz_values, tolerance, units)
      }, error = function(e) {
        NULL
      })
      
      # Check 7: Matching succeeded
      if (is.null(idx) || length(idx) == 0) {
        pixels_failed <- pixels_failed + 1
        next
      }
      
      # Assign intensities to output matrix
      valid <- !is.na(idx)
      
      # Check 8: Valid matches exist
      if (sum(valid) > 0) {
        # Check 9: Boundary check on indices
        valid_idx <- idx[valid]
        valid_rows <- valid_idx > 0 & valid_idx <= nmz_out
        
        if (sum(valid_rows) > 0) {
          dense_mat[valid_idx[valid_rows], pix] <- spec_int[valid][valid_rows]
          pixels_processed <- pixels_processed + 1
        } else {
          pixels_failed <- pixels_failed + 1
        }
      } else {
        pixels_failed <- pixels_failed + 1
      }
    }
    
    # Report processing results
    cat("  Processed pixels:", pixels_processed, "/", npixels, "\n")
    if (pixels_empty > 0)
      cat("  Empty pixels:", pixels_empty, "\n")
    if (pixels_failed > 0)
      cat("  Failed pixels:", pixels_failed, "\n")
    
    if (pixels_processed == 0) {
      stop("No pixels were successfully processed. ",
           "Check that your Cardinal object has valid spectra data.")
    }
    
    expr_matrix <- dense_mat
  }
  
  # ---- transpose to samples (columns) × variables (rows) ------------------
  
  if (nrow(expr_matrix) != length(mz_values)) {
    expr_matrix <- t(expr_matrix)
  }
  
  stopifnot(nrow(expr_matrix) == length(mz_values),
            ncol(expr_matrix) == length(coord_x))
  
  # ---- build sample_info dataframe ----------------------------------------
  
  sample_info_df <- data.frame(
    sample_id = paste0(sample_name, "_", seq_along(coord_x)),
    x         = as.numeric(coord_x),
    y         = as.numeric(coord_y),
    z         = if ("z" %in% colnames(pixelData))
                  as.numeric(pixelData[["z"]]) else rep(1L, length(coord_x)),
    section_id = section_name,
    stringsAsFactors = FALSE
  )
  
  if (!is.na(resolution))
    sample_info_df$resolution <- resolution
  if (!is.na(matrix_name))
    sample_info_df$matrix <- matrix_name
  if ("run" %in% colnames(pixelData))
    sample_info_df$run <- pixelData[["run"]]
  
  # ---- build variable_info dataframe (mz values) --------------------------
  
  variable_info_df <- data.frame(
    variable_id = paste0("mz_", sprintf("%.4f", mz_values)),
    mz          = mz_values,
    ms2         = FALSE,
    stringsAsFactors = FALSE
  )
  
  # ==========================================================================
  # DIMENSION GUARD
  # ==========================================================================
  
  if (nrow(expr_matrix) != nrow(variable_info_df)) {
    min_v <- min(nrow(expr_matrix), nrow(variable_info_df))
    warning("Row mismatch (", nrow(expr_matrix), " vs ",
            nrow(variable_info_df), ") — truncating to ", min_v)
    expr_matrix      <- expr_matrix[seq_len(min_v), , drop = FALSE]
    variable_info_df <- variable_info_df[seq_len(min_v), , drop = FALSE]
  }
  
  if (ncol(expr_matrix) != nrow(sample_info_df)) {
    min_s <- min(ncol(expr_matrix), nrow(sample_info_df))
    warning("Col mismatch (", ncol(expr_matrix), " vs ",
            nrow(sample_info_df), ") — truncating to ", min_s)
    expr_matrix    <- expr_matrix[, seq_len(min_s), drop = FALSE]
    sample_info_df <- sample_info_df[seq_len(min_s), , drop = FALSE]
  }
  
  colnames(expr_matrix) <- sample_info_df$sample_id
  rownames(expr_matrix) <- variable_info_df$variable_id
  
  # ==========================================================================
  # NOTE TABLES
  # ==========================================================================
  
  sinfo_note <- data.frame(
    name    = c("sample_id", "x", "y", "z", "section_id"),
    meaning = c("Pixel ID",
                "Horizontal spatial coordinate",
                "Vertical spatial coordinate",
                "Tissue slice depth",
                "Tissue section identifier"),
    stringsAsFactors = FALSE
  )
  if (!is.na(resolution))
    sinfo_note <- rbind(sinfo_note,
                        data.frame(name = "resolution",
                                   meaning = paste0("Spatial resolution (", resolution, ")"),
                                   stringsAsFactors = FALSE))
  if (!is.na(matrix_name))
    sinfo_note <- rbind(sinfo_note,
                        data.frame(name = "matrix",
                                   meaning = "Chemical matrix for ionization",
                                   stringsAsFactors = FALSE))
  if ("run" %in% colnames(sample_info_df))
    sinfo_note <- rbind(sinfo_note,
                        data.frame(name = "run",
                                   meaning = "Acquisition run identifier",
                                   stringsAsFactors = FALSE))
  
  vinfo_note <- data.frame(
    name    = c("variable_id", "mz", "ms2"),
    meaning = c("Variable ID",
                "Mass-to-charge ratio",
                "MS2 data available (TRUE/FALSE)"),
    stringsAsFactors = FALSE
  )
  
  # ==========================================================================
  # BUILD msi_dataset
  # ==========================================================================
  
  cat("Building msi_dataset...\n")
  msi_obj <- create_msi_dataset(
    expression_data    = expr_matrix,
    sample_info        = sample_info_df,
    variable_info      = variable_info_df,
    sample_info_note   = sinfo_note,
    variable_info_note = vinfo_note
  )
  
  parameter <- new(
    Class         = "tidymass_parameter",
    pacakge_name  = "msidataset",
    function_name = "convert_cardinal2msi_dataset()",
    parameter     = list(
      sample_name    = sample_name,
      section_name   = section_name,
      resolution     = resolution,
      matrix_name    = matrix_name,
      tolerance      = tolerance,
      units          = units,
      cardinal_class = paste(class(cardinal_object), collapse = ", ")
    ),
    time = Sys.time()
  )
  msi_obj@process_info$convert_from_cardinal <- parameter
  
  # ==========================================================================
  # SUMMARY
  # ==========================================================================
  
  cat(paste(rep("=", 55), collapse = ""), "\n")
  cat("Conversion complete!\n")
  cat("  Cardinal class :", class(cardinal_object)[1], "\n")
  cat(sprintf("  Tolerance used : %s %s\n",
              if (use_binning) as.character(tolerance) else "NA (exact)",
              if (use_binning) units else ""))
  cat("  Variables (m/z):", nrow(variable_info_df), "\n")
  cat("  Pixels         :", nrow(sample_info_df), "\n")
  cat(sprintf("  m/z range      : [%.4f - %.4f]\n",
              min(mz_values), max(mz_values)))
  cat(sprintf("  Spatial extent : x[%d-%d]  y[%d-%d]\n",
              min(sample_info_df$x), max(sample_info_df$x),
              min(sample_info_df$y), max(sample_info_df$y)))
  cat(paste(rep("=", 55), collapse = ""), "\n")
  
  invisible(msi_obj)
}


# ===========================================================================
# convert_msi_dataset2cardinal()
# ===========================================================================

#' Convert msi_dataset to Cardinal MSImagingExperiment
#'
#' @param msi_object An \code{msi_dataset} object.
#' @return A \code{Cardinal::MSImagingExperiment} object.
#' @export

convert_msi_dataset2cardinal <- function(msi_object) {
  
  if (!requireNamespace("Cardinal", quietly = TRUE))
    stop("Package 'Cardinal' is required.")
  
  expr  <- expression_data(msi_object)
  sinfo <- sample_info(msi_object)
  vinfo <- variable_info(msi_object)
  
  if (!all(c("x", "y") %in% colnames(sinfo)))
    stop("sample_info must contain 'x' and 'y'.")
  if (!"mz" %in% colnames(vinfo))
    stop("variable_info must contain 'mz'.")
  
  z_vals   <- if ("z" %in% colnames(sinfo)) sinfo$z else rep(1L, nrow(sinfo))
  coord_df <- data.frame(x = sinfo$x, y = sinfo$y, z = z_vals)
  
  Cardinal::MSImagingExperiment(
    imageData   = as.matrix(expr),
    featureData = Cardinal::MassDataFrame(mz = vinfo$mz),
    pixelData   = Cardinal::PositionDataFrame(coord = coord_df)
  )
}

# ===========================================================================
# batch_convert_cardinal2msi_dataset()
# ===========================================================================

#' Batch-convert a list of Cardinal objects to msi_dataset
#'
#' @param cardinal_list Named or unnamed \code{list} of Cardinal objects.
#' @param sample_names  Character vector (auto-generated if \code{NULL}).
#' @param section_names Character vector (auto-generated if \code{NULL}).
#' @param resolution    Recycled to list length.
#' @param matrix_name   Recycled to list length.
#' @param tolerance     m/z binning tolerance passed to each conversion.
#'   Can be \code{10}, \code{50}, \code{100}, or \code{NA}.
#' @param units         \code{"ppm"} or \code{"mz"}.
#' @return Named \code{list} of \code{msi_dataset} objects.
#' @export

batch_convert_cardinal2msi_dataset <- function(cardinal_list,
                                               sample_names  = NULL,
                                               section_names = NULL,
                                               resolution    = NA,
                                               matrix_name   = NA,
                                               tolerance     = 10,
                                               units         = "ppm") {
  n <- length(cardinal_list)
  if (is.null(sample_names))  sample_names  <- paste0("sample_",  seq_len(n))
  if (is.null(section_names)) section_names <- paste0("section_", seq_len(n))
  if (length(resolution)  == 1) resolution  <- rep(resolution,  n)
  if (length(matrix_name) == 1) matrix_name <- rep(matrix_name, n)
  
  msi_list <- lapply(seq_len(n), function(i) {
    cat("\n--- Object", i, "of", n, "---\n")
    convert_cardinal2msi_dataset(
      cardinal_object = cardinal_list[[i]],
      sample_name     = sample_names[i],
      section_name    = section_names[i],
      resolution      = resolution[i],
      matrix_name     = matrix_name[i],
      tolerance       = tolerance,
      units           = units
    )
  })
  
  names(msi_list) <- sample_names
  message("\nDone. ", n, " msi_dataset object(s) created.")
  invisible(msi_list)
}
