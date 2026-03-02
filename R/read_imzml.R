# ============================================================================
# Read imzML data and convert to msi_dataset
# ============================================================================

.find_imzml_file <- function(path_1) {
  if (dir.exists(path_1)) {
    files <- list.files(path_1, pattern = "\\.imzML$", ignore.case = TRUE,
                        full.names = TRUE)
    if (length(files) == 0) {
      stop("No .imzML file found in path_1 directory.")
    }
    return(files[[1]])
  }

  if (!file.exists(path_1)) {
    stop("path_1 does not exist.")
  }

  if (!grepl("\\.imzML$", path_1, ignore.case = TRUE)) {
    stop("path_1 must be an .imzML file or a directory containing one.")
  }

  path_1
}

.parse_mgf <- function(mgf_file) {
  lines <- readLines(mgf_file, warn = FALSE)
  if (length(lines) == 0) {
    return(list())
  }

  blocks <- split(lines, cumsum(lines == "BEGIN IONS"))
  blocks <- blocks[vapply(blocks, function(x) any(x == "BEGIN IONS"), logical(1))]

  spectra <- vector("list", length(blocks))
  for (i in seq_along(blocks)) {
    block <- blocks[[i]]
    block <- block[block != "BEGIN IONS" & block != "END IONS"]

    title <- NA_character_
    pepmass <- NA_real_
    rt <- NA_real_
    sample_id <- NA_character_
    peak_lines <- character()

    for (line in block) {
      if (grepl("^TITLE=", line)) {
        title <- sub("^TITLE=", "", line)
      } else if (grepl("^PEPMASS=", line)) {
        pepmass <- suppressWarnings(
          as.numeric(strsplit(sub("^PEPMASS=", "", line), " ")[[1]][1])
        )
      } else if (grepl("^RTINSECONDS=", line)) {
        rt <- suppressWarnings(as.numeric(sub("^RTINSECONDS=", "", line)))
      } else if (grepl("^SCANS=", line)) {
        sample_id <- sub("^SCANS=", "", line)
      } else if (grepl("^[0-9]", line)) {
        peak_lines <- c(peak_lines, line)
      }
    }

    peak_matrix <- if (length(peak_lines) > 0) {
      x <- do.call(rbind, strsplit(peak_lines, "[[:space:]]+"))
      mode(x) <- "numeric"
      x
    } else {
      matrix(numeric(0), ncol = 2)
    }

    spectra[[i]] <- list(
      title = title,
      precursor_mz = pepmass,
      rt = rt,
      sample_id = sample_id,
      peaks = peak_matrix
    )
  }

  spectra
}

.map_ms2_features_fast <- function(variable_mz, precursor_mz, ppm_tolerance) {
  if (length(variable_mz) == 0) {
    return(character(0))
  }

  precursor_mz <- as.numeric(precursor_mz)
  precursor_mz <- precursor_mz[is.finite(precursor_mz)]
  if (length(precursor_mz) == 0) {
    return(rep("no", length(variable_mz)))
  }

  precursor_mz <- sort(precursor_mz)
  mz <- as.numeric(variable_mz)
  mz[!is.finite(mz)] <- NA_real_

  tol <- mz * ppm_tolerance / 1e6
  lower <- mz - tol
  upper <- mz + tol

  left <- findInterval(lower - 1e-12, precursor_mz)
  right <- findInterval(upper, precursor_mz)

  matched <- right > left
  matched[is.na(mz)] <- FALSE
  ifelse(matched, "Yes", "no")
}

.extract_precursor_mz <- function(spectra) {
  if (length(spectra) == 0) {
    return(numeric(0))
  }
  vapply(spectra, function(x) x$precursor_mz, numeric(1))
}

.coerce_numeric_vector <- function(x) {
  if (is.null(x)) {
    return(numeric(0))
  }

  if (is.numeric(x) || is.integer(x)) {
    return(as.numeric(x))
  }

  if (is.character(x)) {
    out <- suppressWarnings(as.numeric(x))
    return(out[is.finite(out)])
  }

  if (is.matrix(x) || is.data.frame(x)) {
    out <- suppressWarnings(as.numeric(as.matrix(x)))
    return(out[is.finite(out)])
  }

  if (isS4(x)) {
    s4_slots <- methods::slotNames(x)
    for (nm in c("mz", "mzs", "mass", "values", "data")) {
      if (nm %in% s4_slots) {
        out <- .coerce_numeric_vector(methods::slot(x, nm))
        if (length(out) > 0) {
          return(out)
        }
      }
    }

    as_list <- try(as.list(x), silent = TRUE)
    if (!inherits(as_list, "try-error")) {
      out <- .coerce_numeric_vector(as_list)
      if (length(out) > 0) {
        return(out)
      }
    }
  }

  if (is.list(x)) {
    if (!is.null(x$mz)) {
      out <- .coerce_numeric_vector(x$mz)
      if (length(out) > 0) {
        return(out)
      }
    }

    for (i in seq_along(x)) {
      out <- .coerce_numeric_vector(x[[i]])
      if (length(out) > 0) {
        return(out)
      }
    }

    pooled <- suppressWarnings(as.numeric(unlist(x, recursive = TRUE, use.names = FALSE)))
    pooled <- pooled[is.finite(pooled)]
    return(pooled)
  }

  out <- suppressWarnings(as.numeric(as.vector(x)))
  out[is.finite(out)]
}

.extract_mz_list <- function(cardinal_object, n_pixels) {
  raw_mz <- Cardinal::mz(cardinal_object)
  mz_list <- NULL

  if (is.list(raw_mz)) {
    mz_list <- raw_mz
  } else if (methods::is(raw_mz, "List") || isS4(raw_mz) ||
             methods::is(raw_mz, "MSImagingArrays")) {
    as_list_try <- try(as.list(raw_mz), silent = TRUE)
    if (!inherits(as_list_try, "try-error")) {
      mz_list <- as_list_try
    }
  }

  if (is.null(mz_list) || length(mz_list) != n_pixels) {
    return(NULL)
  }

  lapply(mz_list, .coerce_numeric_vector)
}

.parallel_lapply <- function(X,
                             FUN,
                             backend = c("none", "future", "biocparallel"),
                             workers = 1L,
                             verbose = TRUE) {
  backend <- match.arg(backend)
  workers <- max(1L, as.integer(workers))

  if (backend == "future" && requireNamespace("future.apply", quietly = TRUE) &&
      requireNamespace("future", quietly = TRUE)) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    if (workers > 1L) {
      future::plan(future::multisession, workers = workers)
    } else {
      future::plan(future::sequential)
    }
    return(future.apply::future_lapply(X, FUN, future.seed = TRUE))
  }

  if (backend == "biocparallel" && requireNamespace("BiocParallel", quietly = TRUE)) {
    param <- if (.Platform$OS.type == "windows") {
      BiocParallel::SnowParam(workers = workers, progressbar = verbose)
    } else {
      BiocParallel::MulticoreParam(workers = workers, progressbar = verbose)
    }
    return(BiocParallel::bplapply(X, FUN, BPPARAM = param))
  }

  lapply(X, FUN)
}

.cache_key <- function(imzml_file, bin_digits, max_features) {
  ibd_file <- sub("\\.imzML$", ".ibd", imzml_file, ignore.case = TRUE)
  imzml_info <- file.info(imzml_file)
  ibd_info <- if (file.exists(ibd_file)) file.info(ibd_file) else NULL

  key_raw <- paste(
    normalizePath(imzml_file, winslash = "/", mustWork = FALSE),
    imzml_info$size,
    as.numeric(imzml_info$mtime),
    if (!is.null(ibd_info)) ibd_info$size else NA,
    if (!is.null(ibd_info)) as.numeric(ibd_info$mtime) else NA,
    as.integer(bin_digits),
    as.integer(max_features),
    sep = "|"
  )
  key <- gsub("[^A-Za-z0-9]+", "_", key_raw)
  substr(key, 1, 180)
}

.load_alignment_cache <- function(cache_dir, key) {
  cache_file <- file.path(cache_dir, paste0("align_", key, ".rds"))
  if (!file.exists(cache_file)) {
    return(NULL)
  }

  x <- try(readRDS(cache_file), silent = TRUE)
  if (inherits(x, "try-error") || !is.list(x) ||
      is.null(x$expression_data) || is.null(x$mz_values)) {
    return(NULL)
  }
  x
}

.save_alignment_cache <- function(cache_dir, key, expression_res) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_file <- file.path(cache_dir, paste0("align_", key, ".rds"))
  saveRDS(expression_res, cache_file)
  invisible(cache_file)
}

.align_processed_spectra <- function(intensity_list,
                                     mz_list,
                                     bin_digits = 4L,
                                     max_features = 150000L,
                                     verbose = TRUE,
                                     chunk_size = 256L,
                                     parallel_backend = c("none", "future", "biocparallel"),
                                     workers = 1L) {
  parallel_backend <- match.arg(parallel_backend)
  n_pixels <- length(intensity_list)
  if (length(mz_list) != n_pixels) {
    stop("Alignment error: mz_list length does not match intensity_list length.")
  }

  all_mz <- unlist(mz_list, recursive = TRUE, use.names = FALSE)
  all_mz <- suppressWarnings(as.numeric(all_mz))
  all_mz <- all_mz[is.finite(all_mz)]
  if (length(all_mz) == 0) {
    stop("Alignment error: no finite m/z values found in processed spectra.")
  }

  digits <- as.integer(bin_digits)
  if (is.na(digits) || digits < 1L) {
    digits <- 4L
  }
  axis <- sort(unique(round(all_mz, digits)))

  while (length(axis) > as.integer(max_features) && digits > 1L) {
    digits <- digits - 1L
    axis <- sort(unique(round(all_mz, digits)))
  }

  if (verbose) {
    message(
      "Processed alignment: pixels=", n_pixels,
      ", features=", length(axis),
      ", digits=", digits, "."
    )
  }

  chunk_size <- max(1L, as.integer(chunk_size))
  chunks <- split(seq_len(n_pixels), ceiling(seq_len(n_pixels) / chunk_size))

  process_chunk <- function(idx_vec) {
    if (requireNamespace("data.table", quietly = TRUE)) {
      dt_list <- vector("list", length(idx_vec))
      for (k in seq_along(idx_vec)) {
        i <- idx_vec[[k]]
        mzi <- mz_list[[i]]
        inti <- intensity_list[[i]]
        len <- min(length(mzi), length(inti))
        if (len == 0) next
        mzi <- round(mzi[seq_len(len)], digits)
        inti <- as.numeric(inti[seq_len(len)])
        keep <- is.finite(mzi) & is.finite(inti)
        if (!any(keep)) next
        dt_list[[k]] <- data.table::data.table(
          pixel = i,
          mz = mzi[keep],
          intensity = inti[keep]
        )
      }
      dt <- data.table::rbindlist(dt_list, use.names = TRUE, fill = FALSE)
      if (nrow(dt) == 0) return(NULL)
      df <- as.data.frame(dt, stringsAsFactors = FALSE)
      agg <- stats::aggregate(
        x = df[["intensity"]],
        by = list(mz = df[["mz"]], pixel = df[["pixel"]]),
        FUN = function(z) sum(z, na.rm = TRUE)
      )
      names(agg)[names(agg) == "x"] <- "intensity"
      return(agg)
    }

    rows <- vector("list", length(idx_vec))
    for (k in seq_along(idx_vec)) {
      i <- idx_vec[[k]]
      mzi <- mz_list[[i]]
      inti <- intensity_list[[i]]
      len <- min(length(mzi), length(inti))
      if (len == 0) next
      mzi <- round(mzi[seq_len(len)], digits)
      inti <- as.numeric(inti[seq_len(len)])
      keep <- is.finite(mzi) & is.finite(inti)
      if (!any(keep)) next
      agg <- tapply(inti[keep], mzi[keep], sum, na.rm = TRUE)
      rows[[k]] <- data.frame(
        pixel = i,
        mz = as.numeric(names(agg)),
        intensity = as.numeric(agg),
        stringsAsFactors = FALSE
      )
    }
    out <- do.call(rbind, rows)
    if (is.null(out) || nrow(out) == 0) return(NULL)
    out
  }

  chunk_results <- .parallel_lapply(
    X = chunks,
    FUN = process_chunk,
    backend = parallel_backend,
    workers = workers,
    verbose = verbose
  )

  expr <- matrix(0, nrow = length(axis), ncol = n_pixels)
  axis_chr <- as.character(axis)
  for (res in chunk_results) {
    if (is.null(res) || nrow(res) == 0) next
    mz_idx <- match(as.character(res$mz), axis_chr)
    keep <- !is.na(mz_idx)
    if (any(keep)) {
      expr[cbind(mz_idx[keep], as.integer(res$pixel[keep]))] <- as.numeric(res$intensity[keep])
    }
  }

  warning(
    "Processed imzML with inconsistent per-pixel spectrum lengths detected. ",
    "Applied rounded m/z alignment (digits = ", digits, ") to build common feature axis.",
    call. = FALSE
  )

  list(expression_data = expr, mz_values = axis)
}

.extract_expression_matrix <- function(cardinal_object,
                                      n_pixels,
                                      bin_digits = 4L,
                                      max_features = 150000L,
                                      verbose = TRUE,
                                      chunk_size = 256L,
                                      parallel_backend = c("none", "future", "biocparallel"),
                                      workers = 1L) {
  parallel_backend <- match.arg(parallel_backend)
  raw_spec <- Cardinal::spectra(cardinal_object)

  spec_mat <- try(as.matrix(raw_spec), silent = TRUE)
  if (!inherits(spec_mat, "try-error") && is.matrix(spec_mat)) {
    if (is.numeric(spec_mat)) {
      if (ncol(spec_mat) == n_pixels) {
        return(list(expression_data = spec_mat, mz_values = NULL))
      }
      if (nrow(spec_mat) == n_pixels) {
        return(list(expression_data = t(spec_mat), mz_values = NULL))
      }
    }

    if (nrow(spec_mat) == 1 && ncol(spec_mat) == n_pixels) {
      cols <- lapply(seq_len(n_pixels), function(j) .coerce_numeric_vector(spec_mat[1, j]))
      lens <- vapply(cols, length, integer(1))
      if (all(lens > 0) && length(unique(lens)) == 1) {
        return(list(expression_data = do.call(cbind, cols), mz_values = NULL))
      }
      mz_list <- .extract_mz_list(cardinal_object, n_pixels)
      if (!is.null(mz_list)) {
        return(.align_processed_spectra(
          intensity_list = cols,
          mz_list = mz_list,
          bin_digits = bin_digits,
          max_features = max_features,
          verbose = verbose,
          chunk_size = chunk_size,
          parallel_backend = parallel_backend,
          workers = workers
        ))
      }
    }

    if (ncol(spec_mat) == 1 && nrow(spec_mat) == n_pixels) {
      cols <- lapply(seq_len(n_pixels), function(i) .coerce_numeric_vector(spec_mat[i, 1]))
      lens <- vapply(cols, length, integer(1))
      if (all(lens > 0) && length(unique(lens)) == 1) {
        return(list(expression_data = do.call(cbind, cols), mz_values = NULL))
      }
      mz_list <- .extract_mz_list(cardinal_object, n_pixels)
      if (!is.null(mz_list)) {
        return(.align_processed_spectra(
          intensity_list = cols,
          mz_list = mz_list,
          bin_digits = bin_digits,
          max_features = max_features,
          verbose = verbose,
          chunk_size = chunk_size,
          parallel_backend = parallel_backend,
          workers = workers
        ))
      }
    }
  }

  spec_list <- NULL
  if (is.list(raw_spec)) {
    spec_list <- raw_spec
  } else if (methods::is(raw_spec, "List") || isS4(raw_spec)) {
    as_list_try <- try(as.list(raw_spec), silent = TRUE)
    if (!inherits(as_list_try, "try-error")) {
      spec_list <- as_list_try
    }
  }

  if (!is.null(spec_list) && length(spec_list) == n_pixels) {
    cols <- lapply(spec_list, .coerce_numeric_vector)
    lens <- vapply(cols, length, integer(1))
    if (all(lens > 0) && length(unique(lens)) == 1) {
      return(list(expression_data = do.call(cbind, cols), mz_values = NULL))
    }
    mz_list <- .extract_mz_list(cardinal_object, n_pixels)
    if (!is.null(mz_list)) {
      return(.align_processed_spectra(
        intensity_list = cols,
        mz_list = mz_list,
        bin_digits = bin_digits,
        max_features = max_features,
        verbose = verbose,
        chunk_size = chunk_size,
        parallel_backend = parallel_backend,
        workers = workers
      ))
    }
  }

  stop("Failed to extract expression_data matrix from Cardinal spectra container.")
}

.extract_mz_force_numeric <- function(cardinal_object, n_feature) {
  raw_mz <- Cardinal::mz(cardinal_object)
  if (methods::is(raw_mz, "MSImagingArrays") ||
      methods::is(raw_mz, "List") ||
      is.list(raw_mz)) {
    message("Normalizing MSImagingArrays for variable_info...")
  }

  candidates <- list(raw_mz)

  first_axis <- try(raw_mz[[1]], silent = TRUE)
  if (!inherits(first_axis, "try-error")) {
    candidates <- c(candidates, list(first_axis))
  }

  raw_list <- try(as.list(raw_mz), silent = TRUE)
  if (!inherits(raw_list, "try-error") && length(raw_list) > 0) {
    candidates <- c(candidates, list(raw_list[[1]]))
  }

  fdata <- try(Cardinal::fData(cardinal_object), silent = TRUE)
  fdata_df <- try(as.data.frame(fdata), silent = TRUE)
  if (!inherits(fdata_df, "try-error") && "mz" %in% colnames(fdata_df)) {
    candidates <- c(candidates, list(fdata_df$mz))
  }

  numeric_candidates <- lapply(candidates, .coerce_numeric_vector)
  numeric_candidates <- numeric_candidates[vapply(numeric_candidates, length, integer(1)) > 0]

  if (length(numeric_candidates) == 0) {
    stop("Failed to coerce m/z values to numeric vector: no numeric content found in Cardinal mz container.")
  }

  exact_idx <- which(vapply(numeric_candidates, length, integer(1)) == n_feature)
  if (length(exact_idx) > 0) {
    mz_values <- numeric_candidates[[exact_idx[[1]]]]
  } else {
    lens <- vapply(numeric_candidates, length, integer(1))
    stop(
      "m/z axis length mismatch: no candidate matched nrow(expression_data) = ", n_feature,
      ". Candidate lengths: ", paste(lens, collapse = ", "), "."
    )
  }

  if (!is.numeric(mz_values)) {
    stop("Failed to coerce m/z values to numeric vector.")
  }

  if (any(!is.finite(mz_values))) {
    stop("m/z values contain non-finite entries after numeric coercion.")
  }

  if (length(mz_values) != n_feature) {
    stop(
      "m/z axis length mismatch: length(mz_values) = ", length(mz_values),
      ", nrow(expression_data) = ", n_feature, "."
    )
  }

  mz_values
}

.build_msi_dataset_from_components <- function(expression_data,
                                               coords,
                                               mz_values,
                                               sample_id = NULL,
                                               subject_id = "subject_1",
                                               section_id = "section_1",
                                               resolution = NA_real_,
                                               matrix_name = NA_character_) {
  expression_data <- as.matrix(expression_data)
  coords <- as.matrix(coords)
  mz_values <- suppressWarnings(as.numeric(mz_values))
  stopifnot(is.numeric(mz_values))

  if (ncol(coords) < 2) {
    stop("coords must provide at least x and y columns.")
  }
  if (!is.numeric(mz_values) || any(!is.finite(mz_values))) {
    stop("mz_values must be a finite numeric vector.")
  }
  if (length(mz_values) != nrow(expression_data)) {
    stop("Length of mz_values must equal nrow(expression_data).")
  }

  n_samples <- nrow(coords)
  if (is.null(sample_id)) {
    sample_id <- paste0("pixel_", seq_len(n_samples))
  }

  sample_info <- data.frame(
    sample_id = as.character(sample_id),
    x = as.numeric(coords[, 1]),
    y = as.numeric(coords[, 2]),
    z = if (ncol(coords) >= 3) as.numeric(coords[, 3]) else rep(1, n_samples),
    subject_id = as.character(subject_id),
    section_id = as.character(section_id),
    resolution = as.numeric(resolution),
    matrix = as.character(matrix_name),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  variable_info <- data.frame(
    variable_id = generate_variable_id(mz_values),
    mz = as.numeric(mz_values),
    ms2 = rep("no", length(mz_values)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  rownames(expression_data) <- variable_info$variable_id
  colnames(expression_data) <- sample_info$sample_id

  create_msi_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )
}

#' Read imzML files into an msidataset object
#'
#' @description
#' Read MSI data from an `imzML` file (or directory containing one), construct an
#' `msi_dataset` object, parse optional MGF MS2 spectra, and save the resulting
#' object to disk. The function auto-populates MSI-specific `sample_info`
#' including `x`, `y`, and `z` coordinates.
#'
#' Coordinate system convention: pixel coordinates are treated as 1-based image
#' coordinates where `(x = 1, y = 1)` corresponds to the top-left position.
#'
#' `peak_align_tolerance` is applied as a ppm tolerance window when mapping MGF
#' precursor m/z values to MSI features in `variable_info$mz`.
#'
#' @param path_1 Character. Input path. Either:
#'   1) full path to one `.imzML` file, or
#'   2) directory containing at least one `.imzML` file.
#' @param peak_align_tolerance Numeric scalar (ppm). Tolerance used for MS2
#'   precursor-to-feature alignment.
#' @param path_2 Character. Output path. Either:
#'   1) output directory (an `.rds` name will be generated), or
#'   2) full `.rds` output file path.
#' @param bin_digits Integer. Rounded m/z digits used to align processed imzML.
#' @param max_features Integer. Maximum aligned feature count; alignment reduces
#'   precision automatically if exceeded.
#' @param verbose Logical. Print alignment progress.
#' @param chunk_size Integer. Pixels per aggregation chunk for memory control.
#' @param parallel_backend One of `none`, `future`, `biocparallel`.
#' @param workers Integer number of workers when parallel backend is enabled.
#' @param use_cache Logical. Cache aligned processed spectra result for reuse.
#' @param cache_dir Character directory used to store/read alignment cache.
#' @param resolution Numeric. Spatial resolution (for example, micrometers per
#'   pixel) to write into `sample_info$resolution`.
#' @param matrix_name Character. MALDI matrix name to write into
#'   `sample_info$matrix`.
#'
#' @return An object of class `msi_dataset`.
#' @importFrom tools file_path_sans_ext
#' @export
read_imzml <- function(path_1 = "input",
                       peak_align_tolerance = 50,
                       path_2 = "output",
                       bin_digits = 4,
                       max_features = 150000,
                       verbose = TRUE,
                       chunk_size = 256,
                       parallel_backend = c("none", "future", "biocparallel"),
                       workers = 1,
                       use_cache = TRUE,
                       cache_dir = NULL,
                       resolution = NA_real_,
                       matrix_name = NA_character_) {
  parallel_backend <- match.arg(parallel_backend)

  if (!requireNamespace("Cardinal", quietly = TRUE)) {
    stop("please install Cardinal package: BiocManager::install('Cardinal')")
  }

  imzml_file <- .find_imzml_file(path_1)
  cardinal_object <- Cardinal::readMSIData(imzml_file)
  pdata <- as.data.frame(Cardinal::pData(cardinal_object), stringsAsFactors = FALSE)
  coords <- Cardinal::coord(cardinal_object)
  n_pixels <- nrow(pdata)

  if (is.null(cache_dir)) {
    cache_dir <- file.path(dirname(imzml_file), ".msidataset_cache")
  }

  cache_hit <- FALSE
  cache_key <- .cache_key(
    imzml_file = imzml_file,
    bin_digits = bin_digits,
    max_features = max_features
  )

  expression_res <- if (isTRUE(use_cache)) .load_alignment_cache(cache_dir, cache_key) else NULL
  if (!is.null(expression_res) && ncol(expression_res$expression_data) == n_pixels) {
    cache_hit <- TRUE
    if (verbose) {
      message("Loaded aligned spectra from cache: ", cache_key)
    }
  } else {
    expression_res <- .extract_expression_matrix(
      cardinal_object = cardinal_object,
      n_pixels = n_pixels,
      bin_digits = bin_digits,
      max_features = max_features,
      verbose = verbose,
      chunk_size = chunk_size,
      parallel_backend = parallel_backend,
      workers = workers
    )
    if (isTRUE(use_cache) && !is.null(expression_res$mz_values)) {
      .save_alignment_cache(cache_dir, cache_key, expression_res)
    }
  }

  expression_data <- expression_res$expression_data

  if (ncol(expression_data) == n_pixels) {
    # correct orientation
  } else if (nrow(expression_data) == n_pixels) {
    expression_data <- t(expression_data)
  } else {
    stop(
      "Cannot align expression_data with sample_info: n_pixels = ", n_pixels,
      ", dim(expression_data) = ", paste(dim(expression_data), collapse = " x "), "."
    )
  }

  if (ncol(coords) < 2) {
    stop("Could not parse x and y coordinates from imzML.")
  }

  sample_id <- rownames(pdata)
  if (is.null(sample_id) || any(sample_id == "")) {
    sample_id <- paste0("pixel_", seq_len(nrow(pdata)))
  }

  if (!is.null(expression_res$mz_values)) {
    mz_values <- expression_res$mz_values
  } else {
    mz_values <- .extract_mz_force_numeric(
      cardinal_object = cardinal_object,
      n_feature = nrow(expression_data)
    )
  }

  object <- .build_msi_dataset_from_components(
    expression_data = expression_data,
    coords = coords,
    mz_values = mz_values,
    sample_id = sample_id,
    subject_id = file_path_sans_ext(basename(imzml_file)),
    section_id = "section_1",
    resolution = resolution,
    matrix_name = matrix_name
  )

  mgf_candidates <- list.files(dirname(imzml_file), pattern = "\\.mgf$",
                               ignore.case = TRUE, full.names = TRUE)
  if (length(mgf_candidates) > 0) {
    spectra <- .parse_mgf(mgf_candidates[[1]])
    precursor_mz <- .extract_precursor_mz(spectra)

    object@ms2_data <- spectra
    object@variable_info$ms2 <- .map_ms2_features_fast(
      variable_mz = object@variable_info$mz,
      precursor_mz = precursor_mz,
      ppm_tolerance = peak_align_tolerance
    )
  }

  object <- add_process_log(
    object = object,
    function_name = "read_imzml()",
    parameter = list(
      input = imzml_file,
      peak_align_tolerance = peak_align_tolerance,
      n_variable = nrow(object@variable_info),
      n_sample = nrow(object@sample_info),
      n_ms2 = length(object@ms2_data),
      bin_digits = bin_digits,
      max_features = max_features,
      chunk_size = chunk_size,
      parallel_backend = parallel_backend,
      workers = workers,
      cache_hit = cache_hit,
      cache_dir = cache_dir,
      resolution = resolution,
      matrix_name = matrix_name
    )
  )

  output_file <- path_2
  if (dir.exists(path_2) || !grepl("\\.rds$", path_2, ignore.case = TRUE)) {
    dir.create(path_2, recursive = TRUE, showWarnings = FALSE)
    output_file <- file.path(
      path_2,
      paste0(file_path_sans_ext(basename(imzml_file)), "_msidataset.rds")
    )
  }
  saveRDS(object, output_file)

  object
}
