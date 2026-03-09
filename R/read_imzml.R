# ============================================================================
# Read imzML data and convert to msi_dataset
# ============================================================================

.find_imzml_files <- function(path_1) {
  if (length(path_1) == 0) {
    stop("path_1 is empty.")
  }
  all_files <- character()
  for (p in as.character(path_1)) {
    if (dir.exists(p)) {
      files <- list.files(p, pattern = "\\.imzML$", ignore.case = TRUE,
                          full.names = TRUE)
      if (length(files) == 0) {
        stop("No .imzML file found in directory: ", p)
      }
      all_files <- c(all_files, files[[1]])
      next
    }
    if (!file.exists(p)) {
      stop("Path does not exist: ", p)
    }
    if (!grepl("\\.imzML$", p, ignore.case = TRUE)) {
      stop("Each path in path_1 must be an .imzML file or directory.")
    }
    all_files <- c(all_files, p)
  }
  unique(normalizePath(all_files, winslash = "/", mustWork = TRUE))
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

.cli_progress_start <- function(total,
                                name = "MSI import",
                                enable = TRUE,
                                format = "{name}\n{bar} {percent} | ETA: {eta}") {
  if (!isTRUE(enable) || !requireNamespace("cli", quietly = TRUE)) {
    return(NULL)
  }
  tryCatch(
    cli::cli_progress_bar(
      name = name,
      total = as.integer(total),
      format = format,
      clear = FALSE
    ),
    error = function(e) NULL
  )
}

.cli_progress_update <- function(id, msg = NULL, inc = 1L) {
  if (is.null(id) || !requireNamespace("cli", quietly = TRUE)) {
    return(invisible(NULL))
  }
  try(
    cli::cli_progress_update(
      id = id,
      inc = as.integer(inc),
      status = if (is.null(msg)) "" else as.character(msg)
    ),
    silent = TRUE
  )
  invisible(NULL)
}

.cli_progress_done <- function(id) {
  if (is.null(id) || !requireNamespace("cli", quietly = TRUE)) {
    return(invisible(NULL))
  }
  try(cli::cli_progress_done(id = id), silent = TRUE)
  invisible(NULL)
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
                                     workers = 1L,
                                     progress = interactive()) {
  parallel_backend <- match.arg(parallel_backend)
  n_pixels <- length(intensity_list)
  if (length(mz_list) != n_pixels) {
    stop("Alignment error: mz_list length does not match intensity_list length.")
  }

  digits <- as.integer(bin_digits)
  if (is.na(digits) || digits < 1L) {
    digits <- 4L
  }

  chunk_size <- max(1L, as.integer(chunk_size))
  chunks <- split(seq_len(n_pixels), ceiling(seq_len(n_pixels) / chunk_size))
  n_chunks <- length(chunks)

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

  run_chunk_processing <- function(update_fun = NULL) {
    if (parallel_backend == "none") {
      out <- vector("list", n_chunks)
      for (ci in seq_along(chunks)) {
        out[[ci]] <- process_chunk(chunks[[ci]])
        if (!is.null(update_fun)) {
          update_fun(sprintf("Aligning spectra %d / %d", ci, n_chunks))
        }
      }
      return(out)
    }

    out <- .parallel_lapply(
      X = chunks,
      FUN = process_chunk,
      backend = parallel_backend,
      workers = workers,
      verbose = verbose
    )
    if (!is.null(update_fun)) {
      for (ci in seq_along(out)) {
        update_fun(sprintf("Aligning spectra %d / %d", ci, n_chunks))
      }
    }
    out
  }

  read_cli_id <- .cli_progress_start(
    total = n_chunks,
    name = "Reading spectra",
    enable = isTRUE(progress)
  )
  all_mz_chunks <- vector("list", n_chunks)
  if (isTRUE(progress) && requireNamespace("progressr", quietly = TRUE)) {
    progressr::with_progress({
      p_read <- progressr::progressor(steps = n_chunks)
      for (ci in seq_along(chunks)) {
        idx_vec <- chunks[[ci]]
        mz_chunk <- unlist(mz_list[idx_vec], recursive = TRUE, use.names = FALSE)
        mz_chunk <- suppressWarnings(as.numeric(mz_chunk))
        all_mz_chunks[[ci]] <- mz_chunk[is.finite(mz_chunk)]
        msg <- sprintf("Reading spectra %d / %d", ci, n_chunks)
        .cli_progress_update(read_cli_id, msg)
        p_read(msg)
      }
    })
  } else {
    for (ci in seq_along(chunks)) {
      idx_vec <- chunks[[ci]]
      mz_chunk <- unlist(mz_list[idx_vec], recursive = TRUE, use.names = FALSE)
      mz_chunk <- suppressWarnings(as.numeric(mz_chunk))
      all_mz_chunks[[ci]] <- mz_chunk[is.finite(mz_chunk)]
      .cli_progress_update(read_cli_id, sprintf("Reading spectra %d / %d", ci, n_chunks))
    }
  }
  .cli_progress_done(read_cli_id)

  all_mz <- unlist(all_mz_chunks, recursive = TRUE, use.names = FALSE)
  if (length(all_mz) == 0) {
    stop("Alignment error: no finite m/z values found in processed spectra.")
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

  align_cli_id <- .cli_progress_start(
    total = n_chunks,
    name = "Aligning m/z axis",
    enable = isTRUE(progress)
  )
  chunk_results_env <- new.env(parent = emptyenv())
  chunk_results_env$value <- NULL
  if (isTRUE(progress) && requireNamespace("progressr", quietly = TRUE)) {
    progressr::with_progress({
      p_align <- progressr::progressor(steps = n_chunks)
      chunk_results_env$value <- run_chunk_processing(update_fun = function(msg) {
        .cli_progress_update(align_cli_id, msg)
        p_align(msg)
      })
    })
  } else {
    chunk_results_env$value <- run_chunk_processing(update_fun = function(msg) {
      .cli_progress_update(align_cli_id, msg)
    })
  }
  chunk_results <- chunk_results_env$value
  .cli_progress_done(align_cli_id)

  expr <- matrix(0, nrow = length(axis), ncol = n_pixels)
  axis_chr <- as.character(axis)
  fill_cli_id <- .cli_progress_start(
    total = n_chunks,
    name = "Constructing feature matrix",
    enable = isTRUE(progress)
  )
  fill_loop <- function(progress_fun = NULL) {
    for (ci in seq_along(chunk_results)) {
      res <- chunk_results[[ci]]
      if (!is.null(res) && nrow(res) > 0) {
        mz_idx <- match(as.character(res$mz), axis_chr)
        keep <- !is.na(mz_idx)
        if (any(keep)) {
          expr[cbind(mz_idx[keep], as.integer(res$pixel[keep]))] <- as.numeric(res$intensity[keep])
        }
      }
      msg <- sprintf("Constructing feature matrix %d / %d", ci, n_chunks)
      .cli_progress_update(fill_cli_id, msg)
      if (!is.null(progress_fun)) {
        progress_fun(msg)
      }
    }
  }
  if (isTRUE(progress) && requireNamespace("progressr", quietly = TRUE)) {
    progressr::with_progress({
      p_fill <- progressr::progressor(steps = n_chunks)
      fill_loop(progress_fun = p_fill)
    })
  } else {
    fill_loop()
  }
  .cli_progress_done(fill_cli_id)

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
                                      workers = 1L,
                                      progress = interactive()) {
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
          workers = workers,
          progress = progress
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
          workers = workers,
          progress = progress
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
        workers = workers,
        progress = progress
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
    sample_id <- .make_spectrum_ids(n_samples)
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
    has_ms2 = rep("no", length(mz_values)),
    ms2_id = rep(NA_character_, length(mz_values)),
    ms2_type = rep(NA_character_, length(mz_values)),
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

.read_single_imzml_section <- function(imzml_file,
                                       bin_digits,
                                       max_features,
                                       verbose,
                                       chunk_size,
                                       parallel_backend,
                                       workers,
                                       use_cache,
                                       progress,
                                       cache_dir = NULL) {
  cardinal_object <- Cardinal::readMSIData(imzml_file)
  pdata <- as.data.frame(Cardinal::pData(cardinal_object), stringsAsFactors = FALSE)
  coords <- Cardinal::coord(cardinal_object)
  n_pixels <- nrow(pdata)

  if (is.null(cache_dir)) {
    cache_dir <- file.path(dirname(imzml_file), ".msidataset_cache")
  }

  cache_key <- .cache_key(
    imzml_file = imzml_file,
    bin_digits = bin_digits,
    max_features = max_features
  )
  cache_hit <- FALSE
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
      workers = workers,
      progress = progress
    )
    if (isTRUE(use_cache) && !is.null(expression_res$mz_values)) {
      .save_alignment_cache(cache_dir, cache_key, expression_res)
    }
  }

  expression_data <- expression_res$expression_data
  if (ncol(expression_data) == n_pixels) {
    # aligned
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
  if (!is.null(expression_res$mz_values)) {
    mz_values <- expression_res$mz_values
  } else {
    mz_values <- .extract_mz_force_numeric(
      cardinal_object = cardinal_object,
      n_feature = nrow(expression_data)
    )
  }

  list(
    expression_data = as.matrix(expression_data),
    mz_values = as.numeric(mz_values),
    coords = coords,
    cache_hit = cache_hit,
    cache_dir = cache_dir
  )
}

.align_section_matrices <- function(section_data, bin_digits) {
  mz_rounded <- lapply(section_data, function(x) round(as.numeric(x$mz_values), as.integer(bin_digits)))
  all_mz <- sort(unique(unlist(mz_rounded, recursive = TRUE, use.names = FALSE)))
  if (length(all_mz) == 0) {
    stop("No m/z values found across sections.")
  }
  aligned <- lapply(seq_along(section_data), function(i) {
    x <- section_data[[i]]
    expr <- as.matrix(x$expression_data)
    idx <- match(mz_rounded[[i]], all_mz)
    mat <- matrix(0, nrow = length(all_mz), ncol = ncol(expr))
    if (any(duplicated(idx))) {
      agg <- rowsum(expr, group = idx, reorder = FALSE)
      row_idx <- suppressWarnings(as.integer(rownames(agg)))
      mat[row_idx, ] <- agg
    } else {
      mat[idx, ] <- expr
    }
    mat
  })
  list(expression_data = do.call(cbind, aligned), mz_values = all_mz)
}

#' Read imzML files into an msidataset object
#'
#' @description
#' Read MSI data from one or more `imzML` files, construct an
#' `msi_dataset` object, parse optional MGF MS2 spectra, and save the resulting
#' object to disk. The function auto-populates MSI-specific `sample_info`
#' including `x`, `y`, and section-level `z` coordinates.
#'
#' Coordinate system convention: pixel coordinates are treated as 1-based image
#' coordinates where `(x = 1, y = 1)` corresponds to the top-left position.
#'
#' `peak_align_tolerance` is applied as a ppm tolerance window when mapping MGF
#' precursor m/z values to MSI features in `variable_info$mz`.
#'
#' @param path_1 Character vector. Input path(s). Each element can be:
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
#' @param progress Logical. Show real-time progress for major import stages.
#' @param resolution Numeric. Spatial resolution (for example, micrometers per
#'   pixel) to write into `sample_info$resolution`.
#' @param matrix_name Character. MALDI matrix name to write into
#'   `sample_info$matrix`.
#'
#' @return An object of class `msi_dataset`.
#' @importFrom tools file_path_sans_ext
#' @export
read_imzml <- function(path_1 = "input",
                       peak_align_tolerance = 10,
                       path_2 = "output",
                       bin_digits = 4,
                       max_features = 150000,
                       verbose = TRUE,
                       chunk_size = 256,
                       parallel_backend = c("none", "future", "biocparallel"),
                       workers = 1,
                       use_cache = TRUE,
                       cache_dir = NULL,
                       progress = interactive(),
                       resolution = NA_real_,
                       matrix_name = NA_character_) {
  parallel_backend <- match.arg(parallel_backend)

  if (!requireNamespace("Cardinal", quietly = TRUE)) {
    stop("please install Cardinal package: BiocManager::install('Cardinal')")
  }

  imzml_files <- .find_imzml_files(path_1)
  section_data <- lapply(
    imzml_files,
    function(f) {
      .read_single_imzml_section(
        imzml_file = f,
        bin_digits = bin_digits,
        max_features = max_features,
        verbose = verbose,
        chunk_size = chunk_size,
        parallel_backend = parallel_backend,
        workers = workers,
        use_cache = use_cache,
        progress = progress,
        cache_dir = cache_dir
      )
    }
  )

  aligned <- .align_section_matrices(section_data, bin_digits = bin_digits)
  expression_data <- aligned$expression_data
  mz_values <- aligned$mz_values

  sample_info_list <- vector("list", length(section_data))
  offset <- 0L
  for (i in seq_along(section_data)) {
    sec <- section_data[[i]]
    coords <- as.matrix(sec$coords)
    n_pix <- nrow(coords)
    ids <- paste0("spectrum_", offset + seq_len(n_pix))
    sample_info_list[[i]] <- data.frame(
      sample_id = ids,
      x = as.numeric(coords[, 1]),
      y = as.numeric(coords[, 2]),
      z = rep(as.numeric(i), n_pix),
      subject_id = rep(file_path_sans_ext(basename(imzml_files[[i]])), n_pix),
      section_id = rep(paste0("section_", i), n_pix),
      resolution = rep(as.numeric(resolution), n_pix),
      matrix = rep(as.character(matrix_name), n_pix),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    offset <- offset + n_pix
  }
  sample_info <- do.call(rbind, sample_info_list)

  variable_info <- data.frame(
    variable_id = generate_variable_id(mz_values),
    mz = as.numeric(mz_values),
    ms2 = rep("no", length(mz_values)),
    has_ms2 = rep("no", length(mz_values)),
    ms2_id = rep(NA_character_, length(mz_values)),
    ms2_type = rep(NA_character_, length(mz_values)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  rownames(expression_data) <- variable_info$variable_id
  colnames(expression_data) <- sample_info$sample_id

  object <- create_msi_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

  mgf_dirs <- unique(dirname(imzml_files))
  mgf_candidates <- unlist(lapply(mgf_dirs, function(x) {
    list.files(x, pattern = "\\.mgf$", ignore.case = TRUE, full.names = TRUE)
  }))
  if (length(mgf_candidates) > 0) {
    object <- mutate_ms2(
      object = object,
      path = mgf_dirs,
      ms1.ms2.match.mz.tol = peak_align_tolerance
    )
  }

  object <- add_process_log(
    object = object,
    function_name = "read_imzml()",
    parameter = list(
      input = imzml_files,
      n_section = length(imzml_files),
      peak_align_tolerance = peak_align_tolerance,
      n_variable = nrow(object@variable_info),
      n_sample = nrow(object@sample_info),
      n_ms2 = if (!is.null(object@ms2_data$ms2_table)) {
        length(unique(object@ms2_data$ms2_table$ms2_id))
      } else {
        length(object@ms2_data)
      },
      bin_digits = bin_digits,
      max_features = max_features,
      chunk_size = chunk_size,
      parallel_backend = parallel_backend,
      workers = workers,
      progress = progress,
      cache_hit = vapply(section_data, function(x) isTRUE(x$cache_hit), logical(1)),
      cache_dir = unique(vapply(section_data, function(x) x$cache_dir, character(1))),
      resolution = resolution,
      matrix_name = matrix_name
    )
  )

  output_file <- path_2
  if (dir.exists(path_2) || !grepl("\\.rds$", path_2, ignore.case = TRUE)) {
    dir.create(path_2, recursive = TRUE, showWarnings = FALSE)
    out_stub <- if (length(imzml_files) == 1) {
      file_path_sans_ext(basename(imzml_files[[1]]))
    } else {
      paste0("multi_section_", length(imzml_files), "sections")
    }
    output_file <- file.path(
      path_2,
      paste0(out_stub, "_msidataset.rds")
    )
  }
  saveRDS(object, output_file)

  object
}
