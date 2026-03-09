# ============================================================================
# MS2 mutation/filter/extract for msi_dataset
# ============================================================================

.parse_ms2_files <- function(path) {
  path <- as.character(path)
  files <- character()
  for (p in path) {
    if (dir.exists(p)) {
      files <- c(
        files,
        list.files(p, pattern = "\\.mgf$", ignore.case = TRUE,
                   full.names = TRUE, recursive = TRUE)
      )
    } else if (file.exists(p) && grepl("\\.mgf$", p, ignore.case = TRUE)) {
      files <- c(files, p)
    }
  }
  unique(files)
}

#' Match m/z values by ppm tolerance
#'
#' @param msi_mz Numeric vector of MSI feature m/z values.
#' @param precursor_mz Numeric scalar precursor m/z.
#' @param tolerance_ppm Numeric ppm tolerance. Default 10.
#' @return Logical vector indicating matches in `msi_mz`.
#' @export
mz_match_ppm <- function(msi_mz, precursor_mz, tolerance_ppm = 10) {
  msi_mz <- as.numeric(msi_mz)
  precursor_mz <- as.numeric(precursor_mz)
  if (length(precursor_mz) != 1 || !is.finite(precursor_mz) || precursor_mz <= 0) {
    return(rep(FALSE, length(msi_mz)))
  }
  ppm <- abs(msi_mz - precursor_mz) / precursor_mz * 1e6
  is.finite(ppm) & ppm <= as.numeric(tolerance_ppm)
}

.closest_feature_index <- function(feature_mz, precursor_mz, mz_tol_ppm) {
  feature_mz <- as.numeric(feature_mz)
  precursor_mz <- as.numeric(precursor_mz)
  out <- rep(NA_integer_, length(precursor_mz))
  for (i in seq_along(precursor_mz)) {
    pmz <- precursor_mz[[i]]
    if (!is.finite(pmz) || pmz <= 0) {
      next
    }
    matched <- which(mz_match_ppm(feature_mz, pmz, tolerance_ppm = mz_tol_ppm))
    if (length(matched) == 0) {
      next
    }
    ppm <- abs(feature_mz[matched] - pmz) / pmz * 1e6
    out[[i]] <- matched[[which.min(ppm)]]
  }
  out
}

.normalize_ms2_sample_id <- function(raw_id, valid_sample_id) {
  sid <- as.character(raw_id)
  sid[!nzchar(sid)] <- NA_character_
  out <- rep(NA_character_, length(sid))

  direct <- sid %in% valid_sample_id
  out[direct] <- sid[direct]

  pixel_pat <- grepl("^pixel_[0-9]+$", sid)
  if (any(pixel_pat, na.rm = TRUE)) {
    idx <- suppressWarnings(as.integer(sub("^pixel_", "", sid[pixel_pat])))
    ok <- !is.na(idx) & idx >= 1L & idx <= length(valid_sample_id)
    mapped <- rep(NA_character_, length(idx))
    mapped[ok] <- valid_sample_id[idx[ok]]
    out[pixel_pat] <- mapped
  }

  numeric_pat <- grepl("^[0-9]+$", sid)
  if (any(numeric_pat, na.rm = TRUE)) {
    idx <- suppressWarnings(as.integer(sid[numeric_pat]))
    ok <- !is.na(idx) & idx >= 1L & idx <= length(valid_sample_id)
    mapped <- rep(NA_character_, length(idx))
    mapped[ok] <- valid_sample_id[idx[ok]]
    out[numeric_pat] <- mapped
  }

  out
}

.preprocess_ms2_peaks <- function(peaks,
                                  precursor_mz,
                                  precursor_tolerance = 0,
                                  remove_precursor_peak = TRUE) {
  precursor_mz <- as.numeric(precursor_mz)
  precursor_tolerance <- as.numeric(precursor_tolerance)
  if (!is.finite(precursor_tolerance) || precursor_tolerance < 0) {
    precursor_tolerance <- 0
  }

  if (is.null(peaks) || length(peaks) == 0) {
    return(data.frame(
      fragment_mz = numeric(),
      fragment_intensity = numeric(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  peak_mat <- as.matrix(peaks)
  if (ncol(peak_mat) < 2) {
    return(data.frame(
      fragment_mz = numeric(),
      fragment_intensity = numeric(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  frag_mz <- suppressWarnings(as.numeric(peak_mat[, 1]))
  frag_int <- suppressWarnings(as.numeric(peak_mat[, 2]))
  keep <- is.finite(frag_mz) & is.finite(frag_int) & frag_int > 0
  frag_mz <- frag_mz[keep]
  frag_int <- frag_int[keep]

  if (length(frag_mz) == 0) {
    return(data.frame(
      fragment_mz = numeric(),
      fragment_intensity = numeric(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  if (is.finite(precursor_mz)) {
    upper_limit <- precursor_mz + precursor_tolerance
    keep_upper <- frag_mz <= upper_limit
    frag_mz <- frag_mz[keep_upper]
    frag_int <- frag_int[keep_upper]

    if (isTRUE(remove_precursor_peak) && length(frag_mz) > 0) {
      eq_tol <- max(1e-8, precursor_tolerance)
      keep_not_precursor <- abs(frag_mz - precursor_mz) > eq_tol
      frag_mz <- frag_mz[keep_not_precursor]
      frag_int <- frag_int[keep_not_precursor]
    }
  }

  # Normalize after filtering so the retained fragment set is scaled to [0, 1].
  if (length(frag_int) > 0) {
    max_int <- suppressWarnings(max(frag_int, na.rm = TRUE))
    if (is.finite(max_int) && max_int > 0) {
      frag_int <- frag_int / max_int
    }
  }

  data.frame(
    fragment_mz = as.numeric(frag_mz),
    fragment_intensity = as.numeric(frag_int),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.empty_ms2_table <- function() {
  data.frame(
    ms2_id = character(),
    precursor_mz = numeric(),
    precursor_variable_id = character(),
    fragment_mz = numeric(),
    fragment_intensity = numeric(),
    ms2_type = character(),
    source_file = character(),
    sample_id = character(),
    rt = numeric(),
    confidence_score = numeric(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.ensure_ms2_table <- function(x) {
  template <- .empty_ms2_table()
  if (is.null(x) || nrow(as.data.frame(x)) == 0) {
    return(template)
  }

  out <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  needed <- colnames(template)
  missing <- setdiff(needed, colnames(out))
  for (nm in missing) {
    out[[nm]] <- template[[nm]]
  }
  out <- out[, needed, drop = FALSE]

  out$ms2_id <- as.character(out$ms2_id)
  out$precursor_mz <- as.numeric(out$precursor_mz)
  out$precursor_variable_id <- as.character(out$precursor_variable_id)
  out$fragment_mz <- as.numeric(out$fragment_mz)
  out$fragment_intensity <- as.numeric(out$fragment_intensity)
  out$ms2_type <- as.character(out$ms2_type)
  out$source_file <- as.character(out$source_file)
  out$sample_id <- as.character(out$sample_id)
  out$rt <- as.numeric(out$rt)
  out$confidence_score <- as.numeric(out$confidence_score)

  out
}

.ms2_table_from_long <- function(ms2_data) {
  if (!is.list(ms2_data) || length(ms2_data) == 0) {
    return(.empty_ms2_table())
  }
  if (!is.null(ms2_data$ms2_table)) {
    return(.ensure_ms2_table(ms2_data$ms2_table))
  }
  NULL
}

.ms2_table_from_legacy <- function(ms2_data) {
  if (!is.list(ms2_data) || length(ms2_data) == 0) {
    return(.empty_ms2_table())
  }

  nms <- names(ms2_data)
  if (is.null(nms) || any(!nzchar(nms))) {
    nms <- paste0("legacy_", seq_along(ms2_data))
  }

  rec <- lapply(seq_along(ms2_data), function(i) {
    blk <- ms2_data[[i]]
    if (!(is.list(blk) && !is.null(blk$ms2_info) && !is.null(blk$ms2_spectra))) {
      return(NULL)
    }
    info <- as.data.frame(blk$ms2_info, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(info) == 0) {
      return(NULL)
    }
    if (!("ms2_spectrum_id" %in% colnames(info))) {
      info$ms2_spectrum_id <- paste0(nms[[i]], "_", seq_len(nrow(info)))
    }

    do.call(rbind, lapply(seq_len(nrow(info)), function(k) {
      sid <- as.character(info$ms2_spectrum_id[[k]])
      peaks <- blk$ms2_spectra[[sid]]
      if (is.null(peaks)) {
        peaks <- blk$ms2_spectra[[k]]
      }
      peaks <- .preprocess_ms2_peaks(
        peaks = peaks,
        precursor_mz = suppressWarnings(as.numeric(info$precursor_mz[[k]])),
        precursor_tolerance = 0,
        remove_precursor_peak = FALSE
      )
      if (nrow(peaks) == 0) {
        return(NULL)
      }
      data.frame(
        ms2_id = sid,
        precursor_mz = suppressWarnings(as.numeric(info$precursor_mz[[k]])),
        precursor_variable_id = as.character(info$variable_id[[k]]),
        fragment_mz = peaks$fragment_mz,
        fragment_intensity = peaks$fragment_intensity,
        ms2_type = "actual",
        source_file = if ("ms2_file" %in% colnames(info)) as.character(info$ms2_file[[k]]) else nms[[i]],
        sample_id = if ("sample_id" %in% colnames(info)) as.character(info$sample_id[[k]]) else NA_character_,
        rt = if ("rt" %in% colnames(info)) suppressWarnings(as.numeric(info$rt[[k]])) else NA_real_,
        confidence_score = NA_real_,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }))
  })

  out <- do.call(rbind, rec)
  if (is.null(out) || nrow(out) == 0) {
    .empty_ms2_table()
  } else {
    .ensure_ms2_table(out)
  }
}

.get_ms2_table <- function(object) {
  long <- .ms2_table_from_long(object@ms2_data)
  if (!is.null(long)) {
    return(long)
  }
  .ms2_table_from_legacy(object@ms2_data)
}

.set_ms2_table <- function(object, ms2_table) {
  object@ms2_data <- list(ms2_table = .ensure_ms2_table(ms2_table))
  object
}

.sync_variable_ms2 <- function(object) {
  vinfo <- object@variable_info
  if (!all(c("variable_id", "mz") %in% colnames(vinfo))) {
    return(object)
  }

  if (!("ms2" %in% colnames(vinfo))) {
    vinfo$ms2 <- "no"
  }
  if (!("has_ms2" %in% colnames(vinfo))) {
    vinfo$has_ms2 <- "no"
  }
  if (!("ms2_id" %in% colnames(vinfo))) {
    vinfo$ms2_id <- NA_character_
  }
  if (!("ms2_type" %in% colnames(vinfo))) {
    vinfo$ms2_type <- NA_character_
  }

  ms2_table <- .get_ms2_table(object)
  if (nrow(ms2_table) == 0) {
    vinfo$ms2 <- "no"
    vinfo$has_ms2 <- "no"
    vinfo$ms2_id <- NA_character_
    vinfo$ms2_type <- NA_character_
    object@variable_info <- vinfo
    return(object)
  }

  by_var <- split(ms2_table, ms2_table$precursor_variable_id)
  vinfo$ms2 <- "no"
  vinfo$has_ms2 <- "no"
  vinfo$ms2_id <- NA_character_
  vinfo$ms2_type <- NA_character_

  idx <- match(names(by_var), vinfo$variable_id)
  ok <- !is.na(idx)
  for (i in which(ok)) {
    rows <- by_var[[i]]
    cur <- idx[[i]]
    vinfo$ms2[[cur]] <- "Yes"
    vinfo$has_ms2[[cur]] <- "yes"
    vinfo$ms2_id[[cur]] <- paste(unique(rows$ms2_id), collapse = ";")
    vinfo$ms2_type[[cur]] <- paste(unique(rows$ms2_type), collapse = ";")
  }

  object@variable_info <- vinfo
  object
}

.flatten_ms2_data <- function(object) {
  ms2_table <- .get_ms2_table(object)
  if (nrow(ms2_table) == 0) {
    return(list(
      info = data.frame(
        ms2_spectrum_id = character(),
        variable_id = character(),
        precursor_mz = numeric(),
        rt = numeric(),
        sample_id = character(),
        ms2_file = character(),
        source_name = character(),
        ms2_type = character(),
        confidence_score = numeric(),
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      spectra = list(),
      table = ms2_table
    ))
  }

  info <- unique(ms2_table[, c(
    "ms2_id", "precursor_variable_id", "precursor_mz", "rt", "sample_id", "source_file", "ms2_type", "confidence_score"
  ), drop = FALSE])
  colnames(info) <- c(
    "ms2_spectrum_id", "variable_id", "precursor_mz", "rt", "sample_id", "ms2_file", "ms2_type", "confidence_score"
  )
  info$source_name <- "ms2_table"

  spectra_split <- split(ms2_table, ms2_table$ms2_id)
  spectra <- lapply(spectra_split, function(x) {
    as.matrix(x[, c("fragment_mz", "fragment_intensity"), drop = FALSE])
  })

  list(info = info, spectra = spectra, table = ms2_table)
}

#' Mutate MS2 data in msi_dataset object
#'
#' @param object A \code{msi_dataset} object.
#' @param path Directory/file path(s) containing MGF files.
#' @param ms1.ms2.match.mz.tol Numeric ppm tolerance for feature matching.
#' @param precursor_tolerance Numeric precursor tolerance for fragment filtering.
#'   Fragments are kept only when `fragment_mz <= precursor_mz + precursor_tolerance`.
#' @param remove_precursor_peak Logical. Remove precursor peak from fragments.
#' @return Updated \code{msi_dataset} with structured \code{ms2_data}.
#' @export
mutate_ms2 <- function(object, ...) {
  UseMethod("mutate_ms2")
}

#' @export
mutate_ms2.msi_dataset <- function(object,
                                   path = ".",
                                   ms1.ms2.match.mz.tol = 10,
                                   precursor_tolerance = 0,
                                   remove_precursor_peak = TRUE) {
  # Step 1: validate input
  # Step 2: parse and preprocess MGF spectra
  # Step 3: match precursor m/z to MSI variables by ppm
  # Step 4: store long-format ms2_data and update variable_info linkage
  # Step 5: run validObject()
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }

  ms2_files <- .parse_ms2_files(path)
  if (length(ms2_files) == 0) {
    stop("No mgf files found in path.")
  }

  parsed <- lapply(ms2_files, .parse_mgf)
  keep_file <- vapply(parsed, length, integer(1)) > 0
  parsed <- parsed[keep_file]
  ms2_files <- ms2_files[keep_file]
  if (length(parsed) == 0) {
    stop("MGF parsing failed or no spectra found.")
  }

  spectra <- do.call(c, parsed)
  file_labels <- rep(basename(ms2_files), vapply(parsed, length, integer(1)))
  precursor_mz <- vapply(spectra, function(x) as.numeric(x$precursor_mz), numeric(1))

  idx <- .closest_feature_index(
    feature_mz = object@variable_info$mz,
    precursor_mz = precursor_mz,
    mz_tol_ppm = ms1.ms2.match.mz.tol
  )
  keep <- !is.na(idx)
  if (!any(keep)) {
    message(crayon::red("No variable matched with MS2 spectra."))
    object <- .set_ms2_table(object, .empty_ms2_table())
    object <- .sync_variable_ms2(object)
    return(object)
  }

  spectra <- spectra[keep]
  idx <- idx[keep]
  precursor_mz <- precursor_mz[keep]
  file_labels <- file_labels[keep]

  valid_sample_id <- object@sample_info$sample_id
  raw_sample_id <- vapply(spectra, function(x) {
    sid <- x$sample_id
    if (is.null(sid)) "" else as.character(sid)
  }, character(1))
  normalized_sample_id <- .normalize_ms2_sample_id(raw_sample_id, valid_sample_id)

  ms2_rows <- lapply(seq_along(spectra), function(i) {
    peaks <- .preprocess_ms2_peaks(
      peaks = spectra[[i]]$peaks,
      precursor_mz = precursor_mz[[i]],
      precursor_tolerance = precursor_tolerance,
      remove_precursor_peak = remove_precursor_peak
    )
    if (nrow(peaks) == 0) {
      return(NULL)
    }

    data.frame(
      ms2_id = paste0("ms2_", i),
      precursor_mz = precursor_mz[[i]],
      precursor_variable_id = object@variable_info$variable_id[idx[[i]]],
      fragment_mz = peaks$fragment_mz,
      fragment_intensity = peaks$fragment_intensity,
      ms2_type = "actual",
      source_file = file_labels[[i]],
      sample_id = normalized_sample_id[[i]],
      rt = suppressWarnings(as.numeric(spectra[[i]]$rt)),
      confidence_score = NA_real_,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  ms2_table <- do.call(rbind, ms2_rows)
  if (is.null(ms2_table) || nrow(ms2_table) == 0) {
    ms2_table <- .empty_ms2_table()
  }

  object <- .set_ms2_table(object, ms2_table)
  object <- .sync_variable_ms2(object)

  object <- add_process_log(
    object = object,
    function_name = "mutate_ms2()",
    parameter = list(
      path = path,
      n_ms2 = length(unique(ms2_table$ms2_id)),
      n_fragment_rows = nrow(ms2_table),
      ms1.ms2.match.mz.tol = ms1.ms2.match.mz.tol,
      precursor_tolerance = precursor_tolerance,
      remove_precursor_peak = remove_precursor_peak
    )
  )
  methods::validObject(object)
  object
}

#' Extract MS2 information from msi_dataset
#'
#' @param object A \code{msi_dataset} object.
#' @param variable_id Optional precursor variable_id.
#' @param mz Optional precursor m/z. Nearest precursor variable is selected.
#' @param tolerance_ppm PPM tolerance used when selecting by `mz`.
#' @param ms2_type Optional filter (`"actual"`, `"pseudo"`, or
#'   `"pseudo_candidate"`) when multiple
#'   spectra are linked to one precursor.
#' @param what Backward-compatible extractor mode:
#'   one of \code{"ms2_data"}, \code{"ms2_info"}, or \code{"ms2_spectra"}.
#' @return Either legacy extracted content (`what`) or a structured MS2 list.
#' @export
extract_ms2 <- function(object, ...) {
  UseMethod("extract_ms2")
}

#' @export
extract_ms2.msi_dataset <- function(object,
                                    variable_id = NULL,
                                    mz = NULL,
                                    tolerance_ppm = 10,
                                    ms2_type = NULL,
                                    what = c("ms2_data", "ms2_info", "ms2_spectra")) {
  flat <- .flatten_ms2_data(object)

  # Backward-compatible interface
  if (!is.null(what) && is.null(variable_id) && is.null(mz)) {
    what <- match.arg(what)
    if (what == "ms2_data") {
      return(object@ms2_data)
    }
    if (what == "ms2_info") {
      return(flat$info)
    }
    return(flat$spectra)
  }

  ms2_table <- flat$table
  if (nrow(ms2_table) == 0) {
    return(NULL)
  }

  if (!is.null(mz) && is.null(variable_id)) {
    m <- mz_match_ppm(object@variable_info$mz, mz, tolerance_ppm = tolerance_ppm)
    if (!any(m)) {
      d <- abs(object@variable_info$mz - as.numeric(mz))
      variable_id <- object@variable_info$variable_id[[which.min(d)]]
    } else {
      idx <- which(m)
      d <- abs(object@variable_info$mz[idx] - as.numeric(mz))
      variable_id <- object@variable_info$variable_id[idx[[which.min(d)]]]
    }
  }

  if (is.null(variable_id)) {
    stop("Please provide either variable_id or mz.")
  }

  sub <- ms2_table[ms2_table$precursor_variable_id %in% as.character(variable_id), , drop = FALSE]
  if (!is.null(ms2_type)) {
    sub <- sub[sub$ms2_type %in% as.character(ms2_type), , drop = FALSE]
  }
  if (nrow(sub) == 0) {
    return(NULL)
  }

  ms2_id <- unique(sub$ms2_id)[[1]]
  rec <- sub[sub$ms2_id == ms2_id, , drop = FALSE]
  rec <- rec[order(rec$fragment_mz), , drop = FALSE]

  list(
    precursor_mz = unique(rec$precursor_mz)[[1]],
    fragments = as.numeric(rec$fragment_mz),
    intensity = as.numeric(rec$fragment_intensity),
    ms2_type = unique(rec$ms2_type)[[1]],
    metadata = list(
      ms2_id = ms2_id,
      precursor_variable_id = unique(rec$precursor_variable_id)[[1]],
      source_file = unique(rec$source_file)[[1]],
      sample_id = unique(rec$sample_id),
      rt = unique(rec$rt),
      confidence_score = unique(rec$confidence_score)
    )
  )
}

.smooth_spatial_signal <- function(x,
                                   y,
                                   value,
                                   k = 1L,
                                   method = c("gaussian", "median"),
                                   sigma = 1) {
  method <- match.arg(method)
  out <- as.numeric(value)
  k <- max(0L, as.integer(k))
  if (k == 0L || length(out) == 0) {
    return(out)
  }
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(out) || length(y) != length(out)) {
    return(out)
  }

  for (iter in seq_len(k)) {
    prev <- out
    for (i in seq_along(prev)) {
      nb <- (abs(x - x[[i]]) + abs(y - y[[i]])) <= 1
      vals <- prev[nb]
      if (length(vals) > 0 && any(is.finite(vals))) {
        if (method == "median") {
          out[[i]] <- stats::median(vals, na.rm = TRUE)
        } else {
          dx <- x[nb] - x[[i]]
          dy <- y[nb] - y[[i]]
          dist2 <- dx * dx + dy * dy
          w <- exp(-dist2 / (2 * sigma * sigma + 1e-12))
          if (sum(w, na.rm = TRUE) > 0) {
            out[[i]] <- sum(vals * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
          } else {
            out[[i]] <- mean(vals, na.rm = TRUE)
          }
        }
      }
    }
  }
  out
}

#' Infer pseudo-MS2 by spatial correlation
#'
#' @param object A \code{msi_dataset} object.
#' @param variable_id Optional precursor variable_id(s). If NULL, all variables.
#' @param mz Optional precursor m/z value(s), mapped by ppm to variable_id.
#' @param tolerance_ppm Numeric ppm tolerance for precursor mapping.
#' @param correlation_threshold Numeric correlation threshold used when
#'   `adaptive_threshold = FALSE`.
#' @param correlation_method Correlation method passed to \code{stats::cor()}.
#'   One of `"pearson"` or `"spearman"`.
#' @param spatial_smoothing Logical. Apply simple neighborhood smoothing on
#'   precursor/fragment spatial signals before correlation.
#' @param smoothing_k Integer smoothing iterations when
#'   `spatial_smoothing = TRUE`.
#' @param smoothing_method Smoothing method. One of `"gaussian"` or `"median"`.
#' @param smoothing_sigma Sigma used for Gaussian smoothing.
#' @param fragment_intensity_threshold Numeric threshold between 0 and 1 applied to
#'   normalized pseudo-fragment intensity.
#' @param precursor_exclusion_tol Numeric m/z tolerance used to exclude
#'   precursor self-matching from candidate fragments.
#' @param fragment_window Numeric m/z window below precursor used for candidate
#'   selection. Candidates satisfy:
#'   `fragment_mz > precursor_mz - fragment_window`.
#' @param adaptive_threshold Logical. Use adaptive threshold:
#'   `max(0.4, quantile(correlation_values, 0.9))`.
#' @param adaptive_quantile Quantile used when `adaptive_threshold = TRUE`.
#' @param fallback_top_n Integer fallback count. If no fragment passes the
#'   threshold, select top-N scored fragments as `"pseudo_candidate"`.
#' @param replace_existing_pseudo Logical, replace existing pseudo entries.
#' @return Updated \code{msi_dataset} with pseudo MS2 records.
#' @export
infer_pseudo_ms2 <- function(object,
                             variable_id = NULL,
                             mz = NULL,
                             tolerance_ppm = 10,
                             correlation_threshold = 0.8,
                             correlation_method = c("pearson", "spearman"),
                             spatial_smoothing = FALSE,
                             smoothing_k = 1,
                             smoothing_method = c("gaussian", "median"),
                             smoothing_sigma = 1,
                             fragment_intensity_threshold = 0,
                             precursor_exclusion_tol = 0.01,
                             fragment_window = 150,
                             adaptive_threshold = TRUE,
                             adaptive_quantile = 0.9,
                             fallback_top_n = 5,
                             replace_existing_pseudo = TRUE) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }
  correlation_method <- match.arg(correlation_method)
  smoothing_method <- match.arg(smoothing_method)

  vinfo <- object@variable_info
  expr <- as.matrix(object@expression_data)
  if (nrow(expr) != nrow(vinfo)) {
    stop("expression_data and variable_info are not aligned.")
  }

  if (!is.null(mz)) {
    mapped <- vapply(as.numeric(mz), function(pmz) {
      m <- mz_match_ppm(vinfo$mz, pmz, tolerance_ppm = tolerance_ppm)
      if (!any(m)) {
        d <- abs(vinfo$mz - pmz)
        vinfo$variable_id[[which.min(d)]]
      } else {
        idx <- which(m)
        d <- abs(vinfo$mz[idx] - pmz)
        vinfo$variable_id[idx[[which.min(d)]]]
      }
    }, character(1))
    variable_id <- unique(c(variable_id, mapped))
  }

  if (is.null(variable_id)) {
    variable_id <- as.character(vinfo$variable_id)
  }

  precursor_idx <- match(variable_id, vinfo$variable_id)
  precursor_idx <- precursor_idx[!is.na(precursor_idx)]
  if (length(precursor_idx) == 0) {
    warning("No valid precursor variable_id found for pseudo MS2 inference.")
    return(object)
  }

  all_table <- .get_ms2_table(object)
  if (isTRUE(replace_existing_pseudo) && nrow(all_table) > 0) {
    keep <- !(all_table$ms2_type %in% c("pseudo", "pseudo_candidate") &
      all_table$precursor_variable_id %in% vinfo$variable_id[precursor_idx])
    all_table <- all_table[keep, , drop = FALSE]
  }

  pseudo_rows <- list()
  corr_rows <- list()
  for (i in precursor_idx) {
    precursor_id <- as.character(vinfo$variable_id[[i]])
    precursor_mz <- as.numeric(vinfo$mz[[i]])
    precursor_img <- as.numeric(expr[i, ])
    if (isTRUE(spatial_smoothing)) {
      precursor_img <- .smooth_spatial_signal(
        x = object@sample_info$x,
        y = object@sample_info$y,
        value = precursor_img,
        k = smoothing_k,
        method = smoothing_method,
        sigma = smoothing_sigma
      )
    }

    cand_idx <- which(
      vinfo$mz < precursor_mz &
        vinfo$mz > (precursor_mz - as.numeric(fragment_window)) &
        abs(vinfo$mz - precursor_mz) > as.numeric(precursor_exclusion_tol)
    )
    if (length(cand_idx) == 0) {
      next
    }

    cors <- vapply(cand_idx, function(j) {
      frag_img <- as.numeric(expr[j, ])
      if (isTRUE(spatial_smoothing)) {
        frag_img <- .smooth_spatial_signal(
          x = object@sample_info$x,
          y = object@sample_info$y,
          value = frag_img,
          k = smoothing_k,
          method = smoothing_method,
          sigma = smoothing_sigma
        )
      }
      suppressWarnings(stats::cor(
        precursor_img,
        frag_img,
        method = correlation_method,
        use = "pairwise.complete.obs"
      ))
    }, numeric(1))

    frag_mean_int <- rowMeans(expr[cand_idx, , drop = FALSE], na.rm = TRUE)
    frag_norm_int <- as.numeric(frag_mean_int)
    frag_norm_int[!is.finite(frag_norm_int)] <- 0
    if (length(frag_norm_int) > 0 && max(frag_norm_int, na.rm = TRUE) > 0) {
      frag_norm_int <- frag_norm_int / max(frag_norm_int, na.rm = TRUE)
    } else {
      frag_norm_int <- rep(0, length(frag_norm_int))
    }

    scores <- pmax(0, cors) * frag_norm_int
    corr_tbl <- data.frame(
      precursor_variable_id = precursor_id,
      precursor_mz = precursor_mz,
      fragment_mz = as.numeric(vinfo$mz[cand_idx]),
      variable_id = as.character(vinfo$variable_id[cand_idx]),
      correlation = as.numeric(cors),
      normalized_fragment_intensity = as.numeric(frag_norm_int),
      score = as.numeric(scores),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    corr_tbl <- corr_tbl[is.finite(corr_tbl$correlation), , drop = FALSE]
    if (nrow(corr_tbl) == 0) {
      next
    }

    corr_thr <- as.numeric(correlation_threshold)
    if (isTRUE(adaptive_threshold)) {
      corr_thr <- max(0.4, as.numeric(stats::quantile(
        corr_tbl$correlation,
        probs = adaptive_quantile,
        na.rm = TRUE,
        names = FALSE
      )))
    }

    keep <- corr_tbl$correlation > corr_thr &
      corr_tbl$normalized_fragment_intensity >= as.numeric(fragment_intensity_threshold)
    if (!any(keep)) {
      n_pick <- min(max(1L, as.integer(fallback_top_n)), nrow(corr_tbl))
      ord <- order(corr_tbl$score, decreasing = TRUE)
      keep <- rep(FALSE, nrow(corr_tbl))
      keep[ord[seq_len(n_pick)]] <- TRUE
      inferred_type <- "pseudo_candidate"
    } else {
      inferred_type <- "pseudo"
    }

    picked <- corr_tbl[keep, , drop = FALSE]
    picked <- picked[order(picked$score, decreasing = TRUE), , drop = FALSE]

    conf <- pmax(0, pmin(1, (picked$correlation - corr_thr) / (1 - corr_thr + 1e-12)))
    corr_rows[[length(corr_rows) + 1L]] <- corr_tbl[order(corr_tbl$score, decreasing = TRUE), , drop = FALSE]

    pseudo_rows[[length(pseudo_rows) + 1L]] <- data.frame(
      ms2_id = paste0(inferred_type, "_", precursor_id),
      precursor_mz = precursor_mz,
      precursor_variable_id = precursor_id,
      fragment_mz = picked$fragment_mz,
      fragment_intensity = picked$score,
      ms2_type = inferred_type,
      source_file = "pseudo_spatial_correlation",
      sample_id = NA_character_,
      rt = NA_real_,
      confidence_score = as.numeric(conf),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  pseudo_table <- do.call(rbind, pseudo_rows)
  if (!is.null(pseudo_table) && nrow(pseudo_table) > 0) {
    all_table <- rbind(all_table, pseudo_table)
  }

  object <- .set_ms2_table(object, all_table)
  if (length(corr_rows) > 0) {
    object@ms2_data$pseudo_correlation_table <- do.call(rbind, corr_rows)
  } else {
    object@ms2_data$pseudo_correlation_table <- data.frame(
      precursor_variable_id = character(),
      precursor_mz = numeric(),
      fragment_mz = numeric(),
      variable_id = character(),
      correlation = numeric(),
      normalized_fragment_intensity = numeric(),
      score = numeric(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  object@ms2_data$fragment_candidates <- object@ms2_data$pseudo_correlation_table
  object@ms2_data$pseudo_ms2_table <- all_table[
    all_table$ms2_type %in% c("pseudo", "pseudo_candidate"),
    ,
    drop = FALSE
  ]
  object <- .sync_variable_ms2(object)

  object <- add_process_log(
    object = object,
    function_name = "infer_pseudo_ms2()",
    parameter = list(
      n_precursor = length(precursor_idx),
      n_pseudo_ms2 = if (is.null(pseudo_table)) 0L else length(unique(pseudo_table$ms2_id)),
      n_pseudo_fragment = if (is.null(pseudo_table)) 0L else nrow(pseudo_table),
      tolerance_ppm = tolerance_ppm,
      correlation_threshold = correlation_threshold,
      correlation_method = correlation_method,
      spatial_smoothing = spatial_smoothing,
      smoothing_k = smoothing_k,
      smoothing_method = smoothing_method,
      smoothing_sigma = smoothing_sigma,
      fragment_intensity_threshold = fragment_intensity_threshold,
      precursor_exclusion_tol = precursor_exclusion_tol,
      fragment_window = fragment_window,
      adaptive_threshold = adaptive_threshold,
      adaptive_quantile = adaptive_quantile,
      fallback_top_n = fallback_top_n,
      replace_existing_pseudo = replace_existing_pseudo
    )
  )

  methods::validObject(object)
  object
}

#' Plot MS2 spectrum
#'
#' @param object A \code{msi_dataset} object.
#' @param variable_id Optional precursor variable_id.
#' @param mz Optional precursor m/z.
#' @param tolerance_ppm Numeric ppm tolerance when selecting by `mz`.
#' @param ms2_type Optional MS2 type filter (`"actual"`, `"pseudo"`, or
#'   `"pseudo_candidate"`).
#' @return A \code{ggplot2} object.
#' @export
plot_ms2 <- function(object, variable_id = NULL, mz = NULL, tolerance_ppm = 10, ms2_type = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2 for plot_ms2().")
  }

  spec <- extract_ms2(
    object = object,
    variable_id = variable_id,
    mz = mz,
    tolerance_ppm = tolerance_ppm,
    ms2_type = ms2_type,
    what = NULL
  )
  if (is.null(spec)) {
    stop("No MS2 spectrum found for the requested precursor.")
  }

  df <- data.frame(
    fragment_mz = as.numeric(spec$fragments),
    intensity = as.numeric(spec$intensity),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  subtitle <- paste0(
    "precursor_mz=", signif(spec$precursor_mz, 6),
    ", type=", spec$ms2_type,
    ", ms2_id=", spec$metadata$ms2_id
  )

  ggplot2::ggplot(df, ggplot2::aes(x = fragment_mz, y = intensity)) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = fragment_mz, yend = 0),
      linewidth = 0.4,
      color = "#1f77b4"
    ) +
    ggplot2::labs(
      title = "MS2 Spectrum",
      subtitle = subtitle,
      x = "fragment_mz",
      y = "normalized intensity"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

#' Plot spatial maps for MS2 fragments
#'
#' @param object A \code{msi_dataset} object.
#' @param variable_id Optional precursor variable_id.
#' @param mz Optional precursor m/z.
#' @param ms2_type Optional MS2 type filter (`"actual"`, `"pseudo"`, or
#'   `"pseudo_candidate"`).
#' @param section_id Optional section filter.
#' @param roi_id Optional ROI filter.
#' @param overlay_roi Logical. Overlay ROI pixels.
#' @param tolerance_ppm Numeric ppm tolerance for precursor and fragment mapping.
#' @param max_fragments Optional maximum number of fragments to visualize.
#' @param scales Facet scales passed to \code{plot_spatial()}.
#' @param section_scale One of `"global"` or `"section"`.
#' @return A \code{ggplot2} object.
#' @export
plot_spatial_ms2 <- function(object,
                             variable_id = NULL,
                             mz = NULL,
                             ms2_type = NULL,
                             section_id = NULL,
                             roi_id = NULL,
                             overlay_roi = FALSE,
                             tolerance_ppm = 20,
                             max_fragments = NULL,
                             scales = "fixed",
                             section_scale = c("global", "section")) {
  section_scale <- match.arg(section_scale)
  spec <- extract_ms2(
    object = object,
    variable_id = variable_id,
    mz = mz,
    tolerance_ppm = tolerance_ppm,
    ms2_type = ms2_type,
    what = NULL
  )
  if (is.null(spec) || length(spec$fragments) == 0) {
    stop("No MS2 fragments available for spatial plotting.")
  }

  frag_df <- data.frame(
    fragment_mz = as.numeric(spec$fragments),
    intensity = as.numeric(spec$intensity),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  frag_df <- frag_df[order(frag_df$intensity, decreasing = TRUE), , drop = FALSE]
  if (!is.null(max_fragments)) {
    max_fragments <- max(1L, as.integer(max_fragments))
    frag_df <- utils::head(frag_df, max_fragments)
  }

  p <- plot_spatial(
    object = object,
    feature = frag_df$fragment_mz,
    feature_type = "mz",
    section_id = section_id,
    roi_id = roi_id,
    overlay_roi = overlay_roi,
    facet_by = c("feature"),
    scales = scales,
    aggregation = "none",
    section_scale = section_scale,
    tolerance_ppm = tolerance_ppm
  ) +
    ggplot2::labs(
      title = "MS2 Fragment Spatial Maps",
      subtitle = paste0(
        "precursor_mz=", signif(spec$precursor_mz, 6),
        ", ms2_type=", spec$ms2_type
      )
    )
  p
}

#' Build precursor-fragment MS2 network
#'
#' @param object A \code{msi_dataset} object.
#' @param variable_id Optional precursor variable_id.
#' @param mz Optional precursor m/z.
#' @param ms2_type Optional MS2 type filter (`"actual"`, `"pseudo"`, or
#'   `"pseudo_candidate"`).
#' @param edge_weight One of `"spatial_correlation"` or `"fragment_intensity"`.
#' @param correlation_method Correlation method (`"pearson"` or `"spearman"`).
#' @param tolerance_ppm Numeric ppm tolerance used to map fragment m/z to MSI variables.
#' @return An \code{igraph} object.
#' @export
build_ms2_network <- function(object,
                              variable_id = NULL,
                              mz = NULL,
                              ms2_type = NULL,
                              edge_weight = c("spatial_correlation", "fragment_intensity"),
                              correlation_method = c("pearson", "spearman"),
                              tolerance_ppm = 20) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Please install igraph to use build_ms2_network().")
  }
  edge_weight <- match.arg(edge_weight)
  correlation_method <- match.arg(correlation_method)

  spec <- extract_ms2(
    object = object,
    variable_id = variable_id,
    mz = mz,
    tolerance_ppm = tolerance_ppm,
    ms2_type = ms2_type,
    what = NULL
  )
  if (is.null(spec) || length(spec$fragments) == 0) {
    stop("No MS2 fragments available for network construction.")
  }

  vinfo <- object@variable_info
  expr <- as.matrix(object@expression_data)
  precursor_id <- as.character(spec$metadata$precursor_variable_id)
  precursor_idx <- match(precursor_id, vinfo$variable_id)
  if (is.na(precursor_idx)) {
    stop("Precursor variable_id not found in object.")
  }

  frag_map_idx <- vapply(spec$fragments, function(fmz) {
    matched <- which(mz_match_ppm(vinfo$mz, fmz, tolerance_ppm = tolerance_ppm))
    if (length(matched) == 0) {
      which.min(abs(vinfo$mz - fmz))
    } else {
      dd <- abs(vinfo$mz[matched] - fmz)
      matched[[which.min(dd)]]
    }
  }, integer(1))

  edge_df <- data.frame(
    from = precursor_id,
    to = as.character(vinfo$variable_id[frag_map_idx]),
    fragment_mz = as.numeric(spec$fragments),
    fragment_intensity = as.numeric(spec$intensity),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  edge_df <- edge_df[!duplicated(edge_df$to), , drop = FALSE]
  if (nrow(edge_df) == 0) {
    stop("No mapped fragment nodes available for network construction.")
  }

  edge_df$weight <- if (edge_weight == "fragment_intensity") {
    edge_df$fragment_intensity
  } else {
    precursor_vec <- as.numeric(expr[precursor_idx, ])
    frag_idx <- match(edge_df$to, vinfo$variable_id)
    vapply(frag_idx, function(i) {
      suppressWarnings(stats::cor(
        precursor_vec,
        as.numeric(expr[i, ]),
        method = correlation_method,
        use = "pairwise.complete.obs"
      ))
    }, numeric(1))
  }
  edge_df$weight[!is.finite(edge_df$weight)] <- 0

  node_ids <- unique(c(edge_df$from, edge_df$to))
  node_df <- data.frame(
    name = node_ids,
    type = ifelse(node_ids == precursor_id, "precursor", "fragment"),
    mz = vinfo$mz[match(node_ids, vinfo$variable_id)],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  g <- igraph::graph_from_data_frame(
    d = edge_df[, c("from", "to", "weight", "fragment_mz", "fragment_intensity"), drop = FALSE],
    vertices = node_df,
    directed = TRUE
  )
  igraph::graph_attr(g, "precursor_id") <- precursor_id
  igraph::graph_attr(g, "ms2_type") <- spec$ms2_type
  igraph::graph_attr(g, "edge_weight_mode") <- edge_weight
  g
}

#' Plot MS2 fragment network
#'
#' @param graph An \code{igraph} object from \code{build_ms2_network()}.
#' @param layout Graph layout string passed to \code{ggraph::create_layout()}.
#' @return A \code{ggplot2} object.
#' @export
plot_ms2_network <- function(graph, layout = "fr") {
  if (!inherits(graph, "igraph")) {
    stop("graph must be an igraph object.")
  }
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop("Please install ggraph to use plot_ms2_network().")
  }

  ggraph::ggraph(graph, layout = layout) +
    ggraph::geom_edge_link(
      ggplot2::aes(width = weight),
      alpha = 0.7,
      colour = "#636363",
      end_cap = ggraph::circle(2, "mm"),
      arrow = grid::arrow(length = grid::unit(2, "mm"))
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(color = type),
      size = 4
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = signif(mz, 6)),
      repel = TRUE,
      size = 3
    ) +
    ggplot2::scale_color_manual(values = c(precursor = "#d7301f", fragment = "#2b8cbe")) +
    ggplot2::scale_edge_width(range = c(0.2, 1.5)) +
    ggplot2::theme_void() +
    ggplot2::labs(
      title = "MS2 Fragment Network",
      subtitle = paste0(
        "precursor=", igraph::graph_attr(graph, "precursor_id"),
        ", weight=", igraph::graph_attr(graph, "edge_weight_mode")
      ),
      color = "node type",
      edge_width = "weight"
    )
}

.extract_ms2_by_id <- function(object, ms2_id) {
  ms2_table <- .get_ms2_table(object)
  sub <- ms2_table[ms2_table$ms2_id %in% as.character(ms2_id), , drop = FALSE]
  if (nrow(sub) == 0) {
    return(NULL)
  }
  sub[order(sub$fragment_mz), , drop = FALSE]
}

.empty_pseudo_fragment_table <- function() {
  data.frame(
    pseudo_ms2_id = character(),
    precursor_variable_id = character(),
    precursor_mz = numeric(),
    fragment_mz = numeric(),
    fragment_variable_id = character(),
    correlation = numeric(),
    fragment_intensity = numeric(),
    score = numeric(),
    normalized_intensity = numeric(),
    slice_consistency = numeric(),
    ms2_type = character(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.section_correlation <- function(precursor_vec,
                                 fragment_vec,
                                 section_id,
                                 correlation_method) {
  sections <- unique(as.character(section_id))
  sections <- sections[nzchar(sections) & !is.na(sections)]
  if (length(sections) == 0) {
    sections <- "section_1"
    section_id <- rep("section_1", length(precursor_vec))
  }

  cors <- vapply(sections, function(sec) {
    idx <- which(as.character(section_id) %in% sec)
    if (length(idx) < 2L) {
      return(NA_real_)
    }
    suppressWarnings(stats::cor(
      as.numeric(precursor_vec[idx]),
      as.numeric(fragment_vec[idx]),
      method = correlation_method,
      use = "pairwise.complete.obs"
    ))
  }, numeric(1))
  cors[is.nan(cors)] <- NA_real_
  cors
}

.target_precursor_indices <- function(vinfo, precursor_mz, ppm) {
  unique(unlist(lapply(as.numeric(precursor_mz), function(pmz) {
    matched <- which(mz_match_ppm(vinfo$mz, pmz, tolerance_ppm = ppm))
    if (length(matched) == 0) {
      return(integer())
    }
    matched[[which.min(abs(vinfo$mz[matched] - pmz))]]
  })))
}

.joint_fragment_candidates_template <- function() {
  data.frame(
    precursor_variable_id = character(),
    precursor_mz = numeric(),
    fragment_variable_id = character(),
    fragment_mz = numeric(),
    spatial_corr = numeric(),
    spectral_score = numeric(),
    intensity_score = numeric(),
    joint_score = numeric(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.neutral_loss_score <- function(delta_mz,
                                losses = c(18, 28, 44),
                                sigma = 5) {
  delta_mz <- as.numeric(delta_mz)
  if (!is.finite(delta_mz)) {
    return(0)
  }
  max(exp(-abs(delta_mz - losses) / sigma), na.rm = TRUE)
}

#' Generate targeted pseudo-MS2 spectra
#'
#' @param object A \code{msi_dataset} object.
#' @param precursor_mz Numeric precursor m/z target vector.
#' @param ppm Numeric precursor matching tolerance.
#' @param correlation_method Spatial correlation method.
#' @param correlation_threshold Minimum final correlation.
#' @param top_n Maximum fragment count per precursor.
#' @return Updated \code{msi_dataset} with targeted pseudo-MS2 results stored in
#'   \code{ms2_data$pseudo_ms2_table}.
#' @export
generate_pseudo_ms2 <- function(object,
                                precursor_mz,
                                ppm = 10,
                                correlation_method = c("spearman", "pearson"),
                                correlation_threshold = 0.3,
                                top_n = 20) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }
  correlation_method <- match.arg(correlation_method)
  precursor_mz <- as.numeric(precursor_mz)
  precursor_mz <- precursor_mz[is.finite(precursor_mz)]
  if (length(precursor_mz) == 0) {
    stop("precursor_mz must contain at least one finite numeric value.")
  }

  vinfo <- object@variable_info
  expr <- as.matrix(object@expression_data)
  if (nrow(expr) != nrow(vinfo)) {
    stop("expression_data and variable_info are not aligned.")
  }

  precursor_idx <- .target_precursor_indices(vinfo, precursor_mz = precursor_mz, ppm = ppm)
  if (length(precursor_idx) == 0) {
    warning("No precursor features matched the requested precursor_mz values.")
    object@ms2_data$pseudo_ms2_table <- .empty_pseudo_fragment_table()
    object@ms2_data$targeted_pseudo_ms2_table <- .empty_pseudo_fragment_table()
    return(object)
  }

  ms2_table <- .get_ms2_table(object)
  if (nrow(ms2_table) > 0) {
    keep <- !(ms2_table$ms2_type %in% c("pseudo", "pseudo_candidate") &
      ms2_table$precursor_variable_id %in% as.character(vinfo$variable_id[precursor_idx]))
    ms2_table <- ms2_table[keep, , drop = FALSE]
  }

  section_id <- if ("section_id" %in% colnames(object@sample_info)) {
    as.character(object@sample_info$section_id)
  } else {
    rep("section_1", ncol(expr))
  }

  pseudo_rows <- list()
  ms2_rows <- list()
  for (i in precursor_idx) {
    precursor_id <- as.character(vinfo$variable_id[[i]])
    pmz <- as.numeric(vinfo$mz[[i]])
    precursor_vec <- as.numeric(expr[i, ])
    exclusion_tol <- max(pmz * as.numeric(ppm) / 1e6, 1e-8)
    cand_idx <- which(
      as.numeric(vinfo$mz) < pmz &
        abs(as.numeric(vinfo$mz) - pmz) > exclusion_tol
    )
    if (length(cand_idx) == 0) {
      next
    }

    cand_tbl <- lapply(cand_idx, function(j) {
      frag_vec <- as.numeric(expr[j, ])
      corr_slice <- .section_correlation(
        precursor_vec = precursor_vec,
        fragment_vec = frag_vec,
        section_id = section_id,
        correlation_method = correlation_method
      )
      corr_valid <- corr_slice[is.finite(corr_slice)]
      if (length(corr_valid) == 0) {
        return(NULL)
      }
      corr_final <- mean(corr_valid, na.rm = TRUE)
      mean_corr <- mean(corr_valid, na.rm = TRUE)
      slice_consistency <- if (!is.finite(mean_corr) || mean_corr <= 0) {
        0
      } else {
        1 - stats::sd(corr_valid, na.rm = TRUE) / (mean_corr + 1e-12)
      }
      if (!is.finite(slice_consistency)) {
        slice_consistency <- 0
      }
      slice_consistency <- max(0, min(1, slice_consistency))

      frag_int <- mean(frag_vec[is.finite(frag_vec)], na.rm = TRUE)
      if (!is.finite(frag_int)) {
        frag_int <- 0
      }
      score <- as.numeric(corr_final) * log1p(as.numeric(frag_int)) * slice_consistency
      data.frame(
        pseudo_ms2_id = paste0("pseudo_", precursor_id),
        precursor_variable_id = precursor_id,
        precursor_mz = pmz,
        fragment_mz = as.numeric(vinfo$mz[[j]]),
        fragment_variable_id = as.character(vinfo$variable_id[[j]]),
        correlation = as.numeric(corr_final),
        fragment_intensity = as.numeric(frag_int),
        score = as.numeric(score),
        normalized_intensity = NA_real_,
        slice_consistency = as.numeric(slice_consistency),
        ms2_type = "pseudo",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    })

    cand_tbl <- do.call(rbind, cand_tbl)
    if (is.null(cand_tbl) || nrow(cand_tbl) == 0) {
      next
    }
    cand_tbl <- cand_tbl[is.finite(cand_tbl$correlation), , drop = FALSE]
    cand_tbl <- cand_tbl[cand_tbl$correlation >= as.numeric(correlation_threshold), , drop = FALSE]
    if (nrow(cand_tbl) == 0) {
      next
    }
    cand_tbl <- cand_tbl[order(cand_tbl$score, decreasing = TRUE), , drop = FALSE]
    cand_tbl <- utils::head(cand_tbl, max(1L, as.integer(top_n)))
    max_score <- max(cand_tbl$score, na.rm = TRUE)
    cand_tbl$normalized_intensity <- if (is.finite(max_score) && max_score > 0) {
      cand_tbl$score / max_score
    } else {
      rep(0, nrow(cand_tbl))
    }

    pseudo_rows[[length(pseudo_rows) + 1L]] <- cand_tbl
    ms2_rows[[length(ms2_rows) + 1L]] <- data.frame(
      ms2_id = cand_tbl$pseudo_ms2_id,
      precursor_mz = cand_tbl$precursor_mz,
      precursor_variable_id = cand_tbl$precursor_variable_id,
      fragment_mz = cand_tbl$fragment_mz,
      fragment_intensity = cand_tbl$normalized_intensity,
      ms2_type = cand_tbl$ms2_type,
      source_file = "pseudo_targeted_spatial_correlation",
      sample_id = NA_character_,
      rt = NA_real_,
      confidence_score = cand_tbl$correlation,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  pseudo_tbl <- do.call(rbind, pseudo_rows)
  if (is.null(pseudo_tbl)) {
    pseudo_tbl <- .empty_pseudo_fragment_table()
  }
  compat_tbl <- do.call(rbind, ms2_rows)
  if (is.null(compat_tbl)) {
    compat_tbl <- .empty_ms2_table()
  }

  if (nrow(compat_tbl) > 0) {
    ms2_table <- rbind(ms2_table, compat_tbl)
  }
  object <- .set_ms2_table(object, ms2_table)
  object@ms2_data$pseudo_ms2_table <- pseudo_tbl
  object@ms2_data$targeted_pseudo_ms2_table <- pseudo_tbl
  object@ms2_data$pseudo_correlation_table <- pseudo_tbl[, c(
    "precursor_variable_id", "precursor_mz", "fragment_mz",
    "fragment_variable_id", "correlation", "fragment_intensity", "score"
  ), drop = FALSE]
  colnames(object@ms2_data$pseudo_correlation_table)[4] <- "variable_id"
  object@ms2_data$fragment_candidates <- object@ms2_data$pseudo_correlation_table
  object <- .sync_variable_ms2(object)
  object <- add_process_log(
    object = object,
    function_name = "generate_pseudo_ms2()",
    parameter = list(
      precursor_mz = precursor_mz,
      ppm = ppm,
      correlation_method = correlation_method,
      correlation_threshold = correlation_threshold,
      top_n = top_n,
      n_precursor_matched = length(unique(pseudo_tbl$precursor_variable_id)),
      n_fragment_rows = nrow(pseudo_tbl)
    )
  )
  methods::validObject(object)
  object
}

#' Detect spatial-spectral joint fragment candidates
#'
#' @param object A \code{msi_dataset} object.
#' @param precursor_mz Numeric precursor m/z target vector.
#' @param ppm Numeric precursor matching tolerance.
#' @param correlation_method Spatial correlation method.
#' @return Updated \code{msi_dataset} with joint fragment candidates stored in
#'   \code{ms2_data$fragment_candidates}.
#' @export
detect_joint_fragments <- function(object,
                                   precursor_mz,
                                   ppm = 10,
                                   correlation_method = c("spearman", "pearson")) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }
  correlation_method <- match.arg(correlation_method)
  precursor_mz <- as.numeric(precursor_mz)
  precursor_mz <- precursor_mz[is.finite(precursor_mz)]
  if (length(precursor_mz) == 0) {
    stop("precursor_mz must contain at least one finite numeric value.")
  }

  vinfo <- object@variable_info
  expr <- as.matrix(object@expression_data)
  precursor_idx <- .target_precursor_indices(vinfo, precursor_mz = precursor_mz, ppm = ppm)
  if (length(precursor_idx) == 0) {
    warning("No precursor features matched the requested precursor_mz values.")
    object@ms2_data$fragment_candidates <- .joint_fragment_candidates_template()
    return(object)
  }

  candidate_rows <- lapply(precursor_idx, function(i) {
    pmz <- as.numeric(vinfo$mz[[i]])
    precursor_id <- as.character(vinfo$variable_id[[i]])
    precursor_vec <- as.numeric(expr[i, ])
    exclusion_tol <- max(pmz * as.numeric(ppm) / 1e6, 1e-8)
    cand_idx <- which(
      as.numeric(vinfo$mz) < pmz &
        abs(as.numeric(vinfo$mz) - pmz) > exclusion_tol
    )
    if (length(cand_idx) == 0) {
      return(NULL)
    }

    do.call(rbind, lapply(cand_idx, function(j) {
      fragment_vec <- as.numeric(expr[j, ])
      spatial_corr <- suppressWarnings(stats::cor(
        precursor_vec,
        fragment_vec,
        method = correlation_method,
        use = "pairwise.complete.obs"
      ))
      if (!is.finite(spatial_corr)) {
        return(NULL)
      }

      frag_int <- mean(fragment_vec[is.finite(fragment_vec)], na.rm = TRUE)
      if (!is.finite(frag_int)) {
        frag_int <- 0
      }
      delta_mz <- pmz - as.numeric(vinfo$mz[[j]])
      spectral_score <- .neutral_loss_score(delta_mz)
      intensity_score <- log1p(frag_int)
      joint_score <- as.numeric(spatial_corr) * spectral_score * intensity_score

      data.frame(
        precursor_variable_id = precursor_id,
        precursor_mz = pmz,
        fragment_variable_id = as.character(vinfo$variable_id[[j]]),
        fragment_mz = as.numeric(vinfo$mz[[j]]),
        spatial_corr = as.numeric(spatial_corr),
        spectral_score = as.numeric(spectral_score),
        intensity_score = as.numeric(intensity_score),
        joint_score = as.numeric(joint_score),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }))
  })

  candidate_tbl <- do.call(rbind, candidate_rows)
  if (is.null(candidate_tbl)) {
    candidate_tbl <- .joint_fragment_candidates_template()
  } else {
    candidate_tbl <- candidate_tbl[order(candidate_tbl$joint_score, decreasing = TRUE), , drop = FALSE]
  }

  object@ms2_data$fragment_candidates <- candidate_tbl
  object <- add_process_log(
    object = object,
    function_name = "detect_joint_fragments()",
    parameter = list(
      precursor_mz = precursor_mz,
      ppm = ppm,
      correlation_method = correlation_method,
      n_candidate_rows = nrow(candidate_tbl)
    )
  )
  methods::validObject(object)
  object
}

#' Build graph-based pseudo fragmentation network
#'
#' @param object A \code{msi_dataset} object.
#' @param top_k Integer number of top edges kept per precursor.
#' @return Updated \code{msi_dataset} with an \code{igraph} stored in
#'   \code{other_files$pseudo_fragment_network}.
#' @export
build_fragment_network <- function(object, top_k = 10) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Please install igraph to use build_fragment_network().")
  }

  candidate_tbl <- object@ms2_data$fragment_candidates
  if (is.null(candidate_tbl) || nrow(as.data.frame(candidate_tbl)) == 0) {
    stop("No fragment_candidates found. Run detect_joint_fragments() first.")
  }
  candidate_tbl <- as.data.frame(candidate_tbl, stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("precursor_mz", "fragment_mz", "joint_score")
  if (!all(required_cols %in% colnames(candidate_tbl))) {
    stop("fragment_candidates must contain precursor_mz, fragment_mz, and joint_score columns.")
  }

  top_k <- max(1L, as.integer(top_k))
  split_tbl <- split(candidate_tbl, candidate_tbl$precursor_mz)
  edge_df <- do.call(rbind, lapply(split_tbl, function(df) {
    df <- df[order(df$joint_score, decreasing = TRUE), , drop = FALSE]
    utils::head(df, top_k)
  }))
  if (is.null(edge_df) || nrow(edge_df) == 0) {
    stop("No edges available after top_k pruning.")
  }

  edge_df <- data.frame(
    from = paste0("P_", signif(edge_df$precursor_mz, 8)),
    to = paste0("F_", if ("fragment_variable_id" %in% colnames(edge_df)) {
      as.character(edge_df$fragment_variable_id)
    } else {
      signif(edge_df$fragment_mz, 8)
    }),
    weight = as.numeric(edge_df$joint_score),
    precursor_mz = as.numeric(edge_df$precursor_mz),
    fragment_mz = as.numeric(edge_df$fragment_mz),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  node_df <- unique(rbind(
    data.frame(
      name = unique(edge_df$from),
      type = "precursor",
      mz = suppressWarnings(as.numeric(sub("^P_", "", unique(edge_df$from)))),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    data.frame(
      name = unique(edge_df$to),
      type = "fragment",
      mz = vapply(unique(edge_df$to), function(node) {
        row_idx <- match(node, edge_df$to)
        edge_df$fragment_mz[[row_idx]]
      }, numeric(1)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  ))

  graph <- igraph::graph_from_data_frame(
    d = edge_df,
    vertices = node_df,
    directed = TRUE
  )
  igraph::graph_attr(graph, "top_k") <- top_k
  igraph::graph_attr(graph, "graph_type") <- "pseudo_fragment_network"

  object@other_files$pseudo_fragment_network <- graph
  object <- add_process_log(
    object = object,
    function_name = "build_fragment_network()",
    parameter = list(
      top_k = top_k,
      n_edges = nrow(edge_df),
      n_nodes = nrow(node_df)
    )
  )
  methods::validObject(object)
  object
}

#' Plot graph-based pseudo fragmentation network
#'
#' @param object An \code{msi_dataset} object or \code{igraph} graph.
#' @param layout Graph layout string passed to \code{ggraph::create_layout()}.
#' @return A \code{ggplot2} object.
#' @export
plot_fragment_network <- function(object, layout = "fr") {
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop("Please install ggraph to use plot_fragment_network().")
  }

  graph <- if (inherits(object, "igraph")) {
    object
  } else {
    if (!methods::is(object, "msi_dataset")) {
      stop("object must be an msi_dataset or igraph object.")
    }
    object@other_files$pseudo_fragment_network
  }
  if (is.null(graph) || !inherits(graph, "igraph")) {
    stop("No pseudo fragment network found. Run build_fragment_network() first.")
  }

  ggraph::ggraph(graph, layout = layout) +
    ggraph::geom_edge_link(
      ggplot2::aes(width = weight),
      alpha = 0.65,
      colour = "#7f8c8d",
      end_cap = ggraph::circle(2, "mm"),
      arrow = grid::arrow(length = grid::unit(2, "mm"))
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(color = type, size = type)
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = signif(mz, 6)),
      repel = TRUE,
      size = 3
    ) +
    ggplot2::scale_color_manual(values = c(precursor = "#d7301f", fragment = "#3182bd")) +
    ggplot2::scale_size_manual(values = c(precursor = 4, fragment = 3)) +
    ggplot2::scale_edge_width(range = c(0.3, 1.6)) +
    ggplot2::theme_void() +
    ggplot2::labs(
      title = "Pseudo Fragmentation Network",
      subtitle = paste0("top_k=", igraph::graph_attr(graph, "top_k")),
      color = "node type",
      size = "node type",
      edge_width = "joint score"
    )
}

.match_fragments_tol <- function(real_df, pseudo_df, tolerance = 0.01) {
  if (nrow(real_df) == 0 || nrow(pseudo_df) == 0) {
    return(data.frame(
      mz_real = numeric(),
      mz_pseudo = numeric(),
      intensity_real = numeric(),
      intensity_pseudo = numeric(),
      delta_mz = numeric(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  used_pseudo <- rep(FALSE, nrow(pseudo_df))
  matched <- vector("list", nrow(real_df))
  k <- 0L
  for (i in seq_len(nrow(real_df))) {
    deltas <- abs(pseudo_df$fragment_mz - real_df$fragment_mz[[i]])
    idx <- which(deltas <= as.numeric(tolerance) & !used_pseudo)
    if (length(idx) == 0) {
      next
    }
    j <- idx[[which.min(deltas[idx])]]
    used_pseudo[[j]] <- TRUE
    k <- k + 1L
    matched[[k]] <- data.frame(
      mz_real = real_df$fragment_mz[[i]],
      mz_pseudo = pseudo_df$fragment_mz[[j]],
      intensity_real = real_df$fragment_intensity[[i]],
      intensity_pseudo = pseudo_df$fragment_intensity[[j]],
      delta_mz = abs(real_df$fragment_mz[[i]] - pseudo_df$fragment_mz[[j]]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  out <- do.call(rbind, matched[seq_len(k)])
  if (is.null(out)) {
    out <- data.frame(
      mz_real = numeric(),
      mz_pseudo = numeric(),
      intensity_real = numeric(),
      intensity_pseudo = numeric(),
      delta_mz = numeric(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  out
}

#' Compare real and pseudo MS2 spectra
#'
#' @param object A \code{msi_dataset} object.
#' @param ms2_id Real MS2 spectrum ID.
#' @param pseudo_ms2_id Pseudo MS2 spectrum ID.
#' @param tolerance Numeric fragment m/z matching tolerance.
#' @return A list with matched fragments and similarity metrics.
#' @export
compare_spectra <- function(object,
                            ms2_id,
                            pseudo_ms2_id,
                            tolerance = 0.01) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }

  real_df <- .extract_ms2_by_id(object, ms2_id)
  pseudo_df <- .extract_ms2_by_id(object, pseudo_ms2_id)
  if (is.null(real_df)) {
    stop("ms2_id not found: ", ms2_id)
  }
  if (is.null(pseudo_df)) {
    stop("pseudo_ms2_id not found: ", pseudo_ms2_id)
  }

  matched <- .match_fragments_tol(real_df, pseudo_df, tolerance = tolerance)
  matched_n <- nrow(matched)
  union_n <- length(unique(c(real_df$fragment_mz, pseudo_df$fragment_mz)))
  overlap_ratio <- if (union_n > 0) matched_n / union_n else 0

  real_vec <- rep(0, nrow(real_df))
  pseudo_vec <- rep(0, nrow(real_df))
  if (matched_n > 0) {
    for (i in seq_len(matched_n)) {
      r_idx <- which.min(abs(real_df$fragment_mz - matched$mz_real[[i]]))
      real_vec[[r_idx]] <- matched$intensity_real[[i]]
      pseudo_vec[[r_idx]] <- matched$intensity_pseudo[[i]]
    }
  }
  denom <- sqrt(sum(real_vec^2) * sum(pseudo_vec^2))
  cosine_similarity <- if (is.finite(denom) && denom > 0) sum(real_vec * pseudo_vec) / denom else 0

  list(
    matched_fragments = matched,
    cosine_similarity = as.numeric(cosine_similarity),
    overlap_ratio = as.numeric(overlap_ratio),
    matched_fragment_count = as.integer(matched_n)
  )
}

#' Mirror plot of real vs pseudo MS2
#'
#' @param object A \code{msi_dataset} object.
#' @param ms2_id Real MS2 spectrum ID.
#' @param pseudo_ms2_id Pseudo MS2 spectrum ID.
#' @return A \code{ggplot2} object.
#' @export
mirror_plot <- function(object, ms2_id, pseudo_ms2_id) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2 for mirror_plot().")
  }
  real_df <- .extract_ms2_by_id(object, ms2_id)
  pseudo_df <- .extract_ms2_by_id(object, pseudo_ms2_id)
  if (is.null(real_df) || is.null(pseudo_df)) {
    stop("Both ms2_id and pseudo_ms2_id must exist in object.")
  }

  real_plot <- data.frame(
    fragment_mz = real_df$fragment_mz,
    intensity = real_df$fragment_intensity / (max(real_df$fragment_intensity, na.rm = TRUE) + 1e-12),
    spectrum = "real",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  pseudo_plot <- data.frame(
    fragment_mz = pseudo_df$fragment_mz,
    intensity = -pseudo_df$fragment_intensity / (max(pseudo_df$fragment_intensity, na.rm = TRUE) + 1e-12),
    spectrum = "pseudo",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  plot_df <- rbind(real_plot, pseudo_plot)

  ggplot2::ggplot(plot_df, ggplot2::aes(x = fragment_mz, y = intensity, color = spectrum)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
    ggplot2::geom_segment(ggplot2::aes(xend = fragment_mz, yend = 0), linewidth = 0.4) +
    ggplot2::scale_color_manual(values = c(real = "#d7301f", pseudo = "#2b8cbe")) +
    ggplot2::labs(
      title = "Mirror Plot: Real vs Pseudo MS2",
      subtitle = paste0("real=", ms2_id, " | pseudo=", pseudo_ms2_id),
      x = "fragment m/z",
      y = "normalized intensity"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

#' Benchmark pseudo-MS2 reconstruction
#'
#' @param object A \code{msi_dataset} object.
#' @param tolerance Fragment m/z matching tolerance for spectrum comparison.
#' @param ... Additional arguments passed to \code{infer_pseudo_ms2()}.
#' @return Updated \code{msi_dataset} with benchmark table stored in
#'   `other_files$pseudo_ms2_benchmark`.
#' @export
benchmark_pseudo_ms2 <- function(object, tolerance = 0.01, ...) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }

  ms2_table <- .get_ms2_table(object)
  real_tbl <- ms2_table[ms2_table$ms2_type == "actual", , drop = FALSE]
  if (nrow(real_tbl) == 0) {
    warning("No real MS2 spectra found for benchmarking.")
    object@other_files$pseudo_ms2_benchmark <- data.frame()
    return(object)
  }

  precursors <- unique(real_tbl$precursor_variable_id)
  object <- infer_pseudo_ms2(object = object, variable_id = precursors, ...)
  ms2_table <- .get_ms2_table(object)

  bench_rows <- lapply(precursors, function(pid) {
    real_sub <- ms2_table[
      ms2_table$precursor_variable_id == pid & ms2_table$ms2_type == "actual",
      ,
      drop = FALSE
    ]
    pseudo_sub <- ms2_table[
      ms2_table$precursor_variable_id == pid & ms2_table$ms2_type %in% c("pseudo", "pseudo_candidate"),
      ,
      drop = FALSE
    ]
    if (nrow(real_sub) == 0) {
      return(NULL)
    }

    real_id <- unique(real_sub$ms2_id)[[1]]
    pseudo_id <- if (nrow(pseudo_sub) > 0) unique(pseudo_sub$ms2_id)[[1]] else NA_character_
    cmp <- if (!is.na(pseudo_id)) {
      compare_spectra(object, ms2_id = real_id, pseudo_ms2_id = pseudo_id, tolerance = tolerance)
    } else {
      NULL
    }

    data.frame(
      precursor_variable_id = pid,
      precursor_mz = unique(real_sub$precursor_mz)[[1]],
      real_fragment_count = length(unique(real_sub$fragment_mz)),
      pseudo_fragment_count = if (nrow(pseudo_sub) > 0) length(unique(pseudo_sub$fragment_mz)) else 0L,
      matched_fragment_count = if (!is.null(cmp)) cmp$matched_fragment_count else 0L,
      cosine_similarity = if (!is.null(cmp)) cmp$cosine_similarity else NA_real_,
      fragment_overlap_ratio = if (!is.null(cmp)) cmp$overlap_ratio else NA_real_,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  bench_tbl <- do.call(rbind, bench_rows)
  if (is.null(bench_tbl)) {
    bench_tbl <- data.frame(
      precursor_variable_id = character(),
      precursor_mz = numeric(),
      real_fragment_count = integer(),
      pseudo_fragment_count = integer(),
      matched_fragment_count = integer(),
      cosine_similarity = numeric(),
      fragment_overlap_ratio = numeric(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  object@other_files$pseudo_ms2_benchmark <- bench_tbl
  object <- add_process_log(
    object = object,
    function_name = "benchmark_pseudo_ms2()",
    parameter = list(
      n_precursor = nrow(bench_tbl),
      tolerance = tolerance
    )
  )
  object
}

#' Plot pseudo MS2 spectrum
#'
#' @inheritParams plot_ms2
#' @return A \code{ggplot2} object.
#' @export
plot_pseudo_ms2 <- function(object, variable_id = NULL, mz = NULL, tolerance_ppm = 10) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2 for plot_pseudo_ms2().")
  }

  pseudo_tbl <- object@ms2_data$pseudo_ms2_table
  if (!is.null(pseudo_tbl) && nrow(as.data.frame(pseudo_tbl)) > 0 &&
      all(c("precursor_mz", "fragment_mz", "normalized_intensity") %in% colnames(pseudo_tbl))) {
    df <- as.data.frame(pseudo_tbl, stringsAsFactors = FALSE, check.names = FALSE)
    if (!is.null(mz) && is.null(variable_id)) {
      idx <- which(mz_match_ppm(df$precursor_mz, mz, tolerance_ppm = tolerance_ppm))
      if (length(idx) > 0) {
        df <- df[idx, , drop = FALSE]
      }
    }
    if (!is.null(variable_id) && "precursor_variable_id" %in% colnames(df)) {
      df <- df[df$precursor_variable_id %in% as.character(variable_id), , drop = FALSE]
    }
    if (nrow(df) > 0) {
      return(
        ggplot2::ggplot(df, ggplot2::aes(x = fragment_mz, y = normalized_intensity)) +
          ggplot2::geom_segment(
            ggplot2::aes(xend = fragment_mz, yend = 0),
            linewidth = 0.45,
            color = "#d95f02"
          ) +
          ggplot2::labs(
            title = "Pseudo MS2 Spectrum",
            subtitle = paste0("precursor_mz=", paste(unique(signif(df$precursor_mz, 6)), collapse = ", ")),
            x = "fragment m/z",
            y = "normalized intensity"
          ) +
          ggplot2::theme_minimal(base_size = 12)
      )
    }
  }

  spec <- extract_ms2(object, variable_id = variable_id, mz = mz, tolerance_ppm = tolerance_ppm, ms2_type = "pseudo")
  if (is.null(spec)) {
    spec <- extract_ms2(object, variable_id = variable_id, mz = mz, tolerance_ppm = tolerance_ppm, ms2_type = "pseudo_candidate")
  }
  if (is.null(spec)) {
    stop("No pseudo/pseudo_candidate spectrum found for requested precursor.")
  }
  plot_ms2(
    object = object,
    variable_id = spec$metadata$precursor_variable_id,
    tolerance_ppm = tolerance_ppm,
    ms2_type = spec$ms2_type
  )
}

#' Plot spatial precursor-fragment comparison
#'
#' @param object A \code{msi_dataset} object.
#' @param precursor_mz Numeric precursor m/z.
#' @param fragment_mz Optional fragment m/z. If omitted, the top-ranked targeted
#'   pseudo fragment is used.
#' @param ppm Numeric m/z matching tolerance.
#' @param section_id Optional section filter.
#' @return A \code{ggplot2} object.
#' @export
plot_spatial_fragment <- function(object,
                                  precursor_mz,
                                  fragment_mz = NULL,
                                  ppm = 10,
                                  section_id = NULL) {
  pseudo_tbl <- object@ms2_data$pseudo_ms2_table
  if (is.null(pseudo_tbl) || nrow(as.data.frame(pseudo_tbl)) == 0) {
    stop("No pseudo MS2 table found. Run generate_pseudo_ms2() first.")
  }
  pseudo_tbl <- as.data.frame(pseudo_tbl, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("precursor_mz", "fragment_mz") %in% colnames(pseudo_tbl))) {
    stop("pseudo_ms2_table does not contain targeted pseudo-MS2 records.")
  }

  sub <- pseudo_tbl[mz_match_ppm(pseudo_tbl$precursor_mz, precursor_mz, tolerance_ppm = ppm), , drop = FALSE]
  if (nrow(sub) == 0) {
    stop("Requested precursor_mz was not found in pseudo_ms2_table.")
  }
  if (is.null(fragment_mz)) {
    rank_col <- if ("score" %in% colnames(sub)) "score" else "normalized_intensity"
    sub <- sub[order(sub[[rank_col]], decreasing = TRUE), , drop = FALSE]
    fragment_mz <- sub$fragment_mz[[1]]
  }

  plot_spatial(
    object = object,
    feature = c(as.numeric(precursor_mz), as.numeric(fragment_mz)),
    feature_type = "mz",
    section_id = section_id,
    facet_by = c("feature"),
    aggregation = "none",
    section_scale = "section",
    tolerance_ppm = ppm
  ) +
    ggplot2::labs(
      title = "Spatial Precursor vs Fragment",
      subtitle = paste0(
        "precursor_mz=", signif(precursor_mz, 6),
        ", fragment_mz=", signif(fragment_mz, 6)
      )
    )
}

#' Plot pseudo MS2 network
#'
#' @inheritParams build_ms2_network
#' @return A \code{ggplot2} object.
#' @export
plot_pseudo_network <- function(object,
                                variable_id = NULL,
                                mz = NULL,
                                edge_weight = c("spatial_correlation", "fragment_intensity"),
                                correlation_method = c("pearson", "spearman"),
                                tolerance_ppm = 20) {
  g <- build_ms2_network(
    object = object,
    variable_id = variable_id,
    mz = mz,
    ms2_type = c("pseudo", "pseudo_candidate"),
    edge_weight = edge_weight,
    correlation_method = correlation_method,
    tolerance_ppm = tolerance_ppm
  )
  plot_ms2_network(g)
}

#' Filter MS2 entries in msi_dataset
#'
#' @param object A \code{msi_dataset} object.
#' @param ... Filtering expressions evaluated on MS2 info.
#' @param .preserve Passed to \code{dplyr::filter()}.
#' @return Updated \code{msi_dataset}.
#' @importFrom rlang quos
#' @export
filter_ms2 <- function(object, ...) {
  UseMethod("filter_ms2")
}

#' @export
filter_ms2.msi_dataset <- function(object, ..., .preserve = FALSE) {
  dots <- rlang::quos(...)
  flat <- .flatten_ms2_data(object)
  if (nrow(flat$info) == 0) {
    return(object)
  }

  filtered <- dplyr::filter(flat$info, !!!dots, .preserve = .preserve)
  keep_id <- unique(filtered$ms2_spectrum_id)

  ms2_table <- flat$table
  if (length(keep_id) == 0) {
    ms2_table <- .empty_ms2_table()
  } else {
    ms2_table <- ms2_table[ms2_table$ms2_id %in% keep_id, , drop = FALSE]
  }

  object <- .set_ms2_table(object, ms2_table)
  object <- .sync_variable_ms2(object)

  object <- add_process_log(
    object = object,
    function_name = "filter_ms2()",
    parameter = list(n_remain = length(unique(ms2_table$ms2_id)))
  )
  object
}
