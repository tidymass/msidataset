# ============================================================================
# Core user-facing MSI verbs
# ============================================================================

#' Filter an msi_dataset using unified verb interface
#'
#' @param object A \code{msi_dataset} object.
#' @param target One of \code{"sample_info"} or \code{"variable_info"}.
#' @param filter_type Filtering strategy.
#' @param ... Additional arguments passed to \code{filter_msi_dataset()}.
#' @return Filtered \code{msi_dataset}.
#' @export
filter_msi <- function(object, target = "sample_info", filter_type = NULL, ...) {
  # Step 1: validate input
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }
  # Step 2: compute mask / filter strategy
  # Step 3: subset expression_data
  # Step 4: update pixel_info (sample_info)
  # Step 5: run validObject()
  before_dim <- dim(object@expression_data)
  out <- filter_msi_dataset(object = object, target = target, filter_type = filter_type, ...)
  after_dim <- dim(out@expression_data)
  out <- add_process_log(
    object = out,
    function_name = "filter_msi()",
    parameter = list(
      target = target,
      filter_type = filter_type,
      n_variables_before = before_dim[[1]],
      n_samples_before = before_dim[[2]],
      n_variables_after = after_dim[[1]],
      n_samples_after = after_dim[[2]]
    )
  )
  methods::validObject(out)
  out
}

#' Mutate an msi_dataset using unified verb interface
#'
#' @param object A \code{msi_dataset} object.
#' @param method Currently only \code{"ms2"} is supported.
#' @param ... Additional arguments passed to mutate backend.
#' @return Mutated \code{msi_dataset}.
#' @export
mutate_msi <- function(object, method = c("ms2"), ...) {
  # Step 1: validate input
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }
  method <- match.arg(method)
  before_dim <- dim(object@expression_data)
  # Step 2: compute mask / mutation target
  out <- switch(
    method,
    ms2 = mutate_ms2(object = object, ...)
  )
  after_dim <- dim(out@expression_data)
  out <- add_process_log(
    object = out,
    function_name = "mutate_msi()",
    parameter = list(
      method = method,
      n_variables_before = before_dim[[1]],
      n_samples_before = before_dim[[2]],
      n_variables_after = after_dim[[1]],
      n_samples_after = after_dim[[2]]
    )
  )
  # Step 3: subset expression_data (not required for ms2 mutate)
  # Step 4: update pixel_info (sample_info) if needed
  # Step 5: run validObject()
  methods::validObject(out)
  out
}

#' ROI operations via unified verb interface
#'
#' @param object A \code{msi_dataset} object.
#' @param action One of \code{"add"}, \code{"list"}, \code{"subset"},
#'   \code{"remove"}, \code{"rename"}, \code{"data"}, \code{"validate"}.
#' @param ... Additional arguments passed to ROI backend functions.
#' @return Action-dependent ROI result.
#' @export
roi_msi <- function(object,
                    action = c("add", "list", "subset", "remove", "rename", "data", "validate"),
                    ...) {
  # Step 1: validate input
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }
  action <- match.arg(action)
  before_dim <- dim(object@expression_data)
  # Step 2: compute mask / ROI operation
  out <- switch(
    action,
    add = add_roi(object = object, ...),
    list = list_roi(object = object, ...),
    subset = subset_by_roi(object = object, ...),
    remove = remove_roi(object = object, ...),
    rename = rename_roi(object = object, ...),
    data = roi_data(object = object, ...),
    validate = validate_roi_integrity(object = object, ...)
  )
  # Step 3: subset expression_data (for subset action)
  # Step 4: update pixel_info (sample_info) (for subset action)
  # Step 5: run validObject()
  if (methods::is(out, "msi_dataset")) {
    after_dim <- dim(out@expression_data)
    out <- add_process_log(
      object = out,
      function_name = "roi_msi()",
      parameter = list(
        action = action,
        n_variables_before = before_dim[[1]],
        n_samples_before = before_dim[[2]],
        n_variables_after = after_dim[[1]],
        n_samples_after = after_dim[[2]]
      )
    )
    methods::validObject(out)
  }
  out
}

#' Slice/subset an msi_dataset
#'
#' @param object A \code{msi_dataset} object.
#' @param i Variable indices/ids.
#' @param j Sample indices/ids.
#' @param section_id Optional section filter.
#' @return Subsetted \code{msi_dataset}.
#' @export
slice_msi <- function(object, i = NULL, j = NULL, section_id = NULL) {
  # Step 1: validate input
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }
  before_dim <- dim(object@expression_data)
  out <- object
  # Step 2: compute mask
  if (!is.null(section_id)) {
    out <- filter_msi_dataset(
      out,
      target = "sample_info",
      filter_type = "section",
      section_ids = section_id
    )
    if (ncol(out@expression_data) == 0 || nrow(out@sample_info) == 0) {
      stop("No pixels remaining after section filter.")
    }
  }
  # Step 3: subset expression_data
  if (!is.null(i) || !is.null(j)) {
    ii <- if (is.null(i)) seq_len(nrow(out@expression_data)) else i
    jj <- if (is.null(j)) seq_len(ncol(out@expression_data)) else j
    out <- out[ii, jj, drop = FALSE]
  }
  # Step 4: update pixel_info
  out <- update_msi_dataset(out)
  after_dim <- dim(out@expression_data)
  out <- add_process_log(
    object = out,
    function_name = "slice_msi()",
    parameter = list(
      section_id = section_id,
      n_variables_before = before_dim[[1]],
      n_samples_before = before_dim[[2]],
      n_variables_after = after_dim[[1]],
      n_samples_after = after_dim[[2]]
    )
  )
  # Step 5: run validObject()
  methods::validObject(out)
  out
}

#' Summarize key properties of an msi_dataset
#'
#' @param object A \code{msi_dataset} object.
#' @return A named list summary.
#' @export
summarize_msi <- function(object) {
  # Step 1: validate input
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }
  # Step 2: compute summary metrics
  before_dim <- dim(object@expression_data)
  mz_range <- if (nrow(object@variable_info) == 0 ||
                  !("mz" %in% colnames(object@variable_info))) {
    c(NA_real_, NA_real_)
  } else {
    c(
      suppressWarnings(min(object@variable_info$mz, na.rm = TRUE)),
      suppressWarnings(max(object@variable_info$mz, na.rm = TRUE))
    )
  }
  summary_list <- list(
    n_variables = nrow(object@variable_info),
    n_samples = nrow(object@sample_info),
    n_sections = length(unique(object@sample_info$section_id)),
    mz_range = mz_range,
    n_roi = nrow(list_roi(object)),
    n_process_steps = length(object@process_info)
  )
  # Step 3: subset expression_data (not applicable)
  # Step 4: update pixel_info (not applicable)
  # Step 5: run validObject()
  object2 <- add_process_log(
    object = object,
    function_name = "summarize_msi()",
    parameter = list(
      n_variables_before = before_dim[[1]],
      n_samples_before = before_dim[[2]],
      n_variables_after = before_dim[[1]],
      n_samples_after = before_dim[[2]]
    )
  )
  methods::validObject(object2)
  attr(summary_list, "object") <- object2
  summary_list
}

.merge_ms2_tables <- function(objects,
                              mz_axis,
                              variable_ids,
                              sample_maps) {
  out <- lapply(seq_along(objects), function(i) {
    obj <- objects[[i]]
    tbl <- .get_ms2_table(obj)
    if (is.null(tbl) || nrow(tbl) == 0) {
      return(NULL)
    }

    precursor_idx <- vapply(tbl$precursor_mz, function(pmz) {
      diffs <- abs(as.numeric(mz_axis) - as.numeric(pmz))
      diffs[!is.finite(diffs)] <- Inf
      which.min(diffs)
    }, integer(1))
    tbl$precursor_variable_id <- as.character(variable_ids[precursor_idx])

    if (length(sample_maps[[i]]) > 0) {
      mapped_sample <- unname(sample_maps[[i]][tbl$sample_id])
      tbl$sample_id[!is.na(mapped_sample)] <- mapped_sample[!is.na(mapped_sample)]
    }
    tbl$ms2_id <- paste0("merge", i, "_", tbl$ms2_id)
    tbl
  })

  out <- do.call(rbind, out)
  if (is.null(out)) {
    .empty_ms2_table()
  } else {
    .ensure_ms2_table(out)
  }
}

#' Merge multiple msi_dataset objects
#'
#' @param ... \code{msi_dataset} objects or a single list of them.
#' @param bin_digits Integer m/z rounding precision used to align features.
#' @return A merged \code{msi_dataset}.
#' @export
merge_msi_datasets <- function(..., bin_digits = 4) {
  objects <- list(...)
  if (length(objects) == 1L && is.list(objects[[1]]) && !methods::is(objects[[1]], "msi_dataset")) {
    objects <- objects[[1]]
  }
  if (length(objects) < 2L) {
    stop("Please provide at least two msi_dataset objects to merge.")
  }
  if (!all(vapply(objects, methods::is, logical(1), "msi_dataset"))) {
    stop("All inputs must be msi_dataset objects.")
  }

  digits <- max(1L, as.integer(bin_digits))
  mz_rounded <- lapply(objects, function(obj) round(as.numeric(obj@variable_info$mz), digits))
  mz_axis <- sort(unique(unlist(mz_rounded, recursive = TRUE, use.names = FALSE)))
  if (length(mz_axis) == 0) {
    stop("No m/z values available for merging.")
  }

  sample_info_list <- vector("list", length(objects))
  expr_list <- vector("list", length(objects))
  sample_maps <- vector("list", length(objects))
  offset <- 0L

  for (i in seq_along(objects)) {
    obj <- objects[[i]]
    expr <- as.matrix(obj@expression_data)
    idx <- match(mz_rounded[[i]], mz_axis)
    aligned <- matrix(0, nrow = length(mz_axis), ncol = ncol(expr))
    if (anyDuplicated(idx) > 0) {
      agg <- rowsum(expr, group = idx, reorder = FALSE)
      row_idx <- suppressWarnings(as.integer(rownames(agg)))
      aligned[row_idx, ] <- agg
    } else {
      aligned[idx, ] <- expr
    }

    old_ids <- as.character(obj@sample_info$sample_id)
    new_ids <- paste0("spectrum_", offset + seq_len(length(old_ids)))
    sample_maps[[i]] <- stats::setNames(new_ids, old_ids)

    sinfo <- obj@sample_info
    sinfo$sample_id <- new_ids
    sinfo$section_id <- if ("section_id" %in% colnames(sinfo)) {
      paste0(as.character(sinfo$section_id), "_merge", i)
    } else {
      paste0("section_", i)
    }
    sinfo$z <- rep(i, nrow(sinfo))

    colnames(aligned) <- new_ids
    sample_info_list[[i]] <- sinfo
    expr_list[[i]] <- aligned
    offset <- offset + length(old_ids)
  }

  sample_info <- do.call(rbind, sample_info_list)
  variable_info <- data.frame(
    variable_id = generate_variable_id(mz_axis),
    mz = as.numeric(mz_axis),
    ms2 = rep("no", length(mz_axis)),
    has_ms2 = rep("no", length(mz_axis)),
    ms2_id = rep(NA_character_, length(mz_axis)),
    ms2_type = rep(NA_character_, length(mz_axis)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  expression_data <- do.call(cbind, expr_list)
  rownames(expression_data) <- variable_info$variable_id
  colnames(expression_data) <- sample_info$sample_id

  out <- create_msi_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

  merged_ms2 <- .merge_ms2_tables(
    objects = objects,
    mz_axis = mz_axis,
    variable_ids = variable_info$variable_id,
    sample_maps = sample_maps
  )
  out <- .set_ms2_table(out, merged_ms2)
  if (nrow(merged_ms2) > 0) {
    out@ms2_data$pseudo_ms2_table <- merged_ms2[
      merged_ms2$ms2_type %in% c("pseudo", "pseudo_candidate"),
      ,
      drop = FALSE
    ]
  }
  out <- .sync_variable_ms2(out)
  out <- add_process_log(
    object = out,
    function_name = "merge_msi_datasets()",
    parameter = list(
      n_input = length(objects),
      bin_digits = digits,
      n_variables = nrow(variable_info),
      n_samples = nrow(sample_info)
    )
  )
  methods::validObject(out)
  out
}
