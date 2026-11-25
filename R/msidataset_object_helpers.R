# Helpers for preserving a massdataset and its step results.  These functions
# are deliberately internal until the object-facing API is merged.
.msidataset_feature_mz <- function(variables) {
  .msi_resolve_observed_mz(data.frame(
    observed_mz = if ("mz" %in% names(variables)) variables$mz else rep(NA_real_, nrow(variables)),
    observed_mz_full = if ("observed_mz_full" %in% names(variables))
      variables$observed_mz_full else rep(NA_real_, nrow(variables))
  ))
}

.msidataset_validate <- function(object, require_roi = FALSE, require_mz = FALSE) {
  if (!methods::is(object, "mass_dataset")) stop("Expected a mass_dataset", call. = FALSE)
  e <- object@expression_data
  s <- object@sample_info
  v <- object@variable_info
  if (ncol(e) != nrow(s) || nrow(e) != nrow(v) ||
      !all(c("sample_id", "x", "y", "tic") %in% names(s)) ||
      !"variable_id" %in% names(v) || anyDuplicated(s$sample_id) ||
      anyDuplicated(v$variable_id) || anyNA(s$sample_id) || anyNA(v$variable_id) ||
      !identical(as.character(colnames(e)), as.character(s$sample_id)) ||
      !identical(as.character(rownames(e)), as.character(v$variable_id)))
    stop("Inconsistent mass_dataset dimensions or identifiers", call. = FALSE)
  if (require_roi && !"region" %in% names(s)) stop("ROI labels are required", call. = FALSE)
  if (require_mz && (!"mz" %in% names(v) || any(!is.finite(.msidataset_feature_mz(v)$mz))))
    stop("Finite feature m/z values are required", call. = FALSE)
  invisible(object)
}

.msidataset_append_step <- function(object, step, result, info = attr(result, "process_info")) {
  .msidataset_validate(object)
  if (!is.character(step) || length(step) != 1L || is.na(step) || !nzchar(step))
    stop("step must be one nonempty name", call. = FALSE)
  # Never overwrite earlier runs; both result and metadata travel with the object.
  n <- 1L
  key <- step
  while (key %in% names(object@process_info)) {
    n <- n + 1L
    key <- paste0(step, "_", n)
  }
  version <- tryCatch(as.character(utils::packageVersion("msidataset")),
                      error = function(e) NA_character_)
  object@process_info[[key]] <- list(function_name = step, result = result,
    input_summary = list(pixels = nrow(object@sample_info), features = nrow(object@variable_info)),
    parameters = info, package_version = version, time = as.character(Sys.time()))
  object
}

.msidataset_inputs <- function(object, require_roi = FALSE, require_mz = FALSE) {
  .msidataset_validate(object, require_roi, require_mz)
  list(xy = as.matrix(object@sample_info[, c("x", "y"), drop = FALSE]),
       tic = object@sample_info$tic,
       intensity = t(as.matrix(object@expression_data)),
       region = object@sample_info$region,
       mz = object@variable_info$mz,
       feature_id = as.character(object@variable_info$variable_id))
}

.msidataset_input_problems <- function(i, require_region = FALSE, require_mz = FALSE,
                                       integer_xy = FALSE) {
  problems <- list()
  bad_xy <- which(!is.finite(i$xy[, 1]) | !is.finite(i$xy[, 2]))
  if (integer_xy) bad_xy <- union(bad_xy, which(i$xy[, 1] != floor(i$xy[, 1]) |
                                                i$xy[, 2] != floor(i$xy[, 2])))
  if (length(bad_xy)) problems$invalid_coordinates <- bad_xy
  if (integer_xy && anyDuplicated(data.frame(i$xy)))
    problems$duplicate_coordinates <- which(duplicated(data.frame(i$xy)) |
                                             duplicated(data.frame(i$xy), fromLast = TRUE))
  bad_tic <- which(!is.finite(i$tic) | i$tic <= 0)
  if (length(bad_tic)) problems$invalid_tic <- bad_tic
  bad_intensity <- which(!is.finite(i$intensity), arr.ind = TRUE)
  if (nrow(bad_intensity)) problems$invalid_intensity <- bad_intensity
  if (require_region && anyNA(i$region)) problems$missing_region <- which(is.na(i$region))
  if (require_mz && any(!is.finite(i$mz))) problems$invalid_mz <- which(!is.finite(i$mz))
  problems
}
