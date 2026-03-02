# ============================================================================
# Missing value imputation for msi_dataset
# ============================================================================

#' Impute missing values in msi_dataset expression_data
#'
#' @param object A \code{msi_dataset} object.
#' @param method Imputation method: median, mean, half_min, zero.
#' @param by Impute by variable (row-wise) or sample (column-wise).
#'
#' @return Updated \code{msi_dataset} object.
#' @export
impute_mv <- function(object,
                      method = c("median", "mean", "half_min", "zero"),
                      by = c("variable", "sample")) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }

  method <- match.arg(method)
  by <- match.arg(by)

  x <- as.matrix(object@expression_data)
  na_count_before <- sum(!is.finite(x))
  x[!is.finite(x)] <- NA_real_

  fill_fun <- switch(
    method,
    median = function(v) stats::median(v, na.rm = TRUE),
    mean = function(v) mean(v, na.rm = TRUE),
    half_min = function(v) 0.5 * min(v, na.rm = TRUE),
    zero = function(v) 0
  )

  if (by == "variable") {
    for (i in seq_len(nrow(x))) {
      idx <- is.na(x[i, ])
      if (any(idx)) {
        v <- x[i, !idx]
        if (length(v) == 0) v <- 0
        x[i, idx] <- fill_fun(v)
      }
    }
  } else {
    for (j in seq_len(ncol(x))) {
      idx <- is.na(x[, j])
      if (any(idx)) {
        v <- x[!idx, j]
        if (length(v) == 0) v <- 0
        x[idx, j] <- fill_fun(v)
      }
    }
  }

  dimnames(x) <- list(object@variable_info$variable_id, object@sample_info$sample_id)
  object@expression_data <- x

  object <- add_process_log(
    object = object,
    function_name = "impute_mv()",
    parameter = list(
      method = method,
      by = by,
      na_before = na_count_before,
      na_after = sum(!is.finite(as.matrix(object@expression_data)))
    )
  )

  object
}
