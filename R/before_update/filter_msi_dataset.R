# ============================================================================
# dplyr::filter() integration for msi_dataset
# Author: Xiaoyu Xu <XIAOYU021@e.ntu.edu.sg>
# ============================================================================

#' Filter a msi_dataset Object
#'
#' @description
#' Applies \code{dplyr::filter()} to one component of a \code{msi_dataset}.
#'
#' \itemize{
#'   \item \code{target = "sample_info"} — filters pixels (columns of
#'     \code{expression_data}).
#'   \item \code{target = "variable_info"} — filters m/z features (rows of
#'     \code{expression_data}).
#' }
#'
#' @param .data A \code{msi_dataset} object.
#' @param ... Filter expressions passed to \code{dplyr::filter()}.
#' @param .by Passed to \code{dplyr::filter()}.
#' @param .preserve Passed to \code{dplyr::filter()}.
#' @param target Which component to filter.  One of \code{"sample_info"}
#'   (default) or \code{"variable_info"}.  Falls back to
#'   \code{.data@@activated} when \code{target = NULL}.
#'
#' @return A filtered \code{msi_dataset} object with updated
#'   \code{process_info}.
#' @export
#' @method filter msi_dataset
#'
#' @examples
#' \dontrun{
#' # Filter pixels by spatial coordinates
#' obj2 <- dplyr::filter(obj, x > 10, y > 20)
#'
#' # Filter variables by m/z range
#' obj3 <- dplyr::filter(obj, target = "variable_info", mz > 500, mz < 900)
#' }

filter.msi_dataset <- function(.data,
                                ...,
                                .by       = NULL,
                                .preserve = FALSE,
                                target    = NULL) {

  if (!methods::is(.data, "msi_dataset"))
    stop(".data must be a msi_dataset object.")

  ## ---- determine target ----
  valid_targets <- c("sample_info", "variable_info")
  if (is.null(target)) {
    target <- if (.data@activated %in% valid_targets) .data@activated else "sample_info"
  }
  target <- match.arg(target, valid_targets)

  ## ---- capture filter expressions for process_info ----
  filter_exprs      <- as.list(substitute(list(...)))[-1L]
  filter_expr_text  <- if (length(filter_exprs) > 0) {
    vapply(filter_exprs, function(x) paste(deparse(x), collapse = " "),
           character(1))
  } else {
    character(0)
  }

  ## ---- apply filter ----
  if (target == "sample_info") {

    new_sinfo <- dplyr::filter(.data@sample_info, ...,
                                .by = .by, .preserve = .preserve)
    keep_ids  <- new_sinfo$sample_id

    if (length(keep_ids) == 0)
      warning("filter.msi_dataset: no pixels remain after filtering.")

    bad <- keep_ids[!keep_ids %in% colnames(.data@expression_data)]
    if (length(bad) > 0)
      stop("Filtered sample_id not found in expression_data: ",
           paste(bad, collapse = ", "))

    new_expr <- .data@expression_data[, keep_ids, drop = FALSE]

    .data@sample_info    <- new_sinfo
    .data@expression_data <- new_expr

  } else {  # variable_info

    new_vinfo <- dplyr::filter(.data@variable_info, ...,
                                .by = .by, .preserve = .preserve)
    keep_ids  <- new_vinfo$variable_id

    if (length(keep_ids) == 0)
      warning("filter.msi_dataset: no variables remain after filtering.")

    bad <- keep_ids[!keep_ids %in% rownames(.data@expression_data)]
    if (length(bad) > 0)
      stop("Filtered variable_id not found in expression_data: ",
           paste(bad, collapse = ", "))

    new_expr <- .data@expression_data[keep_ids, , drop = FALSE]

    .data@variable_info   <- new_vinfo
    .data@expression_data <- new_expr
  }

  ## ---- record in process_info ----
  process_key <- paste0("filter_", target, "_",
                        format(Sys.time(), "%Y%m%d_%H%M%S"))
  parameter <- new(
    Class         = "tidymass_parameter",
    pacakge_name  = "msidataset",
    function_name = "dplyr::filter()",
    parameter     = list(target      = target,
                         filter_expr = filter_expr_text),
    time          = Sys.time()
  )
  .data@process_info[[process_key]] <- parameter

  .data
}


# ---------------------------------------------------------------------------
# filter_pixels() — convenient wrapper
# ---------------------------------------------------------------------------

#' Filter pixels (samples) in a msi_dataset
#'
#' @description
#' A named wrapper around \code{filter.msi_dataset} that always targets
#' \code{sample_info}.
#'
#' @param object A \code{msi_dataset} object.
#' @param ... Filter expressions (as in \code{dplyr::filter}).
#' @return A filtered \code{msi_dataset} object.
#' @export

filter_pixels <- function(object, ...) {
  if (!methods::is(object, "msi_dataset"))
    stop("object must be a msi_dataset.")
  filter.msi_dataset(object, ..., target = "sample_info")
}


# ---------------------------------------------------------------------------
# filter_variables() — convenient wrapper
# ---------------------------------------------------------------------------

#' Filter variables (m/z features) in a msi_dataset
#'
#' @description
#' A named wrapper around \code{filter.msi_dataset} that always targets
#' \code{variable_info}.
#'
#' @param object A \code{msi_dataset} object.
#' @param ... Filter expressions (as in \code{dplyr::filter}).
#' @return A filtered \code{msi_dataset} object.
#' @export

filter_variables <- function(object, ...) {
  if (!methods::is(object, "msi_dataset"))
    stop("object must be a msi_dataset.")
  filter.msi_dataset(object, ..., target = "variable_info")
}
