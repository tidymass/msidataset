#' Filter a msi_dataset Object
#'
#' @description
#' Apply \code{dplyr::filter()} to one component of a \code{msi_dataset} object.
#' Filtering \code{pixel_info} will subset the columns of \code{expression_data},
#' and filtering \code{variable_info} will subset the rows of
#' \code{expression_data}.
#'
#' @param .data An object of class \code{msi_dataset}.
#' @param ... Filter expressions passed to \code{dplyr::filter()}.
#' @param .by Passed to \code{dplyr::filter()}.
#' @param .preserve Passed to \code{dplyr::filter()}.
#' @param target Which component to filter: \code{"pixel_info"},
#'   \code{"variable_info"}, or \code{"sample_info"}. If \code{NULL}, the
#'   function uses \code{.data@activated} when it matches one of those values;
#'   otherwise it defaults to \code{"pixel_info"}.
#'
#' @return A filtered \code{msi_dataset} object.
#' @export
#' @method filter msi_dataset
#'
#' @examples
#' \dontrun{
#' data("expression_data")
#' data("pixel_info")
#' data("variable_info")
#' data("sample_info")
#' object <-
#'   create_msi_dataset(
#'     expression_data = expression_data,
#'     pixel_info = pixel_info,
#'     variable_info = variable_info,
#'     sample_info = sample_info
#'   )
#'
#' # Filter pixels by spatial coordinates
#' object2 <- dplyr::filter(object, x > 10, y > 10)
#'
#' # Filter variables by m/z
#' object3 <- dplyr::filter(object, target = "variable_info", mz > 500)
#' }
filter.msi_dataset <- function(.data,
                               ...,
                               .by = NULL,
                               .preserve = FALSE,
                               target = NULL) {
  if (!methods::is(.data, "msi_dataset")) {
    stop(".data must be a msi_dataset.")
  }

  if (is.null(target)) {
    if (.data@activated %in% c("pixel_info", "variable_info", "sample_info")) {
      target <- .data@activated
    } else {
      target <- "pixel_info"
    }
  }

  target <- match.arg(target, c("pixel_info", "variable_info", "sample_info"))

  filter_exprs <- as.list(substitute(list(...)))[-1L]
  filter_expr_text <- if (length(filter_exprs) > 0) {
    vapply(filter_exprs, function(x) paste(deparse(x), collapse = " "), character(1))
  } else {
    character(0)
  }

  if (target == "pixel_info") {
    pixel_info_filtered <- dplyr::filter(.data@pixel_info, ..., .by = .by, .preserve = .preserve)
    pixel_ids <- pixel_info_filtered$pixel_id

    if (length(pixel_ids) > 0 && anyNA(match(pixel_ids, colnames(.data@expression_data)))) {
      stop("Filtered pixel_id not found in expression_data column names.")
    }

    expression_data_filtered <- .data@expression_data[, pixel_ids, drop = FALSE]

    .data@pixel_info <- pixel_info_filtered
    .data@expression_data <- expression_data_filtered
  } else if (target == "variable_info") {
    variable_info_filtered <- dplyr::filter(.data@variable_info, ..., .by = .by, .preserve = .preserve)
    variable_ids <- variable_info_filtered$variable_id

    if (length(variable_ids) > 0 && anyNA(match(variable_ids, rownames(.data@expression_data)))) {
      stop("Filtered variable_id not found in expression_data row names.")
    }

    expression_data_filtered <- .data@expression_data[variable_ids, , drop = FALSE]

    .data@variable_info <- variable_info_filtered
    .data@expression_data <- expression_data_filtered
  } else {
    sample_info_filtered <- dplyr::filter(.data@sample_info, ..., .by = .by, .preserve = .preserve)
    .data@sample_info <- sample_info_filtered
  }

  process_key <- paste0("filter_msi_dataset_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "dplyr::filter()",
    parameter = list(
      target = target,
      filter_expr = filter_expr_text
    ),
    time = Sys.time()
  )

  .data@process_info[[process_key]] <- parameter

  .data
}
