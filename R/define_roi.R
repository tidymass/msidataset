#' Attach pixel ROI labels
#'
#' Attach labels without changing or discarding pixels; coordinates must match
#' in the same order. Region labels are recorded as supplied.
#'
#' @param object A `mass_dataset` returned by [read_msi_data()].
#' @param labels Data frame with columns `x`, `y`, and `region` in pixel order.
#' @return The massdataset object with a region column and process history.
#' @examples
#' \dontrun{msi <- define_roi(msi, labels)}
#' @seealso [read_msi_data()]
#' @export
define_roi <- function(object, labels) {
  if (!methods::is(object, "mass_dataset") || !all(c("x", "y", "region") %in% names(labels)) ||
      !all(c("x", "y") %in% names(object@sample_info)) || nrow(labels) != nrow(object@sample_info) ||
      !isTRUE(all(labels$x == object@sample_info$x & labels$y == object@sample_info$y)))
    stop("ROI coordinates do not match the data", call. = FALSE)
  previous_region <- if ("region" %in% names(object@sample_info)) object@sample_info$region else NULL
  object@sample_info$region <- as.character(labels$region)
  note <- object@sample_info_note
  note <- note[note$name != "region", , drop = FALSE]
  object@sample_info_note <- rbind(note, data.frame(name = "region", meaning = "ROI region label"))
  object <- .msidataset_append_step(object, "define_roi",
    list(previous_region = previous_region),
    list(n_pixels = nrow(labels), time = as.character(Sys.time())))
  object
}
