# ============================================================================
# Processing log helpers for msi_dataset
# ============================================================================

#' Append processing information to msi_dataset
#'
#' @param object A \code{msi_dataset} object.
#' @param function_name Function name to record.
#' @param parameter Named list of parameters.
#'
#' @return Updated \code{msi_dataset} object.
#' @export
add_process_log <- function(object, function_name, parameter = list()) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }

  key <- paste0("step_", length(object@process_info) + 1, "_",
                gsub("[^A-Za-z0-9_]", "_", function_name))

  parameter_object <- methods::new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = function_name,
    parameter = parameter,
    time = Sys.time()
  )

  object@process_info[[key]] <- parameter_object
  object
}
