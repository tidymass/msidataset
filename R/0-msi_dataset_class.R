# ============================================================================
# msi_dataset Class Definition and Constructor
# Author: Xiaoyu Xu <XIAOYU021@e.ntu.edu.sg>
# ============================================================================

# ---------------------------------------------------------------------------
# Internal pre-check helper (called BEFORE object construction)
# ---------------------------------------------------------------------------

#' @keywords internal
check_msi_dataset <- function(expression_data,
                               sample_info,
                               variable_info,
                               sample_info_note,
                               variable_info_note) {

  errors <- character()

  if (!("variable_id" %in% colnames(variable_info)))
    errors <- c(errors, "variable_info must have variable_id.")

  if (sum(duplicated(variable_info$variable_id)) > 0)
    errors <- c(errors, "variable_id has duplicated items.")

  if (!("sample_id" %in% colnames(sample_info)))
    errors <- c(errors, "sample_info must have sample_id.")

  if (sum(duplicated(sample_info$sample_id)) > 0)
    errors <- c(errors, "sample_id has duplicated items.")

  if (!("x" %in% colnames(sample_info)))
    errors <- c(errors, "sample_info must have x.")

  if (!("y" %in% colnames(sample_info)))
    errors <- c(errors, "sample_info must have y.")

  if (ncol(expression_data) != nrow(sample_info))
    errors <- c(errors,
                "expression_data column number must equal sample_info row number.")

  if (nrow(expression_data) != nrow(variable_info))
    errors <- c(errors,
                "expression_data row number must equal variable_info row number.")

  if (length(errors) > 0) {
    return(paste("error:", paste(errors, collapse = "; ")))
  }

  return("pass")
}


# ---------------------------------------------------------------------------
# S4 validity checker (called by validObject / new())
# ---------------------------------------------------------------------------

#' Check the Validity of a msi_dataset Object
#'
#' @param object An object of class \code{msi_dataset}.
#' @return \code{TRUE} if valid, otherwise a character vector of error messages.
#' @export
check_msi_dataset_class <- function(object) {

  errors <- character()

  ## variable_info
  if (!("variable_id" %in% colnames(object@variable_info)))
    errors <- c(errors, "variable_info must have variable_id.")

  if (sum(duplicated(object@variable_info$variable_id)) > 0)
    errors <- c(errors, "variable_id has duplicated items.")

  ## sample_info
  if (!("sample_id" %in% colnames(object@sample_info)))
    errors <- c(errors, "sample_info must have sample_id.")

  if (sum(duplicated(object@sample_info$sample_id)) > 0)
    errors <- c(errors, "sample_id has duplicated items.")

  ## MSI spatial columns
  if (!("x" %in% colnames(object@sample_info)))
    errors <- c(errors, "sample_info must have x.")

  if (!("y" %in% colnames(object@sample_info)))
    errors <- c(errors, "sample_info must have y.")

  ## note tables
  if (nrow(object@sample_info_note) > 0) {
    if (!("name" %in% colnames(object@sample_info_note)))
      errors <- c(errors, "sample_info_note must have column: name.")
    if (!("meaning" %in% colnames(object@sample_info_note)))
      errors <- c(errors, "sample_info_note must have column: meaning.")
  }

  if (nrow(object@variable_info_note) > 0) {
    if (!("name" %in% colnames(object@variable_info_note)))
      errors <- c(errors, "variable_info_note must have column: name.")
    if (!("meaning" %in% colnames(object@variable_info_note)))
      errors <- c(errors, "variable_info_note must have column: meaning.")
  }

  ## expression vs sample
  if (ncol(object@expression_data) != nrow(object@sample_info)) {
    errors <- c(errors,
                "expression_data column number must equal sample_info row number.")
  } else if (sum(colnames(object@expression_data) != object@sample_info$sample_id) > 0) {
    errors <- c(errors,
                "expression_data colnames must equal sample_info$sample_id.")
  }

  ## expression vs variable
  if (nrow(object@expression_data) != nrow(object@variable_info)) {
    errors <- c(errors,
                "expression_data row number must equal variable_info row number.")
  } else if (sum(rownames(object@expression_data) != object@variable_info$variable_id) > 0) {
    errors <- c(errors,
                "expression_data rownames must equal variable_info$variable_id.")
  }

  ## note consistency
  if (nrow(object@sample_info_note) > 0) {
    if (ncol(object@sample_info) != nrow(object@sample_info_note)) {
      errors <- c(errors,
                  "sample_info column number must equal sample_info_note row number.")
    } else if (sum(colnames(object@sample_info) != object@sample_info_note$name) > 0) {
      errors <- c(errors,
                  "sample_info colnames must equal sample_info_note$name.")
    }
  }

  if (nrow(object@variable_info_note) > 0) {
    if (ncol(object@variable_info) != nrow(object@variable_info_note)) {
      errors <- c(errors,
                  "variable_info column number must equal variable_info_note row number.")
    } else if (sum(colnames(object@variable_info) != object@variable_info_note$name) > 0) {
      errors <- c(errors,
                  "variable_info colnames must equal variable_info_note$name.")
    }
  }

  if (length(errors) == 0) TRUE else errors
}


# ---------------------------------------------------------------------------
# S4 Class
# ---------------------------------------------------------------------------

#' @title S4 Class Representing a msi_dataset
#'
#' @description
#' Holds all data related to a mass spectrometry imaging (MSI) experiment.
#' The structure mirrors \pkg{massdataset} for full inter-operability:
#' pixels are treated as "samples" and stored in \code{sample_info}.
#'
#' @slot expression_data \code{data.frame} or \code{matrix}.
#'   Rows = variables (m/z features), columns = pixels.
#' @slot ms2_data \code{list} of MS2 spectra objects.
#' @slot annotation_table \code{data.frame} of metabolite annotations.
#' @slot sample_info \code{data.frame} with one row per pixel.
#'   Required columns: \code{sample_id}, \code{x}, \code{y}.
#'   Optional columns: \code{z}, \code{section_id}, \code{resolution},
#'   \code{matrix}.
#' @slot variable_info \code{data.frame} with one row per m/z feature.
#'   Required columns: \code{variable_id}, \code{mz}.
#'   Optional: \code{ms2} (logical).
#' @slot sample_info_note \code{data.frame} with columns \code{name} and
#'   \code{meaning} describing each column of \code{sample_info}.
#' @slot variable_info_note \code{data.frame} with columns \code{name} and
#'   \code{meaning} describing each column of \code{variable_info}.
#' @slot process_info \code{list} of \code{tidymass_parameter} objects
#'   recording processing history.
#' @slot other_files \code{list} for miscellaneous attached files.
#' @slot version \code{character} package version used to create the object.
#' @slot activated \code{character} currently active slot
#'   (\code{"expression_data"}, \code{"sample_info"}, or
#'   \code{"variable_info"}).
#'
#' @name msi_dataset-class
#' @rdname msi_dataset-class
#' @docType class
#' @exportClass msi_dataset

setClass(
  Class = "msi_dataset",
  representation(
    expression_data  = "data.frame",
    ms2_data         = "list",
    annotation_table = "data.frame",
    sample_info      = "data.frame",
    variable_info    = "data.frame",
    sample_info_note     = "data.frame",
    variable_info_note   = "data.frame",
    process_info     = "list",
    other_files      = "list",
    version          = "character",
    activated        = "character"
  ),
  prototype(
    expression_data      = data.frame(),
    ms2_data             = list(),
    annotation_table     = data.frame(),
    sample_info          = data.frame(),
    variable_info        = data.frame(),
    sample_info_note     = data.frame(),
    variable_info_note   = data.frame(),
    process_info         = list(),
    other_files          = list(),
    version              = "0.1.0",
    activated            = "expression_data"
  ),
  validity = check_msi_dataset_class
)


# ---------------------------------------------------------------------------
# Constructor
# ---------------------------------------------------------------------------

#' Create a New msi_dataset Object
#'
#' @description
#' The primary constructor for \code{msi_dataset}.  Performs basic validation
#' before building the object.
#'
#' @param expression_data A \code{data.frame} or \code{matrix}.
#'   Rows are variables, columns are pixels.  Row names must match
#'   \code{variable_info$variable_id}; column names must match
#'   \code{sample_info$sample_id}.
#' @param sample_info A \code{data.frame} with one row per pixel.
#'   Must contain \code{sample_id}, \code{x}, and \code{y}.
#' @param variable_info A \code{data.frame} with one row per m/z feature.
#'   Must contain \code{variable_id} and \code{mz}.
#' @param sample_info_note Optional \code{data.frame} with \code{name} and
#'   \code{meaning} columns describing \code{sample_info}.  Auto-generated
#'   if omitted.
#' @param variable_info_note Optional \code{data.frame} with \code{name} and
#'   \code{meaning} columns describing \code{variable_info}.  Auto-generated
#'   if omitted.
#'
#' @return An object of class \code{msi_dataset}.
#' @export
#'
#' @examples
#' \dontrun{
#' object <- create_msi_dataset(
#'   expression_data = expression_data,
#'   sample_info     = sample_info,
#'   variable_info   = variable_info
#' )
#' object
#' }

create_msi_dataset <- function(expression_data,
                                sample_info,
                                variable_info,
                                sample_info_note,
                                variable_info_note) {

  ## ---- auto-generate notes if missing ----
  if (missing(sample_info_note) || is.null(sample_info_note)) {
    sample_info_note <- data.frame(
      name    = colnames(sample_info),
      meaning = colnames(sample_info),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  if (missing(variable_info_note) || is.null(variable_info_note)) {
    variable_info_note <- data.frame(
      name    = colnames(variable_info),
      meaning = colnames(variable_info),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  ## ---- pre-construction check ----
  check_result <- check_msi_dataset(
    expression_data    = expression_data,
    sample_info        = sample_info,
    variable_info      = variable_info,
    sample_info_note   = sample_info_note,
    variable_info_note = variable_info_note
  )

  if (stringr::str_detect(check_result, "^error")) {
    stop(check_result)
  }

  ## ---- coerce to data.frame (preserving rownames) ----
  old_rn         <- rownames(expression_data)
  expression_data <- as.data.frame(expression_data, check.names = FALSE)
  rownames(expression_data) <- old_rn

  sample_info        <- as.data.frame(sample_info,        check.names = FALSE,
                                       stringsAsFactors = FALSE)
  variable_info      <- as.data.frame(variable_info,      check.names = FALSE,
                                       stringsAsFactors = FALSE)
  sample_info_note   <- as.data.frame(sample_info_note,   check.names = FALSE,
                                       stringsAsFactors = FALSE)
  variable_info_note <- as.data.frame(variable_info_note, check.names = FALSE,
                                       stringsAsFactors = FALSE)

  ## ---- ensure dimnames are consistent ----
  colnames(expression_data) <- sample_info$sample_id
  rownames(expression_data) <- variable_info$variable_id

  ## ---- record creation in process_info ----
  process_info <- list()
  parameter <- new(
    Class         = "tidymass_parameter",
    pacakge_name  = "msidataset",
    function_name = "create_msi_dataset()",
    parameter     = list("no" = "no"),
    time          = Sys.time()
  )
  process_info$create_msi_dataset <- parameter

  ## ---- construct S4 object ----
  object <- new(
    Class              = "msi_dataset",
    expression_data    = expression_data,
    ms2_data           = list(),
    annotation_table   = data.frame(),
    sample_info        = sample_info,
    variable_info      = variable_info,
    sample_info_note   = sample_info_note,
    variable_info_note = variable_info_note,
    process_info       = process_info,
    other_files        = list(),
    version            = tryCatch(
      as.character(utils::packageVersion("msidataset")),
      error = function(e) "0.1.0"
    ),
    activated          = "expression_data"
  )

  invisible(object)
}
