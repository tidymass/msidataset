# ============================================================================
# msi_dataset Class Definition and Constructor - UPDATED
# Updated to match new structure:
# - sample_info: sample_id, x, y, z, subject_id, section_id, resolution, matrix
# - variable_info: variable_id, mz, ms2
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
  
  # Required columns in variable_info
  if (!("variable_id" %in% colnames(variable_info)))
    errors <- c(errors, "variable_info must have variable_id.")
  
  if (sum(duplicated(variable_info$variable_id)) > 0)
    errors <- c(errors, "variable_id has duplicated items.")
  
  if (!("mz" %in% colnames(variable_info)))
    errors <- c(errors, "variable_info must have mz.")
  
  if (!("ms2" %in% colnames(variable_info)))
    errors <- c(errors, "variable_info must have ms2.")
  
  # Required columns in sample_info
  if (!("sample_id" %in% colnames(sample_info)))
    errors <- c(errors, "sample_info must have sample_id.")
  
  if (sum(duplicated(sample_info$sample_id)) > 0)
    errors <- c(errors, "sample_id has duplicated items.")
  
  required_sample_cols <- c(
    "sample_id", "x", "y", "z", "subject_id", "section_id", "resolution", "matrix"
  )
  missing_sample_cols <- setdiff(required_sample_cols, colnames(sample_info))
  if (length(missing_sample_cols) > 0) {
    errors <- c(
      errors,
      paste0("sample_info missing required columns: ",
             paste(missing_sample_cols, collapse = ", "))
    )
  }
  
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
  
  ## variable_info - Required columns
  if (!("variable_id" %in% colnames(object@variable_info)))
    errors <- c(errors, "variable_info must have variable_id.")
  
  if (sum(duplicated(object@variable_info$variable_id)) > 0)
    errors <- c(errors, "variable_id has duplicated items.")
  
  if (!("mz" %in% colnames(object@variable_info)))
    errors <- c(errors, "variable_info must have mz.")
  
  if (!("ms2" %in% colnames(object@variable_info)))
    errors <- c(errors, "variable_info must have ms2.")
  
  ## sample_info - Required columns
  if (!("sample_id" %in% colnames(object@sample_info)))
    errors <- c(errors, "sample_info must have sample_id.")
  
  if (sum(duplicated(object@sample_info$sample_id)) > 0)
    errors <- c(errors, "sample_id has duplicated items.")
  
  required_sample_cols <- c(
    "sample_id", "x", "y", "z", "subject_id", "section_id", "resolution", "matrix"
  )
  missing_sample_cols <- setdiff(required_sample_cols, colnames(object@sample_info))
  if (length(missing_sample_cols) > 0) {
    errors <- c(
      errors,
      paste0("sample_info missing required columns: ",
             paste(missing_sample_cols, collapse = ", "))
    )
  }
  
  ## MSI spatial columns (required: x, y; optional: z)
  ## z is optional so we don't check for it here
  
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
#'   Column names must match \code{sample_info$sample_id}.
#'   Row names must match \code{variable_info$variable_id}.
#'
#' @slot ms2_data \code{list} of MS2 spectra objects.
#'
#' @slot annotation_table \code{data.frame} of metabolite annotations.
#'
#' @slot sample_info \code{data.frame} with one row per pixel.
#'   Required columns: \code{sample_id}, \code{x}, \code{y}.
#'   Optional columns: \code{z}, \code{subject_id}, \code{section_id},
#'   \code{resolution}, \code{matrix}.
#'
#' @slot variable_info \code{data.frame} with one row per m/z feature.
#'   Required columns: \code{variable_id}, \code{mz}.
#'   Optional columns: \code{ms2} (logical, yes/no indicator).
#'
#' @slot sample_info_note \code{data.frame} with columns \code{name} and
#'   \code{meaning} describing each column of \code{sample_info}.
#'
#' @slot variable_info_note \code{data.frame} with columns \code{name} and
#'   \code{meaning} describing each column of \code{variable_info}.
#'
#' @slot process_info \code{list} of \code{tidymass_parameter} objects
#'   recording processing history.
#'
#' @slot other_files \code{list} for miscellaneous attached files.
#'
#' @slot version \code{character} package version used to create the object.
#'
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
    expression_data  = "ANY",
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
    expression_data      = matrix(numeric(0), nrow = 0, ncol = 0),
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
#'   Rows are variables (m/z features), columns are pixels.
#'   Row names must match \code{variable_info$variable_id};
#'   column names must match \code{sample_info$sample_id}.
#'
#' @param sample_info A \code{data.frame} with one row per pixel.
#'   Must contain: \code{sample_id}, \code{x}, \code{y}.
#'   Can optionally contain: \code{z}, \code{subject_id}, \code{section_id},
#'   \code{resolution}, \code{matrix}.
#'
#' @param variable_info A \code{data.frame} with one row per m/z feature.
#'   Must contain: \code{variable_id}, \code{mz}.
#'   Can optionally contain: \code{ms2} (logical, yes/no).
#'
#' @param sample_info_note Optional \code{data.frame} with \code{name} and
#'   \code{meaning} columns describing \code{sample_info}.  Auto-generated
#'   if omitted.
#'
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
  if ("pixel_id" %in% colnames(sample_info) &&
      !("sample_id" %in% colnames(sample_info))) {
    sample_info$sample_id <- sample_info$pixel_id
  }
  
  if (!("z" %in% colnames(sample_info))) {
    sample_info$z <- 1
  }
  if (!("subject_id" %in% colnames(sample_info))) {
    sample_info$subject_id <- "subject_1"
  }
  if (!("section_id" %in% colnames(sample_info))) {
    sample_info$section_id <- "section_1"
  }
  if (!("resolution" %in% colnames(sample_info))) {
    sample_info$resolution <- NA_real_
  }
  if (!("matrix" %in% colnames(sample_info))) {
    sample_info$matrix <- NA_character_
  }
  required_order <- c("sample_id", "x", "y", "z", "subject_id",
                      "section_id", "resolution", "matrix")
  remaining_cols <- setdiff(colnames(sample_info), required_order)
  sample_info <- sample_info[, c(intersect(required_order, colnames(sample_info)),
                                 remaining_cols), drop = FALSE]
  
  if (!("ms2" %in% colnames(variable_info))) {
    variable_info$ms2 <- "no"
  }
  
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
  
  if (grepl("^error", check_result)) {
    stop(check_result)
  }
  
  ## ---- keep matrix-like backend, only coerce when needed ----
  if (is.data.frame(expression_data)) {
    expression_data <- as.matrix(expression_data)
  }
  
  sample_info        <- as.data.frame(sample_info,        check.names = FALSE,
                                      stringsAsFactors = FALSE)
  variable_info      <- as.data.frame(variable_info,      check.names = FALSE,
                                      stringsAsFactors = FALSE)
  sample_info_note   <- as.data.frame(sample_info_note,   check.names = FALSE,
                                      stringsAsFactors = FALSE)
  variable_info_note <- as.data.frame(variable_info_note, check.names = FALSE,
                                      stringsAsFactors = FALSE)
  
  ## ---- ensure dimnames are consistent ----
  dimnames(expression_data) <- list(variable_info$variable_id, sample_info$sample_id)
  
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
