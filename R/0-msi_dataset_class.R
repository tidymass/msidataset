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

.normalize_sample_info <- function(sample_info) {
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
  sample_info[, c(intersect(required_order, colnames(sample_info)),
                  remaining_cols), drop = FALSE]
}

.normalize_variable_info <- function(variable_info) {
  if (!("ms2" %in% colnames(variable_info))) {
    variable_info$ms2 <- "no"
  }
  variable_info
}

#' @keywords internal
check_msi_dataset <- function(expression_data,
                              sample_info,
                              variable_info,
                              sample_info_note,
                              variable_info_note) {
  errors <- character()

  sample_info <- as.data.frame(sample_info, check.names = FALSE, stringsAsFactors = FALSE)
  variable_info <- as.data.frame(variable_info, check.names = FALSE, stringsAsFactors = FALSE)

  if (!is.data.frame(expression_data) && !is.matrix(expression_data)) {
    errors <- c(errors, "expression_data must be a data.frame or matrix.")
  }
  expression_df <- as.data.frame(expression_data, check.names = FALSE, stringsAsFactors = FALSE)

  # Reuse massdataset core checks where possible.
  sample_info_md <- sample_info
  sample_info_note_md <- sample_info_note
  if ("sample_id" %in% colnames(sample_info_md) && !("class" %in% colnames(sample_info_md))) {
    sample_info_md$class <- "MSI_PIXEL"
    if (!missing(sample_info_note) && !is.null(sample_info_note)) {
      sample_info_note_md <- rbind(
        sample_info_note_md,
        data.frame(
          name = "class",
          meaning = "Temporary class column for massdataset-compatible validation.",
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  md_result <- massdataset::check_mass_dataset(
    expression_data = expression_df,
    sample_info = sample_info_md,
    variable_info = variable_info,
    sample_info_note = sample_info_note_md,
    variable_info_note = variable_info_note
  )
  if (!identical(md_result, "all good.")) {
    errors <- c(errors, sub("^error:\\s*", "", md_result))
  }

  if (!("mz" %in% colnames(variable_info))) {
    errors <- c(errors, "variable_info must have mz.")
  }

  required_sample_cols <- c("sample_id", "x", "y")
  missing_sample_cols <- setdiff(required_sample_cols, colnames(sample_info))
  if (length(missing_sample_cols) > 0) {
    errors <- c(errors,
                paste0("sample_info missing required columns: ",
                       paste(missing_sample_cols, collapse = ", ")))
  }

  if ("x" %in% colnames(sample_info) && !is.numeric(sample_info$x)) {
    errors <- c(errors, "sample_info$x must be numeric.")
  }
  if ("y" %in% colnames(sample_info) && !is.numeric(sample_info$y)) {
    errors <- c(errors, "sample_info$y must be numeric.")
  }

  if (length(errors) > 0) {
    return(paste0("error: ", paste(unique(errors), collapse = "; ")))
  }
  "all good."
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
  check_result <- check_msi_dataset(
    expression_data = object@expression_data,
    sample_info = object@sample_info,
    variable_info = object@variable_info,
    sample_info_note = object@sample_info_note,
    variable_info_note = object@variable_info_note
  )
  errors <- character()
  if (!identical(check_result, "all good.")) {
    errors <- c(
      errors,
      unlist(strsplit(sub("^error:\\s*", "", check_result), ";\\s*"),
             use.names = FALSE)
    )
  }

  roi_check <- try(.validate_roi_integrity(object), silent = TRUE)
  if (inherits(roi_check, "try-error")) {
    errors <- c(errors, as.character(roi_check))
  } else if (!isTRUE(roi_check)) {
    errors <- c(errors, roi_check)
  }

  if (length(errors) == 0) {
    TRUE
  } else {
    unique(errors)
  }
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
#'   ROI is stored as \code{other_files$roi} using \code{msi_roi_set}.
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
  if (missing(sample_info_note)) {
    sample_info_note <- NULL
  }
  if (missing(variable_info_note)) {
    variable_info_note <- NULL
  }

  sample_info <- .normalize_sample_info(sample_info)
  variable_info <- .normalize_variable_info(variable_info)
  
  ## ---- auto-generate notes if missing ----
  sample_info_note <- .sync_info_note(sample_info_note, sample_info)
  variable_info_note <- .sync_info_note(variable_info_note, variable_info)
  
  ## ---- pre-construction check ----
  check_result <- check_msi_dataset(
    expression_data    = expression_data,
    sample_info        = sample_info,
    variable_info      = variable_info,
    sample_info_note   = sample_info_note,
    variable_info_note = variable_info_note
  )
  
  if (!identical(check_result, "all good.")) {
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
    process_info       = list(),
    other_files        = list(),
    version            = tryCatch(
      as.character(utils::packageVersion("msidataset")),
      error = function(e) "0.1.0"
    ),
    activated          = "expression_data"
  )
  object <- add_process_log(
    object = object,
    function_name = "create_msi_dataset()",
    parameter = list(no = "no")
  )

  invisible(object)
}
