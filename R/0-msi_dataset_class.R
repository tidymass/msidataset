# ============================================================================
# msi_dataset Class Definition and Constructor
# Structure aligned with massdataset logic and MSI-specific metadata.
# ============================================================================

#' @keywords internal
check_msi_dataset <- function(expression_data,
                              sample_info,
                              variable_info,
                              sample_info_note,
                              variable_info_note) {
  errors <- character()

  required_sample_cols <- c(
    "sample_id", "x", "y", "z", "subject_id", "section_id", "resolution", "matrix"
  )
  required_variable_cols <- c("variable_id", "mz", "ms2")

  miss_sample <- setdiff(required_sample_cols, colnames(sample_info))
  if (length(miss_sample) > 0) {
    errors <- c(errors, paste0("sample_info missing columns: ", paste(miss_sample, collapse = ", "), "."))
  }

  miss_variable <- setdiff(required_variable_cols, colnames(variable_info))
  if (length(miss_variable) > 0) {
    errors <- c(errors, paste0("variable_info missing columns: ", paste(miss_variable, collapse = ", "), "."))
  }

  if ("sample_id" %in% colnames(sample_info) && anyDuplicated(sample_info$sample_id) > 0) {
    errors <- c(errors, "sample_id has duplicated items.")
  }

  if ("variable_id" %in% colnames(variable_info) && anyDuplicated(variable_info$variable_id) > 0) {
    errors <- c(errors, "variable_id has duplicated items.")
  }

  if (ncol(expression_data) != nrow(sample_info)) {
    errors <- c(errors, "expression_data column number must equal sample_info row number.")
  }

  if (nrow(expression_data) != nrow(variable_info)) {
    errors <- c(errors, "expression_data row number must equal variable_info row number.")
  }

  if (length(errors) > 0) {
    return(paste("error:", paste(errors, collapse = "; ")))
  }

  "pass"
}


#' Check the Validity of a msi_dataset Object
#'
#' @param object An object of class \code{msi_dataset}.
#' @return \code{TRUE} if valid, otherwise a character vector of error messages.
#' @export
check_msi_dataset_class <- function(object) {
  errors <- character()

  required_sample_cols <- c(
    "sample_id", "x", "y", "z", "subject_id", "section_id", "resolution", "matrix"
  )
  required_variable_cols <- c("variable_id", "mz", "ms2")

  miss_sample <- setdiff(required_sample_cols, colnames(object@sample_info))
  if (length(miss_sample) > 0) {
    errors <- c(errors, paste0("sample_info missing columns: ", paste(miss_sample, collapse = ", "), "."))
  }

  miss_variable <- setdiff(required_variable_cols, colnames(object@variable_info))
  if (length(miss_variable) > 0) {
    errors <- c(errors, paste0("variable_info missing columns: ", paste(miss_variable, collapse = ", "), "."))
  }

  if ("sample_id" %in% colnames(object@sample_info) && anyDuplicated(object@sample_info$sample_id) > 0) {
    errors <- c(errors, "sample_id has duplicated items.")
  }

  if ("variable_id" %in% colnames(object@variable_info) && anyDuplicated(object@variable_info$variable_id) > 0) {
    errors <- c(errors, "variable_id has duplicated items.")
  }

  if (ncol(object@expression_data) != nrow(object@sample_info)) {
    errors <- c(errors, "expression_data column number must equal sample_info row number.")
  } else if (sum(colnames(object@expression_data) != object@sample_info$sample_id) > 0) {
    errors <- c(errors, "expression_data colnames must equal sample_info$sample_id.")
  }

  if (nrow(object@expression_data) != nrow(object@variable_info)) {
    errors <- c(errors, "expression_data row number must equal variable_info row number.")
  } else if (sum(rownames(object@expression_data) != object@variable_info$variable_id) > 0) {
    errors <- c(errors, "expression_data rownames must equal variable_info$variable_id.")
  }

  if (nrow(object@sample_info_note) > 0) {
    if (!all(c("name", "meaning") %in% colnames(object@sample_info_note))) {
      errors <- c(errors, "sample_info_note must have columns: name, meaning.")
    }
    if (ncol(object@sample_info) != nrow(object@sample_info_note)) {
      errors <- c(errors, "sample_info column number must equal sample_info_note row number.")
    }
  }

  if (nrow(object@variable_info_note) > 0) {
    if (!all(c("name", "meaning") %in% colnames(object@variable_info_note))) {
      errors <- c(errors, "variable_info_note must have columns: name, meaning.")
    }
    if (ncol(object@variable_info) != nrow(object@variable_info_note)) {
      errors <- c(errors, "variable_info column number must equal variable_info_note row number.")
    }
  }

  if (length(errors) == 0) TRUE else errors
}


#' @title S4 Class Representing a msi_dataset
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
    sample_info_note = "data.frame",
    variable_info_note = "data.frame",
    process_info     = "list",
    other_files      = "list",
    version          = "character",
    activated        = "character"
  ),
  prototype(
    expression_data = data.frame(),
    ms2_data = list(),
    annotation_table = data.frame(),
    sample_info = data.frame(),
    variable_info = data.frame(),
    sample_info_note = data.frame(),
    variable_info_note = data.frame(),
    process_info = list(),
    other_files = list(),
    version = "0.1.0",
    activated = "expression_data"
  ),
  validity = check_msi_dataset_class
)


#' Create a New msi_dataset Object
#'
#' @return An object of class \code{msi_dataset}.
#' @export
create_msi_dataset <- function(expression_data,
                               sample_info,
                               variable_info,
                               sample_info_note,
                               variable_info_note) {
  if (missing(sample_info_note) || is.null(sample_info_note)) {
    sample_info_note <- data.frame(
      name = colnames(sample_info),
      meaning = colnames(sample_info),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  if (missing(variable_info_note) || is.null(variable_info_note)) {
    variable_info_note <- data.frame(
      name = colnames(variable_info),
      meaning = colnames(variable_info),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  check_result <- check_msi_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info,
    sample_info_note = sample_info_note,
    variable_info_note = variable_info_note
  )

  if (stringr::str_detect(check_result, "^error")) {
    stop(check_result)
  }

  old_rn <- rownames(expression_data)
  expression_data <- as.data.frame(expression_data, check.names = FALSE)
  rownames(expression_data) <- old_rn

  sample_info <- as.data.frame(sample_info, check.names = FALSE, stringsAsFactors = FALSE)
  variable_info <- as.data.frame(variable_info, check.names = FALSE, stringsAsFactors = FALSE)
  sample_info_note <- as.data.frame(sample_info_note, check.names = FALSE, stringsAsFactors = FALSE)
  variable_info_note <- as.data.frame(variable_info_note, check.names = FALSE, stringsAsFactors = FALSE)

  colnames(expression_data) <- sample_info$sample_id
  rownames(expression_data) <- variable_info$variable_id

  process_info <- list()
  process_info$create_msi_dataset <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "create_msi_dataset",
    parameter = list(n_sample = nrow(sample_info), n_variable = nrow(variable_info)),
    time = Sys.time()
  )

  object <- new(
    Class = "msi_dataset",
    expression_data = expression_data,
    ms2_data = list(),
    annotation_table = data.frame(),
    sample_info = sample_info,
    variable_info = variable_info,
    sample_info_note = sample_info_note,
    variable_info_note = variable_info_note,
    process_info = process_info,
    other_files = list(),
    version = tryCatch(as.character(utils::packageVersion("msidataset")), error = function(e) "0.1.0"),
    activated = "expression_data"
  )

  invisible(object)
}
