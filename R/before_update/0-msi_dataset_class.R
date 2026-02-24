#' Check the Validity of a msi_dataset Object
#' 
#' @description
#' This function checks the validity of a msi_dataset object.
#' It performs various checks on the object's slots to ensure they meet the required criteria.
#' The structure now matches massdataset for better compatibility.
#'
#' @docType methods
#' @author Xiaoyu Xu
#' \email{XIAOYU021@@e.ntu.edu.sg}
#' @param object An object of class \code{msi_dataset}.
#' @return
#' Returns \code{TRUE} if the object is valid, otherwise returns a character vector of error messages.
#'
#' @export
#' @examples
#' data("expression_data")
#' data("sample_info")
#' data("variable_info")
#' object <-
#'   create_msi_dataset(expression_data = expression_data,
#'                       sample_info = sample_info,
#'                       variable_info = variable_info)
#'   check_msi_dataset_class(object)

check_msi_dataset_class <- function(object) {
  
  errors <- character()
  
  ## ---------------- variable_info ----------------
  if (all(colnames(object@variable_info) != "variable_id")) {
    errors <- c(errors, "variable_info must have variable_id.")
  }
  
  if (sum(duplicated(object@variable_info$variable_id)) > 0) {
    errors <- c(errors, "variable_id has duplicated items.")
  }
  
  ## ---------------- sample_info ----------------
  if (all(colnames(object@sample_info) != "sample_id")) {
    errors <- c(errors, "sample_info must have sample_id.")
  }
  
  if (sum(duplicated(object@sample_info$sample_id)) > 0) {
    errors <- c(errors, "sample_id has duplicated items.")
  }
  
  ## MSI spatial info
  if (all(colnames(object@sample_info) != "x")) {
    errors <- c(errors, "sample_info must have x.")
  }
  
  if (all(colnames(object@sample_info) != "y")) {
    errors <- c(errors, "sample_info must have y.")
  }
  
  ## ---------------- note tables ----------------
  if (nrow(object@sample_info_note) > 0) {
    
    if (all(colnames(object@sample_info_note) != "name")) {
      errors <- c(errors, "sample_info_note must have column: name.")
    }
    
    if (all(colnames(object@sample_info_note) != "meaning")) {
      errors <- c(errors, "sample_info_note must have column: meaning")
    }
  }
  
  if (nrow(object@variable_info_note) > 0) {
    
    if (all(colnames(object@variable_info_note) != "name")) {
      errors <- c(errors, "variable_info_note must have column: name.")
    }
    
    if (all(colnames(object@variable_info_note) != "meaning")) {
      errors <- c(errors, "variable_info_note must have column: meaning")
    }
  }
  
  ## ---------------- expression vs sample ----------------
  if (ncol(object@expression_data) != nrow(object@sample_info)) {
    errors <- c(
      errors,
      "expression_data's column number should equal sample_info row number."
    )
  } else {
    if (sum(colnames(object@expression_data) != object@sample_info$sample_id) > 0) {
      errors <- c(
        errors,
        "expression_data colnames must equal sample_info$sample_id."
      )
    }
  }
  
  ## ---------------- expression vs variable ----------------
  if (nrow(object@expression_data) != nrow(object@variable_info)) {
    errors <- c(
      errors,
      "expression_data row number should equal variable_info row number."
    )
  } else {
    if (sum(rownames(object@expression_data) != object@variable_info$variable_id) > 0) {
      errors <- c(
        errors,
        "expression_data rownames must equal variable_info$variable_id."
      )
    }
  }
  
  ## ---------------- note consistency ----------------
  if (nrow(object@sample_info_note) > 0) {
    
    if (ncol(object@sample_info) != nrow(object@sample_info_note)) {
      errors <- c(
        errors,
        "sample_info column number must equal sample_info_note row number."
      )
    } else {
      if (sum(colnames(object@sample_info) != object@sample_info_note$name) > 0) {
        errors <- c(
          errors,
          "sample_info colnames must equal sample_info_note$name."
        )
      }
    }
  }
  
  if (nrow(object@variable_info_note) > 0) {
    
    if (ncol(object@variable_info) != nrow(object@variable_info_note)) {
      errors <- c(
        errors,
        "variable_info column number must equal variable_info_note row number."
      )
    } else {
      if (sum(colnames(object@variable_info) != object@variable_info_note$name) > 0) {
        errors <- c(
          errors,
          "variable_info colnames must equal variable_info_note$name."
        )
      }
    }
  }
  
  if (length(errors) == 0) TRUE else errors
}

#' @title S4 Class Representing a msi_dataset
#' 
#' @description
#' The \code{msi_dataset} class is designed to hold various types of data 
#' related to mass spectrometry imaging (MSI) experiments, including 
#' expression data, MS2 data, annotation tables, spatial pixel information, 
#' and other related information.
#'
#' @slot expression_data A data.frame or matrix containing the expression data.
#' @slot ms2_data A list containing MS2 data.
#' @slot annotation_table A data.frame containing annotation information.
#' @slot sample_info A data.frame containing sample (pixel) information.
#' @slot variable_info A data.frame containing variable information.
#' @slot sample_info_note A data.frame containing notes about the sample information.
#' @slot variable_info_note A data.frame containing notes about the variable information.
#' @slot process_info A list containing information about the data processing steps.
#' @slot other_files A list containing other related files.
#' @slot version A character string indicating the version of the dataset.
#' @slot activated A character string indicating which slot is currently activated for manipulation.
#'
#' 
#' @examples
#' \dontrun{
#' # Create a new msi_dataset object
#' new_data <- new("msi_dataset", 
#'                 expression_data = data.frame(),
#'                 ms2_data = list(),
#'                 annotation_table = data.frame(),
#'                 sample_info = data.frame(),
#'                 variable_info = data.frame(),
#'                 sample_info_note = data.frame(),
#'                 variable_info_note = data.frame(),
#'                 process_info = list(),
#'                 other_files = list(),
#'                 version = "0.1.0",
#'                 activated = "expression_data")
#' }
#'
#' @name msi_dataset-class
#' @rdname msi_dataset-class
#' @docType class
#' @exportClass msi_dataset

setClass(
  Class = "msi_dataset",
  representation(
    expression_data = "data.frame",
    ms2_data = "list",
    annotation_table = "data.frame",
    sample_info = "data.frame",
    variable_info = "data.frame",
    sample_info_note = "data.frame",
    variable_info_note = "data.frame",
    process_info = "list",
    other_files = "list",
    version = "character",
    activated = "character"
  ),
  validity = check_msi_dataset_class
)

#' Create a New msi_dataset Object
#' @description
#' This function creates a new msi_dataset object 
#' by taking in various data frames and performing checks on them.
#'
#' @docType methods
#' @author Xiaoyu Xu
#' \email{XIAOYU021@@e.ntu.edu.sg}
#' @param expression_data MSI peak table name.
#' \url{https://tidymass.github.io/massdataset/articles/data_import_and_export.html}
#' @param sample_info Sample (pixel) information name.
#' \url{https://tidymass.github.io/massdataset/articles/data_import_and_export.html}
#' @param variable_info MSI peak table name.
#' Columns are samples (pixels) and rows are variables.
#' \url{https://tidymass.github.io/massdataset/articles/data_import_and_export.html}
#' @param sample_info_note Sample information name.
#' \url{https://tidymass.github.io/massdataset/articles/data_import_and_export.html}
#' @param variable_info_note Sample information name.
#' \url{https://tidymass.github.io/massdataset/articles/data_import_and_export.html}
#' @return An object of class \code{msi_dataset}.
#' @export
#' @examples
#' data("expression_data")
#' data("sample_info")
#' data("sample_info_note")
#' data("variable_info")
#' data("variable_info_note")
#' object =
#'   create_msi_dataset(
#'     expression_data = expression_data,
#'     sample_info = sample_info,
#'     variable_info = variable_info,
#'     sample_info_note = sample_info_note,
#'     variable_info_note = variable_info_note
#'   )
#' object

create_msi_dataset <-
  function(expression_data,
           sample_info,
           variable_info,
           sample_info_note,
           variable_info_note) {
    
    check_result <-
      check_msi_dataset(
        expression_data = expression_data,
        sample_info = sample_info,
        variable_info = variable_info,
        sample_info_note = sample_info_note,
        variable_info_note = variable_info_note
      )
    
    if (stringr::str_detect(check_result, "error")) {
      stop(check_result)
    }
    
    if (missing(sample_info_note)) {
      sample_info_note <-
        data.frame(
          name = colnames(sample_info),
          meaning = colnames(sample_info),
          check.names = FALSE
        )
    }
    
    if (missing(variable_info_note)) {
      variable_info_note <-
        data.frame(
          name = colnames(variable_info),
          meaning = colnames(variable_info),
          check.names = FALSE
        )
    }
    
    process_info <- list()
    
    parameter <- new(
      Class = "tidymass_parameter",
      pacakge_name = "msidataset",
      function_name = "create_msi_dataset()",
      parameter = list("no" = "no"),
      time = Sys.time()
    )
    
    process_info$create_msi_dataset <- parameter
    
    old_rownames <- rownames(expression_data)
    expression_data <- as.data.frame(expression_data)
    rownames(expression_data) <- old_rownames
    
    sample_info <- as.data.frame(sample_info)
    variable_info <- as.data.frame(variable_info)
    sample_info_note <- as.data.frame(sample_info_note)
    variable_info_note <- as.data.frame(variable_info_note)
    
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
      version = as.character(utils::packageVersion(pkg = "msidataset"))
    )
    
    invisible(object)
  }
