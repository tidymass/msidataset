#' Check the Validity of a msi_dataset Object
#'
#' @description
#' This function checks the validity of a msi_dataset object.
#' It performs various checks on the object's slots to ensure they meet the required criteria.
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
#' data("pixel_info")
#' data("variable_info")
#' data("sample_info")
#' object <-
#'   create_msi_dataset(expression_data = expression_data,
#'                       pixel_info = pixel_info,
#'                       variable_info = variable_info,
#'                       sample_info = sample_info)
#'   check_msi_dataset_class(object)

check_msi_dataset_class <-
  function(object) {
    errors <- character()
    ##check variable_info format
    if (!"variable_id" %in% colnames(object@variable_info)) {
      msg <- "variable_info must have variable_id."
      errors <- c(errors, msg)
    }
    
    if (sum(duplicated(object@variable_info$variable_id)) > 0) {
      msg <- "variable_id has duplicated items."
      errors <- c(errors, msg)
    }
    
    ##check pixel_info format
    if (!"pixel_id" %in% colnames(object@pixel_info)) {
      msg <- "pixel_info must have pixel_id."
      errors <- c(errors, msg)
    }
    
    if (sum(duplicated(object@pixel_info$pixel_id)) > 0) {
      msg <- "pixel_id has duplicated items."
      errors <- c(errors, msg)
    }
    
    if (!"x" %in% colnames(object@pixel_info)) {
      msg <- "pixel_info must have x (horizontal spatial coordinate)."
      errors <- c(errors, msg)
    }
    
    if (!"y" %in% colnames(object@pixel_info)) {
      msg <- "pixel_info must have y (vertical spatial coordinate)."
      errors <- c(errors, msg)
    }
    
    if (!"z" %in% colnames(object@pixel_info)) {
      msg <- "pixel_info must have z (tissue slice depth / relative to a reference slice)."
      errors <- c(errors, msg)
    }
    
    ##check sample_info format
    if (!"sample_id" %in% colnames(object@sample_info)) {
      msg <- "sample_info must have sample_id."
      errors <- c(errors, msg)
    }
    
    if (sum(duplicated(object@sample_info$sample_id)) > 0) {
      msg <- "sample_id has duplicated items."
      errors <- c(errors, msg)
    }
    
    ##check sample_info_note format
    if (nrow(object@sample_info_note) > 0) {
      if (all(colnames(object@sample_info_note) != "name")) {
        msg <- "sample_info_note must have column: name."
        errors <- c(errors, msg)
      }
      
      if (all(colnames(object@sample_info_note) != "meaning")) {
        msg <- "sample_info_note must have column: meaning"
        errors <- c(errors, msg)
      }
    }
    
    ##check pixel_info_note format
    if (nrow(object@pixel_info_note) > 0) {
      if (all(colnames(object@pixel_info_note) != "name")) {
        msg <- "pixel_info_note must have column: name."
        errors <- c(errors, msg)
      }
      
      if (all(colnames(object@pixel_info_note) != "meaning")) {
        msg <- "pixel_info_note must have column: meaning"
        errors <- c(errors, msg)
      }
    }
    
    ##check variable_info_note format
    if (nrow(object@variable_info_note) > 0) {
      if (all(colnames(object@variable_info_note) != "name")) {
        msg <- "variable_info_note must have column: name."
        errors <- c(errors, msg)
      }
      
      if (all(colnames(object@variable_info_note) != "meaning")) {
        msg <- "variable_info_note must have column: meaning"
        errors <- c(errors, msg)
      }
    }
    
    ###expression_data and pixel_info
    if (ncol(object@expression_data) != nrow(object@pixel_info)) {
      msg <-
        "expression_data's column number should be same with pixel_info's row number."
      errors <- c(errors, msg)
    } else{
      if (sum(colnames(object@expression_data) != object@pixel_info$pixel_id) > 0) {
        msg <-
          "expression_data's column names must be identical with pixel_info's pixel_id."
        errors <- c(errors, msg)
      }
    }
    
    ###expression_data and variable_info
    if (nrow(object@expression_data) != nrow(object@variable_info)) {
      msg <-
        "expression_data's row number should be same with variable_info's row number."
      errors <- c(errors, msg)
    } else{
      if (sum(rownames(object@expression_data) != object@variable_info$variable_id) > 0) {
        msg <-
          "expression_data's row names must be identical with variable_info's variable_id"
        errors <- c(errors, msg)
      }
    }
    
    ###pixel_info and pixel_info_note
    if (nrow(object@pixel_info_note) > 0) {
      if (ncol(object@pixel_info) != nrow(object@pixel_info_note)) {
        msg <-
          "pixel_info's column number should be same with pixel_info_note's row number."
        errors <- c(errors, msg)
      } else{
        if (sum(colnames(object@pixel_info) != object@pixel_info_note$name) > 0) {
          msg <-
            "pixel_info's column names must be identical with pixel_info_note's name."
          errors <- c(errors, msg)
        }
      }
    }
    
    ###sample_info and sample_info_note
    if (nrow(object@sample_info_note) > 0) {
      if (ncol(object@sample_info) != nrow(object@sample_info_note)) {
        msg <-
          "sample_info's column number should be same with sample_info_note's row number."
        errors <- c(errors, msg)
      } else{
        if (sum(colnames(object@sample_info) != object@sample_info_note$name) > 0) {
          msg <-
            "sample_info's column names must be identical with sample_info_note's name."
          errors <- c(errors, msg)
        }
      }
    }
    
    ###variable_info and variable_info_note
    if (nrow(object@variable_info_note) > 0) {
      if (ncol(object@variable_info) != nrow(object@variable_info_note)) {
        msg <-
          "variable_info's column number should be same with variable_info_note's row number."
        errors <- c(errors, msg)
      } else{
        if (sum(colnames(object@variable_info) != object@variable_info_note$name) > 0) {
          msg <-
            "variable_info's column names must be identical with variable_info_note's name."
          errors <- c(errors, msg)
        }
      }
    }
    
    if (length(errors) == 0)
      TRUE
    else
      errors
  }

#' @title S4 Class Representing a msi_dataset
#' 
#' @description
#' The \code{msi_dataset} class is designed to hold various types of data 
#' related to mass spectrometry experiments, including expression data, 
#' MS2 data, annotation tables, and other related information.
#'
#' @slot expression_data A data.frame containing the expression data.
#' @slot ms2_data A list containing MS2 data.
#' @slot annotation_table A data.frame containing annotation information.
#' @slot pixel_info A data.frame containing pixel information.
#' @slot sample_info A data.frame containing sample information.
#' @slot variable_info A data.frame containing variable information.
#' @slot sample_info_note A data.frame containing notes about the sample information.
#' @slot pixel_info_note A data.frame containing notes about the pixel information.
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
#'                 pixel_info = data.frame(),
#'                 variable_info = data.frame(),
#'                 sample_info_note = data.frame(),
#'                 pixel_info_note = data.frame(),
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
    pixel_info = "data.frame",
    variable_info = "data.frame",
    sample_info_note = "data.frame",
    pixel_info_note = "data.frame",
    variable_info_note = "data.frame",
    process_info = "list",
    other_files = "list",
    version = "character",
    activated = "character"
  ),
  prototype = list(
    activated = "expression_data"  # 添加默认值
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
#' @param expression_data Expression data matrix (variables × pixels).
#' @param pixel_info Pixel information data.frame with pixel_id, x, y, z columns.
#' @param sample_info Sample information data.frame with sample_id column.
#' @param variable_info Variable information data.frame with variable_id, mz columns.
#' @param pixel_info_note Optional. Descriptions for pixel_info columns.
#' @param sample_info_note Optional. Descriptions for sample_info columns.
#' @param variable_info_note Optional. Descriptions for variable_info columns.
#' @return An object of class \code{msi_dataset}.
#' @export
#' @examples
#' \dontrun{
#' data("expression_data")
#' data("pixel_info")
#' data("pixel_info_note")
#' data("sample_info")
#' data("sample_info_note")
#' data("variable_info")
#' data("variable_info_note")
#' object =
#'   create_msi_dataset(
#'     expression_data = expression_data,
#'     pixel_info = pixel_info,
#'     sample_info = sample_info,
#'     variable_info = variable_info,
#'     pixel_info_note = pixel_info_note,
#'     sample_info_note = sample_info_note,
#'     variable_info_note = variable_info_note
#'   )
#' object
#' }

create_msi_dataset <-
  function(expression_data,
           pixel_info,
           sample_info,
           variable_info,
           pixel_info_note,
           sample_info_note,
           variable_info_note) {
    
    # ========== Pre-creation validation (before creating object) ==========
    # Check variable_info
    if (!"variable_id" %in% colnames(variable_info)) {
      stop("variable_info must have variable_id column.")
    }
    
    if (sum(duplicated(variable_info$variable_id)) > 0) {
      stop("variable_id has duplicated items.")
    }
    
    # Check pixel_info
    if (!"pixel_id" %in% colnames(pixel_info)) {
      stop("pixel_info must have pixel_id column.")
    }
    
    if (sum(duplicated(pixel_info$pixel_id)) > 0) {
      stop("pixel_id has duplicated items.")
    }
    
    if (!"x" %in% colnames(pixel_info)) {
      stop("pixel_info must have x (horizontal spatial coordinate).")
    }
    
    if (!"y" %in% colnames(pixel_info)) {
      stop("pixel_info must have y (vertical spatial coordinate).")
    }
    
    if (!"z" %in% colnames(pixel_info)) {
      stop("pixel_info must have z (tissue slice depth).")
    }
    
    # Check sample_info
    if (!"sample_id" %in% colnames(sample_info)) {
      stop("sample_info must have sample_id column.")
    }
    
    if (sum(duplicated(sample_info$sample_id)) > 0) {
      stop("sample_id has duplicated items.")
    }
    
    # Check dimensions match
    if (ncol(expression_data) != nrow(pixel_info)) {
      stop("expression_data's column number must equal pixel_info's row number.")
    }
    
    if (nrow(expression_data) != nrow(variable_info)) {
      stop("expression_data's row number must equal variable_info's row number.")
    }
    
    # Check IDs match (only if names are already set)
    if (!is.null(colnames(expression_data)) && 
        !all(colnames(expression_data) == "")) {
      if (sum(colnames(expression_data) != pixel_info$pixel_id) > 0) {
        stop("expression_data's column names must match pixel_info's pixel_id.")
      }
    }
    
    if (!is.null(rownames(expression_data)) && 
        !all(rownames(expression_data) == "")) {
      if (sum(rownames(expression_data) != variable_info$variable_id) > 0) {
        stop("expression_data's row names must match variable_info's variable_id.")
      }
    }
    
    # ========== Create info_note if missing ==========
    if (missing(pixel_info_note)) {
      pixel_info_note <-
        data.frame(
          name = colnames(pixel_info),
          meaning = colnames(pixel_info),
          check.names = FALSE
        )
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
    
    # ========== Create process_info ==========
    process_info = list()
    
    parameter <- new(
      Class = "tidymass_parameter",
      pacakge_name = "msidataset",
      function_name = "create_msi_dataset()",
      parameter = list("no" = "no"),
      time = Sys.time()
    )
    
    process_info$create_msi_dataset = parameter
    
    # ========== Convert to data.frame ==========
    old_rownames <- rownames(expression_data)
    expression_data <- as.data.frame(expression_data)
    rownames(expression_data) <- old_rownames
    pixel_info <- as.data.frame(pixel_info)
    sample_info <- as.data.frame(sample_info)
    variable_info <- as.data.frame(variable_info)
    pixel_info_note <- as.data.frame(pixel_info_note)
    sample_info_note <- as.data.frame(sample_info_note)
    variable_info_note <- as.data.frame(variable_info_note)
    
    # ========== Create object ==========
    object <- new(
      Class = "msi_dataset",
      expression_data = expression_data,
      ms2_data = list(),
      annotation_table = data.frame(),
      pixel_info = pixel_info,
      sample_info = sample_info,
      variable_info = variable_info,
      pixel_info_note = pixel_info_note,
      sample_info_note = sample_info_note,
      variable_info_note = variable_info_note,
      process_info = process_info,
      version = as.character(utils::packageVersion(pkg = "msidataset")),
      activated = "expression_data"  # 添加 activated 字段
    )
    invisible(object)
  }