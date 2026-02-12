#' Namespace imports
#' @keywords internal
#' @noRd
#' @importFrom magrittr %>%
NULL

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
    
    ##check sample_info format (now for pixels)
    if (!"sample_id" %in% colnames(object@sample_info)) {
      msg <- "sample_info must have sample_id (pixel_id)."
      errors <- c(errors, msg)
    }
    
    if (sum(duplicated(object@sample_info$sample_id)) > 0) {
      msg <- "sample_id has duplicated items."
      errors <- c(errors, msg)
    }
    
    # Check required spatial columns in sample_info
    if (!"x" %in% colnames(object@sample_info)) {
      msg <- "sample_info must have x (horizontal spatial coordinate)."
      errors <- c(errors, msg)
    }
    
    if (!"y" %in% colnames(object@sample_info)) {
      msg <- "sample_info must have y (vertical spatial coordinate)."
      errors <- c(errors, msg)
    }
    
    if (!"z" %in% colnames(object@sample_info)) {
      msg <- "sample_info must have z (tissue slice depth / relative to a reference slice)."
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
    
    ###expression_data and sample_info (pixels)
    if (ncol(object@expression_data) != nrow(object@sample_info)) {
      msg <-
        "expression_data's column number should be same with sample_info's row number."
      errors <- c(errors, msg)
    } else{
      if (sum(colnames(object@expression_data) != object@sample_info$sample_id) > 0) {
        msg <-
          "expression_data's column names must be identical with sample_info's sample_id."
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
#' related to mass spectrometry imaging experiments. This version has been 
#' restructured to match the massdataset structure for better compatibility.
#' 
#' Key changes from original structure:
#' - pixel_info has been merged into sample_info (pixels are treated as samples)
#' - pixel_info_note has been removed (notes are in sample_info_note)
#' - This allows direct use of massdataset functions with minimal code changes
#'
#' @slot expression_data A data.frame or matrix containing the expression data (variables × samples/pixels).
#' @slot ms2_data A list containing MS2 spectra data.
#' @slot annotation_table A data.frame containing annotation information.
#' @slot sample_info A data.frame containing sample/pixel information with required columns:
#'   sample_id (pixel_id), x, y, z, and optional section_id, resolution, matrix, etc.
#' @slot variable_info A data.frame containing variable information with required column:
#'   variable_id, and typically mz, ms2 (logical).
#' @slot sample_info_note A data.frame containing notes about the sample_info columns.
#' @slot variable_info_note A data.frame containing notes about the variable_info columns.
#' @slot process_info A list containing information about the data processing steps.
#' @slot other_files A list containing other related files.
#' @slot version A character string indicating the version of the dataset.
#' @slot activated A character string indicating which slot is currently activated for manipulation.
#'
#' @examples
#' \dontrun{
#' # Create a new msi_dataset object
#' new_data <- new("msi_dataset", 
#'                 expression_data = matrix(0, nrow = 100, ncol = 1000),
#'                 ms2_data = list(),
#'                 annotation_table = data.frame(),
#'                 sample_info = data.frame(
#'                   sample_id = paste0("pixel_", 1:1000),
#'                   x = rep(1:40, 25),
#'                   y = rep(1:25, each = 40),
#'                   z = 1
#'                 ),
#'                 variable_info = data.frame(
#'                   variable_id = paste0("var_", 1:100),
#'                   mz = seq(100, 1000, length.out = 100)
#'                 ),
#'                 sample_info_note = data.frame(),
#'                 variable_info_note = data.frame(),
#'                 process_info = list(),
#'                 other_files = list(),
#'                 version = "0.2.0",
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
    expression_data = "ANY",  # Can be data.frame, matrix, or DelayedArray
    ms2_data = "list",
    annotation_table = "data.frame",
    sample_info = "data.frame",  # Now contains pixel information (x, y, z, etc.)
    variable_info = "data.frame",
    sample_info_note = "data.frame",  # Now contains notes for pixel/sample info
    variable_info_note = "data.frame",
    process_info = "list",
    other_files = "list",
    version = "character",
    activated = "character"
  ),
  prototype = list(
    activated = "expression_data"
  ),
  validity = check_msi_dataset_class
)


#' Create a New msi_dataset Object
#' @description
#' This function creates a new msi_dataset object with structure matching massdataset.
#' Pixels are now treated as samples, stored in sample_info with spatial coordinates.
#' This allows direct compatibility with massdataset functions.
#'
#' @docType methods
#' @author Xiaoyu Xu
#' \email{XIAOYU021@@e.ntu.edu.sg}
#' @param expression_data Expression data matrix/data.frame (variables × pixels).
#' @param sample_info Sample/pixel information data.frame with required columns:
#'   sample_id (pixel ID), x, y, z. Optional: section_id, resolution, matrix.
#' @param variable_info Variable information data.frame with variable_id, mz columns.
#' @param sample_info_note Optional. Descriptions for sample_info columns.
#' @param variable_info_note Optional. Descriptions for variable_info columns.
#' @return An object of class \code{msi_dataset}.
#' @export
#' @examples
#' \dontrun{
#' # Create sample_info (formerly pixel_info + sample_info combined)
#' sample_info <- data.frame(
#'   sample_id = paste0("pixel_", 1:1000),
#'   x = rep(1:40, 25),
#'   y = rep(1:25, each = 40),
#'   z = 1,
#'   section_id = "section_1",
#'   resolution = 50,
#'   matrix = "DHB"
#' )
#' 
#' # Create variable_info
#' variable_info <- data.frame(
#'   variable_id = paste0("mz_", 1:100),
#'   mz = seq(100, 1000, length.out = 100),
#'   ms2 = FALSE
#' )
#' 
#' # Create expression_data
#' expression_data <- matrix(
#'   rnorm(100 * 1000),
#'   nrow = 100,
#'   ncol = 1000
#' )
#' rownames(expression_data) <- variable_info$variable_id
#' colnames(expression_data) <- sample_info$sample_id
#' 
#' object <- create_msi_dataset(
#'   expression_data = expression_data,
#'   sample_info = sample_info,
#'   variable_info = variable_info
#' )
#' }

create_msi_dataset <-
  function(expression_data,
           sample_info,
           variable_info,
           sample_info_note,
           variable_info_note) {
    
    # ========== Pre-creation validation ==========
    # Check variable_info
    if (!"variable_id" %in% colnames(variable_info)) {
      stop("variable_info must have variable_id column.")
    }
    
    if (sum(duplicated(variable_info$variable_id)) > 0) {
      stop("variable_id has duplicated items.")
    }
    
    # Check sample_info (now contains pixel information)
    if (!"sample_id" %in% colnames(sample_info)) {
      stop("sample_info must have sample_id column (pixel_id).")
    }
    
    if (sum(duplicated(sample_info$sample_id)) > 0) {
      stop("sample_id has duplicated items.")
    }
    
    # Check required spatial columns
    if (!"x" %in% colnames(sample_info)) {
      stop("sample_info must have x (horizontal spatial coordinate).")
    }
    
    if (!"y" %in% colnames(sample_info)) {
      stop("sample_info must have y (vertical spatial coordinate).")
    }
    
    if (!"z" %in% colnames(sample_info)) {
      stop("sample_info must have z (tissue slice depth).")
    }
    
    # Check dimensions match
    if (ncol(expression_data) != nrow(sample_info)) {
      stop("expression_data's column number must equal sample_info's row number.")
    }
    
    if (nrow(expression_data) != nrow(variable_info)) {
      stop("expression_data's row number must equal variable_info's row number.")
    }
    
    # Check IDs match (only if names are already set)
    if (!is.null(colnames(expression_data)) && 
        !all(colnames(expression_data) == "")) {
      if (sum(colnames(expression_data) != sample_info$sample_id) > 0) {
        stop("expression_data's column names must match sample_info's sample_id.")
      }
    }
    
    if (!is.null(rownames(expression_data)) && 
        !all(rownames(expression_data) == "")) {
      if (sum(rownames(expression_data) != variable_info$variable_id) > 0) {
        stop("expression_data's row names must match variable_info's variable_id.")
      }
    }
    
    # ========== Create info_note if missing ==========
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
    
    # ========== Convert to data.frame (keep expression_data flexible) ==========
    # expression_data can be matrix, data.frame, or DelayedArray
    # Only convert if it's not already in a suitable format
    if (!is.matrix(expression_data) && 
        !is.data.frame(expression_data) &&
        !methods::is(expression_data, "DelayedArray")) {
      expression_data <- as.data.frame(expression_data)
    }
    
    # Ensure rownames and colnames are set
    if (is.null(rownames(expression_data))) {
      rownames(expression_data) <- variable_info$variable_id
    }
    if (is.null(colnames(expression_data))) {
      colnames(expression_data) <- sample_info$sample_id
    }
    
    sample_info <- as.data.frame(sample_info)
    variable_info <- as.data.frame(variable_info)
    sample_info_note <- as.data.frame(sample_info_note)
    variable_info_note <- as.data.frame(variable_info_note)
    
    # ========== Create object ==========
    object <- new(
      Class = "msi_dataset",
      expression_data = expression_data,
      ms2_data = list(),
      annotation_table = data.frame(),
      sample_info = sample_info,  # Now contains pixel info (x, y, z, etc.)
      variable_info = variable_info,
      sample_info_note = sample_info_note,  # Now contains pixel info notes
      variable_info_note = variable_info_note,
      process_info = process_info,
      version = as.character(utils::packageVersion(pkg = "msidataset")),
      activated = "expression_data"
    )
    invisible(object)
  }


#' Convert Old msi_dataset to New Structure
#' @description
#' Helper function to convert old msi_dataset objects (with separate pixel_info)
#' to new structure (pixel_info merged into sample_info).
#'
#' @param old_object An old msi_dataset object with pixel_info slot
#' @return A new msi_dataset object with restructured data
#' @export
#' @examples
#' \dontrun{
#' # Convert old object to new structure
#' new_object <- convert_msi_dataset_structure(old_object)
#' }

convert_msi_dataset_structure <- function(old_object) {
  
  # Check if this is an old-style object
  if (!.hasSlot(old_object, "pixel_info")) {
    message("Object already has new structure. Returning as is.")
    return(old_object)
  }
  
  # Merge pixel_info into sample_info structure
  pixel_info <- old_object@pixel_info
  old_sample_info <- old_object@sample_info
  
  # Create new sample_info by adding pixel info to each pixel
  # Assuming pixel_id links to sample_id
  new_sample_info <- pixel_info
  colnames(new_sample_info)[colnames(new_sample_info) == "pixel_id"] <- "sample_id"
  
  # If there was additional sample-level info, add it
  if (nrow(old_sample_info) > 0 && "sample_id" %in% colnames(pixel_info)) {
    # Merge sample-level metadata with pixel info
    sample_columns <- setdiff(colnames(old_sample_info), colnames(new_sample_info))
    if (length(sample_columns) > 0) {
      merge_cols <- old_sample_info[, c("sample_id", sample_columns), drop = FALSE]
      new_sample_info <- merge(new_sample_info, merge_cols, 
                               by.x = "section_id", by.y = "sample_id", 
                               all.x = TRUE, sort = FALSE)
    }
  }
  
  # Create new sample_info_note
  new_sample_info_note <- old_object@pixel_info_note
  if (nrow(old_object@sample_info_note) > 0) {
    # Combine notes
    new_sample_info_note <- rbind(
      old_object@pixel_info_note,
      old_object@sample_info_note[
        !old_object@sample_info_note$name %in% old_object@pixel_info_note$name,
      ]
    )
  }
  
  # Create new object
  create_msi_dataset(
    expression_data = old_object@expression_data,
    sample_info = new_sample_info,
    variable_info = old_object@variable_info,
    sample_info_note = new_sample_info_note,
    variable_info_note = old_object@variable_info_note
  )
}