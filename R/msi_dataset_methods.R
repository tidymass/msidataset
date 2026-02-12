#' Show method for msi_dataset
#'
#' @description
#' Display method for msi_dataset objects. Shows a summary of the dataset
#' including dimensions, sample information, variable information, and 
#' processing history. This version is compatible with massdataset structure.
#'
#' @param object A msi_dataset object
#' 
#' @return Invisibly returns the object
#' 
#' @docType methods
#' @name show,msi_dataset-method
#' @rdname show-methods
#' @aliases show,msi_dataset-method
#'
#' @export
setMethod(
  f = "show",
  signature = "msi_dataset",
  definition = function(object) {
    ###check again
    check_result = check_msi_dataset_class(object)
    
    if (!isTRUE(check_result)) {
      cat(crayon::red("Validation errors found:\n"))
      cat(crayon::red(paste(check_result, collapse = "\n"), "\n"))
      cat(crayon::red(
        "You may have changed the slots. Try to use update_msi_dataset().\n"
      ))
    }
    
    cat(crayon::yellow(paste(rep("-", 20), collapse = ""), "\n"))
    cat(crayon::green("msidataset version:", object@version, "\n"))
    cat(crayon::yellow("Structure: massdataset-compatible\n"))
    cat(crayon::yellow(paste(rep("-", 20), collapse = ""), "\n"))
    
    # 1. expression_data
    cat(crayon::green("1. expression_data:"))
    if (methods::is(object@expression_data, "DelayedArray")) {
      cat(
        "[",
        nrow(object@expression_data),
        "x",
        ncol(object@expression_data),
        "DelayedArray]\n"
      )
    } else if (is.matrix(object@expression_data)) {
      cat(
        "[",
        nrow(object@expression_data),
        "x",
        ncol(object@expression_data),
        "matrix]\n"
      )
    } else {
      cat(
        "[",
        nrow(object@expression_data),
        "x",
        ncol(object@expression_data),
        "data.frame]\n"
      )
    }
    
    # 2. variable_info
    cat(crayon::green("2. variable_info:"))
    cat(
      "[",
      nrow(object@variable_info),
      "x",
      ncol(object@variable_info),
      "data.frame]\n"
    )
    if (nrow(object@variable_info) <= 5) {
      variable_id <- object@variable_info$variable_id
    } else{
      variable_id <-
        c(head(object@variable_info$variable_id, 3), "...",
          tail(object@variable_info$variable_id, 2))
    }
    cat(nrow(object@variable_info), "variables:")
    cat(variable_id)
    cat("\n")
    
    # 3. sample_info (now contains pixel data)
    cat(crayon::green("3. sample_info:"))
    cat("[",
        nrow(object@sample_info),
        "x",
        ncol(object@sample_info),
        "data.frame]\n")
    if (nrow(object@sample_info) <= 5) {
      sample_id <- object@sample_info$sample_id
    } else{
      sample_id <-
        c(head(object@sample_info$sample_id, 3), "...",
          tail(object@sample_info$sample_id, 2))
    }
    cat(nrow(object@sample_info), "samples (pixels):")
    cat(sample_id)
    cat("\n")
    
    # Show spatial extent if available
    if (all(c("x", "y", "z") %in% colnames(object@sample_info))) {
      cat("  Spatial extent: x[", 
          min(object@sample_info$x), "-", max(object@sample_info$x), "], ",
          "y[", min(object@sample_info$y), "-", max(object@sample_info$y), "], ",
          "z[", min(object@sample_info$z), "-", max(object@sample_info$z), "]\n",
          sep = "")
    }
    
    # Show sections if available
    if ("section_id" %in% colnames(object@sample_info)) {
      sections <- unique(object@sample_info$section_id)
      if (length(sections) <= 5) {
        cat("  Sections:", paste(sections, collapse = ", "), "\n")
      } else {
        cat("  Sections:", length(sections), "sections\n")
      }
    }
    
    # 4. sample_info_note
    cat(crayon::green("4. sample_info_note:"))
    cat(
      "[",
      nrow(object@sample_info_note),
      "x",
      ncol(object@sample_info_note),
      "data.frame]\n"
    )
    
    # 5. variable_info_note
    cat(crayon::green("5. variable_info_note:"))
    cat(
      "[",
      nrow(object@variable_info_note),
      "x",
      ncol(object@variable_info_note),
      "data.frame]\n"
    )
    
    # 6. ms2_data
    cat(crayon::green("6. ms2_data:"))
    if (length(object@ms2_data) > 0) {
      cat("[",
          sum(unlist(lapply(object@ms2_data, function(x) {
            length(unique(x@variable_id))
          }))),
          "variables x",
          sum(unlist(lapply(object@ms2_data, function(x) {
            length(unique(x@ms2_spectra))
          }))),
          "MS2 spectra]\n")
    } else {
      cat("[0 MS2 spectra]\n")
    }
    
    cat(crayon::yellow(paste(rep("-", 20), collapse = ""), "\n"))
    cat(crayon::green("Processing information\n"))
    if (.hasSlot(object = object, name = "process_info") &
        length(object@process_info) != 0) {
      process_info <- object@process_info
      cat(crayon::green(length(process_info), "processings in total\n"))
      if (length(process_info) > 5) {
        cat(crayon::green("Latest 3 processings shown\n"))
        process_info <- tail(process_info, 3)
      }
      for (idx in seq_along(process_info)) {
        cat(crayon::green(names(process_info)[idx], paste(rep("-", 10), collapse = ""), "\n"))
        if (length(process_info[[idx]]) == 1) {
          data.frame(
            "Package" = process_info[[idx]]@pacakge_name,
            "Function used" = process_info[[idx]]@function_name,
            "Time" = process_info[[idx]]@time
          ) %>%
            print()
        } else{
          data.frame(
            "Package" = process_info[[idx]] %>% lapply(function(x)
              x@pacakge_name) %>% unlist(),
            "Function used" = process_info[[idx]] %>% lapply(function(x)
              x@function_name) %>% unlist(),
            "Time" = process_info[[idx]] %>% lapply(function(x)
              as.character(x@time)) %>% unlist()
          ) %>%
            print()
        }
      }
    } else{
      cat(crayon::red("No processing information available.\n"))
    }
    
    cat(crayon::yellow(paste(rep("-", 20), collapse = ""), "\n"))
    cat(crayon::blue("Note: This msi_dataset has massdataset-compatible structure.\n"))
    cat(crayon::blue("Pixels are stored as 'samples' in sample_info.\n"))
  }
)


#' Extract/subset method for msi_dataset
#'
#' @description
#' Subset a msi_dataset object by variables (rows) and samples/pixels (columns).
#'
#' @param x A msi_dataset object
#' @param i Row indices (variables). Can be numeric, logical, or character vector
#' @param j Column indices (samples/pixels). Can be numeric, logical, or character vector
#' @param drop Logical, whether to drop dimensions. Default is FALSE
#' 
#' @return A subsetted msi_dataset object
#' 
#' @docType methods
#' @rdname extract-methods
#' @aliases [,msi_dataset,ANY,ANY,ANY-method
#' 
#' @export
setMethod(
  "[",
  signature = "msi_dataset",
  function(x, i, j, drop = FALSE) {
    
    # Handle missing indices
    if (missing(i)) {
      i <- seq_len(nrow(x@expression_data))
    }
    if (missing(j)) {
      j <- seq_len(ncol(x@expression_data))
    }
    
    # Subset expression_data
    new_expression_data <- x@expression_data[i, j, drop = drop]
    
    # Subset variable_info
    new_variable_info <- x@variable_info[i, , drop = FALSE]
    
    # Subset sample_info (pixels)
    new_sample_info <- x@sample_info[j, , drop = FALSE]
    
    # Update notes if needed
    new_variable_info_note <- x@variable_info_note
    new_sample_info_note <- x@sample_info_note
    
    # Create new object
    new_object <- new(
      Class = "msi_dataset",
      expression_data = new_expression_data,
      ms2_data = x@ms2_data,  # Keep all MS2 data for now
      annotation_table = x@annotation_table,
      sample_info = new_sample_info,
      variable_info = new_variable_info,
      sample_info_note = new_sample_info_note,
      variable_info_note = new_variable_info_note,
      process_info = x@process_info,
      other_files = x@other_files,
      version = x@version,
      activated = x@activated
    )
    
    return(new_object)
  }
)


#' Get dimensions of msi_dataset
#'
#' @description
#' Returns the dimensions of the expression data in a msi_dataset object.
#'
#' @param x A msi_dataset object
#' 
#' @return Integer vector of length 2: c(n_variables, n_samples)
#' 
#' @docType methods
#' @rdname dim-methods
#' @aliases dim,msi_dataset-method
#' 
#' @export
setMethod(
  "dim",
  "msi_dataset",
  function(x) {
    c(nrow(x@expression_data), ncol(x@expression_data))
  }
)


#' Get dimension names of msi_dataset
#'
#' @description
#' Returns the row and column names of the expression data.
#'
#' @param x A msi_dataset object
#' 
#' @return List with two elements: variable IDs (rows) and sample IDs (columns)
#' 
#' @docType methods
#' @rdname dimnames-methods
#' @aliases dimnames,msi_dataset-method
#' 
#' @export
setMethod(
  "dimnames",
  "msi_dataset",
  function(x) {
    list(
      rownames(x@expression_data),
      colnames(x@expression_data)
    )
  }
)


#' Update msi_dataset object
#'
#' @description
#' Updates an msi_dataset object to ensure consistency between slots.
#' Useful after manual modifications.
#'
#' @param object A msi_dataset object
#' @return Updated msi_dataset object
#' @export
update_msi_dataset <- function(object) {
  
  # Ensure expression_data dimnames match info tables
  if (!is.null(object@variable_info) && nrow(object@variable_info) > 0) {
    rownames(object@expression_data) <- object@variable_info$variable_id
  }
  
  if (!is.null(object@sample_info) && nrow(object@sample_info) > 0) {
    colnames(object@expression_data) <- object@sample_info$sample_id
  }
  
  # Update info_notes if they don't match
  if (nrow(object@sample_info_note) == 0 || 
      !all(object@sample_info_note$name == colnames(object@sample_info))) {
    object@sample_info_note <- data.frame(
      name = colnames(object@sample_info),
      meaning = colnames(object@sample_info),
      check.names = FALSE
    )
  }
  
  if (nrow(object@variable_info_note) == 0 || 
      !all(object@variable_info_note$name == colnames(object@variable_info))) {
    object@variable_info_note <- data.frame(
      name = colnames(object@variable_info),
      meaning = colnames(object@variable_info),
      check.names = FALSE
    )
  }
  
  # Add update to process_info
  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "update_msi_dataset()",
    parameter = list("no" = "no"),
    time = Sys.time()
  )
  
  object@process_info$update_msi_dataset <- parameter
  
  # Validate
  check_result <- check_msi_dataset_class(object)
  if (!isTRUE(check_result)) {
    warning("Validation errors after update:\n", 
            paste(check_result, collapse = "\n"))
  }
  
  return(object)
}