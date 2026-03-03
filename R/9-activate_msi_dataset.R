#' @title Activate a Specific Dataset in msi_dataset
#'
#' @description Activates a dataset within the msi_dataset object for further analysis.
#'
#' @param .data A `msi_dataset` object.
#' @param what A character string specifying which dataset to activate.
#'   The available options are:
#'   \itemize{
#'     \item "sample_info": Activates the sample information dataset.
#'     \item "variable_info": Activates the variable information dataset.
#'     \item "expression_data": Activates the expression data (intensity matrix).
#'     \item "annotation_table": Activates the annotation table.
#'   }
#'   By default, it activates the "sample_info" dataset.
#'
#' @details 
#' This function allows users to activate a specific dataset within the 
#' `msi_dataset` object. The activated dataset can be then used in subsequent
#' tidyverse-like analysis steps (e.g., mutate, filter).
#'
#' @return 
#' Returns the `msi_dataset` object with the specified dataset activated.
#'
#' @author Xiaoyu Xu \email{XIAOYU021@@e.ntu.edu.sg}
#' @export
#' @examples
#' # Assume object is a msi_dataset
#' # object = activate_msi_dataset(.data = object, what = "sample_info")
#' # object@@activated

activate_msi_dataset <-
  function(.data,
           what = c(
             "sample_info",
             "variable_info",
             "expression_data",
             "annotation_table"
           )) {
    if (!methods::is(.data, "msi_dataset")) {
      stop("object must be a msi_dataset.")
    }

    what <- match.arg(what)

    slot(object = .data, name = "activated") <- what

    return(.data)
  }
