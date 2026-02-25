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
#'     \item "pixel_info": Activates the pixel/spatial information (x, y coordinates).
#'     \item "annotation_table": Activates the annotation table.
#'   }
#'   By default, it activates the "pixel_info" dataset.
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
#' # object = activate_msi_dataset(.data = object, what = "pixel_info")
#' # object@@activated

activate_msi_dataset <-
  function(.data,
           what = c(
             "pixel_info",
             "sample_info",
             "variable_info",
             "expression_data",
             "annotation_table"
           )) {
    check_object_class(object = .data, class = "msi_dataset")
    
    # 匹配参数
    what <- match.arg(what)
    
    # 修改对象槽位
    slot(object = .data, name = "activated") <- what
    
    return(.data)
  }