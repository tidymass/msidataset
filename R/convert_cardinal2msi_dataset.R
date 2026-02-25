#' @title Convert Cardinal Object to msi_dataset
#' @description 
#' Convert Cardinal object to msi_dataset
#'
#' @param cardinal_object Cardinal
#' @param sample_name sample name
#' @param section_name section name
#' @param resolution um
#' @param matrix_name CHCA|DHB|NEDC
#' @export

convert_cardinal2msi_dataset <- function(cardinal_object,
                                         sample_name = "sample_1",
                                         section_name = "section_1",
                                         resolution = NA,
                                         matrix_name = NA) {
  
  # 1. check cardinal
  if (!requireNamespace("Cardinal", quietly = TRUE)) {
    stop("please install Cardinal package: BiocManager::install('Cardinal')")
  }
  
  # 2. Expression Data
  cat(crayon::blue("Extracting intensity matrix...\n"))
  expression_data <- as.matrix(Cardinal::spectra(cardinal_object))
  
  # 3. Sample Info
  sample_info <- as.data.frame(Cardinal::pData(cardinal_object))
  
  # Unify ID naming
  sample_info$sample_id <- paste0(sample_name, "_pixel_", seq_len(nrow(sample_info)))
  sample_info$section_id <- section_name
  if (!is.na(resolution)) sample_info$resolution <- resolution
  if (!is.na(matrix_name)) sample_info$matrix <- matrix_name
  
  # 4. Variable Info
  variable_info <- as.data.frame(Cardinal::fData(cardinal_object))
  
  # Ensure presence of mz column
  if (!"mz" %in% colnames(variable_info)) {
    variable_info$mz <- Cardinal::mz(cardinal_object)
  }
  
  if (!exists("generate_variable_id")) {
    stop("Function generate_variable_id not found. Please ensure source('generate_variable_id.R') has been run.")
  }
  
  variable_info$variable_id <- generate_variable_id(variable_info$mz)
  
  # 5. Match row and column names
  rownames(expression_data) <- as.character(variable_info$variable_id)
  colnames(expression_data) <- as.character(sample_info$sample_id)
  
  # 6. Create msi_dataset object
  # create_msi_dataset defined in 0-msi_dataset_class.R
  object <- create_msi_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )
  
  # Log processing activities
  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "convert_cardinal2msi_dataset()",
    parameter = list(
      sample_name = sample_name,
      section_name = section_name,
      cardinal_class = class(cardinal_object)[1]
    ),
    time = Sys.time()
  )
  object@process_info$convert_cardinal2msi_dataset <- parameter
  
  cat(crayon::green("Successfully converted!\n"))
  return(object)
}