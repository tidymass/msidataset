#' @title Convert Cardinal object to msi_dataset
#' @param cardinal_object Cardinal MSI object.
#' @param sample_id Sample identifier.
#' @param subject_id Subject identifier.
#' @param section_id Section identifier.
#' @param resolution Resolution (um/pixel).
#' @param matrix Matrix name.
#' @export
convert_cardinal2msi_dataset <- function(cardinal_object,
                                         sample_id = "sample_1",
                                         subject_id = "subject_1",
                                         section_id = "section_1",
                                         resolution = NA_real_,
                                         matrix = NA_character_) {
  if (!requireNamespace("Cardinal", quietly = TRUE)) {
    stop("please install Cardinal package: BiocManager::install('Cardinal')")
  }

  expression_data <- as.matrix(Cardinal::spectra(cardinal_object))
  sample_info_raw <- as.data.frame(Cardinal::pData(cardinal_object), stringsAsFactors = FALSE)
  variable_info_raw <- as.data.frame(Cardinal::fData(cardinal_object), stringsAsFactors = FALSE)

  coords <- Cardinal::coord(cardinal_object)
  sample_n <- nrow(sample_info_raw)
  sample_info <- data.frame(
    sample_id = if (length(sample_id) == 1) paste0(sample_id, "_", seq_len(sample_n)) else sample_id,
    x = coords[, 1],
    y = coords[, 2],
    z = if (ncol(coords) >= 3) coords[, 3] else 1,
    subject_id = subject_id,
    section_id = section_id,
    resolution = resolution,
    matrix = matrix,
    stringsAsFactors = FALSE
  )

  mz_vec <- if ("mz" %in% colnames(variable_info_raw)) variable_info_raw$mz else Cardinal::mz(cardinal_object)
  variable_info <- data.frame(
    variable_id = generate_variable_id(mz_vec),
    mz = mz_vec,
    ms2 = "no",
    stringsAsFactors = FALSE
  )

  rownames(expression_data) <- variable_info$variable_id
  colnames(expression_data) <- sample_info$sample_id

  object <- create_msi_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

  object@process_info$convert_cardinal2msi_dataset <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "convert_cardinal2msi_dataset",
    parameter = list(sample_id = sample_id, subject_id = subject_id, section_id = section_id),
    time = Sys.time()
  )

  object
}
