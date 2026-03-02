#' @title Read imzML as msi_dataset
#' @description Encapsulates imzML reading + conversion workflow.
#' @param path_1 Input imzML path.
#' @param peak_align_tolerance Numeric ppm tolerance placeholder for compatibility.
#' @param path_2 Optional output .rds path.
#' @param sample_id_prefix Prefix for generated sample IDs.
#' @param subject_id Subject identifier.
#' @param section_id Section identifier.
#' @param resolution Resolution value.
#' @param matrix Matrix name.
#' @return msi_dataset object.
#' @export
read_imzml <- function(path_1 = "input",
                       peak_align_tolerance = 50,
                       path_2 = "output",
                       sample_id_prefix = "sample",
                       subject_id = "subject_1",
                       section_id = "section_1",
                       resolution = NA_real_,
                       matrix = NA_character_) {
  if (!file.exists(path_1)) {
    stop("path_1 does not exist.")
  }
  if (!requireNamespace("Cardinal", quietly = TRUE)) {
    stop("please install Cardinal package: BiocManager::install('Cardinal')")
  }

  cardinal_object <- Cardinal::readImzML(path_1)
  object <- convert_cardinal2msi_dataset(
    cardinal_object = cardinal_object,
    sample_id = sample_id_prefix,
    subject_id = subject_id,
    section_id = section_id,
    resolution = resolution,
    matrix = matrix
  )

  object@process_info$read_imzml <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "read_imzml",
    parameter = list(path_1 = path_1, peak_align_tolerance = peak_align_tolerance, path_2 = path_2),
    time = Sys.time()
  )

  if (!identical(path_2, "output")) {
    saveRDS(object, file = path_2)
  }

  object
}
