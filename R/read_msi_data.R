#' Read pixel-table CSV or imzML MSI data
#'
#' Convert a headerless pixel table (x, y, TIC, then features) or an imzML file
#' into a massdataset object. Continuous imzML retains its native axis.
#' For processed imzML, the default axis sums native intensities in rounded
#' m/z bins within each pixel and keeps bins that pass the presence filter.
#' For that axis, feature metadata stores bin labels in `mz` and
#' intensity-weighted native centroids in `observed_mz_full`.
#' Sample TIC is the sum across all aligned bins before filtering.
#'
#' @param path Path to the CSV or imzML file.
#' @param format Either `"csv"` or `"imzml"`.
#' @param mz Optional numeric vector of feature m/z values for CSV, in order.
#' @param feature_id Optional CSV feature identifiers; defaults to F003 onward.
#' @param axis Either `"aligned"` (native processed spectra) or
#'   `"cardinal_grid"` (Cardinal's common grid). Continuous imzML uses its
#'   native axis in either mode.
#' @param bin_digits Initial decimal precision of processed m/z bins.
#' @param max_features If the initial axis exceeds this many bins, decrease
#'   precision one decimal place at a time, stopping at one decimal place.
#' @param presence_fraction Fraction of pixels required to have a positive
#'   summed intensity in a bin.
#' @param min_pixels Absolute floor on the number of positive pixels.
#' @param max_pixels_required Cap on the required positive pixel count. The
#'   requirement is `min(max_pixels_required, max(min_pixels,
#'   ceiling(presence_fraction * n_pixels)))`.
#' @param noise_intensity_fraction Optional legacy intensity QC. A retained
#'   bin's maximum TIC-normalized intensity must reach this fraction of the
#'   global maximum TIC-normalized intensity. Zero disables intensity QC.
#'   The returned expression values remain on their original scale.
#' @param chunk_size Maximum number of imzML pixels read as a block.
#' @return A `mass_dataset` with pixels as samples and process history.
#' @examples
#' \dontrun{msi <- read_msi_data("roi.csv", format = "csv")}
#' @seealso [massdataset::create_mass_dataset()]
#' @export
read_msi_data <- function(path, format = c("csv", "imzml"), mz = NULL, feature_id = NULL,
                          axis = c("aligned", "cardinal_grid"), bin_digits = 4L,
                          max_features = 150000L, presence_fraction = 0.30,
                          min_pixels = 8L, max_pixels_required = 1000L,
                          noise_intensity_fraction = 0, chunk_size = 256L) {
  format <- match.arg(format)
  if (length(path) != 1L || !file.exists(path)) stop("Input file not found", call. = FALSE)
  imzml_info <- NULL
  feature_metadata <- NULL
  if (format == "csv") {
    raw <- data.table::fread(path, header = FALSE, data.table = FALSE, showProgress = FALSE)
    if (ncol(raw) < 4L) stop("CSV must contain coordinates, TIC and features", call. = FALSE)
    xy <- as.matrix(raw[, 1:2]); tic <- raw[[3]]
    expression <- t(as.matrix(raw[, -(1:3), drop = FALSE]))
    if (is.null(feature_id)) feature_id <- sprintf("F%03d", seq_len(nrow(expression)) + 2L)
  } else {
    axis <- match.arg(axis)
    .imzml_check_number(chunk_size, "chunk_size", 1, .Machine$integer.max, whole = TRUE)
    if (axis == "aligned") {
      .imzml_check_number(bin_digits, "bin_digits", 1, 15, whole = TRUE)
      .imzml_check_number(max_features, "max_features", 1, .Machine$integer.max, whole = TRUE)
      .imzml_check_number(presence_fraction, "presence_fraction", 0, 1)
      .imzml_check_number(min_pixels, "min_pixels", 0, .Machine$integer.max, whole = TRUE)
      .imzml_check_number(max_pixels_required, "max_pixels_required", min_pixels,
                          .Machine$integer.max, whole = TRUE)
      .imzml_check_number(noise_intensity_fraction, "noise_intensity_fraction", 0, 1)
    }
    loaded <- .imzml_read_pair(path, axis, chunk_size, bin_digits = bin_digits,
      max_features = max_features, presence_fraction = presence_fraction, min_pixels = min_pixels,
      max_pixels_required = max_pixels_required, noise_intensity_fraction = noise_intensity_fraction)
    xy <- loaded$xy
    tic <- loaded$tic
    mz <- loaded$mz
    feature_metadata <- loaded$variables
    imzml_info <- c(loaded$info, list(requested_axis = axis, bin_digits = bin_digits,
      max_features = max_features, presence_fraction = presence_fraction, min_pixels = min_pixels,
      max_pixels_required = max_pixels_required, noise_intensity_fraction = noise_intensity_fraction,
      chunk_size = chunk_size))
    if (is.null(feature_id)) feature_id <- sprintf("F%03d", seq_along(mz))
  }
  n_features <- if (format == "csv") nrow(expression) else length(mz)
  n_samples <- if (format == "csv") ncol(expression) else length(loaded$columns)
  if (n_samples != nrow(xy) || length(feature_id) != n_features ||
      anyNA(feature_id) || anyDuplicated(feature_id) || (!is.null(mz) && length(mz) != n_features))
    stop("Dimensions or feature IDs are invalid", call. = FALSE)
  sample_id <- sprintf("pixel_%d", seq_len(nrow(xy)))
  if (format == "csv") {
    colnames(expression) <- sample_id
    rownames(expression) <- feature_id
    expr <- as.data.frame(expression, check.names = FALSE)
  } else {
    expr <- structure(loaded$columns, names = sample_id, row.names = as.character(feature_id),
                      class = "data.frame")
  }
  samples <- data.frame(sample_id = sample_id, class = "unassigned", x = xy[, 1], y = xy[, 2], tic = tic)
  variables <- data.frame(variable_id = feature_id, mz = if (is.null(mz)) rep(NA_real_, length(feature_id)) else mz)
  if (!is.null(feature_metadata)) variables <- cbind(variables, feature_metadata)
  result <- massdataset::create_mass_dataset(expr, samples, variables)
  result@process_info$read_msi_data <- c(list(function_name = "read_msi_data", format = format,
    input = .msi_file_record(path), time = as.character(Sys.time())), imzml_info)
  result
}
