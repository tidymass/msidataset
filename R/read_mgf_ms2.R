#' @title Read MGF and attach MS2 to msi_dataset
#' @param object msi_dataset object.
#' @param mgf_path MGF file path.
#' @param mz_tolerance matching tolerance (Da).
#' @return updated msi_dataset object.
#' @export
read_mgf_ms2 <- function(object, mgf_path, mz_tolerance = 0.01) {
  if (!methods::is(object, "msi_dataset")) stop("object must be a msi_dataset.")
  if (!file.exists(mgf_path)) stop("mgf_path does not exist.")

  lines <- readLines(mgf_path, warn = FALSE)
  idx_begin <- which(lines == "BEGIN IONS")
  idx_end <- which(lines == "END IONS")
  if (length(idx_begin) != length(idx_end)) stop("Invalid MGF format.")

  ms2_list <- vector("list", length(idx_begin))
  names(ms2_list) <- paste0("ms2_", seq_along(ms2_list))

  for (i in seq_along(idx_begin)) {
    block <- lines[(idx_begin[i] + 1):(idx_end[i] - 1)]
    pepmass_line <- block[grepl("^PEPMASS=", block)][1]
    precursor_mz <- as.numeric(strsplit(sub("^PEPMASS=", "", pepmass_line), "\\s+")[[1]][1])
    peak_lines <- block[grepl("^[0-9]", block)]
    peak_mat <- do.call(rbind, strsplit(peak_lines, "\\s+"))
    ms2_list[[i]] <- list(
      precursor_mz = precursor_mz,
      mz = as.numeric(peak_mat[, 1]),
      intensity = as.numeric(peak_mat[, 2])
    )
  }

  object@ms2_data <- ms2_list
  vi <- object@variable_info
  vi$ms2 <- "no"
  for (i in seq_len(nrow(vi))) {
    precursors <- vapply(ms2_list, function(x) x$precursor_mz, numeric(1))
    if (any(abs(precursors - vi$mz[i]) <= mz_tolerance)) {
      vi$ms2[i] <- "yes"
    }
  }
  object@variable_info <- vi

  object@process_info$read_mgf_ms2 <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "read_mgf_ms2",
    parameter = list(mgf_path = mgf_path, mz_tolerance = mz_tolerance, n_ms2 = length(ms2_list)),
    time = Sys.time()
  )

  methods::validObject(object)
  object
}
