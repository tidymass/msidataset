# Cardinal canonicalizes paths before pairing the XML and IBD files.
.imzml_read_pair <- function(path, axis, chunk_size, ...) {
  path <- normalizePath(path, mustWork = TRUE)
  stem <- tools::file_path_sans_ext(path)
  binary <- paste0(stem, c(".ibd", ".IBD"))
  binary <- binary[file.exists(binary)][1L]
  if (is.na(binary)) stop("Paired imzML .ibd file not found", call. = FALSE)

  work <- tempfile("msidataset-imzml-")
  if (!dir.create(work)) stop("Cannot create imzML work directory", call. = FALSE)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)

  read <- function(file) {
    if (axis == "aligned") {
      parsed <- Cardinal::readMSIData(file, parse.only = TRUE)
      content <- names(parsed[["fileDescription"]][["fileContent"]])
      if ("IMS:1000031" %in% content) {
        spectra <- parsed[["run"]][["spectrumList"]]
        positions <- spectra[["positions"]]
        xy <- cbind(x = as.numeric(positions[["position x"]]),
                    y = as.numeric(positions[["position y"]]))
        arrays <- parsed[["ibd"]]
        return(.imzml_align(arrays[["mz"]], arrays[["intensity"]], xy,
                            work, chunk_size = chunk_size, ...))
      }
    }
    .imzml_cardinal(file, chunk_size)
  }

  canonical <- identical(tools::file_ext(path), "imzML") && endsWith(binary, ".ibd")
  if (canonical) {
    result <- tryCatch(read(path), error = function(e) {
      if (!.imzml_open_error(e)) stop(e)
      NULL
    })
    if (!is.null(result)) {
      result$info$staging <- "original"
      return(result)
    }
  }

  source <- c(path, binary)
  target <- paste0(file.path(work, basename(stem)), c(".imzML", ".ibd"))
  for (method in c("symlink", "hardlink", "copy")) {
    staged <- switch(method,
      symlink = suppressWarnings(file.symlink(source, target)),
      hardlink = suppressWarnings(file.link(source, target)),
      copy = file.copy(source, target)
    )
    # A symlink to a differently cased/unpaired XML resolves to the old name.
    usable <- all(staged) && (method != "symlink" || canonical)
    if (usable) {
      result <- tryCatch(read(target[1L]), error = function(e) {
        if (method == "copy" || !.imzml_open_error(e)) stop(e)
        NULL
      })
      if (!is.null(result)) {
        result$info$staging <- method
        return(result)
      }
    }
    unlink(target)
  }
  stop("Cannot stage imzML file pair", call. = FALSE)
}

.imzml_open_error <- function(error) {
  grepl("could not open file", conditionMessage(error), fixed = TRUE)
}

.imzml_pixel_blocks <- function(n, chunk_size) {
  if (n == 0L) stop("imzML contains no pixel spectra", call. = FALSE)
  starts <- seq.int(1, n, by = chunk_size)
  lapply(starts, function(start) seq.int(start, min(n, start + chunk_size - 1)))
}

.imzml_check_output_size <- function(n_features, n_pixels) {
  # mass_dataset stores one dense double vector per pixel.
  bytes <- 8 * as.double(n_features) * n_pixels
  if (bytes > 8e9) {
    stop("The retained imzML expression table exceeds 8 GB. ",
         "Use a stricter presence filter or a smaller input.", call. = FALSE)
  }
}

.imzml_cardinal <- function(path, chunk_size) {
  image <- Cardinal::readMSIData(path)
  spectra <- Cardinal::spectra(image)
  xy <- as.matrix(Cardinal::coord(image))[, 1:2, drop = FALSE]
  mz <- as.numeric(Cardinal::mz(image))
  if (length(mz) != nrow(spectra)) stop("A common m/z axis is required", call. = FALSE)
  .imzml_check_output_size(length(mz), nrow(xy))
  # The two-dimensional extraction materializes matter sparse_mat correctly.
  chunk_size <- min(chunk_size, max(1, floor(1e7 / max(1, length(mz)))))
  columns <- vector("list", nrow(xy))
  tic <- numeric(nrow(xy))
  for (ids in .imzml_pixel_blocks(nrow(xy), chunk_size)) {
    block <- as.matrix(spectra[, ids, drop = FALSE])
    tic[ids] <- colSums(block)
    columns[ids] <- lapply(seq_along(ids), function(i) block[, i])
  }
  list(columns = columns, xy = xy, tic = tic, mz = mz, variables = NULL,
       info = list(axis = if (inherits(image, "MSProcessedImagingExperiment"))
         "Cardinal_processed_grid" else "native"))
}
