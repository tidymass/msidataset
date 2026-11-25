.imzml_check_number <- function(value, name, lower, upper, whole = FALSE) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < lower || value > upper || (whole && value != floor(value))) {
    stop("Invalid imzML argument: ", name, call. = FALSE)
  }
}

.imzml_align <- function(masses, intensities, xy, work, bin_digits, max_features,
                         presence_fraction, min_pixels, max_pixels_required,
                         noise_intensity_fraction, chunk_size) {
  n <- nrow(xy)
  if (length(masses) != n || length(intensities) != n) {
    stop("imzML coordinates and spectrum counts differ", call. = FALSE)
  }
  blocks <- .imzml_pixel_blocks(n, chunk_size)

  # Match the single-section reader: build the union at the requested
  # precision, then round the existing axis again if it exceeds the limit.
  mz <- numeric()
  for (ids in blocks) {
    block_mz <- unlist(masses[ids], use.names = FALSE)
    mz <- union(mz, unique(round(block_mz[is.finite(block_mz)], bin_digits)))
  }
  if (!length(mz)) stop("No finite m/z values found in processed imzML", call. = FALSE)
  digits <- bin_digits
  bin_counts <- stats::setNames(length(mz), digits)
  while (length(mz) > max_features && digits > 1L) {
    digits <- digits - 1L
    mz <- unique(round(mz, digits))
    bin_counts[as.character(digits)] <- length(mz)
  }
  mz <- sort(mz)

  # Count detection after summing duplicate bins within each pixel. Only
  # grouped pixel blocks are spooled to disk, so F1 never becomes dense.
  n_features <- length(mz)
  presence <- centroid_count <- numeric(n_features)
  invalid_by_pixel <- unmatched_by_pixel <- tic <- numeric(n)
  numerator <- denominator <- feature_max <- numeric(n_features)
  block_files <- file.path(work, paste0("pixels-", seq_along(blocks), ".rds"))
  for (i in seq_along(blocks)) {
    ids <- blocks[[i]]
    block_mz <- masses[ids]
    block_intensity <- intensities[ids]
    sizes <- lengths(block_mz)
    if (any(sizes != lengths(block_intensity))) {
      pixel <- ids[which(sizes != lengths(block_intensity))[1L]]
      stop("Unequal imzML m/z and intensity array lengths at pixel ", pixel, call. = FALSE)
    }
    block_index <- match(round(unlist(block_mz, use.names = FALSE), digits), mz)
    ends <- cumsum(sizes)
    grouped <- vector("list", length(ids))
    for (j in seq_along(ids)) {
      if (!sizes[j]) next
      pixel <- ids[j]
      raw_mz <- block_mz[[j]]
      intensity <- block_intensity[[j]]
      index <- block_index[seq.int(ends[j] - sizes[j] + 1, ends[j])]
      valid <- is.finite(raw_mz) & is.finite(intensity)
      invalid_by_pixel[pixel] <- sum(!valid)
      assigned <- valid & !is.na(index)
      unmatched_by_pixel[pixel] <- sum(valid & !assigned)
      index <- index[assigned]
      if (!length(index)) next
      signal <- intensity[assigned]
      values <- cbind(signal, raw_mz[assigned] * signal, 1)

      if (anyDuplicated(index)) {
        values <- rowsum(values, index, reorder = FALSE)
        index <- as.integer(rownames(values))
      }
      signal <- values[, 1L]
      tic[pixel] <- sum(signal)
      presence[index] <- presence[index] + (signal > 0)
      centroid_count[index] <- centroid_count[index] + values[, 3L]
      denominator[index] <- denominator[index] + signal
      numerator[index] <- numerator[index] + values[, 2L]
      if (noise_intensity_fraction > 0) {
        scale <- if (is.finite(tic[pixel]) && tic[pixel] > 0) tic[pixel] else 1
        feature_max[index] <- pmax(feature_max[index], signal / scale)
      }
      grouped[[j]] <- list(index = index, signal = signal)
    }
    saveRDS(grouped, block_files[i], compress = FALSE)
  }

  required <- min(max_pixels_required, max(min_pixels, ceiling(presence_fraction * n)))
  pass_presence <- presence >= required
  noise_threshold <- noise_intensity_fraction * max(feature_max)
  keep <- which(pass_presence & feature_max >= noise_threshold)
  .imzml_check_output_size(length(keep), n)
  lookup <- integer(n_features)
  lookup[keep] <- seq_along(keep)
  columns <- vector("list", n)
  for (i in seq_along(blocks)) {
    grouped <- readRDS(block_files[i])
    columns[blocks[[i]]] <- lapply(grouped, function(pixel) {
      column <- numeric(length(keep))
      index <- lookup[pixel$index]
      retained <- index > 0L
      column[index[retained]] <- pixel$signal[retained]
      column
    })
    unlink(block_files[i])
  }

  observed <- rep(NA_real_, length(keep))
  valid <- is.finite(numerator[keep]) & is.finite(denominator[keep]) & denominator[keep] != 0
  observed[valid] <- numerator[keep[valid]] / denominator[keep[valid]]
  list(columns = columns, xy = xy, tic = tic, mz = mz[keep],
       variables = data.frame(observed_mz_full = observed,
                              n_contributing_centroids = centroid_count[keep],
                              n_positive_pixels = presence[keep]),
       info = list(axis = "aligned", effective_bin_digits = digits, bin_counts = bin_counts,
                   n_features_f1 = n_features, n_features_after_presence = sum(pass_presence),
                   n_features_f2 = length(keep), presence_required_pixels = required,
                   noise_intensity_threshold = noise_threshold,
                   n_invalid_peaks = sum(invalid_by_pixel),
                   n_unmatched_peaks = sum(unmatched_by_pixel)))
}
