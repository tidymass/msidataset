#' @keywords internal
.build_spatial_matrix <- function(object, values) {
  si <- object@sample_info
  ux <- sort(unique(si$x))
  uy <- sort(unique(si$y))
  m <- matrix(NA_real_, nrow = length(uy), ncol = length(ux),
              dimnames = list(rev(uy), ux))
  for (i in seq_along(values)) {
    xi <- match(si$x[i], ux)
    yi <- length(uy) - match(si$y[i], uy) + 1
    m[yi, xi] <- values[i]
  }
  m
}

#' @export
plot_ion_image <- function(object, mz, tolerance = 0.01, col = hcl.colors(64, "YlOrRd", rev = TRUE)) {
  idx <- which(abs(object@variable_info$mz - mz) <= tolerance)
  if (length(idx) == 0) stop("No feature found for input mz.")
  vals <- colMeans(as.matrix(object@expression_data[idx, , drop = FALSE]))
  mat <- .build_spatial_matrix(object, vals)
  graphics::image(x = seq_len(ncol(mat)), y = seq_len(nrow(mat)), z = t(mat[nrow(mat):1, ]), col = col,
                  xlab = "x", ylab = "y", main = paste0("Ion image m/z=", signif(mz, 6)))
  invisible(mat)
}

#' @export
plot_multi_ion_image <- function(object, mz, tolerance = 0.01, alpha = 0.5) {
  if (length(mz) < 2) stop("Please provide >=2 mz values.")
  mats <- lapply(mz, function(m) .build_spatial_matrix(object, colMeans(as.matrix(object@expression_data[which(abs(object@variable_info$mz - m) <= tolerance), , drop = FALSE]))))
  base <- mats[[1]] / max(mats[[1]], na.rm = TRUE)
  rgb_arr <- array(0, dim = c(nrow(base), ncol(base), 3))
  for (k in seq_along(mats)) {
    cur <- mats[[k]] / max(mats[[k]], na.rm = TRUE)
    ch <- ((k - 1) %% 3) + 1
    rgb_arr[, , ch] <- pmin(1, rgb_arr[, , ch] + alpha * cur)
  }
  img <- rgb(rgb_arr[, , 1], rgb_arr[, , 2], rgb_arr[, , 3])
  graphics::plot.new(); graphics::rasterImage(as.raster(matrix(img, nrow(base), ncol(base))), 0, 0, 1, 1)
  title(main = "Multi m/z overlay")
  invisible(rgb_arr)
}

#' @export
plot_mean_spectrum <- function(object) {
  mean_int <- rowMeans(as.matrix(object@expression_data))
  graphics::plot(object@variable_info$mz, mean_int, type = "l", xlab = "m/z", ylab = "Mean intensity", main = "Mean Spectrum")
  invisible(data.frame(mz = object@variable_info$mz, mean_intensity = mean_int))
}

#' @export
plot_tic_image <- function(object) {
  tic <- colSums(as.matrix(object@expression_data), na.rm = TRUE)
  mat <- .build_spatial_matrix(object, tic)
  graphics::image(t(mat[nrow(mat):1, ]), col = hcl.colors(64, "Inferno"), xlab = "x", ylab = "y", main = "TIC spatial map")
  invisible(mat)
}

#' @export
plot_pixel_spectrum <- function(object, sample_id) {
  idx <- match(sample_id, object@sample_info$sample_id)
  if (is.na(idx)) stop("sample_id not found.")
  graphics::plot(object@variable_info$mz, object@expression_data[, idx], type = "l", xlab = "m/z", ylab = "Intensity",
                 main = paste("Pixel spectrum", sample_id))
}

#' @export
plot_roi <- function(object, roi_samples, border = "cyan") {
  graphics::plot(object@sample_info$x, object@sample_info$y, pch = 15, col = "grey80", asp = 1, xlab = "x", ylab = "y")
  idx <- object@sample_info$sample_id %in% roi_samples
  graphics::points(object@sample_info$x[idx], object@sample_info$y[idx], pch = 15, col = border)
}

#' @export
plot_similarity_map <- function(object, sample_id, method = "pearson") {
  idx <- match(sample_id, object@sample_info$sample_id)
  if (is.na(idx)) stop("sample_id not found.")
  ref <- as.numeric(object@expression_data[, idx])
  sims <- apply(object@expression_data, 2, function(x) stats::cor(ref, x, method = method, use = "pairwise.complete.obs"))
  mat <- .build_spatial_matrix(object, sims)
  graphics::image(t(mat[nrow(mat):1, ]), col = hcl.colors(64, "Blue-Red 3"), main = "Similarity map")
  invisible(mat)
}

#' @export
plot_spectral_coverage <- function(object, threshold = 0) {
  cov <- colMeans(as.matrix(object@expression_data) > threshold, na.rm = TRUE)
  mat <- .build_spatial_matrix(object, cov)
  graphics::image(t(mat[nrow(mat):1, ]), col = hcl.colors(64, "Viridis"), main = "Spectral coverage")
  invisible(mat)
}

#' @export
plot_spatial_cluster <- function(object, k = 4) {
  km <- stats::kmeans(t(scale(as.matrix(object@expression_data))), centers = k)
  cols <- grDevices::rainbow(k)[km$cluster]
  graphics::plot(object@sample_info$x, object@sample_info$y, pch = 15, col = cols, asp = 1,
                 xlab = "x", ylab = "y", main = "Spatial clustering")
  invisible(km)
}

#' @export
plot_embedding <- function(object, method = c("pca", "tsne")) {
  method <- match.arg(method)
  x <- t(scale(as.matrix(object@expression_data)))
  if (method == "pca") {
    emb <- stats::prcomp(x)$x[, 1:2, drop = FALSE]
  } else {
    if (!requireNamespace("Rtsne", quietly = TRUE)) stop("Install Rtsne for tsne embedding.")
    emb <- Rtsne::Rtsne(x, check_duplicates = FALSE)$Y
  }
  graphics::plot(emb[, 1], emb[, 2], pch = 16, xlab = "Dim1", ylab = "Dim2", main = paste("Embedding:", method))
  invisible(emb)
}

#' @export
plot_ms2_spatial <- function(object, precursor_mz, tolerance = 0.01) {
  idx <- which(abs(object@variable_info$mz - precursor_mz) <= tolerance)
  if (length(idx) == 0) stop("No precursor found in variable_info.")
  vals <- colMeans(as.matrix(object@expression_data[idx, , drop = FALSE]))
  mat <- .build_spatial_matrix(object, vals)
  graphics::image(t(mat[nrow(mat):1, ]), col = hcl.colors(64, "Magma"), main = "MS2 spatial visualization")
  invisible(mat)
}
