#' Generate spatially correlated pseudo-MS2 fragments
#'
#' Apply masked Gaussian smoothing to mean-TIC normalized pixel intensities,
#' then select fragment/precursor pairs by Spearman correlation and m/z gap.
#' Accept a mass_dataset or matrices.
#' For a mass_dataset, use measured centroids when available and otherwise
#' use the feature axis.
#'
#' @param xy A mass_dataset or a two-column matrix of integer pixel coordinates.
#' @param tic Positive TIC for each pixel.
#' @param intensity Pixel-by-feature numeric intensity matrix.
#' @param mz Feature m/z values in column order for matrix input.
#' @param feature_id Feature identifiers in column order.
#' @return A tibble of qualifying pairs and their Spearman and spatial association scores;
#'   parameters are recorded in its `process_info` attribute. For a mass_dataset
#'   input, return the object with results and parameters in process_info.
#' @examples
#' \dontrun{generate_pseudo_ms2(xy, tic, intensity, mz, feature_id)}
#' @seealso [detect_de_features()], [stats::cor()]
#' @export
generate_pseudo_ms2 <- function(xy, tic = NULL, intensity = NULL, mz = NULL, feature_id = NULL) {
  if (methods::is(xy, "mass_dataset")) {
    object <- xy
    i <- .msidataset_inputs(object)
    if (is.null(i$mz)) stop("Feature m/z values are required", call. = FALSE)
    i$mz <- .msidataset_feature_mz(object@variable_info)$mz
    bad <- .msidataset_input_problems(i, require_mz = TRUE, integer_xy = TRUE)
    if (length(bad)) return(.msidataset_append_step(object, "generate_pseudo_ms2",
      list(status = "UNRESOLVED_INPUT", problems = bad), list(status = "UNRESOLVED_INPUT")))
    result <- generate_pseudo_ms2(i$xy, i$tic, i$intensity, i$mz, i$feature_id)
    return(.msidataset_append_step(object, "generate_pseudo_ms2", result))
  }
  xy <- as.matrix(xy); intensity <- as.matrix(intensity); n <- nrow(xy)
  if (ncol(xy) != 2L || nrow(intensity) != n || length(tic) != n ||
      ncol(intensity) != length(mz) || length(mz) != length(feature_id) ||
      any(!is.finite(xy)) || any(xy != floor(xy)) || anyDuplicated(data.frame(xy)) ||
      any(!is.finite(tic) | tic <= 0) || any(!is.finite(intensity)) || any(!is.finite(mz)))
    stop("Invalid pixel inputs", call. = FALSE)
  xy <- sweep(xy, 2L, apply(xy, 2L, min), "-") + 1L
  dims <- apply(xy, 2L, max); ind <- as.integer(xy[, 1] + (xy[, 2] - 1) * dims[1])
  offsets <- -4:4; kernel <- exp(-offsets^2 / 2); kernel <- kernel / sum(kernel)
  reflect <- function(idx, size) {
    q <- (idx - 1L) %% (2L * size)
    ifelse(q >= size, 2L * size - q - 1L, q) + 1L
  }
  smooth <- function(img) {
    a <- matrix(0, dims[1], dims[2]); b <- a
    for (k in seq_along(offsets)) a <- a + kernel[k] * img[reflect(seq_len(dims[1]) + offsets[k], dims[1]), , drop = FALSE]
    for (k in seq_along(offsets)) b <- b + kernel[k] * a[, reflect(seq_len(dims[2]) + offsets[k], dims[2]), drop = FALSE]
    b
  }
  mask <- matrix(0, dims[1], dims[2]); mask[ind] <- 1
  den <- smooth(mask)[ind]
  values <- matrix(0, n, ncol(intensity)); ranked <- values; constant <- character()
  for (j in seq_len(ncol(intensity))) {
    img <- matrix(0, dims[1], dims[2]); img[ind] <- intensity[, j] * (mean(tic) / tic)
    v <- smooth(img)[ind] / den; values[, j] <- v
    r <- rank(v, ties.method = "average"); r <- r - mean(r); norm <- sqrt(sum(r * r))
    if (norm > 0) ranked[, j] <- r / norm else constant <- c(constant, feature_id[j])
  }
  corr <- crossprod(ranked); pairs <- which(outer(mz, mz, "-") >= 10 & is.finite(corr) & corr >= 0.8, arr.ind = TRUE)
  if (nrow(pairs)) pairs <- pairs[order(pairs[,1], pairs[,2]), , drop = FALSE]
  sas <- vapply(seq_len(nrow(pairs)), function(i) {
    x <- values[, pairs[i, 1]]; y <- values[, pairs[i, 2]]; u <- x > 0 | y > 0
    if (!any(u)) return(NA_real_)
    rho <- suppressWarnings(stats::cor(x[u], y[u], method = "spearman"))
    if (!is.finite(rho)) return(NA_real_)
    0.5 * rho + 0.5 * sum(x > 0 & y > 0) / sum(u)
  }, numeric(1))
  a <- pairs[, 1]; b <- pairs[, 2]
  ans <- tibble::tibble(sas_d7 = sas, precursor_feature_id = feature_id[a], precursor_mz = mz[a],
    fragment_feature_id = feature_id[b], fragment_mz = mz[b], spearman_rho = corr[pairs])
  attr(ans, "process_info") <- list(function_name = "generate_pseudo_ms2", sigma_pixels = 1,
    min_rho = 0.8, mz_gap_da = 10, constant_features = constant, time = as.character(Sys.time()))
  ans
}
