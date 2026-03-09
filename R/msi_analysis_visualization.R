# ============================================================================
# Core MSI analysis and visualization functions
# ============================================================================

.get_mz_index <- function(object, mz, tolerance_ppm = 20) {
  deltas <- abs(object@variable_info$mz - mz)
  tol <- mz * tolerance_ppm / 1e6
  idx <- which(deltas <= tol)
  if (length(idx) == 0) {
    idx <- which.min(deltas)
  }
  idx[[1]]
}

.intensity_df <- function(object, values) {
  data.frame(
    sample_id = object@sample_info$sample_id,
    x = as.numeric(object@sample_info$x),
    y = as.numeric(object@sample_info$y),
    value = as.numeric(values),
    stringsAsFactors = FALSE
  )
}

.as_grid_matrix <- function(df, value_col = "value") {
  x_vals <- sort(unique(df$x))
  y_vals <- sort(unique(df$y))
  mat <- matrix(NA_real_, nrow = length(y_vals), ncol = length(x_vals))
  colnames(mat) <- as.character(x_vals)
  rownames(mat) <- as.character(y_vals)

  for (i in seq_len(nrow(df))) {
    xi <- match(df$x[i], x_vals)
    yi <- match(df$y[i], y_vals)
    mat[yi, xi] <- df[[value_col]][i]
  }

  list(x = x_vals, y = y_vals, z = mat)
}

.point_in_polygon <- function(x, y, poly_x, poly_y) {
  n <- length(poly_x)
  inside <- rep(FALSE, length(x))
  j <- n
  for (i in seq_len(n)) {
    yi <- poly_y[i]
    yj <- poly_y[j]
    xi <- poly_x[i]
    xj <- poly_x[j]
    intersect <- ((yi > y) != (yj > y)) &
      (x < (xj - xi) * (y - yi) / (yj - yi + 1e-12) + xi)
    inside <- xor(inside, intersect)
    j <- i
  }
  inside
}

.plot_spatial_map <- function(df,
                              title = "Spatial map",
                              value_col = "value",
                              palette = "YlOrRd",
                              alpha = 1,
                              xlab = "x",
                              ylab = "y") {
  g <- .as_grid_matrix(df, value_col = value_col)
  z <- g$z
  z[is.na(z)] <- min(z, na.rm = TRUE)
  cols <- grDevices::hcl.colors(128, palette, rev = TRUE)
  if (!isTRUE(all.equal(alpha, 1))) {
    cols <- grDevices::adjustcolor(cols, alpha.f = alpha)
  }
  graphics::image(
    x = g$x,
    y = g$y,
    z = t(z),
    xlab = xlab,
    ylab = ylab,
    col = cols,
    main = title
  )
}

#' Ion image for one or multiple m/z values
#'
#' @param object A \code{msi_dataset} object.
#' @param mz Numeric vector of target m/z values.
#' @param tolerance_ppm Numeric tolerance in ppm.
#' @param overlays Number of overlays (<= length(mz)); if >1, averaged overlay
#'   image is shown.
#' @param color_scale Character palette name passed to
#'   \code{grDevices::hcl.colors()}.
#' @param alpha Numeric scalar between 0 and 1 controlling map transparency.
#' @param xlab Character x-axis label.
#' @param ylab Character y-axis label.
#' @param plot Logical.
#'
#' @return List with intensity data.
#' @importFrom methods is
#' @importFrom graphics image
#' @importFrom grDevices hcl.colors adjustcolor
#' @export
ion_image <- function(object,
                      mz,
                      tolerance_ppm = 20,
                      overlays = 1,
                      color_scale = "YlOrRd",
                      alpha = 1,
                      xlab = "x",
                      ylab = "y",
                      plot = TRUE) {
  stopifnot(methods::is(object, "msi_dataset"))
  stopifnot(length(mz) >= 1)

  overlays <- min(overlays, length(mz))
  idx <- vapply(mz, function(x) .get_mz_index(object, x, tolerance_ppm), integer(1))
  intensities <- object@expression_data[idx, , drop = FALSE]

  out <- vector("list", nrow(intensities))
  for (i in seq_len(nrow(intensities))) {
    out[[i]] <- .intensity_df(object, intensities[i, ])
  }
  names(out) <- paste0("mz_", signif(object@variable_info$mz[idx], 6))

  if (plot) {
    if (overlays == 1) {
      .plot_spatial_map(
        out[[1]],
        title = names(out)[1],
        palette = color_scale,
        alpha = alpha,
        xlab = xlab,
        ylab = ylab
      )
    } else {
      merged <- out[[1]]
      merged$value <- 0
      for (i in seq_len(overlays)) {
        v <- out[[i]]$value
        v <- v / (max(v, na.rm = TRUE) + 1e-12)
        merged$value <- merged$value + v
      }
      .plot_spatial_map(
        merged,
        title = paste0("Overlay(", overlays, " ions)"),
        palette = color_scale,
        alpha = alpha,
        xlab = xlab,
        ylab = ylab
      )
    }
  }

  object <- add_process_log(
    object,
    "ion_image()",
    list(mz = mz, tolerance_ppm = tolerance_ppm, overlays = overlays)
  )

  invisible(list(object = object, ion_data = out))
}

#' Mean spectrum across all pixels
#' @param object A \code{msi_dataset} object.
#' @param plot Logical.
#' @return Data frame with mz and mean intensity.
#' @importFrom methods is
#' @importFrom graphics plot
#' @export
mean_spectrum <- function(object, plot = TRUE) {
  stopifnot(methods::is(object, "msi_dataset"))
  mean_intensity <- rowMeans(object@expression_data, na.rm = TRUE)
  out <- data.frame(mz = object@variable_info$mz, intensity = mean_intensity)

  if (plot) {
    graphics::plot(out$mz, out$intensity, type = "l", col = "#0B4F6C",
                   xlab = "m/z", ylab = "Mean intensity",
                   main = "Mean spectrum")
  }

  object <- add_process_log(object, "mean_spectrum()", list(n_variables = nrow(out)))
  attr(out, "object") <- object
  out
}

#' Total ion current spatial map
#' @param object A \code{msi_dataset} object.
#' @param plot Logical.
#' @return Data frame with TIC values by pixel.
#' @importFrom methods is
#' @export
tic_map <- function(object, plot = TRUE) {
  stopifnot(methods::is(object, "msi_dataset"))
  tic <- colSums(object@expression_data, na.rm = TRUE)
  out <- .intensity_df(object, tic)
  if (plot) {
    .plot_spatial_map(out, title = "TIC map")
  }
  object <- add_process_log(object, "tic_map()", list())
  attr(out, "object") <- object
  out
}

#' Spectrum for one pixel
#' @param object A \code{msi_dataset} object.
#' @param sample_id Pixel/sample id.
#' @param plot Logical.
#' @return Data frame with mz and intensity.
#' @importFrom methods is
#' @importFrom graphics plot
#' @export
pixel_spectrum <- function(object, sample_id, plot = TRUE) {
  stopifnot(methods::is(object, "msi_dataset"))
  idx <- match(sample_id, object@sample_info$sample_id)
  if (is.na(idx)) {
    stop("sample_id not found.")
  }

  out <- data.frame(
    mz = object@variable_info$mz,
    intensity = as.numeric(object@expression_data[, idx]),
    stringsAsFactors = FALSE
  )

  if (plot) {
    graphics::plot(out$mz, out$intensity, type = "l", col = "#1B998B",
                   xlab = "m/z", ylab = "Intensity",
                   main = paste("Pixel spectrum:", sample_id))
  }

  object <- add_process_log(object, "pixel_spectrum()", list(sample_id = sample_id))
  attr(out, "object") <- object
  out
}

#' Visualize region of interest
#' @param object A \code{msi_dataset} object.
#' @param type ROI type: rectangle, polygon, freehand.
#' @param roi List describing ROI coordinates.
#' @param plot Logical.
#' @param save_roi Save selected ROI into object ROI subsystem.
#' @param roi_id ROI identifier used when \code{save_roi = TRUE}.
#' @param roi_name ROI display name.
#' @param overwrite Overwrite an existing ROI with the same id.
#' @return Selected sample_info rows.
#' @importFrom methods is
#' @importFrom graphics points
#' @export
roi_visualize <- function(object,
                          type = c("rectangle", "polygon", "freehand"),
                          roi,
                          plot = TRUE,
                          save_roi = FALSE,
                          roi_id = NULL,
                          roi_name = roi_id,
                          overwrite = FALSE) {
  stopifnot(methods::is(object, "msi_dataset"))
  type <- match.arg(type)

  coords <- object@sample_info[, c("sample_id", "x", "y")]
  inside <- rep(FALSE, nrow(coords))

  if (type == "rectangle") {
    stopifnot(all(c("xmin", "xmax", "ymin", "ymax") %in% names(roi)))
    inside <- coords$x >= roi$xmin & coords$x <= roi$xmax &
      coords$y >= roi$ymin & coords$y <= roi$ymax
  } else {
    stopifnot(all(c("x", "y") %in% names(roi)))
    inside <- .point_in_polygon(coords$x, coords$y, roi$x, roi$y)
  }

  selected <- coords[inside, , drop = FALSE]

  if (isTRUE(save_roi)) {
    if (is.null(roi_id)) {
      stop("roi_id must be provided when save_roi = TRUE.")
    }
    object <- add_roi(
      object = object,
      roi_id = roi_id,
      roi_name = roi_name,
      sample_id = selected$sample_id,
      overwrite = overwrite
    )
  }

  if (plot) {
    tic <- colSums(object@expression_data, na.rm = TRUE)
    df <- .intensity_df(object, tic)
    .plot_spatial_map(df, title = paste("ROI:", type))
    graphics::points(selected$x, selected$y, pch = 16, col = "cyan")
  }

  object <- add_process_log(
    object,
    "roi_visualize()",
    list(
      type = type,
      n_selected = nrow(selected),
      save_roi = save_roi,
      roi_id = roi_id
    )
  )
  attr(selected, "object") <- object
  selected
}

#' Multi-m/z overlay image
#' @param object A \code{msi_dataset} object.
#' @param mz Numeric vector (2-3 recommended).
#' @param alpha Numeric alpha between 0 and 1.
#' @param tolerance_ppm Numeric ppm tolerance.
#' @return RGB array invisibly.
#' @importFrom methods is
#' @importFrom grDevices rgb
#' @importFrom graphics plot rasterImage
#' @export
overlay_image <- function(object, mz, alpha = 0.6, tolerance_ppm = 20) {
  stopifnot(methods::is(object, "msi_dataset"))
  stopifnot(length(mz) >= 2)

  idx <- vapply(mz, function(x) .get_mz_index(object, x, tolerance_ppm), integer(1))
  channels <- lapply(seq_along(idx), function(i) {
    .as_grid_matrix(.intensity_df(object, object@expression_data[idx[i], ]))$z
  })

  nr <- nrow(channels[[1]])
  nc <- ncol(channels[[1]])
  rgb_arr <- array(0, dim = c(nr, nc, 3))
  for (i in seq_len(min(3, length(channels)))) {
    v <- channels[[i]]
    v <- v / (max(v, na.rm = TRUE) + 1e-12)
    rgb_arr[, , i] <- v
  }

  img <- grDevices::rgb(rgb_arr[, , 1], rgb_arr[, , 2], rgb_arr[, , 3], alpha = alpha)
  graphics::plot(1, type = "n", xlab = "", ylab = "", axes = FALSE,
                 xlim = c(0, 1), ylim = c(0, 1), main = "overlay_image")
  graphics::rasterImage(as.raster(matrix(img, nrow = nr, ncol = nc)), 0, 0, 1, 1)

  object <- add_process_log(object, "overlay_image()", list(mz = mz, alpha = alpha))
  invisible(list(object = object, rgb = rgb_arr))
}

#' Pixel similarity map
#' @param object A \code{msi_dataset} object.
#' @param sample_id Reference sample id.
#' @param method Correlation method.
#' @param plot Logical.
#' @return Data frame with similarity scores.
#' @importFrom methods is
#' @importFrom stats cor
#' @export
pixel_similarity_map <- function(object,
                                 sample_id,
                                 method = c("pearson", "spearman"),
                                 plot = TRUE) {
  stopifnot(methods::is(object, "msi_dataset"))
  method <- match.arg(method)
  idx <- match(sample_id, object@sample_info$sample_id)
  if (is.na(idx)) {
    stop("sample_id not found.")
  }

  ref <- as.numeric(object@expression_data[, idx])
  sim <- apply(object@expression_data, 2, function(v) {
    suppressWarnings(stats::cor(ref, as.numeric(v), method = method, use = "pairwise.complete.obs"))
  })
  out <- .intensity_df(object, sim)

  if (plot) {
    .plot_spatial_map(out, title = paste("Similarity to", sample_id))
  }
  object <- add_process_log(object, "pixel_similarity_map()", list(sample_id = sample_id, method = method))
  attr(out, "object") <- object
  out
}

#' Mass spectrum coverage map
#' @param object A \code{msi_dataset} object.
#' @param threshold Intensity threshold for detected peaks.
#' @param plot Logical.
#' @return Data frame of coverage ratio per pixel.
#' @importFrom methods is
#' @export
coverage_map <- function(object, threshold = 0, plot = TRUE) {
  stopifnot(methods::is(object, "msi_dataset"))
  coverage <- colMeans(object@expression_data > threshold, na.rm = TRUE)
  out <- .intensity_df(object, coverage)
  if (plot) {
    .plot_spatial_map(out, title = "Coverage map")
  }
  object <- add_process_log(object, "coverage_map()", list(threshold = threshold))
  attr(out, "object") <- object
  out
}

#' Spatial clustering visualization
#' @param object A \code{msi_dataset} object.
#' @param k Number of clusters.
#' @param n_pcs Number of PCs used for clustering.
#' @param plot Logical.
#' @return Data frame with cluster labels.
#' @importFrom methods is
#' @importFrom stats prcomp kmeans
#' @export
spatial_clustering <- function(object, k = 4, n_pcs = 10, plot = TRUE) {
  stopifnot(methods::is(object, "msi_dataset"))
  x <- t(object@expression_data)
  pca <- stats::prcomp(x, center = TRUE, scale. = TRUE)
  emb <- pca$x[, seq_len(min(n_pcs, ncol(pca$x))), drop = FALSE]
  cl <- stats::kmeans(emb, centers = k)
  out <- .intensity_df(object, cl$cluster)
  names(out)[names(out) == "value"] <- "cluster"

  if (plot) {
    .plot_spatial_map(transform(out, value = cluster), title = paste("Spatial clustering (k=", k, ")", sep = ""))
  }
  object <- add_process_log(object, "spatial_clustering()", list(k = k, n_pcs = n_pcs))
  attr(out, "object") <- object
  out
}

#' Embedding visualization (PCA/t-SNE/UMAP)
#' @param object A \code{msi_dataset} object.
#' @param method Embedding method.
#' @param dims Number of output dimensions.
#' @param color_by One of: tic, section_id, subject_id.
#' @return Embedding data frame.
#' @importFrom methods is
#' @importFrom stats prcomp
#' @importFrom graphics plot
#' @export
embedding_visualize <- function(object,
                                method = c("PCA", "TSNE", "UMAP"),
                                dims = 2,
                                color_by = c("tic", "section_id", "subject_id")) {
  stopifnot(methods::is(object, "msi_dataset"))
  method <- match.arg(method)
  color_by <- match.arg(color_by)

  x <- t(object@expression_data)

  if (method == "PCA") {
    fit <- stats::prcomp(x, center = TRUE, scale. = TRUE)
    emb <- fit$x[, seq_len(dims), drop = FALSE]
  } else if (method == "TSNE") {
    if (!requireNamespace("Rtsne", quietly = TRUE)) {
      stop("Please install Rtsne for method = 'TSNE'.")
    }
    emb <- Rtsne::Rtsne(x, dims = dims)$Y
  } else {
    if (!requireNamespace("uwot", quietly = TRUE)) {
      stop("Please install uwot for method = 'UMAP'.")
    }
    emb <- uwot::umap(x, n_components = dims)
  }

  emb <- as.data.frame(emb, stringsAsFactors = FALSE)
  colnames(emb) <- paste0("dim", seq_len(ncol(emb)))
  emb$sample_id <- object@sample_info$sample_id

  if (color_by == "tic") {
    emb$color <- colSums(object@expression_data, na.rm = TRUE)
  } else {
    emb$color <- object@sample_info[[color_by]]
  }

  graphics::plot(emb$dim1, emb$dim2,
                 col = as.factor(emb$color), pch = 16,
                 xlab = "dim1", ylab = "dim2",
                 main = paste(method, "embedding"))

  object <- add_process_log(object, "embedding_visualize()", list(method = method, dims = dims, color_by = color_by))
  attr(emb, "object") <- object
  emb
}

#' MS2 spatial map
#' @param object A \code{msi_dataset} object.
#' @param mode One of binary or count.
#' @param plot Logical.
#' @return Data frame with MS2-derived values by pixel.
#' @importFrom methods is
#' @export
ms2_spatial_map <- function(object, mode = c("binary", "count"), plot = TRUE) {
  stopifnot(methods::is(object, "msi_dataset"))
  mode <- match.arg(mode)

  ms2_flags <- if ("has_ms2" %in% colnames(object@variable_info)) {
    object@variable_info$has_ms2 %in% c("Yes", "yes", TRUE, 1)
  } else {
    object@variable_info$ms2 %in% c("Yes", "yes", TRUE, 1)
  }
  if (!any(ms2_flags)) {
    out <- .intensity_df(object, rep(0, ncol(object@expression_data)))
  } else {
    ms2_expr <- object@expression_data[ms2_flags, , drop = FALSE]
    values <- if (mode == "count") {
      colSums(ms2_expr > 0, na.rm = TRUE)
    } else {
      as.numeric(colSums(ms2_expr > 0, na.rm = TRUE) > 0)
    }
    out <- .intensity_df(object, values)
  }

  if (plot) {
    .plot_spatial_map(out, title = paste("MS2 spatial map:", mode))
  }

  object <- add_process_log(object, "ms2_spatial_map()", list(mode = mode))
  attr(out, "object") <- object
  out
}
