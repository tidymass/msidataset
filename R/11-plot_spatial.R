# ============================================================================
# Unified spatial plotting API (ggplot2-friendly)
# ============================================================================

.resolve_feature_index <- function(object,
                                   feature = NULL,
                                   feature_type = c("mz", "variable_id", "tic"),
                                   tolerance_ppm = 20) {
  feature_type <- match.arg(feature_type)

  if (feature_type == "tic") {
    return(list(
      type = "tic",
      index = NA_integer_,
      label = "TIC"
    ))
  }

  if (is.null(feature) || length(feature) == 0) {
    stop("feature must be provided for feature_type = 'mz' or 'variable_id'.")
  }

  if (feature_type == "mz") {
    mz_vals <- as.numeric(feature)
    idx <- vapply(mz_vals, function(x) .get_mz_index(object, x, tolerance_ppm), integer(1))
    label <- paste0("mz_", signif(object@variable_info$mz[idx], 6))
    return(list(type = "mz", index = idx, label = label))
  }

  feature_id <- as.character(feature)
  idx <- match(feature_id, object@variable_info$variable_id)
  if (any(is.na(idx))) {
    miss <- feature_id[is.na(idx)]
    stop("Unknown variable_id: ", paste(miss, collapse = ", "))
  }
  list(type = "variable_id", index = idx, label = feature_id)
}

.spatial_value_long <- function(object, feature_info) {
  sinfo <- object@sample_info
  base <- data.frame(
    sample_id = sinfo$sample_id,
    x = as.numeric(sinfo$x),
    y = as.numeric(sinfo$y),
    section_id = if ("section_id" %in% colnames(sinfo)) as.character(sinfo$section_id) else "section_1",
    stringsAsFactors = FALSE
  )

  if (identical(feature_info$type, "tic")) {
    base$feature <- "TIC"
    base$value <- as.numeric(colSums(object@expression_data, na.rm = TRUE))
    return(base)
  }

  expr <- as.matrix(object@expression_data[feature_info$index, , drop = FALSE])
  out <- lapply(seq_len(nrow(expr)), function(i) {
    df <- base
    df$feature <- feature_info$label[[i]]
    df$value <- as.numeric(expr[i, ])
    df
  })
  do.call(rbind, out)
}

.attach_roi_membership <- function(object, df, roi_id = NULL) {
  membership <- .get_roi_membership(object, roi_id = roi_id)
  if (nrow(membership) == 0) {
    if (!is.null(roi_id)) {
      stop("Requested roi_id not found in object.")
    }
    df$roi_id <- "all"
    df$roi_name <- "all"
    return(df)
  }
  merged <- merge(
    df,
    unique(membership[, c("sample_id", "roi_id", "roi_name"), drop = FALSE]),
    by = "sample_id",
    all = FALSE
  )
  if (nrow(merged) == 0) {
    stop("No pixels matched requested ROI.")
  }
  merged
}

.apply_value_transform <- function(x, transform) {
  switch(
    transform,
    none = x,
    log10 = log10(pmax(x, 1e-12)),
    sqrt = sqrt(pmax(x, 0))
  )
}

.aggregate_spatial_values <- function(df, aggregation = c("none", "mean", "median", "sum")) {
  aggregation <- match.arg(aggregation)
  if (aggregation == "none") {
    return(df)
  }

  stat_fun <- switch(
    aggregation,
    mean = function(x) mean(x, na.rm = TRUE),
    median = function(x) median(x, na.rm = TRUE),
    sum = function(x) sum(x, na.rm = TRUE)
  )

  stats::aggregate(
    value ~ feature + section_id + roi_id + roi_name + x + y,
    data = df,
    FUN = stat_fun
  )
}

#' Unified spatial plotting entry for msi_dataset
#'
#' @description
#' A coherent \code{ggplot2}-based interface for MSI spatial maps. Supports
#' single/multiple features, ROI filtering/overlay, section filtering, and
#' faceting with predictable semantics.
#'
#' @param object A \code{msi_dataset} object.
#' @param feature Feature specification. For \code{feature_type = "mz"} use
#'   numeric m/z values. For \code{feature_type = "variable_id"} use IDs.
#' @param feature_type One of \code{"mz"}, \code{"variable_id"}, \code{"tic"}.
#' @param section_id Optional section filter.
#' @param roi_id Optional ROI filter.
#' @param overlay_roi If \code{TRUE}, draw ROI pixel markers on top.
#' @param facet_by Facet dimensions; subset of \code{c("feature", "section", "roi")}.
#' @param scales Passed to \code{ggplot2::facet_wrap(scales = ...)}.
#' @param aggregation Aggregation mode for duplicated points:
#'   \code{"none"}, \code{"mean"}, \code{"median"}, \code{"sum"}.
#' @param section_scale Optional per-section scaling:
#'   \code{"global"} or \code{"section"}.
#' @param tolerance_ppm ppm tolerance when \code{feature_type = "mz"}.
#' @param value_transform Value transformation: \code{"none"}, \code{"log10"},
#'   \code{"sqrt"}.
#' @param low_color Low fill color.
#' @param high_color High fill color.
#'
#' @return A \code{ggplot} object.
#' @export
plot_spatial <- function(object,
                         feature = NULL,
                         feature_type = c("mz", "variable_id", "tic"),
                         section_id = NULL,
                         roi_id = NULL,
                         overlay_roi = FALSE,
                         facet_by = c("feature"),
                         scales = "fixed",
                         aggregation = c("none", "mean", "median", "sum"),
                         section_scale = c("global", "section"),
                         tolerance_ppm = 20,
                         value_transform = c("none", "log10", "sqrt"),
                         low_color = "#f7fbff",
                         high_color = "#08306b") {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2 to use plot_spatial().")
  }

  feature_type <- match.arg(feature_type)
  aggregation <- match.arg(aggregation)
  section_scale <- match.arg(section_scale)
  value_transform <- match.arg(value_transform)

  feature_info <- .resolve_feature_index(
    object = object,
    feature = feature,
    feature_type = feature_type,
    tolerance_ppm = tolerance_ppm
  )
  df <- .spatial_value_long(object, feature_info)

  if (!is.null(section_id)) {
    df <- df[df$section_id %in% as.character(section_id), , drop = FALSE]
  }
  if (nrow(df) == 0) {
    stop("No pixels remaining after section filtering.")
  }

  df <- .attach_roi_membership(object, df, roi_id = roi_id)
  if (nrow(df) == 0) {
    stop("No pixels remaining after ROI filtering.")
  }

  overlay_df <- df
  df <- .aggregate_spatial_values(df, aggregation = aggregation)

  if (section_scale == "section") {
    split_key <- interaction(df$section_id, df$feature, df$roi_id, drop = TRUE)
    df$value <- stats::ave(df$value, split_key, FUN = function(v) {
        rng <- range(v, na.rm = TRUE)
        if (is.finite(rng[2] - rng[1]) && (rng[2] - rng[1]) > 0) {
          (v - rng[1]) / (rng[2] - rng[1])
        } else {
          rep(0, length(v))
        }
      })
  }
  df$value <- .apply_value_transform(df$value, transform = value_transform)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = low_color, high = high_color, na.value = "grey90") +
    ggplot2::labs(
      x = "x",
      y = "y",
      fill = "intensity",
      title = "MSI spatial map"
    ) +
    ggplot2::theme_bw()

  # facet free scales are incompatible with fixed-ratio coordinates
  if (identical(scales, "fixed")) {
    p <- p + ggplot2::coord_fixed()
  } else {
    p <- p + ggplot2::coord_cartesian()
  }

  if (isTRUE(overlay_roi) && !all(df$roi_id == "all")) {
    overlay_pts <- unique(overlay_df[, c("sample_id", "roi_id", "x", "y"), drop = FALSE])
    p <- p + ggplot2::geom_point(
      data = overlay_pts,
      mapping = ggplot2::aes(x = x, y = y),
      inherit.aes = FALSE,
      color = "black",
      size = 0.2,
      alpha = 0.6
    )
  }

  facet_by <- unique(facet_by)
  facet_by <- facet_by[facet_by %in% c("feature", "section", "roi")]
  if (length(facet_by) > 0) {
    facet_cols <- c(
      feature = "feature",
      section = "section_id",
      roi = "roi_name"
    )[facet_by]
    p <- p + ggplot2::facet_wrap(
      facets = stats::as.formula(paste("~", paste(facet_cols, collapse = " + "))),
      scales = scales
    )
  }

  object <- add_process_log(
    object = object,
    function_name = "plot_spatial()",
    parameter = list(
      feature_type = feature_type,
      n_features = length(unique(df$feature)),
      section_id = section_id,
      roi_id = roi_id,
      facet_by = facet_by,
      aggregation = aggregation,
      section_scale = section_scale,
      value_transform = value_transform
    )
  )
  attr(p, "object") <- object
  p
}
