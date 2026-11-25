#' Detect spatial-block differential features
#'
#' Compare paired periportal and pericentral tile means after mean-TIC
#' normalization. Tiles in one section are spatial blocks, not biological
#' replicates. Accept a mass_dataset or matrices.
#'
#' @param xy A mass_dataset or a two-column numeric matrix of pixel coordinates.
#' @param tic Positive TIC for each pixel.
#' @param intensity Pixel-by-feature numeric matrix of unnormalised intensity.
#' @param region Pixel region labels, including `periportal`, `pericentral`,
#'   and `VESSEL`.
#' @param feature_id Feature identifiers in the same order as columns.
#' @return A tibble of all features with log2 fold change, p value, FDR,
#'   direction, selection flag and tile count; `process_info` attribute
#'   records the frozen parameters and epsilon. For a mass_dataset input,
#'   return the object with results and parameters in process_info.
#' @examples
#' \dontrun{detect_de_features(xy, tic, intensity, region, feature_id)}
#' @seealso [stats::wilcox.test()], [stats::p.adjust()]
#' @export
detect_de_features <- function(xy, tic = NULL, intensity = NULL, region = NULL, feature_id = NULL) {
  if (methods::is(xy, "mass_dataset")) {
    object <- xy
    i <- .msidataset_inputs(object, require_roi = TRUE)
    bad <- .msidataset_input_problems(i, require_region = TRUE)
    if (length(bad)) return(.msidataset_append_step(object, "detect_de_features",
      list(status = "UNRESOLVED_INPUT", problems = bad), list(status = "UNRESOLVED_INPUT")))
    result <- detect_de_features(i$xy, i$tic, i$intensity, i$region, i$feature_id)
    return(.msidataset_append_step(object, "detect_de_features", result))
  }
  xy <- as.matrix(xy); intensity <- as.matrix(intensity)
  n <- nrow(xy)
  if (ncol(xy) != 2L || nrow(intensity) != n || length(tic) != n ||
      length(region) != n || length(feature_id) != ncol(intensity) ||
      anyNA(xy) || anyNA(region) || any(!is.finite(xy)) ||
      any(!is.finite(tic) | tic <= 0) || any(!is.finite(intensity)))
    stop("Invalid pixel inputs", call. = FALSE)
  values <- intensity * (mean(tic) / tic)
  pos <- values[region != "VESSEL", , drop = FALSE]
  pos <- pos[pos > 0]
  if (!length(pos)) stop("No positive non-vessel intensities", call. = FALSE)
  epsilon <- min(pos) / 2
  rm(pos)
  tile <- floor(sweep(xy, 2L, apply(xy, 2L, min), "-") / 20)
  tile_id <- paste(tile[, 1], tile[, 2], sep = ":")
  pp <- region == "periportal"; pc <- region == "pericentral"
  pp_count <- table(tile_id[pp]); pc_count <- table(tile_id[pc])
  valid <- intersect(names(pp_count)[pp_count >= 30], names(pc_count)[pc_count >= 30])
  if (!length(valid)) stop("No valid paired tiles", call. = FALSE)
  ppsum <- rowsum(values[pp & tile_id %in% valid, , drop = FALSE], tile_id[pp & tile_id %in% valid], reorder = FALSE)
  pcsum <- rowsum(values[pc & tile_id %in% valid, , drop = FALSE], tile_id[pc & tile_id %in% valid], reorder = FALSE)
  ppmean <- ppsum[valid, , drop = FALSE] / as.numeric(pp_count[valid])
  pcmean <- pcsum[valid, , drop = FALSE] / as.numeric(pc_count[valid])
  fc <- log2((pcmean + epsilon) / (ppmean + epsilon))
  p <- vapply(seq_len(ncol(fc)), function(i) {
    z <- fc[, i]
    if (any(!is.finite(z))) return(NA_real_)
    if (all(z == 0)) return(1)
    if (length(z) <= 50L && all(z != 0)) {
      w <- sum(rank(abs(z), ties.method = "average")[z > 0])
      return(min(1, 2 * min(stats::psignrank(floor(w), length(z)),
                               stats::psignrank(ceiling(w) - 1, length(z), lower.tail = FALSE))))
    }
    suppressWarnings(stats::wilcox.test(z, alternative = "two.sided", exact = FALSE, correct = FALSE)$p.value)
  }, numeric(1))
  avg <- unname(colMeans(fc)); fdr <- stats::p.adjust(p, method = "BH")
  ans <- tibble::tibble(feature_id = as.character(feature_id), log2FC = avg,
    p_value = p, fdr = fdr, direction = ifelse(avg > 0, "pericentral", ifelse(avg < 0, "periportal", "undetermined")),
    is_de = fdr < 0.05 & abs(avg) > 0.25, n_tiles = length(valid))
  attr(ans, "process_info") <- list(function_name = "detect_de_features", tile_pixels = 20L,
    region_min_pixels = 30L, epsilon = epsilon, method = "paired Wilcoxon, two-sided, auto exact/normal, no continuity correction; BH",
    fdr_lt = 0.05, absolute_log2fc_gt = 0.25, time = as.character(Sys.time()))
  ans
}
