#' Choose the number of feature clusters
#'
#' Run centered PCA, k-means and a cosine silhouette permutation baseline.
#' Ties in the eligible observed-minus-baseline criterion select the smaller k.
#' @param vectors Feature unit matrix or result of [embed_features()].
#' @return List containing the selection table, selected k, PCA scores and fits.
#'   Status is `"INSUFFICIENT_FEATURES"` when there are too few vectors or
#'   dimensions, or `"NO_ELIGIBLE_K"` when no eligible count can be selected.
#' @examples
#' \dontrun{selection <- select_k(embedding)}
#' @export
select_k <- function(vectors) {
  x <- .msi_unit(vectors)
  if (nrow(x) < 23L || min(dim(x)) < 20L) {
    return(list(status = "INSUFFICIENT_FEATURES", selected_k = NA_integer_,
                scores = matrix(numeric(), nrow(x), 0L), table = data.frame(), fits = list()))
  }
  restore_rng <- .msi_rng_restorer()
  on.exit(restore_rng(), add = TRUE)
  p <- stats::prcomp(x, center = TRUE, scale. = FALSE, rank. = 20)$x[, 1:20, drop = FALSE]
  distinct <- nrow(unique(p))
  if (any(rowSums(p^2) == 0) || distinct < 2L) {
    return(list(status = "NO_ELIGIBLE_K", selected_k = NA_integer_,
                scores = p, table = data.frame(), fits = list()))
  }
  rows <- vector("list", 7); fits <- vector("list", 7); names(fits) <- as.character(2:8)
  for (k in 2:8) {
    if (k > distinct) {
      rows[[k-1L]] <- data.frame(k = k, silhouette = NA_real_, null_mean = NA_real_,
                                difference = NA_real_, min_cluster_size = NA_integer_, eligible = FALSE)
      next
    }
    set.seed(20261111); fit <- .msi_fit(p, k); fits[[as.character(k)]] <- fit
    obs <- mean(.msi_silhouette(p, fit$cluster)); null <- numeric(20)
    for (s in 0:19) {
      set.seed(s); perm <- apply(p, 2, sample, replace = FALSE)
      set.seed(20261111); pred <- .msi_fit(perm, k, 10L)$cluster
      null[s + 1L] <- mean(.msi_silhouette(perm, pred))
    }
    rows[[k-1L]] <- data.frame(k = k, silhouette = obs, null_mean = mean(null), difference = obs - mean(null), min_cluster_size = min(fit$size), eligible = min(fit$size) >= 5)
  }
  tab <- do.call(rbind, rows); eligible <- which(tab$eligible); selected <- if (length(eligible)) tab$k[eligible[order(-tab$difference[eligible], tab$k[eligible])][1]] else NA_integer_
  tab$selected <- tab$k == selected
  list(status = if (is.na(selected)) "NO_ELIGIBLE_K" else "OK", table = tab, selected_k = selected, scores = p, fits = fits, process_info = list(function_name = "select_k", seed = 20261111L, null_permutations = 20L, time = as.character(Sys.time())))
}
