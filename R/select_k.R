#' Choose the number of feature clusters
#'
#' Run centered PCA, k-means and a cosine silhouette permutation baseline.
#' Ties in the eligible observed-minus-baseline criterion select the smaller k.
#' @param vectors Feature unit matrix or result of [embed_features()].
#' @return List containing the selection table, selected k, PCA scores and fits.
#' @examples
#' \dontrun{selection <- select_k(embedding)}
#' @export
select_k <- function(vectors) {
  x <- .msi_unit(vectors); if (nrow(x) < 23L || min(dim(x)) < 20L) stop("At least 23 vectors and 20 dimensions required")
  p <- stats::prcomp(x, center = TRUE, scale. = FALSE, rank. = 20)$x[, 1:20, drop = FALSE]
  rows <- vector("list", 7); fits <- vector("list", 7); names(fits) <- as.character(2:8)
  for (k in 2:8) {
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
  list(table = tab, selected_k = selected, scores = p, fits = fits, process_info = list(function_name = "select_k", seed = 20261111L, null_permutations = 20L, time = as.character(Sys.time())))
}
