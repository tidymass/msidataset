#' Cluster features and estimate subsample stability
#'
#' Uses the selected PCA fit and 50 independent 90% subsamples, each refitting
#' PCA and k-means. Jaccard scores compare best overlapping subsample groups.
#' @param vectors Feature matrix or [embed_features()] result.
#' @param selection Result of [select_k()].
#' @return Members, stability table, co-clustering matrices, and process history.
#'   An unsuccessful selection yields an empty membership table and its status.
#' @examples
#' \dontrun{clusters <- cluster_features(embedding, selection)}
#' @export
cluster_features <- function(vectors, selection) {
  x <- .msi_unit(vectors); k <- selection$selected_k
  if (is.null(selection$scores) || !identical(nrow(selection$scores), nrow(x))) {
    stop("Selection score rows must match the number of input vectors")
  }
  status <- .msi_stage_status(selection)
  if (!is.null(status) || length(k) != 1L || is.na(k)) {
    if (is.null(status)) status <- "NO_ELIGIBLE_K"
    members <- data.frame(feature_id = character(), cluster_id = integer(),
                          silhouette_sample = numeric(), PC1 = numeric(), PC2 = numeric())
    attr(members, "status") <- status
    return(list(status = status, members = members))
  }
  fit <- selection$fits[[as.character(k)]]; labels <- fit$cluster - 1L
  ids <- if (is.list(vectors)) vectors$feature_ids else rownames(x)
  if (is.null(ids)) stop("Feature IDs are required")
  n <- nrow(x); js <- matrix(0, 50, k); co <- matrix(0, n, n); exposure <- co
  restore_rng <- .msi_rng_restorer()
  on.exit(restore_rng(), add = TRUE)
  set.seed(20261110); samples <- replicate(50, sort(sample.int(n, floor(.9*n))), simplify = FALSE)
  for (s in seq_len(50)) {
    ix <- samples[[s]]; p <- stats::prcomp(x[ix, , drop = FALSE], rank. = 20)$x[, 1:20, drop = FALSE]
    set.seed(20261111); pred <- .msi_fit(p, k)$cluster
    for (j in seq_len(k)) {
      base <- which(labels[ix] == j - 1L)
      js[s,j] <- max(vapply(seq_len(k), function(g) { v <- which(pred == g); u <- union(base,v); if (length(u)) length(intersect(base,v))/length(u) else 0 }, numeric(1)))
    }
    exposure[ix,ix] <- exposure[ix,ix] + 1; co[ix,ix] <- co[ix,ix] + outer(pred,pred,`==`)
  }
  med <- apply(js, 2, stats::median)
  members <- data.frame(feature_id = ids, cluster_id = labels, silhouette_sample = .msi_silhouette(selection$scores, labels), PC1 = selection$scores[,1], PC2 = selection$scores[,2])
  stability <- data.frame(cluster_id = 0:(k-1L), n = as.integer(tabulate(labels + 1L, k)), median_jaccard = med, status = ifelse(med > .5, "OK", "UNSTABLE"))
  list(members = members, stability = stability, coclustering = ifelse(exposure > 0, co/exposure, 0), exposure = exposure, process_info = list(function_name = "cluster_features", stability_seed = 20261110L, samples = 50L, time = as.character(Sys.time())))
}
