#' Select representative actual features for each cluster
#'
#' Choose the member nearest the centroid by cosine, then its nearest members;
#' ties use the lexicographic feature ID. Small and zero-centroid clusters are
#' marked in the `issues` attribute rather than silently omitted.
#' @param vectors Unit matrix or [embed_features()] result.
#' @param members Data frame with feature_id and cluster_id.
#' @return Data frame of representatives with an `issues` attribute.
#' @examples
#' \dontrun{representatives <- select_representatives(embedding, clusters$members)}
#' @export
select_representatives <- function(vectors, members) {
  x <- .msi_unit(vectors); ids <- if (is.list(vectors)) vectors$feature_ids else rownames(x)
  if (is.null(ids) || anyDuplicated(ids) || anyDuplicated(members$feature_id) || anyNA(match(members$feature_id,ids))) stop("Invalid feature mapping")
  x <- x / sqrt(rowSums(x^2)); result <- list(); issues <- list()
  for (cid in sort(unique(members$cluster_id))) {
    names <- as.character(members$feature_id[members$cluster_id == cid]); a <- x[match(names,ids),,drop=FALSE]; n <- length(names); centroid <- colMeans(a); den <- sqrt(sum(centroid^2))
    if (n < 5L || !is.finite(den) || den == 0) { issues[[length(issues)+1L]] <- data.frame(cluster_id=cid, status=if (n < 5L) "SMALL_CLUSTER_LT5" else "ZERO_CENTROID"); next }
    center <- order(1-as.vector(a %*% centroid)/den, names)[1]; rest <- setdiff(seq_len(n), center); rest <- rest[order(1-as.vector(a[rest,,drop=FALSE] %*% a[center,]), names[rest])]
    chosen <- c(center, utils::head(rest, min(n,max(5,ceiling(.1*n)))-1L))
    result[[length(result)+1L]] <- data.frame(cluster_id=cid, feature_id=names[chosen], rank=seq_along(chosen), role=c("center",rep("nearest",length(chosen)-1L)), distance_to_center=1-as.vector(a[chosen,,drop=FALSE] %*% a[center,]), n_members=n)
  }
  out <- if(length(result)) do.call(rbind,result) else data.frame(cluster_id=integer(),feature_id=character(),rank=integer(),role=character(),distance_to_center=numeric(),n_members=integer())
  attr(out,"issues") <- issues; attr(out,"process_info") <- list(function_name="select_representatives", time=as.character(Sys.time())); out
}
