.msi_unit <- function(x) {
  if (is.list(x) && !is.null(x$unit)) x <- x$unit
  x <- as.matrix(x)
  if (!is.numeric(x) || any(!is.finite(x)) || any(rowSums(x^2) == 0)) stop("Vectors must be finite and nonzero")
  x
}
.msi_silhouette <- function(p, labels) {
  norms <- sqrt(rowSums(p^2)); if (any(norms == 0)) stop("Zero PCA score has undefined cosine distance")
  u <- p / norms; d <- 1 - tcrossprod(u); d[d < 0] <- 0; diag(d) <- 0
  groups <- split(seq_along(labels), labels)
  vapply(seq_along(labels), function(i) {
    own <- groups[[as.character(labels[i])]]
    if (length(own) == 1L) return(0)
    a <- sum(d[i, own]) / (length(own) - 1L)
    b <- min(vapply(groups[names(groups) != as.character(labels[i])], function(g) mean(d[i, g]), numeric(1)))
    if (max(a,b) == 0) 0 else (b-a)/max(a,b)
  }, numeric(1))
}
.msi_fit <- function(p, k, nstart = 50L) stats::kmeans(p, centers = k, nstart = nstart, iter.max = 300L)

.msi_stage_status <- function(result) {
  status <- if (is.data.frame(result)) attr(result, "status") else if (is.list(result)) result$status
  if (identical(status, "NO_ELIGIBLE_K") || identical(status, "INSUFFICIENT_FEATURES")) status
}

.msi_rng_restorer <- function() {
  seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  function() {
    if (!is.null(seed)) {
      assign(".Random.seed", seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }
}
