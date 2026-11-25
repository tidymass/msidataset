#' Build feature embeddings from MFC candidates
#'
#' Normalize each available library vector, form the MFC-priority-weighted
#' mean, then L2 normalize. Missing or unusable vectors are flagged in the
#' coverage table and never cause the input candidate table to be modified.
#' @param candidates Data frame containing feature_id, HMDB_ID and mfc_priority_score.
#' @param library_path Directory understood by [read_vector_library()].
#' @param feature_ids Feature order, including features without candidates.
#' @param weight_rule Only `"mfc"` is supported.
#' @return A list containing `unit`, `feature_ids`, `coverage`, and `process_info`. A result with no embedded features has
#'   status `"INSUFFICIENT_FEATURES"`.
#' @examples
#' \dontrun{embed_features(candidates, library_path, feature_ids)}
#' @export
embed_features <- function(candidates, library_path, feature_ids = unique(as.character(candidates$feature_id)), weight_rule = "mfc") {
  if (!identical(weight_rule, "mfc")) stop("Only the frozen mfc weighting rule is supported")
  if (!all(c("feature_id", "HMDB_ID", "mfc_priority_score") %in% names(candidates)) || anyDuplicated(feature_ids) || anyNA(feature_ids)) stop("Invalid candidate columns or feature IDs")
  info <- list(function_name = "embed_features", weight_rule = weight_rule,
               time = as.character(Sys.time()))
  if (!length(feature_ids) || !nrow(candidates)) {
    coverage <- data.frame(feature_id = as.character(feature_ids),
                           coverage = rep(0, length(feature_ids)),
                           n_candidates = rep(0L, length(feature_ids)),
                           status = rep("NO_CANDIDATE", length(feature_ids)),
                           vector_available = rep(FALSE, length(feature_ids)))
    return(list(status = "INSUFFICIENT_FEATURES", unit = matrix(numeric(), 0L, 0L),
                feature_ids = character(), coverage = coverage, process_info = info))
  }
  ids <- utils::read.delim(file.path(library_path, "ids.tsv"), check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (anyDuplicated(ids$HMDB_ID)) stop("Duplicate library IDs")
  idx <- match(candidates$HMDB_ID, ids$HMDB_ID)
  unique_rows <- sort(unique(idx[!is.na(idx)]))
  lib <- read_vector_library(library_path, unique_rows)
  norm <- sqrt(rowSums(lib$vectors^2))
  valid <- is.finite(norm) & norm > 0 & rowSums(!is.finite(lib$vectors)) == 0
  if (any(valid)) lib$vectors[valid, ] <- lib$vectors[valid, , drop = FALSE] / norm[valid]
  loc <- match(idx, unique_rows)
  dim <- ncol(lib$vectors); result <- list(); coverage <- vector("list", length(feature_ids))
  for (i in seq_along(feature_ids)) {
    ix <- which(candidates$feature_id == feature_ids[i]); w <- candidates$mfc_priority_score[ix]
    ok <- !is.na(loc[ix]); ok[ok] <- valid[loc[ix][ok]]
    cov <- 0; v <- NULL; status <- "LOW_COVERAGE"
    if (length(ix) && all(is.finite(w) & w >= 0) && sum(w) > 0) {
      cov <- sum(w[ok]) / sum(w)
      if (any(ok) && sum(w[ok]) > 0) {
        v <- colSums(lib$vectors[loc[ix[ok]], , drop = FALSE] * (w[ok] / sum(w[ok])))
        n <- sqrt(sum(v^2)); if (is.finite(n) && n > 0) { v <- v / n; status <- "OK" } else v <- NULL
      }
    }
    if (!length(ix)) status <- "NO_CANDIDATE"
    if (!is.null(v)) result[[feature_ids[i]]] <- v
    coverage[[i]] <- data.frame(feature_id = feature_ids[i], coverage = cov, n_candidates = length(ix), status = status, vector_available = !is.null(v))
  }
  mat <- if (length(result)) do.call(rbind, result) else matrix(numeric(), 0, dim)
  info$library <- list(basename = basename(library_path),
                       ids = .msi_file_record(file.path(library_path, "ids.tsv")))
  list(status = if (nrow(mat)) "OK" else "INSUFFICIENT_FEATURES",
       unit = mat, feature_ids = if (length(result)) names(result) else character(),
       coverage = do.call(rbind, coverage), process_info = info)
}
