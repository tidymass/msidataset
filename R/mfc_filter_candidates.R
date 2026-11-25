# Post-score embedding-library filter.  Kept and excluded rows are both returned;
# the scorer and its input evidence are not modified here.
.mfc_filter_candidates <- function(candidates, library_ids) {
  required <- c("feature_id", "HMDB_ID", "mfc_priority_score")
  if (!is.data.frame(candidates) || !all(required %in% names(candidates)) ||
      !is.numeric(candidates$mfc_priority_score) || !is.character(library_ids))
    stop("Candidate columns or library IDs are invalid", call. = FALSE)
  in_library <- !is.na(candidates$HMDB_ID) & candidates$HMDB_ID %in% library_ids
  kept <- candidates[in_library, , drop = FALSE]
  removed <- candidates[!in_library, , drop = FALSE]
  removed$filter_reason <- ifelse(is.na(removed$HMDB_ID) | !nzchar(removed$HMDB_ID),
                                  "missing_hmdb_id", "no_embedding")
  if (nrow(kept)) {
    idx <- order(kept$feature_id, -kept$mfc_priority_score, kept$HMDB_ID,
                 na.last = TRUE, method = "radix")
    kept <- kept[idx, , drop = FALSE]
    group <- split(seq_len(nrow(kept)), kept$feature_id, drop = TRUE)
    kept$rank_within_feature <- NA_integer_
    for (i in group) {
      kept$rank_within_feature[i] <- seq_along(i)
      vals <- kept$mfc_priority_score[i]
      best <- if (all(is.na(vals))) NA_real_ else max(vals, na.rm = TRUE)
      if ("score_gap" %in% names(kept)) kept$score_gap[i] <- best - vals
      if ("relative_support" %in% names(kept))
        kept$relative_support[i] <- if (is.na(best)) NA_real_ else if (best > 0) vals / best else 0
    }
  }
  rownames(kept) <- rownames(removed) <- NULL
  list(kept = kept, removed = removed,
       unrankable_rows = which(is.na(kept$feature_id) | !is.finite(kept$mfc_priority_score)))
}
