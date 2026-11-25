# Project the already merged retrieval table using an explicit allowlist. This
# helper does not create a complete evidence view.
.msi_evidence_neighbors <- function(merged, cards = NULL) {
  required <- c("cluster_id", "feature_id", "rank", "HMDB_ID", "name", "taxonomy_path",
                "card_completeness", "cosine", "n_reps_hit", "max_cosine")
  if (!is.data.frame(merged) || !all(required %in% names(merged))) stop("Invalid merged neighbors")
  if (anyNA(merged$HMDB_ID) || any(!nzchar(as.character(merged$HMDB_ID))) ||
      anyDuplicated(merged[c("cluster_id", "HMDB_ID")])) stop("Missing or duplicated neighbor ID")
  if (!is.null(cards) && (!is.data.frame(cards) || !all(c("HMDB_ID", "description") %in% names(cards)) ||
                          anyDuplicated(cards$HMDB_ID))) stop("Invalid cards")
  lapply(seq_len(nrow(merged)), function(i) {
    x <- merged[i, , drop = FALSE]
    nullable <- function(y) if (length(y) == 0L || is.na(y) || !nzchar(as.character(y))) NULL else as.character(y)
    cid <- as.character(x$cluster_id)
    item <- list(evidence_id = .msi_evidence_id("semantic_neighbor", cid, x$feature_id,
                       as.character(x$rank), x$HMDB_ID),
      feature_id = as.character(x$feature_id), rank = as.integer(x$rank),
      HMDB_ID = as.character(x$HMDB_ID), name = as.character(x$name),
      taxonomy_path = as.character(x$taxonomy_path), cosine = as.numeric(x$cosine),
      interpretation = "semantic neighbor, not identity", n_reps_hit = as.integer(x$n_reps_hit),
      max_cosine = as.numeric(x$max_cosine), cluster_id = cid, query_type = "feature",
      card_completeness = nullable(x$card_completeness))
    # The evidence view omits the description. Keep the source value separately so
    # that an eventual view builder can resolve the explicit task conflict.
    if (!is.null(cards)) item["description_source"] <- list(nullable(cards$description[match(x$HMDB_ID, cards$HMDB_ID)]))
    item
  })
}
