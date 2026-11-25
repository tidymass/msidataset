test_that("neighbor evidence uses an allowlist and retains supplied descriptions", {
  x <- data.frame(cluster_id = 0, feature_id = "F1", rank = 1, HMDB_ID = "H1",
                  name = "candidate", taxonomy_path = "class", card_completeness = NA_character_,
                  cosine = .9, n_reps_hit = 2, max_cosine = .9,
                  author_name = "must not be projected")
  cards <- data.frame(HMDB_ID = "H1", description = "card description")
  result <- .msi_evidence_neighbors(x, cards)[[1]]
  expect_equal(result$description_source, "card description")
  expect_null(result$card_completeness)
  expect_false("author_name" %in% names(result))
  expect_equal(result$evidence_id, .msi_evidence_id("semantic_neighbor", "0", "F1", "1", "H1"))
  expect_error(.msi_evidence_neighbors(rbind(x, x)), "duplicated")
})
