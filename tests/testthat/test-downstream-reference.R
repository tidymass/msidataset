test_that("Downstream results match the locked reference", {
  skip_if_not(identical(Sys.getenv("MSIDATASET_RUN_REFERENCE_DOWNSTREAM"), "true"),
              "Set MSIDATASET_RUN_REFERENCE_DOWNSTREAM=true for the external reference regression")
  root <- reference_fixture_root()
  data <- reference_fixture_root()
  read_tsv <- function(p) data.table::fread(p, data.table = FALSE, showProgress = FALSE)
  file <- function(p) file.path(root, p)
  reference_require_files(c(file(c("reference_runs/liver/mfc/candidates.tsv",
    "reference_runs/liver/de/features.tsv", "reference_runs/liver/clustering/members.tsv",
    "reference_runs/liver/clustering/representatives.tsv")),
    file(paste0("reference_runs/liver/evidence/base/C", 0:2, ".json")),
    file.path(data, c("reference_runs/liver/embedding/unit.tsv", "reference_runs/liver/embedding/ids.txt",
      "vector_library/meta.json", "vector_library/ids.tsv", "vector_library/embeddings.f32"))))
  cands <- read_tsv(file("reference_runs/liver/mfc/candidates.tsv"))
  de <- read_tsv(file("reference_runs/liver/de/features.tsv"))
  emb <- embed_features(cands, file.path(data, "vector_library"), de$feature_id)
  # The converted reference has no header.
  ref <- as.matrix(data.table::fread(file.path(data, "reference_runs/liver/embedding/unit.tsv"), header = FALSE, showProgress = FALSE))
  expect_identical(emb$feature_ids, readLines(file.path(data, "reference_runs/liver/embedding/ids.txt")))
  expect_lt(max(abs(emb$unit - ref)), 1e-6)
  selection <- select_k(emb)
  expect_equal(selection$selected_k, 3)
  clusters <- cluster_features(emb, selection)
  old <- read_tsv(file("reference_runs/liver/clustering/members.tsv"))
  aligned <- merge(clusters$members, old, by = "feature_id")
  expect_equal(nrow(aligned), nrow(old))
  tab <- table(aligned$cluster_id.x, aligned$cluster_id.y)
  pairs <- function(x) sum(x * (x - 1) / 2)
  a <- pairs(rowSums(tab)); b <- pairs(colSums(tab)); total <- choose(sum(tab), 2)
  ari <- (pairs(tab) - a * b / total) / ((a + b) / 2 - a * b / total)
  expect_equal(ari, 1)
  expect_true(all(clusters$stability$median_jaccard > .5))
  mapping <- apply(tab, 1, function(x) colnames(tab)[which.max(x)])
  reps <- select_representatives(emb, clusters$members)
  reps$reference_cluster <- unname(mapping[as.character(reps$cluster_id)])
  expected <- read_tsv(file("reference_runs/liver/clustering/representatives.tsv"))
  pairs <- merge(reps, expected, by.x = c("reference_cluster", "rank"),
                 by.y = c("cluster_id", "rank"))
  expect_equal(nrow(pairs), nrow(expected))
  expect_identical(pairs$feature_id.x, pairs$feature_id.y)
  expect_identical(pairs$role.x, pairs$role.y)
  near <- retrieve_neighbors(emb, reps, file.path(data, "vector_library"))$merged
  for (cid in unique(near$cluster_id)) {
    ref_id <- unname(mapping[as.character(cid)])
    view <- jsonlite::fromJSON(file(paste0("reference_runs/liver/evidence/base/C", ref_id, ".json")))
    expect_identical(as.character(near$HMDB_ID[near$cluster_id == cid]), as.character(view$semantic_neighbors$HMDB_ID))
  }
})
