test_that("D15 downstream results match the locked references", {
  skip_if_not(identical(Sys.getenv("MSIDATASET_RUN_D15_DOWNSTREAM"), "true"),
              "Set MSIDATASET_RUN_D15_DOWNSTREAM=true for the external D15 regression")
  root <- d15_fixture_root("MSIDATASET_D15_STUDY_ROOT")
  data <- d15_fixture_root("MSIDATASET_D15_DATA_ROOT")
  read_tsv <- function(p) data.table::fread(p, data.table = FALSE, showProgress = FALSE)
  file <- function(p) file.path(root, p)
  d15_require_files(c(file(c("runs/d15_resume/mfc_final_d15/mfc_candidates_filtered.tsv",
    "runs/d15_roi2/de_d15/de_features.tsv", "runs/d15_resume/cluster_d15/cluster_members.tsv",
    "runs/d15_resume/cluster_d15/representatives.tsv")),
    file(paste0("runs/d15_neighbors_v2/neighbors_v2_d15/views/C", 0:2, ".json")),
    file.path(data, c("pk2/d15_unit.tsv", "pk2/d15_ids.txt",
      "hmdb_full_bge_v2/meta.json", "hmdb_full_bge_v2/ids.tsv", "hmdb_full_bge_v2/embeddings.f32"))))
  cands <- read_tsv(file("runs/d15_resume/mfc_final_d15/mfc_candidates_filtered.tsv"))
  de <- read_tsv(file("runs/d15_roi2/de_d15/de_features.tsv"))
  emb <- embed_features(cands, file.path(data, "hmdb_full_bge_v2"), de$feature_id)
  # The converted reference has no header.
  ref <- as.matrix(data.table::fread(file.path(data, "pk2/d15_unit.tsv"), header = FALSE, showProgress = FALSE))
  expect_identical(emb$feature_ids, readLines(file.path(data, "pk2/d15_ids.txt")))
  expect_lt(max(abs(emb$unit - ref)), 1e-6)
  selection <- select_k(emb)
  expect_equal(selection$selected_k, 3)
  clusters <- cluster_features(emb, selection)
  old <- read_tsv(file("runs/d15_resume/cluster_d15/cluster_members.tsv"))
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
  expected <- read_tsv(file("runs/d15_resume/cluster_d15/representatives.tsv"))
  pairs <- merge(reps, expected, by.x = c("reference_cluster", "rank"),
                 by.y = c("cluster_id", "rank"))
  expect_equal(nrow(pairs), nrow(expected))
  expect_identical(pairs$feature_id.x, pairs$feature_id.y)
  expect_identical(pairs$role.x, pairs$role.y)
  near <- retrieve_neighbors(emb, reps, file.path(data, "hmdb_full_bge_v2"))$merged
  for (cid in unique(near$cluster_id)) {
    ref_id <- unname(mapping[as.character(cid)])
    view <- jsonlite::fromJSON(file(paste0("runs/d15_neighbors_v2/neighbors_v2_d15/views/C", ref_id, ".json")))
    expect_identical(as.character(near$HMDB_ID[near$cluster_id == cid]), as.character(view$semantic_neighbors$HMDB_ID))
  }
})
