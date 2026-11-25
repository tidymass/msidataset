test_that("Reference differential features reproduce the locked results", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("yaml")
  p <- reference_lock()
  files <- vapply(p[c("csv", "labels", "de_all_reference")], `[[`, character(1), "path")
  skip_if_not(all(file.exists(files)), "reference data unavailable")
  raw <- data.table::fread(files[[1]], header = FALSE, showProgress = FALSE)
  labels <- data.table::fread(files[[2]], showProgress = FALSE)
  expect_true(all(as.matrix(raw[, 1:2]) == as.matrix(labels[, c("x", "y")])))
  observed <- detect_de_features(as.matrix(raw[, 1:2]), raw[[3]], as.matrix(raw[, 4:159]),
                                 labels$region, sprintf("F%03d", 3:158))
  expected <- data.table::fread(files[[3]])
  expect_identical(observed$feature_id, expected$feature_id)
  expect_identical(observed$direction, expected$direction)
  expect_identical(observed$is_de, expected$is_de)
  expect_equal(sum(observed$is_de), 58)
  expect_equal(sum(observed$is_de & observed$direction == "pericentral"), 38)
  expect_equal(sum(observed$is_de & observed$direction == "periportal"), 20)
  for (name in c("log2FC", "p_value", "fdr")) {
    error <- abs(observed[[name]] - expected[[name]]) / pmax(abs(expected[[name]]), .Machine$double.xmin)
    expect_lt(max(error), 1e-6)
  }
})

test_that("Reference pseudo-MS2 pairs reproduce the locked results", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("yaml")
  p <- reference_lock()
  files <- vapply(p[c("csv", "meta", "pseudo_reference")], `[[`, character(1), "path")
  skip_if_not(all(file.exists(files)), "reference data unavailable")
  raw <- data.table::fread(files[[1]], header = FALSE, showProgress = FALSE)
  meta <- data.table::fread(files[[2]], showProgress = FALSE)
  observed <- generate_pseudo_ms2(as.matrix(raw[, 1:2]), raw[[3]], as.matrix(raw[, 4:159]),
                                  meta[["m/z"]] - 1.007276, sprintf("F%03d", meta$col_index))
  expected <- data.table::fread(files[[3]])
  key <- function(pairs) paste(pairs$precursor_feature_id, pairs$fragment_feature_id)
  expect_equal(nrow(observed), 144L)
  expect_setequal(key(observed), key(expected))
  archived <- expected[match(key(observed), key(expected)), ]
  expect_equal(observed$sas, archived$sas_d7, tolerance = 1e-6)
  expect_equal(observed$spearman_rho, archived$spearman_rho, tolerance = 1e-6)
})
