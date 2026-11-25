test_that("reference library filter preserves excluded candidates and resolves ties", {
  x <- data.frame(feature_id = c("F1", "F1", "F1", "F2"),
                  HMDB_ID = c("B", "A", "X", NA), mfc_priority_score = c(.3, .3, .9, .2),
                  score_gap = 0, relative_support = 1)
  z <- .mfc_filter_candidates(x, c("A", "B"))
  expect_equal(z$kept$HMDB_ID, c("A", "B"))
  expect_equal(z$kept$rank_within_feature, 1:2)
  expect_equal(z$kept$score_gap, c(0, 0))
  expect_equal(z$removed$filter_reason, c("no_embedding", "missing_hmdb_id"))
  expect_equal(nrow(z$kept) + nrow(z$removed), nrow(x))
})

test_that("reference query preparation retains and flags unresolved fragments", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("1,1,100,2,4", "2,1,100,4,8"), path)
  x <- read_msi_data(path, mz = c(100, 50))
  pairs <- data.frame(precursor_feature_id = c("F003", "F003"), precursor_mz = 100,
                      fragment_feature_id = c("F004", "missing"), fragment_mz = 50,
                      sas_d7 = c(.9, .8))
  z <- .mfc_prepare_pseudo(x, pairs)
  expect_equal(nrow(z), 2)
  expect_equal(z$intensity[1], 100)
  expect_equal(z$preparation_status, c("ok", "missing_fragment_feature"))
  expect_equal(z$cluster_id, c("P_F003", "P_F003"))
})

test_that("MFC requires a compatible step6 object", {
  path <- tempfile(fileext = ".csv")
  writeLines(c("1,1,100,2,4", "2,1,100,4,8"), path)
  x <- read_msi_data(path, mz = c(100, 50))
  x <- define_roi(x, data.frame(x = c(1, 2), y = c(1, 1), region = c("periportal", "pericentral")))
  expect_error(calculate_mfc(x, list(), tempdir()), "step6")
})
