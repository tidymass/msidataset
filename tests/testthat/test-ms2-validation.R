test_that("compare_spectra returns similarity metrics", {
  object <- make_test_msi_object()
  mgf <- make_temp_mgf(precursor_mz = c(191.02), scans = c("spectrum_1"))
  object <- mutate_ms2(object, path = mgf, ms1.ms2.match.mz.tol = 20)
  object <- infer_pseudo_ms2(
    object = object,
    variable_id = object@variable_info$variable_id[which.min(abs(object@variable_info$mz - 191.02))],
    adaptive_threshold = TRUE,
    fallback_top_n = 3,
    fragment_window = 200
  )

  tbl <- object@ms2_data$ms2_table
  real_id <- unique(tbl$ms2_id[tbl$ms2_type == "actual"])[1]
  pseudo_id <- unique(tbl$ms2_id[tbl$ms2_type %in% c("pseudo", "pseudo_candidate")])[1]

  cmp <- compare_spectra(object, ms2_id = real_id, pseudo_ms2_id = pseudo_id, tolerance = 0.01)
  expect_true(is.list(cmp))
  expect_true(all(c("matched_fragments", "cosine_similarity", "overlap_ratio", "matched_fragment_count") %in% names(cmp)))
  expect_true(is.numeric(cmp$cosine_similarity))
  expect_true(is.numeric(cmp$overlap_ratio))
})

test_that("benchmark_pseudo_ms2 stores benchmark table", {
  object <- make_test_msi_object()
  mgf <- make_temp_mgf(precursor_mz = c(191.02, 306.07), scans = c("spectrum_1", "spectrum_2"))
  object <- mutate_ms2(object, path = mgf, ms1.ms2.match.mz.tol = 20)

  object <- benchmark_pseudo_ms2(
    object = object,
    tolerance = 0.01,
    adaptive_threshold = TRUE,
    fallback_top_n = 3,
    fragment_window = 200
  )

  expect_true("pseudo_ms2_benchmark" %in% names(object@other_files))
  bench <- object@other_files$pseudo_ms2_benchmark
  expect_true(is.data.frame(bench))
  expect_true(all(c(
    "precursor_variable_id", "precursor_mz", "real_fragment_count", "pseudo_fragment_count",
    "matched_fragment_count", "cosine_similarity", "fragment_overlap_ratio"
  ) %in% colnames(bench)))
})

test_that("mirror_plot returns ggplot", {
  testthat::skip_if_not_installed("ggplot2")

  object <- make_test_msi_object()
  mgf <- make_temp_mgf(precursor_mz = c(191.02), scans = c("spectrum_1"))
  object <- mutate_ms2(object, path = mgf, ms1.ms2.match.mz.tol = 20)
  object <- infer_pseudo_ms2(
    object = object,
    variable_id = object@variable_info$variable_id[which.min(abs(object@variable_info$mz - 191.02))],
    adaptive_threshold = TRUE,
    fallback_top_n = 3,
    fragment_window = 200
  )

  tbl <- object@ms2_data$ms2_table
  real_id <- unique(tbl$ms2_id[tbl$ms2_type == "actual"])[1]
  pseudo_id <- unique(tbl$ms2_id[tbl$ms2_type %in% c("pseudo", "pseudo_candidate")])[1]

  p <- mirror_plot(object, ms2_id = real_id, pseudo_ms2_id = pseudo_id)
  expect_s3_class(p, "ggplot")
})
