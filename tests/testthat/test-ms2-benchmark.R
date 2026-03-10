test_that("benchmarking workflow produces metrics and plot", {
  testthat::skip_if_not_installed("ggplot2")

  object <- make_targeted_pseudo_test_object()
  object <- generate_pseudo_ms2(
    object = object,
    precursor_mz = c(151.02),
    top_n = 3,
    method = "spearman"
  )

  pseudo_tbl <- object@ms2_data$pseudo_ms2_table
  expect_true(nrow(pseudo_tbl) > 0)

  experimental_ms2 <- data.frame(
    precursor_mz = pseudo_tbl$precursor_mz,
    fragment_mz = pseudo_tbl$fragment_mz,
    intensity = pseudo_tbl$pseudo_intensity,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  out <- benchmark_pseudo_ms2(
    object = object,
    experimental_ms2 = experimental_ms2,
    ppm = 10
  )

  expect_true("pseudo_ms2_benchmark" %in% names(out@ms2_data))
  bench_tbl <- out@ms2_data$pseudo_ms2_benchmark
  expect_s3_class(bench_tbl, "data.frame")
  expect_true(all(c(
    "precursor_mz", "cosine_similarity",
    "fragment_recall", "fragment_precision"
  ) %in% colnames(bench_tbl)))
  expect_true(all(bench_tbl$fragment_recall <= 1, na.rm = TRUE))
  expect_true(all(bench_tbl$fragment_precision <= 1, na.rm = TRUE))

  plot_obj <- plot_pseudo_ms2_benchmark(out)
  expect_s3_class(plot_obj, "ggplot")
})
