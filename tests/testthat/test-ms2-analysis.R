test_that("infer_pseudo_ms2 supports spearman smoothing and threshold", {
  object <- make_test_msi_object(n_features = 5, n_samples = 8, n_sections = 2)
  out <- infer_pseudo_ms2(
    object = object,
    variable_id = object@variable_info$variable_id[3],
    correlation_threshold = -1,
    correlation_method = "spearman",
    spatial_smoothing = TRUE,
    smoothing_k = 1,
    fragment_intensity_threshold = 0
  )

  expect_valid_msi(out)
  expect_true("ms2_table" %in% names(out@ms2_data))
  expect_true(any(out@ms2_data$ms2_table$ms2_type %in% c("pseudo", "pseudo_candidate")))
  expect_true("confidence_score" %in% colnames(out@ms2_data$ms2_table))
  expect_true("pseudo_correlation_table" %in% names(out@ms2_data))
})

test_that("plot_spatial_ms2 returns ggplot object", {
  testthat::skip_if_not_installed("ggplot2")

  object <- make_test_msi_object()
  mgf <- make_temp_mgf(precursor_mz = c(191.02), scans = c("spectrum_1"))
  object <- mutate_ms2(object, path = mgf, ms1.ms2.match.mz.tol = 20)

  p <- plot_spatial_ms2(object, mz = 191.02, tolerance_ppm = 20, max_fragments = 2)
  expect_s3_class(p, "ggplot")
})

test_that("build_ms2_network returns igraph object", {
  testthat::skip_if_not_installed("igraph")

  object <- make_test_msi_object()
  mgf <- make_temp_mgf(precursor_mz = c(191.02), scans = c("spectrum_1"))
  object <- mutate_ms2(object, path = mgf, ms1.ms2.match.mz.tol = 20)

  g <- build_ms2_network(object, mz = 191.02, edge_weight = "fragment_intensity", tolerance_ppm = 20)
  expect_true(inherits(g, "igraph"))
})

test_that("generate_pseudo_ms2 stores targeted pseudo fragments across sections", {
  object <- make_targeted_pseudo_test_object()

  out <- generate_pseudo_ms2(
    object = object,
    precursor_mz = c(151.02, 289.13),
    ppm = 10,
    correlation_method = "spearman",
    correlation_threshold = 0.3,
    top_n = 2
  )

  expect_valid_msi(out)
  expect_true("pseudo_ms2_table" %in% names(out@ms2_data))
  tbl <- out@ms2_data$pseudo_ms2_table
  expect_true(all(c(
    "pseudo_ms2_id", "precursor_mz", "fragment_mz", "fragment_variable_id",
    "correlation", "fragment_intensity", "score", "normalized_intensity",
    "slice_consistency", "ms2_type"
  ) %in% colnames(tbl)))
  expect_true(all(tbl$ms2_type == "pseudo"))
  expect_true(all(tbl$fragment_mz < tbl$precursor_mz))
  expect_true(any(tbl$precursor_mz == 151.02))
  expect_true(any(tbl$precursor_mz == 289.13))
  expect_true(any(out@ms2_data$ms2_table$ms2_type == "pseudo"))
})

test_that("plot_pseudo_ms2 and plot_spatial_fragment return ggplot for targeted pseudo data", {
  testthat::skip_if_not_installed("ggplot2")

  object <- make_targeted_pseudo_test_object()
  object <- generate_pseudo_ms2(
    object = object,
    precursor_mz = c(151.02),
    correlation_threshold = 0.3,
    top_n = 2
  )

  p1 <- plot_pseudo_ms2(object, mz = 151.02, tolerance_ppm = 10)
  p2 <- plot_spatial_fragment(object, precursor_mz = 151.02, ppm = 10)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("merge_msi_datasets combines independently loaded objects", {
  object_a <- make_test_msi_object(n_features = 4, n_samples = 4, n_sections = 1)
  object_b <- make_test_msi_object(n_features = 4, n_samples = 4, n_sections = 1)

  merged <- merge_msi_datasets(object_a, object_b, bin_digits = 4)

  expect_valid_msi(merged)
  expect_equal(ncol(merged@expression_data), 8)
  expect_equal(length(unique(merged@sample_info$section_id)), 2)
  expect_true(last_process_fun(merged) == "merge_msi_datasets()")
})

test_that("detect_joint_fragments stores spatial-spectral candidate table", {
  object <- make_targeted_pseudo_test_object()

  out <- detect_joint_fragments(
    object = object,
    precursor_mz = c(151.02, 289.13),
    ppm = 10,
    correlation_method = "spearman"
  )

  expect_valid_msi(out)
  expect_true("fragment_candidates" %in% names(out@ms2_data))
  tbl <- out@ms2_data$fragment_candidates
  expect_true(all(c(
    "precursor_mz", "fragment_mz", "spatial_corr",
    "spectral_score", "intensity_score", "joint_score"
  ) %in% colnames(tbl)))
  expect_true(any(tbl$precursor_mz == 151.02))
  expect_true(any(tbl$precursor_mz == 289.13))
})

test_that("build_fragment_network stores igraph and plot_fragment_network returns ggplot", {
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("ggraph")
  testthat::skip_if_not_installed("ggplot2")

  object <- make_targeted_pseudo_test_object()
  object <- detect_joint_fragments(object, precursor_mz = c(151.02, 289.13))
  object <- build_fragment_network(object, top_k = 2)

  expect_true(inherits(object@other_files$pseudo_fragment_network, "igraph"))
  p <- plot_fragment_network(object)
  expect_s3_class(p, "ggplot")
})
