test_that("plot_spatial returns ggplot and supports section filtering", {
  testthat::skip_if_not_installed("ggplot2")
  object <- make_test_msi_object()
  sec <- unique(object@sample_info$section_id)[1]

  dev_before <- grDevices::dev.cur()
  p <- plot_spatial(
    object = object,
    feature_type = "tic",
    section_id = sec,
    facet_by = c("section")
  )
  dev_after <- grDevices::dev.cur()

  expect_s3_class(p, "ggplot")
  expect_equal(dev_after, dev_before)
  expect_true(all(p$data$section_id == sec))
})

test_that("plot_spatial handles invalid feature id gracefully", {
  testthat::skip_if_not_installed("ggplot2")
  object <- make_test_msi_object()

  expect_error(
    plot_spatial(object, feature = "not_exist", feature_type = "variable_id"),
    "Unknown variable_id"
  )
})

test_that("mean_spectrum and pixel_spectrum return deterministic outputs", {
  object <- make_test_msi_object()

  ms <- mean_spectrum(object, plot = FALSE)
  px <- pixel_spectrum(object, sample_id = "spectrum_1", plot = FALSE)

  expect_s3_class(ms, "data.frame")
  expect_s3_class(px, "data.frame")
  expect_equal(nrow(ms), nrow(object@variable_info))
  expect_equal(nrow(px), nrow(object@variable_info))
  expect_true(all(c("mz", "intensity") %in% colnames(ms)))
  expect_true(all(c("mz", "intensity") %in% colnames(px)))
})

test_that("pixel_spectrum throws error for invalid sample id", {
  object <- make_test_msi_object()
  expect_error(pixel_spectrum(object, sample_id = "spectrum_999", plot = FALSE), "not found")
})
