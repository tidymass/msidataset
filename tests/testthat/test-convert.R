demo_data_dir <- function() {
  file.path(testthat::test_path("..", ".."), "data", "demo data")
}

test_that("convert_cardinal2msi_dataset basic functionality", {
  testthat::skip_if_not_installed("Cardinal")
  testthat::skip_on_ci()
  file1 <- file.path(demo_data_dir(), "03-total ion count.imzML")
  testthat::skip_if_not(file.exists(file1), "demo imzML file not found")

  cardinal_object <- Cardinal::readMSIData(file1)
  testthat::skip_if_not(
    methods::is(cardinal_object, "MSImagingExperiment"),
    "Current Cardinal backend returned non-MSImagingExperiment object for this demo file."
  )
  out <- convert_cardinal2msi_dataset(
    cardinal_object = cardinal_object,
    sample_name = "subject_demo",
    section_name = "section_demo"
  )

  expect_valid_msi(out)
  expect_true(all(grepl("^spectrum_[0-9]+$", out@sample_info$sample_id)))
  expect_equal(colnames(out@expression_data), out@sample_info$sample_id)
  expect_equal(rownames(out@expression_data), out@variable_info$variable_id)
})

test_that("convert_cardinal2msi_dataset input validation", {
  testthat::skip_if_not_installed("Cardinal")
  expect_error(convert_cardinal2msi_dataset(cardinal_object = "bad"), "Cardinal")
})

test_that("convert_cardinal2msi_dataset edge case single feature/pixel integrity", {
  object <- make_test_msi_object(n_features = 1, n_samples = 2, n_sections = 2)
  expect_valid_msi(object)
  expect_equal(nrow(object@expression_data), 1)
  expect_equal(ncol(object@expression_data), 2)
})
