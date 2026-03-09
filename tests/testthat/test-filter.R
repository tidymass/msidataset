test_that("filter_msi_dataset spatial filtering updates object and process log", {
  object <- make_test_msi_object()
  n0 <- ncol(object@expression_data)
  p0 <- length(object@process_info)

  object2 <- filter_msi_dataset(
    object,
    target = "sample_info",
    filter_type = "spatial",
    x_min = 1,
    x_max = 2,
    y_min = 1,
    y_max = 2
  )

  expect_valid_msi(object2)
  expect_lt(ncol(object2@expression_data), n0)
  expect_equal(ncol(object2@expression_data), nrow(object2@sample_info))
  expect_equal(length(object2@process_info), p0 + 1)
  expect_equal(last_process_fun(object2), "filter_msi_dataset(spatial)")
})

test_that("filter_msi_dataset mz range filtering updates variable dimension", {
  object <- make_test_msi_object()
  n0 <- nrow(object@expression_data)

  object2 <- filter_msi_dataset(
    object,
    target = "variable_info",
    filter_type = "mz_range",
    mz_min = 180,
    mz_max = 320
  )

  expect_valid_msi(object2)
  expect_lt(nrow(object2@expression_data), n0)
  expect_equal(nrow(object2@expression_data), nrow(object2@variable_info))
  expect_equal(colnames(object2@expression_data), object2@sample_info$sample_id)
})

test_that("filter_msi works on mz range", {
  object <- make_test_msi_object()
  out <- filter_msi(
    object,
    target = "variable_info",
    filter_type = "mz_range",
    mz_min = 190,
    mz_max = 310
  )
  expect_valid_msi(out)
  expect_lt(nrow(out@expression_data), nrow(object@expression_data))
  lg <- last_process_log(out)
  expect_equal(lg@function_name, "filter_msi()")
  expect_true(all(c(
    "n_variables_before", "n_samples_before", "n_variables_after", "n_samples_after"
  ) %in% names(lg@parameter)))
})

test_that("filter_msi_dataset validates invalid inputs", {
  object <- make_test_msi_object()

  expect_error(
    filter_msi_dataset(object, target = "sample_info", filter_type = "unknown"),
    "Unknown filter_type"
  )

  expect_error(
    filter_msi_dataset(object, target = "sample_info", filter_type = "section"),
    "section_ids must be provided"
  )

  expect_error(
    filter_msi("bad_object", target = "variable_info", filter_type = "mz_range"),
    "msi_dataset"
  )
})

test_that("filter_msi handles edge empty results and preserves integrity", {
  object <- make_test_msi_object()
  out <- filter_msi(
    object,
    target = "variable_info",
    filter_type = "mz_range",
    mz_min = 9999,
    mz_max = 10000
  )
  expect_valid_msi(out)
  expect_equal(nrow(out@expression_data), 0)
  expect_equal(nrow(out@variable_info), 0)
  expect_equal(ncol(out@expression_data), nrow(out@sample_info))
})
