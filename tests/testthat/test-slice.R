test_that("bracket subsetting keeps metadata integrity and z values", {
  object <- make_test_msi_object()
  sub <- object[1:3, 1:4]

  expect_valid_msi(sub)
  expect_equal(dim(sub@expression_data), c(3, 4))
  expect_equal(rownames(sub@expression_data), sub@variable_info$variable_id)
  expect_equal(colnames(sub@expression_data), sub@sample_info$sample_id)
  expect_true(all(sub@sample_info$z %in% object@sample_info$z))
})

test_that("slice_msi works with indices", {
  object <- make_test_msi_object()
  out <- slice_msi(object, i = 1:2, j = 1:3)
  expect_valid_msi(out)
  expect_equal(dim(out@expression_data), c(2, 3))
  lg <- last_process_log(out)
  expect_equal(lg@function_name, "slice_msi()")
})

test_that("section slicing preserves sample_id consistency and z dimension", {
  object <- make_test_msi_object()
  section <- unique(object@sample_info$section_id)[1]

  sub <- filter_msi_dataset(
    object,
    target = "sample_info",
    filter_type = "section",
    section_ids = section
  )

  expect_valid_msi(sub)
  expect_true(all(sub@sample_info$section_id == section))
  expect_equal(length(unique(sub@sample_info$z)), 1)
  expect_equal(colnames(sub@expression_data), sub@sample_info$sample_id)
})

test_that("slice_msi supports section slicing and keeps z", {
  object <- make_test_msi_object()
  sec <- unique(object@sample_info$section_id)[1]
  out <- slice_msi(object, section_id = sec)
  expect_valid_msi(out)
  expect_true(all(out@sample_info$section_id == sec))
  expect_equal(length(unique(out@sample_info$z)), 1)
})

test_that("subset_by_roi returns coherent object without orphan metadata", {
  object <- add_test_rois(make_test_msi_object())
  sub <- subset_by_roi(object, roi_id = "roi_1", keep_roi = TRUE)

  expect_valid_msi(sub)
  expect_true(all(sub@sample_info$sample_id %in% roi_data(object, "roi_1")$sample_id))
  expect_equal(colnames(sub@expression_data), sub@sample_info$sample_id)
  expect_equal(rownames(sub@expression_data), sub@variable_info$variable_id)
})

test_that("slice_msi validates bad inputs", {
  object <- make_test_msi_object()
  expect_error(slice_msi("bad"), "msi_dataset")
  expect_error(slice_msi(object, section_id = "not_exist"), "No pixels")
})
