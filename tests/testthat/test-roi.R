test_that("add_roi and list_roi register memberships correctly", {
  object <- make_test_msi_object()
  ids <- object@sample_info$sample_id[1:3]

  object <- add_roi(
    object = object,
    roi_id = "roi_a",
    roi_name = "ROI A",
    sample_id = ids,
    section_id = object@sample_info$section_id[1],
    overwrite = TRUE
  )

  summary_df <- list_roi(object)
  membership <- roi_data(object, roi_id = "roi_a")

  expect_equal(nrow(summary_df), 1)
  expect_equal(summary_df$n_pixels[[1]], 3)
  expect_true(all(membership$sample_id %in% ids))
  expect_true(isTRUE(validate_roi_integrity(object)))
})

test_that("roi_msi subsets spatial region", {
  object <- add_test_rois(make_test_msi_object())
  sub <- roi_msi(object, action = "subset", roi_id = "roi_1", keep_roi = TRUE)
  expect_valid_msi(sub)
  expect_true(all(sub@sample_info$sample_id %in% roi_data(object, "roi_1")$sample_id))
  lg <- last_process_log(sub)
  expect_equal(lg@function_name, "roi_msi()")
  expect_equal(lg@parameter$action, "subset")
})

test_that("rename_roi and remove_roi behave as expected", {
  object <- add_test_rois(make_test_msi_object())

  object <- rename_roi(object, roi_id = "roi_1", new_name = "ROI 1 renamed")
  expect_equal(list_roi(object)$roi_name[list_roi(object)$roi_id == "roi_1"], "ROI 1 renamed")

  object <- remove_roi(object, roi_id = "roi_2")
  expect_false("roi_2" %in% list_roi(object)$roi_id)
})

test_that("ROI interactions with filtering preserve consistency", {
  object <- add_test_rois(make_test_msi_object())
  filtered <- filter_msi_dataset(
    object,
    target = "sample_info",
    filter_type = "roi",
    roi_id = "roi_1"
  )

  roi_ids <- roi_data(object, "roi_1")$sample_id
  expect_valid_msi(filtered)
  expect_true(all(filtered@sample_info$sample_id %in% roi_ids))
})

test_that("add_roi throws error for unknown sample ids", {
  object <- make_test_msi_object()
  expect_error(
    add_roi(object, roi_id = "bad", sample_id = "spectrum_999"),
    "unknown sample_id"
  )
  expect_error(roi_msi("bad_object", action = "list"), "msi_dataset")
  expect_error(roi_msi(object, action = "bad"), "arg")
})
