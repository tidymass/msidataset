test_that("compare_msi reports key state deltas", {
  before <- add_test_rois(make_test_msi_object())
  after <- filter_msi(
    before,
    target = "variable_info",
    filter_type = "mz_range",
    mz_min = 180,
    mz_max = 320
  )

  cmp <- compare_msi(before, after)

  expect_s3_class(cmp, "compare_msi")
  expect_true(all(c("metric", "before", "after", "delta") %in% colnames(cmp)))
  expect_true("feature_count" %in% cmp$metric)
  expect_true("pixel_count" %in% cmp$metric)
  expect_true("roi_membership_count" %in% cmp$metric)
  expect_true("ms2_attachment_count" %in% cmp$metric)
  expect_true("process_log_steps" %in% cmp$metric)
})

test_that("compare_msi validates inputs", {
  object <- make_test_msi_object()
  expect_error(compare_msi("bad", object), "before must be")
  expect_error(compare_msi(object, "bad"), "after must be")
})
