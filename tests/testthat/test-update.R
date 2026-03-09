test_that("update_msi_dataset migrates legacy IDs and remaps ROI membership", {
  object <- make_test_msi_object()
  object <- add_roi(
    object = object,
    roi_id = "r1",
    sample_id = "spectrum_1",
    overwrite = TRUE
  )

  object@sample_info$sample_id <- c("pixel_1", "pixel_2", "pixel_3", "pixel_4",
                                    "pixel_5", "pixel_6", "pixel_7", "pixel_8")
  colnames(object@expression_data) <- object@sample_info$sample_id
  object@other_files$roi@rois[[1]]@sample_id <- "pixel_1"

  out <- update_msi_dataset(object)

  expect_valid_msi(out)
  expect_true(all(grepl("^spectrum_[0-9]+$", out@sample_info$sample_id)))
  expect_equal(colnames(out@expression_data), out@sample_info$sample_id)
  expect_identical(out@other_files$roi@rois[[1]]@sample_id, "spectrum_1")
})

test_that("update_msi_dataset repairs dimnames without data loss when valid", {
  object <- make_test_msi_object()
  expr0 <- as.matrix(object@expression_data)
  info0 <- object@sample_info

  # break column names only
  colnames(object@expression_data) <- paste0("pixel_", seq_len(ncol(object@expression_data)))
  out <- update_msi_dataset(object)

  expect_valid_msi(out)
  expect_equal(dim(out@expression_data), dim(expr0))
  expect_equal(as.numeric(out@expression_data), as.numeric(expr0))
  expect_equal(out@sample_info$sample_id, info0$sample_id)
})

test_that("update_msi_dataset keeps already valid objects stable", {
  object <- make_test_msi_object()
  dim0 <- dim(object@expression_data)
  p0 <- length(object@process_info)

  out <- update_msi_dataset(object)

  expect_valid_msi(out)
  expect_equal(dim(out@expression_data), dim0)
  expect_equal(length(out@process_info), p0 + 1)
  expect_equal(last_process_fun(out), "update_msi_dataset()")
})
