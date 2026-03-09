test_that("summarize_msi returns expected structure", {
  object <- add_test_rois(make_test_msi_object())
  s <- summarize_msi(object)

  expect_type(s, "list")
  expect_true(all(c("n_variables", "n_samples", "n_sections", "mz_range", "n_roi", "n_process_steps") %in% names(s)))
  expect_equal(s$n_variables, nrow(object@variable_info))
  expect_equal(s$n_samples, nrow(object@sample_info))
  obj2 <- attr(s, "object")
  expect_s4_class(obj2, "msi_dataset")
  expect_equal(last_process_fun(obj2), "summarize_msi()")
  lg <- last_process_log(obj2)
  expect_true(all(c(
    "n_variables_before", "n_samples_before", "n_variables_after", "n_samples_after"
  ) %in% names(lg@parameter)))
})

test_that("summarize_msi input validation", {
  expect_error(summarize_msi("bad"), "msi_dataset")
})

test_that("summarize_msi handles edge empty feature object", {
  object <- make_test_msi_object()
  object2 <- filter_msi(
    object,
    target = "variable_info",
    filter_type = "mz_range",
    mz_min = 9999,
    mz_max = 10000
  )
  s <- summarize_msi(object2)
  expect_equal(s$n_variables, 0)
  expect_equal(s$n_samples, nrow(object2@sample_info))
})
