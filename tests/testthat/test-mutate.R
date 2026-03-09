test_that("mutate_ms2 populates ms2_data and updates variable_info ms2 flags", {
  object <- make_test_msi_object()
  mgf <- make_temp_mgf()

  out <- mutate_ms2(object, path = mgf, ms1.ms2.match.mz.tol = 20)
  info <- extract_ms2(out, what = "ms2_info")
  tbl <- out@ms2_data$ms2_table

  expect_valid_msi(out)
  expect_true(nrow(info) >= 1)
  expect_true(any(out@variable_info$ms2 == "Yes"))
  expect_true(any(out@variable_info$has_ms2 == "yes"))
  expect_true(any(!is.na(out@variable_info$ms2_id)))
  expect_true(any(!is.na(out@variable_info$ms2_type)))
  expect_true(all(c(
    "ms2_id", "precursor_mz", "precursor_variable_id",
    "fragment_mz", "fragment_intensity", "ms2_type", "source_file"
  ) %in% colnames(tbl)))
  expect_true(all(tbl$fragment_intensity > 0 & tbl$fragment_intensity <= 1))
  expect_true(all(tbl$fragment_mz <= tbl$precursor_mz))
  expect_equal(last_process_fun(out), "mutate_ms2()")
})

test_that("mutate_msi adds ms2 metadata correctly", {
  object <- make_test_msi_object()
  mgf <- make_temp_mgf()
  out <- mutate_msi(object, method = "ms2", path = mgf, ms1.ms2.match.mz.tol = 20)
  expect_valid_msi(out)
  expect_true(any(out@variable_info$ms2 == "Yes"))
  lg <- last_process_log(out)
  expect_equal(lg@function_name, "mutate_msi()")
  expect_true(all(c(
    "n_variables_before", "n_samples_before", "n_variables_after", "n_samples_after"
  ) %in% names(lg@parameter)))
})

test_that("mutate_ms2 repeated run is well-defined", {
  object <- make_test_msi_object()
  mgf <- make_temp_mgf()

  out1 <- mutate_ms2(object, path = mgf, ms1.ms2.match.mz.tol = 20)
  out2 <- mutate_ms2(out1, path = mgf, ms1.ms2.match.mz.tol = 20)

  info1 <- extract_ms2(out1, what = "ms2_info")
  info2 <- extract_ms2(out2, what = "ms2_info")
  expect_equal(nrow(info1), nrow(info2))
  expect_equal(out1@variable_info$ms2, out2@variable_info$ms2)
})

test_that("mutate_msi validates bad inputs", {
  object <- make_test_msi_object()
  expect_error(mutate_msi("bad", method = "ms2"), "msi_dataset")
  expect_error(mutate_msi(object, method = "bad"), "arg")
})

test_that("mz_match_ppm works as expected", {
  mz <- c(191.019, 191.020, 191.050)
  res <- mz_match_ppm(mz, precursor_mz = 191.020, tolerance_ppm = 10)
  expect_equal(length(res), length(mz))
  expect_true(any(res))
  expect_false(res[[3]])
})

test_that("mutate_ms2 preprocessing removes invalid fragments and normalizes", {
  object <- make_test_msi_object()
  mgf <- tempfile(fileext = ".mgf")
  writeLines(c(
    "BEGIN IONS",
    "TITLE=test_preprocess",
    "PEPMASS=191.0200",
    "SCANS=spectrum_1",
    "100.0 0",
    "120.0 50",
    "191.0200 1000",
    "192.0 10",
    "END IONS"
  ), mgf)

  out <- mutate_ms2(object, path = mgf, ms1.ms2.match.mz.tol = 20, precursor_tolerance = 0)
  tbl <- out@ms2_data$ms2_table

  expect_true(nrow(tbl) >= 1)
  expect_false(any(tbl$fragment_intensity == 0))
  expect_false(any(abs(tbl$fragment_mz - tbl$precursor_mz) < 1e-8))
  expect_false(any(tbl$fragment_mz > tbl$precursor_mz))
  expect_equal(max(tbl$fragment_intensity), 1)
})

test_that("mutate_ms2 tolerance matching no-hit behavior is stable", {
  object <- make_test_msi_object()
  mgf <- make_temp_mgf(precursor_mz = c(9999.0), scans = c("spectrum_1"))

  out <- mutate_ms2(object, path = mgf, ms1.ms2.match.mz.tol = 1)

  expect_true("ms2_table" %in% names(out@ms2_data))
  expect_equal(nrow(out@ms2_data$ms2_table), 0)
  expect_true(all(out@variable_info$ms2 == "no"))
  expect_true(all(out@variable_info$has_ms2 == "no"))
})

test_that("filter_ms2 and extract_ms2 work together", {
  object <- make_test_msi_object()
  mgf <- make_temp_mgf()
  object <- mutate_ms2(object, path = mgf, ms1.ms2.match.mz.tol = 20)

  out <- filter_ms2(object, precursor_mz > 190, precursor_mz < 200)
  info <- extract_ms2(out, what = "ms2_info")

  expect_valid_msi(out)
  expect_true(nrow(info) >= 1)
  expect_true(all(info$precursor_mz > 190 & info$precursor_mz < 200))
  expect_equal(last_process_fun(out), "filter_ms2()")
})

test_that("extract_ms2 supports variable_id and mz interfaces", {
  object <- make_test_msi_object()
  mgf <- make_temp_mgf()
  out <- mutate_ms2(object, path = mgf, ms1.ms2.match.mz.tol = 20)

  var_id <- out@variable_info$variable_id[which.min(abs(out@variable_info$mz - 191.02))]
  spec_by_id <- extract_ms2(out, variable_id = var_id)
  spec_by_mz <- extract_ms2(out, mz = 191.02, tolerance_ppm = 20)

  expect_true(is.list(spec_by_id))
  expect_true(all(c("precursor_mz", "fragments", "intensity", "ms2_type", "metadata") %in% names(spec_by_id)))
  expect_true(length(spec_by_id$fragments) == length(spec_by_id$intensity))
  expect_true(is.list(spec_by_mz))
})

test_that("infer_pseudo_ms2 adds pseudo entries", {
  object <- make_test_msi_object(n_features = 5, n_samples = 8, n_sections = 2)
  out <- infer_pseudo_ms2(
    object = object,
    variable_id = object@variable_info$variable_id[3],
    correlation_threshold = -1
  )

  tbl <- out@ms2_data$ms2_table
  expect_true(nrow(tbl) > 0)
  expect_true(any(tbl$ms2_type %in% c("pseudo", "pseudo_candidate")))
  expect_true(any(grepl("pseudo", out@variable_info$ms2_type), na.rm = TRUE))
  expect_equal(last_process_fun(out), "infer_pseudo_ms2()")
})

test_that("mutate_msi handles single pixel object", {
  object <- make_test_msi_object(n_features = 3, n_samples = 2, n_sections = 2)
  one <- object[, 1, drop = FALSE]
  mgf <- make_temp_mgf(precursor_mz = c(191.02), scans = c(colnames(one@expression_data)[1]))
  out <- mutate_msi(one, method = "ms2", path = mgf, ms1.ms2.match.mz.tol = 20)
  expect_valid_msi(out)
})
