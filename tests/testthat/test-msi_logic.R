library(testthat)

# ----------------------------------------------------------------------------
# read_imzml mock-import logic
# ----------------------------------------------------------------------------

test_that("mock read pipeline populates subject_id and section_id", {
  expr <- matrix(
    c(10, 0, 5,
      2,  3, 1),
    nrow = 2,
    byrow = TRUE
  )

  coords <- cbind(
    x = c(1, 2, 3),
    y = c(1, 1, 2),
    z = c(1, 1, 1)
  )

  mz <- c(100.1, 200.2)

  object <- msidataset:::.build_msi_dataset_from_components(
    expression_data = expr,
    coords = coords,
    mz_values = mz,
    sample_id = c("pixel_1", "pixel_2", "pixel_3"),
    subject_id = "subject_A",
    section_id = "section_01"
  )

  expect_true(methods::is(object, "msi_dataset"))
  expect_true(all(c("subject_id", "section_id") %in% colnames(object@sample_info)))
  expect_equal(unique(object@sample_info$subject_id), "subject_A")
  expect_equal(unique(object@sample_info$section_id), "section_01")
})

# ----------------------------------------------------------------------------
# impute_mv logic + process_info logging
# ----------------------------------------------------------------------------

test_that("impute_mv logs method and preserves matrix dimensions", {
  expr <- matrix(
    c(1, NA, 3,
      4,  5, NA,
      7,  8,  9),
    nrow = 3,
    byrow = TRUE
  )

  sample_info <- data.frame(
    sample_id = c("p1", "p2", "p3"),
    x = c(1, 2, 3),
    y = c(1, 1, 2),
    z = c(1, 1, 1),
    subject_id = c("s1", "s1", "s1"),
    section_id = c("sec1", "sec1", "sec1"),
    resolution = c(50, 50, 50),
    matrix = c("CHCA", "CHCA", "CHCA"),
    stringsAsFactors = FALSE
  )

  variable_info <- data.frame(
    variable_id = c("m1", "m2", "m3"),
    mz = c(100.1, 101.1, 102.1),
    ms2 = c("no", "no", "no"),
    stringsAsFactors = FALSE
  )

  colnames(expr) <- sample_info$sample_id
  rownames(expr) <- variable_info$variable_id

  object <- create_msi_dataset(
    expression_data = expr,
    sample_info = sample_info,
    variable_info = variable_info
  )

  original_dim <- dim(object@expression_data)

  object2 <- impute_mv(object, method = "median", by = "variable")

  expect_equal(dim(object2@expression_data), original_dim)

  last_log <- object2@process_info[[length(object2@process_info)]]
  expect_s4_class(last_log, "tidymass_parameter")
  expect_equal(last_log@function_name, "impute_mv()")
  expect_equal(last_log@parameter$method, "median")
  expect_true(inherits(last_log@time, c("POSIXct", "POSIXt")))
})

# ----------------------------------------------------------------------------
# S4 data validation: missing spatial columns
# ----------------------------------------------------------------------------

test_that("create_msi_dataset fails when x or y is missing", {
  expr <- matrix(c(1, 2, 3, 4), nrow = 2)

  variable_info <- data.frame(
    variable_id = c("m1", "m2"),
    mz = c(100.1, 200.2),
    ms2 = c("no", "no"),
    stringsAsFactors = FALSE
  )

  sample_info_missing_x <- data.frame(
    sample_id = c("p1", "p2"),
    y = c(1, 2),
    z = c(1, 1),
    subject_id = c("s1", "s1"),
    section_id = c("sec1", "sec1"),
    resolution = c(50, 50),
    matrix = c("CHCA", "CHCA"),
    stringsAsFactors = FALSE
  )

  sample_info_missing_y <- data.frame(
    sample_id = c("p1", "p2"),
    x = c(1, 2),
    z = c(1, 1),
    subject_id = c("s1", "s1"),
    section_id = c("sec1", "sec1"),
    resolution = c(50, 50),
    matrix = c("CHCA", "CHCA"),
    stringsAsFactors = FALSE
  )

  expect_error(
    create_msi_dataset(expr, sample_info_missing_x, variable_info),
    "sample_info missing required columns: x"
  )

  expect_error(
    create_msi_dataset(expr, sample_info_missing_y, variable_info),
    "sample_info missing required columns: y"
  )
})
