test_that("processed bins sum native peaks before counting positive pixels", {
  work <- tempfile()
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  xy <- cbind(x = 1:3, y = 1)
  mz <- list(c(100.00001, 100.00002, 200), c(100, 200, 300), c(100, 300))
  signal <- list(c(2, 3, 4), c(0, 8, 1), c(9, -1))
  run <- function(chunk_size, min_pixels = 2L, noise = 0) {
    .imzml_align(mz, signal, xy, work, bin_digits = 4L, max_features = 150000L,
      presence_fraction = 0, min_pixels = min_pixels, max_pixels_required = 1000L,
      noise_intensity_fraction = noise, chunk_size = chunk_size)
  }
  result <- run(1L)
  expect_identical(result, run(3L))
  expect_equal(result$mz, c(100, 200))
  expect_equal(result$columns, list(c(5, 4), c(0, 8), c(9, 0)))
  expect_equal(result$tic, c(9, 9, 8))
  expect_equal(result$variables$n_positive_pixels, c(2, 2))
  expect_equal(result$variables$n_contributing_centroids, c(4, 2))
  expect_equal(result$variables$observed_mz_full[1], (2 * mz[[1]][1] + 3 * mz[[1]][2] + 9 * 100) / 14)
  expect_equal(result$info$n_features_f1, 3)
  expect_equal(result$info$n_features_f2, 2)
  expect_equal(run(2L, min_pixels = 3L)$columns, rep(list(numeric()), 3))
  expect_equal(run(2L, noise = 1)$mz, 100)
})

test_that("legacy axis coarsens in successive steps and reports unassigned peaks", {
  work <- tempfile()
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  result <- .imzml_align(list(c(100.0046, 100.0051, 100.0012)), list(c(1, 2, 3)),
    cbind(x = 1, y = 1), work, bin_digits = 4L, max_features = 1L,
    presence_fraction = 0, min_pixels = 0L, max_pixels_required = 1000L,
    noise_intensity_fraction = 0, chunk_size = 1L)
  expect_equal(result$info$effective_bin_digits, 2)
  expect_equal(unname(result$info$bin_counts), c(3, 2, 1))
  expect_equal(result$info$n_unmatched_peaks, 1)
  expect_equal(result$mz, 100)
  expect_equal(result$columns, list(4))
})

test_that("continuous imzML keeps its measured axis and supports case variants", {
  work <- tempfile()
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  path <- file.path(work, "continuous.imzML")
  mz <- c(100.125, 200.375)
  values <- matrix(c(2, 3, 5, 7, 11, 13), nrow = 2L)
  image <- Cardinal::MSImagingExperiment(values, featureData = Cardinal::MassDataFrame(mz))
  suppressWarnings(Cardinal::writeMSIData(image, path))
  case_path <- file.path(work, "continuous.imzml")
  stopifnot(file.rename(path, case_path),
            file.rename(file.path(work, "continuous.ibd"), file.path(work, "continuous.IBD")))
  object <- read_msi_data(case_path, "imzml", chunk_size = 1L)
  expect_identical(object@process_info$read_msi_data$axis, "native")
  expect_equal(object@variable_info$mz, mz)
  expect_equal(unname(as.matrix(object@expression_data)), values)
  expect_equal(object@sample_info$tic, colSums(values))
  expect_equal(object@sample_info$x, 1:3)
  expect_equal(object@sample_info$y, rep(1L, 3))
  expect_error(read_msi_data(case_path, "imzml", chunk_size = 0), "chunk_size")
  unlink(file.path(work, "continuous.IBD"))
  expect_error(read_msi_data(case_path, "imzml"), "Paired imzML")
})

test_that("processed imzML defaults to aligned spectra with numeric coordinates", {
  work <- tempfile()
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  path <- file.path(work, "processed.imzML")
  mz <- c(100.125, 200.375)
  values <- matrix(c(2, 3, 5, 7, 11, 13), nrow = 2L)
  image <- Cardinal::MSImagingExperiment(values, featureData = Cardinal::MassDataFrame(mz),
                                         centroided = TRUE)
  processed <- methods::as(image, "MSProcessedImagingExperiment")
  suppressWarnings(Cardinal::writeMSIData(processed, path))
  object <- read_msi_data(path, "imzml", min_pixels = 1L, chunk_size = 1L)
  expect_identical(object@process_info$read_msi_data$axis, "aligned")
  expect_equal(object@variable_info$mz, mz)
  expect_equal(unname(as.matrix(object@expression_data)), values)
  expect_equal(object@sample_info$tic, colSums(values))
  expect_type(object@sample_info$x, "double")
  expect_equal(object@sample_info$x, 1:3)
})
