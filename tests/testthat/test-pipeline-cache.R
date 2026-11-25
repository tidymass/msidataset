test_that("pipeline checks input content before reusing a checkpoint", {
  dir <- tempfile(); dir.create(dir); on.exit(unlink(dir, recursive = TRUE))
  input <- file.path(dir, "pixels.csv")
  cache <- file.path(dir, "cache")
  writeLines(c("1,1,10,2", "2,1,11,3"), input)
  cfg <- list(input = list(path = input, format = "csv", mz = 10, feature_id = "F003"),
              labels = data.frame(x = c(1, 2), y = c(1, 1), region = c("periportal", "pericentral")))
  x <- run_pipeline(cfg, to = "read_msi_data", cache_dir = cache)
  y <- run_pipeline(cfg, to = "read_msi_data", cache_dir = cache)
  expect_identical(x$cache_keys, y$cache_keys)
  resumed <- run_pipeline(cfg, from = "define_roi", to = "define_roi", cache_dir = cache)
  expect_identical(resumed$cache_keys$read_msi_data, x$cache_keys$read_msi_data)
  expect_equal(as.character(resumed$object@sample_info$region), cfg$labels$region)
  writeLines(c("1,1,10,4", "2,1,11,3"), input)
  expect_error(run_pipeline(cfg, from = "define_roi", to = "define_roi", cache_dir = cache), "Missing or stale")
  z <- run_pipeline(cfg, to = "read_msi_data", cache_dir = cache)
  expect_false(identical(x$cache_keys, z$cache_keys))
  expect_equal(z$object@expression_data[1, 1], 4)
  expect_error(run_pipeline(cfg, from = "define_roi", to = "read_msi_data", cache_dir = cache), "Invalid")
})

test_that("pipeline refuses to resume from missing earlier checkpoints", {
  dir <- tempfile(); dir.create(dir); on.exit(unlink(dir, recursive = TRUE))
  input <- file.path(dir, "pixels.csv"); writeLines("1,1,10,2", input)
  cfg <- list(input = list(path = input, format = "csv", mz = 10))
  expect_error(run_pipeline(cfg, from = "define_roi", to = "define_roi",
                            cache_dir = file.path(dir, "cache")), "Missing or stale")
})
