test_that("msi_show_flow returns deterministic stage flow", {
  flow <- msi_show_flow(verb = "filter")
  expect_s3_class(flow, "msi_flow")
  expect_equal(nrow(flow$nodes), 6)
  expect_equal(nrow(flow$edges), 5)
  expect_equal(
    tolower(gsub("^[0-9]+\\.\\s*", "", flow$nodes$label)),
    c(
      "input validation",
      "data alignment",
      "core computation",
      "slot updates",
      "process log update",
      "output object"
    )
  )
})

test_that("msi_show_flow validates input object class", {
  expect_error(msi_show_flow(x = "bad", verb = "filter"), "msi_dataset")
})
