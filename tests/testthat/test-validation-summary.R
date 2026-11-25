validation_dir <- function() {
  local <- testthat::test_path("..", "..", "inst", "validation")
  if (dir.exists(local)) return(local)
  system.file("validation", package = "msidataset", mustWork = TRUE)
}

validation_tools <- function() {
  env <- new.env(parent = globalenv())
  sys.source(file.path(validation_dir(), "recompute_validation.R"), envir = env)
  env
}

test_that("validation summary follows its TSV row schema", {
  tool <- validation_tools()
  schema <- jsonlite::read_json(file.path(validation_dir(), "validation_summary.tsv.schema.json"),
                                simplifyVector = TRUE)
  table <- tool$validation_read_summary(file.path(validation_dir(), "validation_summary.tsv"))
  expect_identical(names(table), schema$items$required)
  expect_gte(nrow(table), schema$minItems)
  expect_false(anyDuplicated(table$validation_id) > 0L)
  for (key in names(table)) {
    rule <- schema$items$properties[[key]]
    expect_identical(rule$type, "string")
    expect_type(table[[key]], "character")
    if (!is.null(rule$minLength)) expect_true(all(nchar(table[[key]]) >= rule$minLength), info = key)
    if (!is.null(rule$enum)) expect_true(all(table[[key]] %in% rule$enum), info = key)
    if (!is.null(rule$pattern)) expect_true(all(grepl(rule$pattern, table[[key]], perl = TRUE)), info = key)
  }
  expect_identical(lengths(strsplit(table$source_file, "; ", fixed = TRUE)),
                   lengths(strsplit(table$source_sha256, "; ", fixed = TRUE)))
})

test_that("validation rows recompute exactly from available source files", {
  tool <- validation_tools()
  root <- Sys.getenv("MSIDATASET_VALIDATION_ROOT")
  skip_if_not(nzchar(root) && dir.exists(root), "Set MSIDATASET_VALIDATION_ROOT to rerun external validation")
  roots <- tool$validation_roots()
  table <- tool$validation_read_summary(file.path(validation_dir(), "validation_summary.tsv"))
  sources <- strsplit(table$source_file, "; ", fixed = TRUE)
  paths <- unique(unlist(sources))
  missing <- paths[!file.exists(tool$validation_paths(paths, roots))]
  optional <- intersect(missing, tool$validation_optional_sources())
  required <- setdiff(missing, optional)
  skip_if_not(!length(required), "The validation archive is missing required sources")
  affected <- vapply(sources, function(s) any(s %in% optional), TRUE)
  omitted <- table$validation_id[affected]
  actual <- tool$recompute_validation(roots, omit = omitted)
  expect_equal(actual$summary, table[!affected, ], ignore_attr = "row.names",
               tolerance = 0)
  old_differences <- tool$validation_read_summary(file.path(validation_dir(), "discrepancies.tsv"))
  expect_equal(actual$discrepancies, old_differences[!old_differences$validation_id %in% omitted, ],
               ignore_attr = "row.names", tolerance = 0)
  if (length(optional)) skip(paste("Optional sources unavailable. Other rows were verified:", paste(optional, collapse = ", ")))
})
