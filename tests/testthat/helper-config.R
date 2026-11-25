locked_d15_path <- function() {
  p <- testthat::test_path("..", "..", "inst", "config", "locked_d15.yaml")
  if (file.exists(p)) return(p)
  system.file("config", "locked_d15.yaml", package = "msidataset")
}

locked_d15 <- function() {
  .msi_locked_resources(locked_d15_path())
}

d15_fixture_root <- function(variable) {
  path <- Sys.getenv(variable)
  skip_if_not(nzchar(path) && dir.exists(path), paste("Set", variable, "to the fixture directory"))
  path
}

d15_require_files <- function(paths) {
  skip_if_not(all(nzchar(paths) & file.exists(paths)), "External D15 fixtures are unavailable")
}
