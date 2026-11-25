reference_lock_path <- function() {
  p <- testthat::test_path("..", "..", "inst", "config", "reference_lock.yaml")
  if (file.exists(p)) return(p)
  system.file("config", "reference_lock.yaml", package = "msidataset")
}

reference_lock <- function() {
  .msi_locked_resources(reference_lock_path())
}

reference_fixture_root <- function() {
  path <- .msi_reference_root()
  skip_if_not(nzchar(path) && dir.exists(path), "Set MSIDATASET_REFERENCE_ROOT to the fixture archive")
  path
}

reference_require_files <- function(paths) {
  skip_if_not(all(nzchar(paths) & file.exists(paths)), "External reference fixtures are unavailable")
}
