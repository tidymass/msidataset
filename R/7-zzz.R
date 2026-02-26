# Version and date information - extracted from DESCRIPTION file during package build
msidataset_version <- utils::packageVersion("msidataset")
update_date <- Sys.Date()

.onAttach <- function(libname, pkgname) {
  needed <- msidataset_core[!is_attached(msidataset_core)]
  if (length(needed) == 0)
    return()
  
  crayon::num_colors(TRUE)
  msidataset_attach()
  
  # if (!"package:conflicted" %in% search()) {
  #   x <- msidataset_conflicts()
  #   msg(msidataset_conflict_message(x), startup = TRUE)
  # }
  packageStartupMessage(paste0("msidataset ", msidataset_version, " (", update_date, ')'))
}

is_attached <- function(x) {
  paste0("package:", x) %in% search()
}

convert_ms2_mzxml2mgf <- function(x) {
  purrr::map(x, function(y) {
    y$info <-
      y$info[, c("mz", "rt")]
    y
  })
}