.msi_locked_resources <- function(lock_path) {
  resources <- yaml::read_yaml(lock_path)$resources
  lapply(resources, function(entry) {
    if (is.null(entry$root)) return(entry)
    root <- Sys.getenv(entry$root)
    entry$path <- if (!nzchar(root)) "" else if (entry$path == ".") root else file.path(root, entry$path)
    entry
  })
}

.msi_check_locked_ms1_path <- function(configured_path, entry) {
  if (identical(normalizePath(configured_path, mustWork = FALSE),
                normalizePath(entry$path))) return(invisible(NULL))
  # A hash-verified configuration can still name the original library location.
  root <- if (is.null(entry$root)) "" else Sys.getenv(entry$root)
  if (nzchar(root) && !is.null(configured_path) && length(configured_path) == 1L) {
    relative_path <- substring(entry$path, nchar(root) + 1L)
    relative_path <- paste0("/", sub("^/+", "", relative_path))
    if (endsWith(gsub("\\\\", "/", configured_path), relative_path)) return(invisible(NULL))
  }
  stop("ADP MS1 library does not match locked resource", call. = FALSE)
}
