.msi_reference_lock <- function() {
  system.file("config", "reference_lock.yaml", package = "msidataset", mustWork = TRUE)
}

.msi_reference_root <- function() {
  root <- Sys.getenv("MSIDATASET_REFERENCE_ROOT")
  if (nzchar(root)) return(root)
  .msi_deprecated_env("MSIDATASET_VALIDATION_ROOT")
}

.msi_deprecated_env <- function(name, replacement = "MSIDATASET_REFERENCE_ROOT") {
  value <- Sys.getenv(name)
  if (nzchar(value)) message(name, " is deprecated. Use ", replacement, ".")
  value
}

.msi_locked_resources <- function(lock_path = .msi_reference_lock()) {
  resources <- yaml::read_yaml(lock_path)$resources
  root <- Sys.getenv("MSIDATASET_REFERENCE_ROOT")
  if (nzchar(root)) return(lapply(resources, function(entry) {
    entry$path <- file.path(root, entry$path)
    entry
  }))
  .msi_legacy_reference_paths(resources, .msi_reference_root())
}

.msi_file_record <- function(path, sha256 = digest::digest(file = path, algo = "sha256")) {
  list(basename = basename(path), sha256 = sha256)
}

.msi_resource_entries <- function(resources) {
  if (is.null(resources)) return(.msi_locked_resources())
  root <- "."
  if (is.character(resources) && length(resources) == 1L && !is.na(resources)) {
    if (!file.exists(resources) || dir.exists(resources)) {
      stop("Resource YAML file not found", call. = FALSE)
    }
    root <- dirname(resources)
    resources <- yaml::read_yaml(resources)
  }
  if (!is.list(resources) || is.null(names(resources))) {
    stop("resources must be a named list or a YAML file", call. = FALSE)
  }
  .msi_check_ion_mode(resources$ion_mode)
  .msi_check_ion_mode(resources$polarity)
  if (!is.null(resources$resources)) {
    resources <- resources$resources
    if (!is.list(resources) || is.null(names(resources))) {
      stop("The resource YAML must contain named entries", call. = FALSE)
    }
    .msi_check_ion_mode(resources$ion_mode)
    .msi_check_ion_mode(resources$polarity)
  }
  resources$ion_mode <- NULL
  resources$polarity <- NULL
  entries <- lapply(resources, function(entry) {
    if (is.character(entry)) entry <- list(path = entry)
    if (!is.list(entry) || !is.character(entry$path) ||
        length(entry$path) != 1L || is.na(entry$path) || !nzchar(entry$path)) {
      stop("Each resource requires a file or directory path", call. = FALSE)
    }
    .msi_check_ion_mode(entry$ion_mode)
    .msi_check_ion_mode(entry$polarity)
    if (!grepl("^(/|[A-Za-z]:|~)", entry$path)) entry$path <- file.path(root, entry$path)
    entry$path <- path.expand(entry$path)
    entry
  })
  if (is.null(entries$ms1_project_root) && !is.null(entries$ms1_library)) {
    entries$ms1_project_root <- list(path = dirname(entries$ms1_library$path))
  }
  entries
}

.msi_check_ion_mode <- function(mode) {
  if (!is.null(mode) && !identical(tolower(as.character(mode)), "negative")) {
    stop("Only negative ion mode is supported ([M-H]-, [M+Cl]-)", call. = FALSE)
  }
}

.msi_verify_resources <- function(entries, keys, locked) {
  records <- lapply(keys, function(key) {
    entry <- entries[[key]]
    if (is.null(entry$path) || !nzchar(entry$path) || !file.exists(entry$path) ||
        dir.exists(entry$path) != identical(key, "ms1_project_root")) {
      hint <- if (locked) ". Set MSIDATASET_REFERENCE_ROOT to the reference archive." else "."
      stop("Missing resource: ", key, hint, call. = FALSE)
    }
    if (dir.exists(entry$path)) return(list(basename = basename(entry$path)))
    sha256 <- digest::digest(file = entry$path, algo = "sha256")
    if ((locked || !is.null(entry$sha256)) && !identical(sha256, entry$sha256)) {
      stop("Resource SHA-256 mismatch: ", key, call. = FALSE)
    }
    .msi_file_record(entry$path, sha256)
  })
  stats::setNames(records, keys)
}

.msi_check_locked_ms1_path <- function(configured_path, entry) {
  if (identical(normalizePath(configured_path, mustWork = FALSE), normalizePath(entry$path)) ||
      endsWith(gsub("\\\\", "/", configured_path), paste0("/derived_library/", basename(entry$path)))) {
    return(invisible(NULL))
  }
  stop("MS1 library does not match the pinned resource", call. = FALSE)
}
