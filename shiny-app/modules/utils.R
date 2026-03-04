app_project_root <- function() {
  wd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

  if (basename(wd) == "shiny-app") {
    return(dirname(wd))
  }

  if (dir.exists(file.path(wd, "shiny-app"))) {
    return(wd)
  }

  wd
}

format_bytes <- function(size_bytes) {
  if (is.na(size_bytes) || !is.finite(size_bytes)) {
    return("unknown size")
  }
  units <- c("B", "KB", "MB", "GB", "TB")
  value <- as.numeric(size_bytes)
  idx <- 1
  while (value >= 1024 && idx < length(units)) {
    value <- value / 1024
    idx <- idx + 1
  }
  sprintf("%.2f %s", value, units[[idx]])
}

find_demo_rds <- function() {
  root <- app_project_root()
  candidates <- c(
    file.path(root, "data", "demo data", "output1", "03-total ion count_msidataset.rds"),
    file.path(root, "data", "output", "my_msi_dataset_full.rds"),
    file.path(root, "data", "output", "my_msi_dataset_03.rds")
  )

  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    return(NULL)
  }
  hit[[1]]
}

find_demo_imzml_dir <- function() {
  root <- app_project_root()
  candidate <- file.path(root, "data", "demo data")
  if (dir.exists(candidate)) candidate else NULL
}

read_table_flexible <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv")) {
    return(utils::read.csv(path, row.names = 1, check.names = FALSE))
  }
  utils::read.delim(path, row.names = 1, check.names = FALSE)
}

validate_msi_object <- function(obj) {
  methods::is(obj, "msi_dataset")
}

ensure_tidymass_parameter_class <- function() {
  if (!requireNamespace("massdataset", quietly = TRUE)) {
    return(FALSE)
  }

  ns <- tryCatch(loadNamespace("massdataset"), error = function(e) NULL)
  if (is.null(ns)) {
    return(FALSE)
  }

  class_def <- tryCatch(
    methods::getClassDef("tidymass_parameter", where = ns),
    error = function(e) NULL
  )

  !is.null(class_def)
}

initialize_msidataset_runtime <- function() {
  if (!requireNamespace("methods", quietly = TRUE)) {
    return(list(ok = FALSE, message = "methods package is required."))
  }

  if (!requireNamespace("massdataset", quietly = TRUE)) {
    return(list(
      ok = FALSE,
      message = "Package 'massdataset' is required. Please install/load massdataset."
    ))
  }

  class_ok <- ensure_tidymass_parameter_class()
  if (!class_ok) {
    return(list(
      ok = FALSE,
      message = paste(
        "Class 'tidymass_parameter' is not available from massdataset namespace.",
        "Please reinstall/update massdataset."
      )
    ))
  }

  list(ok = TRUE, message = "msidataset runtime initialized.")
}

detect_imzml_pair <- function(folder_path) {
  if (is.null(folder_path) || !nzchar(trimws(folder_path))) {
    return(list(ok = FALSE, message = "Please provide a folder path."))
  }

  path <- normalizePath(folder_path, winslash = "/", mustWork = FALSE)
  if (!dir.exists(path)) {
    return(list(ok = FALSE, message = paste0("Folder does not exist: ", folder_path)))
  }

  imzml_files <- list.files(path, pattern = "\\.imzML$", ignore.case = TRUE, full.names = TRUE)
  if (length(imzml_files) == 0) {
    return(list(ok = FALSE, message = "No .imzML file found in the folder."))
  }

  matched <- NULL
  for (imzml_file in imzml_files) {
    ibd_file <- sub("\\.imzML$", ".ibd", imzml_file, ignore.case = TRUE)
    if (file.exists(ibd_file)) {
      matched <- list(imzml = imzml_file, ibd = ibd_file)
      break
    }
  }

  if (is.null(matched)) {
    return(list(
      ok = FALSE,
      message = "No matching .ibd found. The .imzML and .ibd filenames must share the same base name."
    ))
  }

  ibd_info <- file.info(matched$ibd)
  size_label <- format_bytes(ibd_info$size)

  list(
    ok = TRUE,
    folder = path,
    imzml = matched$imzml,
    ibd = matched$ibd,
    ibd_size = ibd_info$size,
    message = paste0("Detected pair: ", basename(matched$imzml), " + ", basename(matched$ibd),
                     " (", size_label, ")")
  )
}

safe_section_choices <- function(obj) {
  if (is.null(obj)) {
    return("All")
  }
  sec <- try(msidataset::get_section_ids(obj), silent = TRUE)
  if (inherits(sec, "try-error") || is.null(sec) || length(sec) == 0) {
    return("All")
  }
  c("All", as.character(sec))
}

safe_variable_ids <- function(obj, max_n = 5000) {
  if (is.null(obj)) {
    return(character())
  }
  vid <- try(msidataset::get_variable_id(obj), silent = TRUE)
  if (inherits(vid, "try-error") || length(vid) == 0) {
    return(character())
  }
  as.character(utils::head(vid, max_n))
}

safe_sample_ids <- function(obj, max_n = 5000) {
  if (is.null(obj)) {
    return(character())
  }
  sid <- try(msidataset::get_sample_id(obj), silent = TRUE)
  if (inherits(sid, "try-error") || length(sid) == 0) {
    return(character())
  }
  as.character(utils::head(sid, max_n))
}
