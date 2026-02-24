msidataset_core <- c(
  "masstools",
  "ggplot2",
  "dplyr",
  "magrittr"
)

#' Check Unloaded Core Packages
#'
#' Identify which core dependencies of msidataset
#' are not yet attached.
#'
#' @return Character vector of unloaded packages.
#' @keywords internal
msidataset_core_unloaded <- function() {
  search_path <- paste0("package:", msidataset_core)
  msidataset_core[!search_path %in% search()]
}

# =====================================================
# attach helper
# =====================================================

#' Attach Package from Same Library
#'
#' Ensures dependent packages are loaded from the
#' same library location as msidataset.
#'
#' @param pkg Package name.
#' @keywords internal
same_library <- function(pkg) {
  loc <-
    if (pkg %in% loadedNamespaces())
      dirname(getNamespaceInfo(pkg, "path"))
  
  do.call(
    "library",
    list(
      pkg,
      lib.loc = loc,
      character.only = TRUE,
      warn.conflicts = FALSE
    )
  )
}

# =====================================================
# main attach
# =====================================================

#' Attach msidataset Core Dependencies
#'
#' Automatically attaches required core packages
#' when msidataset is loaded.
#'
#' @return Invisible NULL.
#' @keywords internal
#' @importFrom utils packageVersion
#' @importFrom cli symbol
#' @importFrom crayon green blue col_align col_nchar red
msidataset_attach <- function() {
  to_load <- msidataset_core_unloaded()
  
  if (length(to_load) == 0)
    return(invisible())
  
  versions <-
    vapply(to_load, msidataset_package_version, character(1))
  
  packages <- paste0(
    crayon::green(cli::symbol$tick),
    " ",
    crayon::blue(format(to_load)),
    " ",
    crayon::col_align(
      versions,
      max(crayon::col_nchar(versions))
    )
  )
  
  # 两列排版（对齐 massdataset）
  if (length(packages) %% 2 == 1) {
    packages <- append(packages, "")
  }
  
  col1 <- seq_len(length(packages) / 2)
  info <- paste0(packages[col1], "     ", packages[-col1])
  
  suppressPackageStartupMessages(
    lapply(to_load, same_library)
  )
  
  invisible()
}

# =====================================================
# version helper
# =====================================================

#' Format Package Version
#'
#' Pretty formatting of package version numbers.
#'
#' @param x Package name.
#' @return Character string.
#' @keywords internal
msidataset_package_version <- function(x) {
  version <- as.character(unclass(utils::packageVersion(x))[[1]])
  
  if (length(version) > 3) {
    version[4:length(version)] <-
      crayon::red(as.character(version[4:length(version)]))
  }
  
  paste0(version, collapse = ".")
}