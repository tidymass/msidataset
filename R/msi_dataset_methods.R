# ============================================================================
# Core S4 Methods for msi_dataset (show, [, dim, dimnames, update)
# Author: Xiaoyu Xu <XIAOYU021@e.ntu.edu.sg>
# ============================================================================

# ---------------------------------------------------------------------------
# show()
# ---------------------------------------------------------------------------

#' Display a msi_dataset object
#'
#' @description
#' Prints a formatted summary of the object including dimensions, spatial
#' extent, section information, and processing history.
#'
#' @param object A \code{msi_dataset} object.
#' @return Invisibly returns the object.
#'
#' @docType methods
#' @name show,msi_dataset-method
#' @rdname show-methods
#' @aliases show,msi_dataset-method
#' @export

setMethod(
  f          = "show",
  signature  = "msi_dataset",
  definition = function(object) {

    ## ---- re-validate ----
    check_result <- check_msi_dataset_class(object)
    if (!isTRUE(check_result)) {
      cat(crayon::red("Validation errors:\n"))
      cat(crayon::red(paste(check_result, collapse = "\n"), "\n"))
      cat(crayon::red(
        "Slots may have been manually modified. ",
        "Try update_msi_dataset().\n"
      ))
    }

    sep <- paste(rep("-", 20), collapse = "")
    cat(crayon::yellow(sep, "\n"))
    cat(crayon::green("msidataset version:", object@version, "\n"))
    cat(crayon::yellow("Structure: massdataset-compatible\n"))
    cat(crayon::yellow("Activated slot:", object@activated, "\n"))
    cat(crayon::yellow(sep, "\n"))

    ## 1. expression_data
    cat(crayon::green("1. expression_data: "))
    data_type <- if (methods::is(object@expression_data, "DelayedArray")) {
      "DelayedArray"
    } else if (is.matrix(object@expression_data)) {
      "matrix"
    } else {
      "data.frame"
    }
    cat("[", nrow(object@expression_data), "x",
        ncol(object@expression_data), data_type, "]\n")

    ## 2. variable_info
    cat(crayon::green("2. variable_info: "))
    cat("[", nrow(object@variable_info), "x",
        ncol(object@variable_info), "data.frame]\n")
    nv <- nrow(object@variable_info)
    vid <- if (nv <= 5) {
      object@variable_info$variable_id
    } else {
      c(head(object@variable_info$variable_id, 3), "...",
        tail(object@variable_info$variable_id, 2))
    }
    cat(nv, "variables:", paste(vid, collapse = " "), "\n")

    ## m/z range
    if ("mz" %in% colnames(object@variable_info)) {
      mz_vals <- object@variable_info$mz
      cat("  m/z range: [",
          round(min(mz_vals), 4), " - ", round(max(mz_vals), 4), "]\n",
          sep = "")
    }

    ## 3. sample_info
    cat(crayon::green("3. sample_info: "))
    cat("[", nrow(object@sample_info), "x",
        ncol(object@sample_info), "data.frame]\n")
    np <- nrow(object@sample_info)
    sid <- if (np <= 5) {
      object@sample_info$sample_id
    } else {
      c(head(object@sample_info$sample_id, 3), "...",
        tail(object@sample_info$sample_id, 2))
    }
    cat(np, "pixels:", paste(sid, collapse = " "), "\n")

    ## spatial extent
    if (all(c("x", "y") %in% colnames(object@sample_info))) {
      cat("  Spatial extent: x[",
          min(object@sample_info$x), "-", max(object@sample_info$x), "], ",
          "y[", min(object@sample_info$y), "-", max(object@sample_info$y), "]",
          sep = "")
      if ("z" %in% colnames(object@sample_info)) {
        cat(", z[",
            min(object@sample_info$z), "-", max(object@sample_info$z), "]",
            sep = "")
      }
      cat("\n")
    }

    ## sections
    if ("section_id" %in% colnames(object@sample_info)) {
      secs <- unique(object@sample_info$section_id)
      if (length(secs) <= 5) {
        cat("  Sections:", paste(secs, collapse = ", "), "\n")
      } else {
        cat("  Sections:", length(secs), "sections\n")
      }
    }

    ## 4. sample_info_note
    cat(crayon::green("4. sample_info_note: "))
    cat("[", nrow(object@sample_info_note), "x",
        ncol(object@sample_info_note), "data.frame]\n")

    ## 5. variable_info_note
    cat(crayon::green("5. variable_info_note: "))
    cat("[", nrow(object@variable_info_note), "x",
        ncol(object@variable_info_note), "data.frame]\n")

    ## 6. ms2_data
    cat(crayon::green("6. ms2_data: "))
    if (length(object@ms2_data) > 0) {
      cat("[", length(object@ms2_data), "MS2 spectrum/spectra]\n")
    } else {
      cat("[0 MS2 spectra]\n")
    }

    ## 7. process_info
    cat(crayon::yellow(sep, "\n"))
    cat(crayon::green("7. Processing information\n"))
    pi_list <- object@process_info
    if (.hasSlot(object, "process_info") && length(pi_list) > 0) {
      cat(crayon::green(length(pi_list), "processing step(s) in total\n"))
      show_pi <- if (length(pi_list) > 5) {
        cat(crayon::green("Latest 3 steps shown:\n"))
        tail(pi_list, 3)
      } else {
        pi_list
      }
      for (nm in names(show_pi)) {
        cat(crayon::green(nm, paste(rep("-", 10), collapse = ""), "\n"))
        entry <- show_pi[[nm]]
        if (methods::is(entry, "tidymass_parameter")) {
          data.frame(
            Package       = entry@pacakge_name,
            `Function`    = entry@function_name,
            Time          = as.character(entry@time),
            check.names   = FALSE
          ) |> print()
        } else if (is.list(entry)) {
          data.frame(
            Package     = vapply(entry, function(x) x@pacakge_name,  character(1)),
            `Function`  = vapply(entry, function(x) x@function_name, character(1)),
            Time        = vapply(entry, function(x) as.character(x@time), character(1)),
            check.names = FALSE
          ) |> print()
        }
      }
    } else {
      cat(crayon::red("No processing information.\n"))
    }

    cat(crayon::yellow(sep, "\n"))
    invisible(object)
  }
)


# ---------------------------------------------------------------------------
# [ (subset)
# ---------------------------------------------------------------------------

#' Subset a msi_dataset object
#'
#' @description
#' Subset by variables (rows) and/or pixels (columns).
#'
#' @param x A \code{msi_dataset} object.
#' @param i Row indices (variables): numeric, logical, or character.
#' @param j Column indices (pixels): numeric, logical, or character.
#' @param drop Logical; default \code{FALSE}.
#' @return A subsetted \code{msi_dataset} object.
#'
#' @docType methods
#' @rdname extract-methods
#' @aliases [,msi_dataset,ANY,ANY,ANY-method
#' @export

setMethod(
  "[",
  signature = "msi_dataset",
  function(x, i, j, drop = FALSE) {

    if (missing(i)) i <- seq_len(nrow(x@expression_data))
    if (missing(j)) j <- seq_len(ncol(x@expression_data))

    new_expr  <- x@expression_data[i, j, drop = FALSE]
    new_vinfo <- x@variable_info[i, , drop = FALSE]
    new_sinfo <- x@sample_info[j, , drop = FALSE]

    ## Rebuild notes if needed (keep existing columns only)
    new_snote <- x@sample_info_note
    new_vnote <- x@variable_info_note

    new_obj <- new(
      Class              = "msi_dataset",
      expression_data    = as.data.frame(new_expr, check.names = FALSE),
      ms2_data           = x@ms2_data,
      annotation_table   = x@annotation_table,
      sample_info        = new_sinfo,
      variable_info      = new_vinfo,
      sample_info_note   = new_snote,
      variable_info_note = new_vnote,
      process_info       = x@process_info,
      other_files        = x@other_files,
      version            = x@version,
      activated          = x@activated
    )
    new_obj
  }
)


# ---------------------------------------------------------------------------
# dim() and dimnames()
# ---------------------------------------------------------------------------

#' Dimensions of a msi_dataset
#' @param x A \code{msi_dataset} object.
#' @return Integer vector \code{c(n_variables, n_pixels)}.
#' @docType methods
#' @rdname dim-methods
#' @aliases dim,msi_dataset-method
#' @export
setMethod("dim", "msi_dataset", function(x) {
  c(nrow(x@expression_data), ncol(x@expression_data))
})

#' Dimension names of a msi_dataset
#' @param x A \code{msi_dataset} object.
#' @return A list: variable IDs (rows) and pixel IDs (columns).
#' @docType methods
#' @rdname dimnames-methods
#' @aliases dimnames,msi_dataset-method
#' @export
setMethod("dimnames", "msi_dataset", function(x) {
  list(rownames(x@expression_data), colnames(x@expression_data))
})


# ---------------------------------------------------------------------------
# update_msi_dataset()
# ---------------------------------------------------------------------------

#' Update a msi_dataset Object
#'
#' @description
#' Re-synchronises \code{expression_data} dimnames with the id columns in
#' \code{variable_info} / \code{sample_info} and regenerates note tables
#' when they are missing or have the wrong number of rows.
#'
#' @param object A \code{msi_dataset} object.
#' @return The updated \code{msi_dataset} object (invisibly).
#' @export

update_msi_dataset <- function(object) {

  ## sync dimnames
  if (nrow(object@variable_info) > 0)
    rownames(object@expression_data) <- object@variable_info$variable_id
  if (nrow(object@sample_info) > 0)
    colnames(object@expression_data) <- object@sample_info$sample_id

  ## rebuild notes when stale
  if (nrow(object@sample_info_note) != ncol(object@sample_info) ||
      !all(object@sample_info_note$name == colnames(object@sample_info))) {
    object@sample_info_note <- data.frame(
      name    = colnames(object@sample_info),
      meaning = colnames(object@sample_info),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }

  if (nrow(object@variable_info_note) != ncol(object@variable_info) ||
      !all(object@variable_info_note$name == colnames(object@variable_info))) {
    object@variable_info_note <- data.frame(
      name    = colnames(object@variable_info),
      meaning = colnames(object@variable_info),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }

  ## record
  parameter <- new(
    Class         = "tidymass_parameter",
    pacakge_name  = "msidataset",
    function_name = "update_msi_dataset()",
    parameter     = list("no" = "no"),
    time          = Sys.time()
  )
  object@process_info$update_msi_dataset <- parameter

  ## validate
  check_result <- check_msi_dataset_class(object)
  if (!isTRUE(check_result))
    warning("Validation errors after update:\n",
            paste(check_result, collapse = "\n"))

  invisible(object)
}


# ---------------------------------------------------------------------------
# activate() — switch the "active" slot (mirrors massdataset)
# ---------------------------------------------------------------------------

#' Activate a slot for default operations
#'
#' @description
#' Sets \code{object@activated} to one of \code{"expression_data"},
#' \code{"sample_info"}, or \code{"variable_info"}, controlling the
#' default target of \code{dplyr::filter()}, \code{mutate()}, etc.
#'
#' @param object A \code{msi_dataset} object.
#' @param slot_name One of \code{"expression_data"}, \code{"sample_info"},
#'   \code{"variable_info"}.
#' @return The modified \code{msi_dataset} object (invisibly).
#' @export

activate_msi_dataset <- function(object,
                                  slot_name = c("expression_data",
                                                "sample_info",
                                                "variable_info")) {
  slot_name <- match.arg(slot_name)
  object@activated <- slot_name
  invisible(object)
}
