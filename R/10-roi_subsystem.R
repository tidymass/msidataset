# ============================================================================
# ROI subsystem for msi_dataset
# ============================================================================

#' @title S4 Class for a Single MSI ROI
#'
#' @slot roi_id Scalar character ROI identifier.
#' @slot roi_name Scalar character ROI display name.
#' @slot sample_id Character vector of member pixel/sample IDs.
#' @slot section_id Scalar character section identifier (or NA).
#' @slot metadata Named list for optional attributes.
#'
#' @rdname msi_roi-class
#' @name msi_roi-class
#' @docType class
#' @exportClass msi_roi
setClass(
  Class = "msi_roi",
  representation(
    roi_id = "character",
    roi_name = "character",
    sample_id = "character",
    section_id = "character",
    metadata = "list"
  ),
  prototype(
    roi_id = character(),
    roi_name = character(),
    sample_id = character(),
    section_id = NA_character_,
    metadata = list()
  ),
  validity = function(object) {
    errors <- character()
    if (length(object@roi_id) != 1L || !nzchar(object@roi_id)) {
      errors <- c(errors, "roi_id must be one non-empty string.")
    }
    if (length(object@roi_name) != 1L || !nzchar(object@roi_name)) {
      errors <- c(errors, "roi_name must be one non-empty string.")
    }
    if (anyNA(object@sample_id) || any(!nzchar(object@sample_id))) {
      errors <- c(errors, "sample_id cannot contain missing/empty values.")
    }
    if (length(unique(object@sample_id)) != length(object@sample_id)) {
      errors <- c(errors, "sample_id membership must be unique within one ROI.")
    }
    if (length(object@section_id) != 1L) {
      errors <- c(errors, "section_id must be length 1 (use NA if unknown).")
    }
    if (length(errors) == 0) TRUE else errors
  }
)

#' @title S4 Class for MSI ROI Collection
#'
#' @slot rois Named list of \code{msi_roi} objects.
#' @slot version Character ROI schema version.
#'
#' @rdname msi_roi_set-class
#' @name msi_roi_set-class
#' @docType class
#' @exportClass msi_roi_set
setClass(
  Class = "msi_roi_set",
  representation(
    rois = "list",
    version = "character"
  ),
  prototype(
    rois = list(),
    version = "0.1.0"
  ),
  validity = function(object) {
    errors <- character()
    rois <- object@rois
    if (length(rois) > 0) {
      if (is.null(names(rois)) || any(names(rois) == "")) {
        errors <- c(errors, "ROI set list must be named by roi_id.")
      }
      is_roi <- vapply(rois, function(x) methods::is(x, "msi_roi"), logical(1))
      if (!all(is_roi)) {
        errors <- c(errors, "All items in rois must be msi_roi objects.")
      } else {
        ids <- vapply(rois, function(x) x@roi_id, character(1))
        if (length(unique(ids)) != length(ids)) {
          errors <- c(errors, "roi_id must be unique in ROI set.")
        }
      }
    }
    if (length(errors) == 0) TRUE else errors
  }
)

.empty_roi_set <- function() {
  methods::new("msi_roi_set", rois = list(), version = "0.1.0")
}

.get_roi_set <- function(object, create = TRUE) {
  roi_obj <- object@other_files$roi
  if (is.null(roi_obj)) {
    if (create) return(.empty_roi_set())
    return(NULL)
  }
  if (methods::is(roi_obj, "msi_roi_set")) {
    return(roi_obj)
  }
  stop("other_files$roi exists but is not a valid msi_roi_set object.")
}

.set_roi_set <- function(object, roi_set) {
  methods::validObject(roi_set)
  object@other_files$roi <- roi_set
  object
}

.get_roi_membership <- function(object, roi_id = NULL) {
  roi_set <- .get_roi_set(object, create = FALSE)
  if (is.null(roi_set) || length(roi_set@rois) == 0) {
    return(data.frame(
      roi_id = character(),
      roi_name = character(),
      sample_id = character(),
      section_id = character(),
      stringsAsFactors = FALSE
    ))
  }

  rois <- roi_set@rois
  if (!is.null(roi_id)) {
    keep <- names(rois) %in% roi_id
    rois <- rois[keep]
  }
  if (length(rois) == 0) {
    return(data.frame(
      roi_id = character(),
      roi_name = character(),
      sample_id = character(),
      section_id = character(),
      stringsAsFactors = FALSE
    ))
  }

  do.call(
    rbind,
    lapply(rois, function(r) {
      data.frame(
        roi_id = r@roi_id,
        roi_name = r@roi_name,
        sample_id = r@sample_id,
        section_id = rep(r@section_id, length(r@sample_id)),
        stringsAsFactors = FALSE
      )
    })
  )
}

.validate_roi_integrity <- function(object, roi_set = NULL) {
  if (is.null(roi_set)) {
    roi_set <- .get_roi_set(object, create = FALSE)
  }
  if (is.null(roi_set)) {
    return(TRUE)
  }

  errors <- character()
  all_samples <- object@sample_info$sample_id
  rois <- roi_set@rois

  for (i in seq_along(rois)) {
    r <- rois[[i]]
    missing_samples <- setdiff(r@sample_id, all_samples)
    if (length(missing_samples) > 0) {
      errors <- c(
        errors,
        paste0("ROI '", r@roi_id, "' contains unknown sample_id: ",
               paste(head(missing_samples, 10), collapse = ", "))
      )
    }
    if (!is.na(r@section_id) && ("section_id" %in% colnames(object@sample_info))) {
      sec <- unique(object@sample_info$section_id[object@sample_info$sample_id %in% r@sample_id])
      if (length(sec) > 1) {
        errors <- c(errors, paste0("ROI '", r@roi_id, "' spans multiple sections."))
      } else if (length(sec) == 1 && sec != r@section_id) {
        errors <- c(
          errors,
          paste0("ROI '", r@roi_id, "' section mismatch: declared=", r@section_id,
                 ", observed=", sec, ".")
        )
      }
    }
  }

  if (length(errors) == 0) TRUE else errors
}

#' Create an ROI object
#'
#' @param roi_id ROI identifier.
#' @param sample_id Pixel/sample IDs in this ROI.
#' @param roi_name ROI display name. Default is \code{roi_id}.
#' @param section_id Section linkage (optional).
#' @param metadata Named list of optional metadata.
#'
#' @return A \code{msi_roi} object.
#' @export
create_roi <- function(roi_id,
                       sample_id,
                       roi_name = roi_id,
                       section_id = NA_character_,
                       metadata = list()) {
  methods::new(
    "msi_roi",
    roi_id = as.character(roi_id),
    roi_name = as.character(roi_name),
    sample_id = as.character(unique(sample_id)),
    section_id = as.character(section_id),
    metadata = metadata
  )
}

#' List ROI summary from msi_dataset
#'
#' @param object A \code{msi_dataset} object.
#'
#' @return A data.frame ROI summary.
#' @export
list_roi <- function(object) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }
  roi_set <- .get_roi_set(object, create = FALSE)
  if (is.null(roi_set) || length(roi_set@rois) == 0) {
    return(data.frame(
      roi_id = character(),
      roi_name = character(),
      n_pixels = integer(),
      section_id = character(),
      stringsAsFactors = FALSE
    ))
  }

  do.call(
    rbind,
    lapply(roi_set@rois, function(r) {
      data.frame(
        roi_id = r@roi_id,
        roi_name = r@roi_name,
        n_pixels = length(r@sample_id),
        section_id = r@section_id,
        stringsAsFactors = FALSE
      )
    })
  )
}

#' Extract ROI membership table
#'
#' @param object A \code{msi_dataset} object.
#' @param roi_id Optional ROI id filter.
#'
#' @return A long-format data.frame with columns
#'   \code{roi_id, roi_name, sample_id, section_id}.
#' @export
roi_data <- function(object, roi_id = NULL) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }
  .get_roi_membership(object, roi_id = roi_id)
}

#' Validate ROI integrity against msi_dataset sample_info
#'
#' @param object A \code{msi_dataset} object.
#'
#' @return \code{TRUE} if valid, otherwise character vector.
#' @export
validate_roi_integrity <- function(object) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }
  .validate_roi_integrity(object)
}

#' Add ROI into msi_dataset
#'
#' @param object A \code{msi_dataset} object.
#' @param roi A \code{msi_roi} object. If \code{NULL}, construct from arguments.
#' @param roi_id ROI identifier when \code{roi} is \code{NULL}.
#' @param sample_id ROI member sample IDs when \code{roi} is \code{NULL}.
#' @param roi_name ROI name when \code{roi} is \code{NULL}.
#' @param section_id ROI section linkage.
#' @param metadata Optional named list metadata.
#' @param overwrite Overwrite existing ROI with same id.
#'
#' @return Updated \code{msi_dataset} object.
#' @export
add_roi <- function(object,
                    roi = NULL,
                    roi_id = NULL,
                    sample_id = NULL,
                    roi_name = roi_id,
                    section_id = NA_character_,
                    metadata = list(),
                    overwrite = FALSE) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }

  if (is.null(roi)) {
    if (is.null(roi_id) || is.null(sample_id)) {
      stop("Provide roi or both roi_id and sample_id.")
    }
    roi <- create_roi(
      roi_id = roi_id,
      sample_id = sample_id,
      roi_name = roi_name,
      section_id = section_id,
      metadata = metadata
    )
  }
  if (!methods::is(roi, "msi_roi")) {
    stop("roi must be a msi_roi object.")
  }

  unknown <- setdiff(roi@sample_id, object@sample_info$sample_id)
  if (length(unknown) > 0) {
    stop("ROI contains unknown sample_id: ", paste(head(unknown, 10), collapse = ", "))
  }

  roi_set <- .get_roi_set(object, create = TRUE)
  exists_same <- roi@roi_id %in% names(roi_set@rois)
  if (exists_same && !isTRUE(overwrite)) {
    stop("ROI id already exists: ", roi@roi_id, ". Use overwrite = TRUE to replace.")
  }
  roi_set@rois[[roi@roi_id]] <- roi
  methods::validObject(roi_set)

  valid_roi <- .validate_roi_integrity(object, roi_set = roi_set)
  if (!isTRUE(valid_roi)) {
    stop(paste(valid_roi, collapse = "\n"))
  }

  object <- .set_roi_set(object, roi_set)
  object <- add_process_log(
    object = object,
    function_name = "add_roi()",
    parameter = list(
      roi_id = roi@roi_id,
      roi_name = roi@roi_name,
      n_pixels = length(roi@sample_id),
      overwrite = overwrite
    )
  )
  object
}

#' Remove ROI from msi_dataset
#'
#' @param object A \code{msi_dataset} object.
#' @param roi_id ROI identifier to remove.
#'
#' @return Updated \code{msi_dataset} object.
#' @export
remove_roi <- function(object, roi_id) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }
  roi_set <- .get_roi_set(object, create = TRUE)
  if (!(roi_id %in% names(roi_set@rois))) {
    stop("ROI id not found: ", roi_id)
  }
  roi_set@rois[[roi_id]] <- NULL
  object <- .set_roi_set(object, roi_set)
  object <- add_process_log(
    object = object,
    function_name = "remove_roi()",
    parameter = list(roi_id = roi_id)
  )
  object
}

#' Rename ROI in msi_dataset
#'
#' @param object A \code{msi_dataset} object.
#' @param roi_id ROI identifier.
#' @param new_name New ROI name.
#'
#' @return Updated \code{msi_dataset} object.
#' @export
rename_roi <- function(object, roi_id, new_name) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }
  roi_set <- .get_roi_set(object, create = TRUE)
  if (!(roi_id %in% names(roi_set@rois))) {
    stop("ROI id not found: ", roi_id)
  }
  roi_set@rois[[roi_id]]@roi_name <- as.character(new_name)
  object <- .set_roi_set(object, roi_set)
  object <- add_process_log(
    object = object,
    function_name = "rename_roi()",
    parameter = list(roi_id = roi_id, new_name = new_name)
  )
  object
}

#' Subset msi_dataset by ROI membership
#'
#' @param object A \code{msi_dataset} object.
#' @param roi_id Character ROI ids.
#' @param keep_roi Keep only selected ROI definitions in object.
#'
#' @return Subsetted \code{msi_dataset} object.
#' @export
subset_by_roi <- function(object, roi_id, keep_roi = TRUE) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }
  membership <- .get_roi_membership(object, roi_id = roi_id)
  if (nrow(membership) == 0) {
    stop("No ROI membership found for roi_id.")
  }
  keep_samples <- unique(membership$sample_id)
  idx <- match(keep_samples, object@sample_info$sample_id)
  idx <- idx[!is.na(idx)]

  # Avoid transient invalid state during subsetting: ROI memberships can
  # temporarily point to dropped sample_id values.
  original_roi_set <- .get_roi_set(object, create = TRUE)
  object <- .set_roi_set(object, .empty_roi_set())
  object <- object[, idx, drop = FALSE]
  object <- update_msi_dataset(object)

  if (keep_roi) {
    keep_ids <- unique(object@sample_info$sample_id)
    selected <- original_roi_set@rois[names(original_roi_set@rois) %in% roi_id]
    if (length(selected) > 0) {
      for (i in seq_along(selected)) {
        selected[[i]]@sample_id <- intersect(selected[[i]]@sample_id, keep_ids)
      }
      selected <- selected[vapply(selected, function(r) length(r@sample_id) > 0, logical(1))]
    }
    original_roi_set@rois <- selected
    object <- .set_roi_set(object, original_roi_set)
  }

  object <- add_process_log(
    object = object,
    function_name = "subset_by_roi()",
    parameter = list(roi_id = roi_id, keep_roi = keep_roi)
  )
  object
}

#' Aggregate intensity by ROI
#'
#' @param object A \code{msi_dataset} object.
#' @param feature Feature specification (\code{NULL} = all).
#' @param feature_type One of \code{"variable_id"} or \code{"mz"}.
#' @param fun Aggregation function: \code{"mean"}, \code{"median"}, \code{"sum"}.
#' @param tolerance_ppm ppm tolerance when \code{feature_type = "mz"}.
#'
#' @return Long-format data.frame with ROI-level summary.
#' @export
aggregate_by_roi <- function(object,
                             feature = NULL,
                             feature_type = c("variable_id", "mz"),
                             fun = c("mean", "median", "sum"),
                             tolerance_ppm = 20) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }
  feature_type <- match.arg(feature_type)
  fun <- match.arg(fun)

  membership <- .get_roi_membership(object)
  if (nrow(membership) == 0) {
    stop("No ROI found. Add ROI first.")
  }

  if (is.null(feature)) {
    feature_id <- object@variable_info$variable_id
    idx <- seq_len(nrow(object@variable_info))
  } else if (feature_type == "variable_id") {
    feature_id <- as.character(feature)
    idx <- match(feature_id, object@variable_info$variable_id)
  } else {
    mz_vals <- as.numeric(feature)
    idx <- vapply(mz_vals, function(x) .get_mz_index(object, x, tolerance_ppm), integer(1))
    feature_id <- object@variable_info$variable_id[idx]
  }

  if (any(is.na(idx))) {
    stop("Some requested feature cannot be matched.")
  }

  expr <- as.matrix(object@expression_data[idx, , drop = FALSE])
  rownames(expr) <- feature_id

  agg_fun <- switch(
    fun,
    mean = function(x) mean(x, na.rm = TRUE),
    median = function(x) median(x, na.rm = TRUE),
    sum = function(x) sum(x, na.rm = TRUE)
  )

  split_member <- split(membership$sample_id, membership$roi_id)
  out <- lapply(names(split_member), function(id) {
    cols <- match(unique(split_member[[id]]), colnames(expr))
    cols <- cols[!is.na(cols)]
    vals <- apply(expr[, cols, drop = FALSE], 1, agg_fun)
    data.frame(
      roi_id = id,
      feature = names(vals),
      value = as.numeric(vals),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

#' Create and add ROI from spatial polygon or rectangle
#'
#' @param object A \code{msi_dataset} object.
#' @param roi_id ROI identifier.
#' @param roi_name ROI display name.
#' @param type ROI shape type: \code{"rectangle"} or \code{"polygon"}.
#' @param roi Coordinate list. Rectangle needs \code{xmin/xmax/ymin/ymax};
#'   polygon needs \code{x/y}.
#' @param section_id Optional section linkage.
#' @param metadata Optional named list metadata.
#' @param overwrite Overwrite existing ROI with same id.
#'
#' @return Updated \code{msi_dataset} object.
#' @export
add_roi_from_spatial <- function(object,
                                 roi_id,
                                 roi_name = roi_id,
                                 type = c("rectangle", "polygon"),
                                 roi,
                                 section_id = NA_character_,
                                 metadata = list(),
                                 overwrite = FALSE) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }
  type <- match.arg(type)
  coords <- object@sample_info[, c("sample_id", "x", "y"), drop = FALSE]
  inside <- rep(FALSE, nrow(coords))

  if (type == "rectangle") {
    if (!all(c("xmin", "xmax", "ymin", "ymax") %in% names(roi))) {
      stop("rectangle ROI requires xmin, xmax, ymin, ymax.")
    }
    inside <- coords$x >= roi$xmin & coords$x <= roi$xmax &
      coords$y >= roi$ymin & coords$y <= roi$ymax
  } else {
    if (!all(c("x", "y") %in% names(roi))) {
      stop("polygon ROI requires x and y coordinates.")
    }
    inside <- .point_in_polygon(coords$x, coords$y, roi$x, roi$y)
  }

  sample_id <- coords$sample_id[inside]
  if (length(sample_id) == 0) {
    stop("No pixels selected for ROI.")
  }

  add_roi(
    object = object,
    roi_id = roi_id,
    roi_name = roi_name,
    sample_id = sample_id,
    section_id = section_id,
    metadata = metadata,
    overwrite = overwrite
  )
}
