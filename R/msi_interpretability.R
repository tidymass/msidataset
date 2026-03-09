# ============================================================================
# Flow visualization and object comparison utilities
# ============================================================================

.msi_flow_stages <- function() {
  c(
    "input validation",
    "data alignment",
    "core computation",
    "slot updates",
    "process log update",
    "output object"
  )
}

#' Show standardized core-verb flow for msidataset
#'
#' @param x Optional msi_dataset object for context.
#' @param verb One of convert, filter, mutate, roi, slice, summarize.
#' @return A flow object (DiagrammeR graph if available, otherwise text-flow object).
#' @export
msi_show_flow <- function(x = NULL,
                          verb = c("convert", "filter", "mutate", "roi", "slice", "summarize")) {
  if (!is.null(x) && !methods::is(x, "msi_dataset")) {
    stop("x must be NULL or a msi_dataset object.")
  }
  verb <- match.arg(verb)
  stages <- .msi_flow_stages()
  nodes <- data.frame(
    id = seq_along(stages),
    label = paste(seq_along(stages), stages, sep = ". "),
    stringsAsFactors = FALSE
  )
  edges <- data.frame(
    from = seq_len(length(stages) - 1L),
    to = seq(2L, length(stages)),
    stringsAsFactors = FALSE
  )

  if (requireNamespace("DiagrammeR", quietly = TRUE)) {
    g <- DiagrammeR::create_graph() |>
      DiagrammeR::add_global_graph_attrs(
        attr = "layout",
        value = "dot",
        attr_type = "graph"
      ) |>
      DiagrammeR::add_global_graph_attrs(
        attr = "rankdir",
        value = "LR",
        attr_type = "graph"
      )

    for (i in seq_len(nrow(nodes))) {
      g <- DiagrammeR::add_node(
        graph = g,
        label = nodes$label[[i]],
        type = verb
      )
    }
    for (i in seq_len(nrow(edges))) {
      g <- DiagrammeR::add_edge(g, from = edges$from[[i]], to = edges$to[[i]])
    }
    out <- list(verb = verb, nodes = nodes, edges = edges, graph = g)
    class(out) <- c("msi_flow", "msi_flow_diagrammer")
    return(out)
  }

  out <- list(
    verb = verb,
    nodes = nodes,
    edges = edges,
    text = paste0(
      "Flow(", verb, "): ",
      paste(stages, collapse = " -> ")
    )
  )
  class(out) <- c("msi_flow", "msi_flow_text")
  out
}

#' @export
print.msi_flow <- function(x, ...) {
  if (inherits(x, "msi_flow_diagrammer")) {
    print(x$graph)
    return(invisible(x))
  }
  cat(x$text, "\n")
  invisible(x)
}

#' @keywords internal
autoplot.msi_flow <- function(object, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2 for autoplot.msi_flow().")
  }
  nodes <- object$nodes
  nodes$y <- rev(seq_len(nrow(nodes)))
  nodes$x <- 1
  ggplot2::ggplot(nodes, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(size = 2, color = "#2C7FB8") +
    ggplot2::geom_segment(
      data = nodes[-nrow(nodes), , drop = FALSE],
      ggplot2::aes(x = x, y = y, xend = x, yend = y - 1),
      arrow = ggplot2::arrow(length = grid::unit(0.15, "inches"))
    ) +
    ggplot2::geom_text(ggplot2::aes(label = label), hjust = -0.05, size = 3.5) +
    ggplot2::coord_cartesian(xlim = c(0.95, 2), clip = "off") +
    ggplot2::theme_void() +
    ggplot2::ggtitle(paste("msidataset flow:", object$verb))
}

.ms2_attachment_count <- function(object) {
  info <- try(extract_ms2(object, what = "ms2_info"), silent = TRUE)
  if (inherits(info, "try-error") || is.null(info)) {
    return(0L)
  }
  as.integer(nrow(info))
}

.roi_membership_count <- function(object) {
  membership <- try(roi_data(object), silent = TRUE)
  if (inherits(membership, "try-error") || is.null(membership)) {
    return(0L)
  }
  as.integer(nrow(membership))
}

#' Compare two msi_dataset object states
#'
#' @param before msi_dataset object before operation.
#' @param after msi_dataset object after operation.
#' @return A tibble with metric-level before/after/delta.
#' @export
compare_msi <- function(before, after) {
  if (!methods::is(before, "msi_dataset")) {
    stop("before must be a msi_dataset object.")
  }
  if (!methods::is(after, "msi_dataset")) {
    stop("after must be a msi_dataset object.")
  }

  b_mz <- if (nrow(before@variable_info) > 0) range(before@variable_info$mz, na.rm = TRUE) else c(NA_real_, NA_real_)
  a_mz <- if (nrow(after@variable_info) > 0) range(after@variable_info$mz, na.rm = TRUE) else c(NA_real_, NA_real_)

  res <- tibble::tibble(
    metric = c(
      "feature_count",
      "pixel_count",
      "mz_min",
      "mz_max",
      "roi_membership_count",
      "ms2_attachment_count",
      "ms2_feature_count",
      "process_log_steps"
    ),
    before = c(
      nrow(before@variable_info),
      nrow(before@sample_info),
      b_mz[[1]],
      b_mz[[2]],
      .roi_membership_count(before),
      .ms2_attachment_count(before),
      sum(
        if ("has_ms2" %in% colnames(before@variable_info)) {
          before@variable_info$has_ms2 %in% c("Yes", "yes", TRUE, 1)
        } else {
          before@variable_info$ms2 %in% c("Yes", "yes", TRUE, 1)
        },
        na.rm = TRUE
      ),
      length(before@process_info)
    ),
    after = c(
      nrow(after@variable_info),
      nrow(after@sample_info),
      a_mz[[1]],
      a_mz[[2]],
      .roi_membership_count(after),
      .ms2_attachment_count(after),
      sum(
        if ("has_ms2" %in% colnames(after@variable_info)) {
          after@variable_info$has_ms2 %in% c("Yes", "yes", TRUE, 1)
        } else {
          after@variable_info$ms2 %in% c("Yes", "yes", TRUE, 1)
        },
        na.rm = TRUE
      ),
      length(after@process_info)
    )
  )
  res$delta <- suppressWarnings(as.numeric(res$after) - as.numeric(res$before))

  roi_before <- try(unique(roi_data(before)$roi_id), silent = TRUE)
  roi_after <- try(unique(roi_data(after)$roi_id), silent = TRUE)
  roi_before <- if (inherits(roi_before, "try-error")) character() else as.character(roi_before)
  roi_after <- if (inherits(roi_after, "try-error")) character() else as.character(roi_after)

  attr(res, "details") <- list(
    roi_added = setdiff(roi_after, roi_before),
    roi_removed = setdiff(roi_before, roi_after),
    process_log_delta = length(after@process_info) - length(before@process_info)
  )
  class(res) <- c("compare_msi", class(res))
  res
}

#' @export
print.compare_msi <- function(x, ...) {
  print(as.data.frame(x), row.names = FALSE)
  d <- attr(x, "details")
  if (!is.null(d)) {
    cat("ROI added:", if (length(d$roi_added) == 0) "none" else paste(d$roi_added, collapse = ", "), "\n")
    cat("ROI removed:", if (length(d$roi_removed) == 0) "none" else paste(d$roi_removed, collapse = ", "), "\n")
    cat("Process log delta:", d$process_log_delta, "\n")
  }
  invisible(x)
}
