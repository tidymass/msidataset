cosine_similarity <- function(x, y) {
  if (length(x) != length(y)) {
    stop("x and y must have the same length.")
  }
  x <- as.numeric(x)
  y <- as.numeric(y)
  num <- sum(x * y, na.rm = TRUE)
  den <- sqrt(sum(x * x, na.rm = TRUE)) * sqrt(sum(y * y, na.rm = TRUE))
  if (!is.finite(den) || den == 0) {
    return(0)
  }
  num / den
}

#' Benchmark pseudo-MS2 inference
#'
#' @param object A \code{msi_dataset} object.
#' @param experimental_ms2 Optional data frame with \code{precursor_mz}, \code{fragment_mz},
#'   and \code{intensity}. If omitted, the internal MS2 table is used.
#' @param tolerance Fragment matching tolerance (ppm).
#' @param ms2_library Optional subset of experimental MS2 entries (used when
#'   \code{experimental_ms2} is NULL).
#' @param ... Additional arguments passed to the internal benchmark helper.
#' @return Updated \code{msi_dataset} with benchmark results stored in
#'   \code{ms2_data$pseudo_ms2_benchmark} and
#'   \code{other_files$pseudo_ms2_benchmark}.
#' @export
benchmark_pseudo_ms2 <- function(object,
                                 experimental_ms2 = NULL,
                                 ppm = 10,
                                 tolerance = 0.01,
                                 ms2_library = NULL,
                                 ...) {
  if (!methods::is(object, "msi_dataset")) {
    stop("object must be a msi_dataset.")
  }

  bench_tbl <- NULL
  if (!is.null(experimental_ms2)) {
    bench_df <- as.data.frame(experimental_ms2, stringsAsFactors = FALSE, check.names = FALSE)
    required_cols <- c("precursor_mz", "fragment_mz", "intensity")
    if (!all(required_cols %in% colnames(bench_df))) {
      stop("experimental_ms2 must contain precursor_mz, fragment_mz, and intensity columns.")
    }
    pseudo_tbl <- object@ms2_data$pseudo_ms2_table
    if (is.null(pseudo_tbl) || nrow(as.data.frame(pseudo_tbl)) == 0) {
      stop("No pseudo-MS2 table found. Run generate_pseudo_ms2() first.")
    }
    pseudo_tbl <- as.data.frame(pseudo_tbl, stringsAsFactors = FALSE, check.names = FALSE)

    precursors <- sort(unique(bench_df$precursor_mz))
    bench_rows <- lapply(precursors, function(pmz) {
      exp_sub <- bench_df[mz_match_ppm(bench_df$precursor_mz, pmz, tolerance_ppm = ppm), , drop = FALSE]
      if (nrow(exp_sub) == 0) {
        return(NULL)
      }
      pseudo_sub <- pseudo_tbl[mz_match_ppm(pseudo_tbl$precursor_mz, pmz, tolerance_ppm = ppm), , drop = FALSE]
      total_exp <- nrow(exp_sub)
      total_pseudo <- nrow(pseudo_sub)
      cosine_sim <- 0
      matched_exp <- 0
      matched_pseudo <- 0
      if (total_pseudo > 0 && total_exp > 0) {
        matched_exp_idx <- logical(total_exp)
        matched_pseudo_idx <- logical(total_pseudo)
        for (i in seq_len(total_pseudo)) {
          ppm_diff <- abs(pseudo_sub$fragment_mz[[i]] - exp_sub$fragment_mz) /
            exp_sub$fragment_mz * 1e6
          matches <- which(ppm_diff <= ppm)
          if (length(matches) > 0) {
            matched_pseudo_idx[i] <- TRUE
            matched_exp_idx[matches] <- TRUE
          }
        }
        matched_exp <- sum(matched_exp_idx)
        matched_pseudo <- sum(matched_pseudo_idx)
        all_mz <- sort(unique(c(pseudo_sub$fragment_mz, exp_sub$fragment_mz)))
        if (length(all_mz) > 0) {
          pseudo_vec <- numeric(length(all_mz))
          exp_vec <- numeric(length(all_mz))
          pseudo_idx <- match(pseudo_sub$fragment_mz, all_mz)
          pseudo_vec[pseudo_idx] <- pseudo_sub$pseudo_intensity
          exp_idx <- match(exp_sub$fragment_mz, all_mz)
          exp_vec[exp_idx] <- exp_sub$intensity
          cosine_sim <- cosine_similarity(pseudo_vec, exp_vec)
        }
      }
      recall <- if (total_exp > 0) matched_exp / total_exp else NA_real_
      precision <- if (total_pseudo > 0) matched_pseudo / total_pseudo else NA_real_
      data.frame(
        precursor_mz = pmz,
        cosine_similarity = as.numeric(cosine_sim),
        fragment_recall = recall,
        fragment_precision = precision,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    })
    valid_rows <- bench_rows[!vapply(bench_rows, is.null, logical(1))]
    bench_tbl <- if (length(valid_rows) == 0) {
      NULL
    } else {
      do.call(rbind, valid_rows)
    }
  } else {
    object <- .benchmark_pseudo_ms2_actual(
      object = object,
      ms2_library = ms2_library,
      ppm = ppm,
      tolerance = tolerance,
      ...
    )
    bench_tbl <- object@other_files$pseudo_ms2_benchmark
  }
  if (is.null(bench_tbl)) {
    bench_tbl <- data.frame(
      precursor_mz = numeric(),
      cosine_similarity = numeric(),
      fragment_recall = numeric(),
      fragment_precision = numeric(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  object@ms2_data$pseudo_ms2_benchmark <- bench_tbl
  object@other_files$pseudo_ms2_benchmark <- bench_tbl
  object <- add_process_log(
    object = object,
    function_name = "benchmark_pseudo_ms2()",
    parameter = list(
      ppm = ppm,
      experimental = !is.null(experimental_ms2),
      n_precursors = nrow(bench_tbl)
    )
  )
  object
}

#' Plot pseudo-MS2 benchmark
#'
#' @param object A \code{msi_dataset} object with benchmark results.
#' @return A \code{ggplot2} object.
#' @export
plot_pseudo_ms2_benchmark <- function(object) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2 for plot_pseudo_ms2_benchmark().")
  }
  bench_tbl <- object@ms2_data$pseudo_ms2_benchmark
  if (is.null(bench_tbl) || nrow(as.data.frame(bench_tbl)) == 0) {
    bench_tbl <- object@other_files$pseudo_ms2_benchmark
  }
  if (is.null(bench_tbl) || nrow(as.data.frame(bench_tbl)) == 0) {
    stop("No pseudo-MS2 benchmark data found. Run benchmark_pseudo_ms2() first.")
  }
  bench_df <- as.data.frame(bench_tbl, stringsAsFactors = FALSE, check.names = FALSE)
  ggplot2::ggplot(bench_df, ggplot2::aes(x = precursor_mz, y = cosine_similarity)) +
    ggplot2::geom_point(color = "#1f78b4") +
    ggplot2::labs(
      title = "Pseudo-MS2 Benchmark",
      x = "Precursor m/z",
      y = "Cosine similarity"
    ) +
    ggplot2::theme_minimal()
}
