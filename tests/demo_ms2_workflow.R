# Professional MSI-MS2 demonstration workflow for msidataset
#
# This script is designed to run end-to-end with:
#   source("demo_ms2_workflow.R")
#
# It is compatible with development usage after:
#   devtools::load_all()

suppressPackageStartupMessages({
  if (!"package:msidataset" %in% search()) {
    library(msidataset)
  }
})

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Please install ggplot2 to run this workflow.")
}

if (requireNamespace("progressr", quietly = TRUE)) {
  progressr::handlers(global = TRUE)
  progressr::handlers("txtprogressbar")
}

.step_header <- function(i, title) {
  cat("\n", paste0(rep("=", 72), collapse = ""), "\n", sep = "")
  cat(sprintf("Step %d. %s\n", i, title))
  cat(paste0(rep("=", 72), collapse = ""), "\n", sep = "")
}

.resolve_demo_dir <- function() {
  # Primary package-internal lookup.
  p <- system.file("data/demo data", package = "msidataset")
  if (nzchar(p) && dir.exists(p)) {
    return(p)
  }

  # Fallback for devtools::load_all() environments where package root is exposed.
  pkg_root <- system.file(package = "msidataset")
  p2 <- file.path(pkg_root, "data", "demo data")
  if (nzchar(pkg_root) && dir.exists(p2)) {
    return(p2)
  }

  stop("Could not locate demo data directory via package paths.")
}

# -----------------------------------------------------------------------------
# Step 1 — Load demo MSI data
# -----------------------------------------------------------------------------
.step_header(1, "Load Demo MSI Data")
demo_dir <- .resolve_demo_dir()
imzml_file <- file.path(demo_dir, "03-total ion count.imzML")
mgf_files <- c(
  file.path(demo_dir, "CHCA-neg-ms2-10ev_MSMS191.0200.mgf"),
  file.path(demo_dir, "NEDC-neg-ms2-10ev_MSMS130.0600.mgf"),
  file.path(demo_dir, "CHCA-neg-ms2-10ev_MSMS306.0700.mgf")
)

cat("Demo directory:", demo_dir, "\n")
cat("imzML file:", basename(imzml_file), "\n")
cat("MGF files:", paste(basename(mgf_files), collapse = ", "), "\n")

if (!file.exists(imzml_file)) {
  stop("imzML file not found: ", imzml_file)
}
if (!all(file.exists(mgf_files))) {
  stop("One or more MGF files were not found in demo directory.")
}

out_dir <- file.path("results", "demo_ms2_workflow")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("Output directory:", normalizePath(out_dir, winslash = "/", mustWork = FALSE), "\n")

# -----------------------------------------------------------------------------
# Auto strategy — detect large files and downgrade read_imzml settings
# -----------------------------------------------------------------------------
imzml_info <- file.info(imzml_file)
ibd_file <- sub("\\.imzML$", ".ibd", imzml_file, ignore.case = TRUE)
ibd_size <- if (file.exists(ibd_file)) file.info(ibd_file)$size else 0
total_size_mb <- (imzml_info$size + ibd_size) / (1024^2)

size_threshold_mb <- 300
is_large_input <- is.finite(total_size_mb) && total_size_mb >= size_threshold_mb

read_params <- if (is_large_input) {
  list(
    peak_align_tolerance = 10,
    bin_digits = 1,
    max_features = 20000,
    chunk_size = 128,
    parallel_backend = "none",
    workers = 1,
    use_cache = TRUE,
    verbose = TRUE
  )
} else {
  list(
    peak_align_tolerance = 10,
    bin_digits = 4,
    max_features = 150000,
    chunk_size = 256,
    parallel_backend = "none",
    workers = 1,
    use_cache = TRUE,
    verbose = TRUE
  )
}

cat("Input size (imzML + ibd):", sprintf("%.1f MB", total_size_mb), "\n")
cat("Auto strategy:", if (is_large_input) "lightweight (large input detected)" else "standard", "\n")
cat(
  "read_imzml params:",
  paste(
    paste(names(read_params), unlist(read_params), sep = "="),
    collapse = ", "
  ),
  "\n"
)

# -----------------------------------------------------------------------------
# Step 2 — Convert to msidataset object
# -----------------------------------------------------------------------------
.step_header(2, "Convert to msidataset Object")
msi_obj <- do.call(
  read_imzml,
  c(
    list(
      path_1 = imzml_file,
      path_2 = out_dir
    ),
    read_params
  )
)

if (!methods::is(msi_obj, "msi_dataset")) {
  stop("read_imzml() did not return an msi_dataset object.")
}
cat("Object class:", class(msi_obj)[1], "\n")
cat("Dimensions (features x pixels):", nrow(expression_data(msi_obj)), "x", ncol(expression_data(msi_obj)), "\n")

# -----------------------------------------------------------------------------
# Step 3 — Import MS2 spectra from MGF files
# -----------------------------------------------------------------------------
.step_header(3, "Import MS2 Spectra from MGF Files")
cat("Importing MS2 spectra from", length(mgf_files), "MGF files.\n")

# -----------------------------------------------------------------------------
# Step 4 — Perform MS2 preprocessing
# -----------------------------------------------------------------------------
.step_header(4, "Perform MS2 Preprocessing")
cat("Preprocessing rules:\n")
cat("  - remove zero intensity peaks\n")
cat("  - normalize intensities to max = 1\n")
cat("  - remove precursor peak\n")
cat("  - remove fragments above precursor m/z\n")

msi_obj <- mutate_ms2(
  object = msi_obj,
  path = mgf_files,
  ms1.ms2.match.mz.tol = 10,
  precursor_tolerance = 0,
  remove_precursor_peak = TRUE
)

ms2_tbl <- msi_obj@ms2_data$ms2_table
cat("Integrated MS2 spectra:", length(unique(ms2_tbl$ms2_id)), "\n")
cat("Integrated fragment rows:", nrow(ms2_tbl), "\n")

# -----------------------------------------------------------------------------
# Step 5 — Perform MS1–MS2 matching using ppm tolerance
# -----------------------------------------------------------------------------
.step_header(5, "Perform MS1-MS2 Matching Using PPM Tolerance")
cat("Matching tolerance (ppm): 10\n")
cat("Matched precursor variables:", length(unique(ms2_tbl$precursor_variable_id)), "\n")

# -----------------------------------------------------------------------------
# Step 6 — Inspect variable_info MS2 annotations
# -----------------------------------------------------------------------------
.step_header(6, "Inspect variable_info MS2 Annotations")
vinfo <- variable_info(msi_obj)
ms2_annot <- vinfo[vinfo$has_ms2 == "yes", c("variable_id", "mz", "has_ms2", "ms2_id", "ms2_type"), drop = FALSE]
cat("Variables with MS2 annotation:", nrow(ms2_annot), "\n")
print(utils::head(ms2_annot, 10))

# Choose a target precursor near 191.02 m/z for the tutorial.
target_var <- vinfo$variable_id[which.min(abs(vinfo$mz - 191.02))]
target_mz <- vinfo$mz[match(target_var, vinfo$variable_id)]
cat("Selected precursor for demonstration:", target_var, "(m/z", signif(target_mz, 6), ")\n")

# -----------------------------------------------------------------------------
# Step 7 — Extract MS2 spectra
# -----------------------------------------------------------------------------
.step_header(7, "Extract MS2 Spectra")
ms2_actual <- extract_ms2(msi_obj, variable_id = target_var, ms2_type = "actual")
if (is.null(ms2_actual)) {
  stop("No actual MS2 spectrum found for selected precursor.")
}
cat("Actual MS2 ID:", ms2_actual$metadata$ms2_id, "\n")
cat("Actual fragment count:", length(ms2_actual$fragments), "\n")

# -----------------------------------------------------------------------------
# Step 8 — Build pseudo MS2 using spatial correlation
# -----------------------------------------------------------------------------
.step_header(8, "Build Pseudo MS2 Using Spatial Correlation")
threshold <- 0.8
cat("Pseudo-MS2 settings: method = spearman, threshold =", threshold, "\n")

msi_obj <- infer_pseudo_ms2(
  object = msi_obj,
  variable_id = target_var,
  correlation_threshold = threshold,
  correlation_method = "spearman",
  spatial_smoothing = TRUE,
  smoothing_k = 1,
  smoothing_method = "gaussian",
  precursor_exclusion_tol = 0.01,
  fragment_window = 150,
  adaptive_threshold = TRUE,
  fallback_top_n = 5,
  fragment_intensity_threshold = 0.1,
  replace_existing_pseudo = TRUE
)

ms2_pseudo <- extract_ms2(msi_obj, variable_id = target_var, ms2_type = "pseudo")
if (is.null(ms2_pseudo)) {
  ms2_pseudo <- extract_ms2(msi_obj, variable_id = target_var, ms2_type = "pseudo_candidate")
}
if (is.null(ms2_pseudo)) {
  cat("No pseudo fragments passed threshold; rerunning with threshold = 0.7 for demonstration.\n")
  msi_obj <- infer_pseudo_ms2(
    object = msi_obj,
    variable_id = target_var,
    correlation_threshold = 0.7,
    correlation_method = "spearman",
    spatial_smoothing = TRUE,
    smoothing_k = 1,
    smoothing_method = "gaussian",
    precursor_exclusion_tol = 0.01,
    fragment_window = 150,
    adaptive_threshold = TRUE,
    fallback_top_n = 5,
    fragment_intensity_threshold = 0.05,
    replace_existing_pseudo = TRUE
  )
  ms2_pseudo <- extract_ms2(msi_obj, variable_id = target_var, ms2_type = "pseudo")
  if (is.null(ms2_pseudo)) {
    ms2_pseudo <- extract_ms2(msi_obj, variable_id = target_var, ms2_type = "pseudo_candidate")
  }
}

if (!is.null(ms2_pseudo)) {
  pseudo_tbl <- msi_obj@ms2_data$ms2_table
  pseudo_tbl <- pseudo_tbl[
    pseudo_tbl$precursor_variable_id == target_var & pseudo_tbl$ms2_type %in% c("pseudo", "pseudo_candidate"),
    c("fragment_mz", "fragment_intensity", "confidence_score"),
    drop = FALSE
  ]
  pseudo_tbl <- pseudo_tbl[order(pseudo_tbl$fragment_intensity, decreasing = TRUE), , drop = FALSE]
  cat("Pseudo fragment count:", nrow(pseudo_tbl), "\n")
  print(utils::head(pseudo_tbl, 10))
} else {
  cat("Pseudo fragment count: 0\n")
}

# -----------------------------------------------------------------------------
# Step 9 — Visualize MS2 spectra
# -----------------------------------------------------------------------------
.step_header(9, "Visualize MS2 Spectra")
p_ms2_actual <- plot_ms2(msi_obj, variable_id = target_var, ms2_type = "actual")

if (!is.null(ms2_pseudo)) {
  p_ms2_pseudo <- plot_ms2(msi_obj, variable_id = target_var, ms2_type = ms2_pseudo$ms2_type)
} else {
  p_ms2_pseudo <- NULL
}

# Select one fragment for spatial comparison (highest pseudo intensity if available).
if (!is.null(ms2_pseudo) && length(ms2_pseudo$fragments) > 0) {
  frag_mz_selected <- ms2_pseudo$fragments[which.max(ms2_pseudo$intensity)]
} else {
  frag_mz_selected <- ms2_actual$fragments[which.max(ms2_actual$intensity)]
}
cat("Selected fragment for spatial visualization: m/z", signif(frag_mz_selected, 6), "\n")

p_spatial_precursor <- plot_spatial(
  object = msi_obj,
  feature = target_var,
  feature_type = "variable_id",
  facet_by = c("feature"),
  section_scale = "global",
  scales = "fixed"
) + ggplot2::labs(title = paste0("Precursor Spatial Map: ", target_var))

p_spatial_fragment <- plot_spatial(
  object = msi_obj,
  feature = frag_mz_selected,
  feature_type = "mz",
  tolerance_ppm = 20,
  facet_by = c("feature"),
  section_scale = "global",
  scales = "fixed"
) + ggplot2::labs(title = paste0("Fragment Spatial Map: m/z ", signif(frag_mz_selected, 6)))

# -----------------------------------------------------------------------------
# Step 10 — Visualize spatial correlation between precursor and fragments
# -----------------------------------------------------------------------------
.step_header(10, "Visualize Spatial Correlation Between Precursor and Fragments")
precursor_exclusion_tol <- 0.01
fragment_window <- 150

corr_df <- msi_obj@ms2_data$pseudo_correlation_table
if (is.null(corr_df) || nrow(corr_df) == 0) {
  corr_df <- data.frame(
    precursor_variable_id = character(),
    precursor_mz = numeric(),
    fragment_mz = numeric(),
    variable_id = character(),
    correlation = numeric(),
    normalized_fragment_intensity = numeric(),
    score = numeric(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
corr_df <- corr_df[corr_df$precursor_variable_id == target_var, , drop = FALSE]
corr_df <- corr_df[
  abs(corr_df$fragment_mz - target_mz) > precursor_exclusion_tol &
    corr_df$fragment_mz < target_mz &
    corr_df$fragment_mz > target_mz - fragment_window,
  ,
  drop = FALSE
]
corr_df <- corr_df[order(corr_df$score, decreasing = TRUE), , drop = FALSE]

if (nrow(corr_df) > 0) {
  frag_mz_top <- corr_df$fragment_mz[[1]]
} else {
  frag_mz_top <- frag_mz_selected
}

# Explicit Step-10 maps: precursor and top-ranked correlated fragment.
p_step10_precursor <- plot_spatial(
  object = msi_obj,
  feature = target_var,
  feature_type = "variable_id",
  facet_by = c("feature"),
  section_scale = "global",
  scales = "fixed"
) + ggplot2::labs(title = paste0("Step 10: Precursor Map (", target_var, ")"))

p_step10_fragment <- plot_spatial(
  object = msi_obj,
  feature = frag_mz_top,
  feature_type = "mz",
  tolerance_ppm = 20,
  facet_by = c("feature"),
  section_scale = "global",
  scales = "fixed"
) + ggplot2::labs(title = paste0("Step 10: Top Correlated Fragment Map (m/z ", signif(frag_mz_top, 6), ")"))

p_corr <- ggplot2::ggplot(corr_df, ggplot2::aes(x = reorder(variable_id, score), y = score)) +
  ggplot2::geom_col(fill = "#2c7fb8") +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("m/z %.4f | r=%.2f", fragment_mz, correlation)),
    hjust = -0.02,
    size = 3
  ) +
  ggplot2::coord_flip(clip = "off") +
  ggplot2::labs(
    title = paste0("Step 10: Fragment Correlation Ranking (", target_var, ")"),
    subtitle = "score = correlation * normalized fragment intensity",
    x = "fragment variable_id",
    y = "score"
  ) +
  ggplot2::theme_minimal(base_size = 12)

cat("Top correlated fragments:\n")
print(utils::head(corr_df, 10))

# -----------------------------------------------------------------------------
# Step 11 — Compare real MS2 and pseudo MS2
# -----------------------------------------------------------------------------
.step_header(11, "Compare Real MS2 and Pseudo MS2")
ms2_tbl_current <- msi_obj@ms2_data$ms2_table
real_id <- unique(ms2_tbl_current$ms2_id[ms2_tbl_current$precursor_variable_id == target_var & ms2_tbl_current$ms2_type == "actual"])[1]
pseudo_id <- unique(ms2_tbl_current$ms2_id[ms2_tbl_current$precursor_variable_id == target_var &
  ms2_tbl_current$ms2_type %in% c("pseudo", "pseudo_candidate")])[1]

if (!is.na(real_id) && !is.na(pseudo_id)) {
  cmp_res <- compare_spectra(
    object = msi_obj,
    ms2_id = real_id,
    pseudo_ms2_id = pseudo_id,
    tolerance = 0.01
  )
  cat("Real MS2 ID:", real_id, "\n")
  cat("Pseudo MS2 ID:", pseudo_id, "\n")
  cat("Matched fragment count:", cmp_res$matched_fragment_count, "\n")
  cat("Cosine similarity:", signif(cmp_res$cosine_similarity, 4), "\n")
  cat("Fragment overlap ratio:", signif(cmp_res$overlap_ratio, 4), "\n")
} else {
  cmp_res <- NULL
  cat("Spectrum comparison skipped (missing real or pseudo spectrum ID).\\n")
}

# -----------------------------------------------------------------------------
# Step 12 — Generate mirror plot
# -----------------------------------------------------------------------------
.step_header(12, "Generate Mirror Plot")
if (!is.null(cmp_res)) {
  p_mirror <- mirror_plot(msi_obj, ms2_id = real_id, pseudo_ms2_id = pseudo_id)
} else {
  p_mirror <- NULL
}

# -----------------------------------------------------------------------------
# Step 13 — Benchmark pseudo MS2 reconstruction
# -----------------------------------------------------------------------------
.step_header(13, "Benchmark Pseudo MS2 Reconstruction")
msi_obj <- benchmark_pseudo_ms2(
  object = msi_obj,
  tolerance = 0.01,
  correlation_method = "spearman",
  spatial_smoothing = TRUE,
  smoothing_k = 1,
  smoothing_method = "gaussian",
  precursor_exclusion_tol = 0.01,
  fragment_window = 150,
  adaptive_threshold = TRUE,
  fallback_top_n = 5,
  fragment_intensity_threshold = 0.05
)
benchmark_tbl <- msi_obj@other_files$pseudo_ms2_benchmark
cat("Benchmark rows:", nrow(benchmark_tbl), "\n")
print(utils::head(benchmark_tbl, 10))

# Required figure outputs for pseudo-MS2 validation module
p_plot_pseudo_ms2 <- plot_pseudo_ms2(msi_obj, variable_id = target_var, tolerance_ppm = 20)
spcorr_obj <- plot_spatial_correlation(msi_obj, n_vars = 20, method = "pearson", show_labels = FALSE)
p_plot_spatial_correlation <- spcorr_obj$plot_x
if (requireNamespace("ggraph", quietly = TRUE)) {
  p_plot_pseudo_network <- plot_pseudo_network(
    object = msi_obj,
    variable_id = target_var,
    edge_weight = "spatial_correlation",
    correlation_method = "spearman",
    tolerance_ppm = 20
  )
} else {
  p_plot_pseudo_network <- NULL
  cat("Skipping pseudo network plot: package 'ggraph' is not installed.\n")
}

# Save outputs
plot_list <- list(
  p_ms2_actual = p_ms2_actual,
  p_ms2_pseudo = p_ms2_pseudo,
  p_spatial_precursor = p_spatial_precursor,
  p_spatial_fragment = p_spatial_fragment,
  p_step10_precursor = p_step10_precursor,
  p_step10_fragment = p_step10_fragment,
  p_corr = p_corr,
  p_plot_pseudo_ms2 = p_plot_pseudo_ms2,
  p_plot_spatial_correlation = p_plot_spatial_correlation,
  p_mirror = p_mirror,
  p_plot_pseudo_network = p_plot_pseudo_network
)

if (!is.null(plot_list$p_ms2_actual)) {
  ggplot2::ggsave(file.path(out_dir, "plot_ms2_actual.png"), plot_list$p_ms2_actual, width = 7, height = 4, dpi = 150)
}
if (!is.null(plot_list$p_ms2_pseudo)) {
  ggplot2::ggsave(file.path(out_dir, "plot_ms2_pseudo.png"), plot_list$p_ms2_pseudo, width = 7, height = 4, dpi = 150)
}
ggplot2::ggsave(file.path(out_dir, "plot_spatial_precursor.png"), plot_list$p_spatial_precursor, width = 7, height = 5, dpi = 150)
ggplot2::ggsave(file.path(out_dir, "plot_spatial_fragment.png"), plot_list$p_spatial_fragment, width = 7, height = 5, dpi = 150)
ggplot2::ggsave(file.path(out_dir, "plot_step10_precursor.png"), plot_list$p_step10_precursor, width = 7, height = 5, dpi = 150)
ggplot2::ggsave(file.path(out_dir, "plot_step10_top_fragment.png"), plot_list$p_step10_fragment, width = 7, height = 5, dpi = 150)
ggplot2::ggsave(file.path(out_dir, "plot_fragment_spatial_correlation.png"), plot_list$p_corr, width = 8, height = 4.5, dpi = 150)
ggplot2::ggsave(file.path(out_dir, "plot_pseudo_ms2.png"), plot_list$p_plot_pseudo_ms2, width = 7, height = 4, dpi = 150)
ggplot2::ggsave(file.path(out_dir, "plot_spatial_correlation.png"), plot_list$p_plot_spatial_correlation, width = 8, height = 5, dpi = 150)
if (!is.null(plot_list$p_mirror)) {
  ggplot2::ggsave(file.path(out_dir, "plot_mirror_real_vs_pseudo.png"), plot_list$p_mirror, width = 8, height = 4, dpi = 150)
}
if (!is.null(plot_list$p_plot_pseudo_network)) {
  ggplot2::ggsave(file.path(out_dir, "plot_pseudo_network.png"), plot_list$p_plot_pseudo_network, width = 7, height = 6, dpi = 150)
}

saveRDS(msi_obj, file.path(out_dir, "demo_ms2_workflow_final.rds"))
saveRDS(plot_list, file.path(out_dir, "demo_ms2_workflow_plots.rds"))
saveRDS(benchmark_tbl, file.path(out_dir, "pseudo_ms2_benchmark.rds"))

cat("\nWorkflow complete. Results written to:\n")
cat(normalizePath(out_dir, winslash = "/", mustWork = FALSE), "\n")

invisible(list(
  object = msi_obj,
  plots = plot_list,
  correlation = corr_df,
  benchmark = benchmark_tbl,
  spectrum_comparison = cmp_res
))
