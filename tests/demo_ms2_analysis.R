# End-to-end MSI-MS2 analysis demo for msidataset

suppressPackageStartupMessages({
  library(msidataset)
})

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Please install ggplot2 to run this demo.")
}

demo_dir <- file.path("data", "demo data")
imzml_file <- file.path(demo_dir, "03-total ion count.imzML")
mgf_files <- c(
  file.path(demo_dir, "CHCA-neg-ms2-10ev_MSMS191.0200.mgf"),
  file.path(demo_dir, "NEDC-neg-ms2-10ev_MSMS130.0600.mgf"),
  file.path(demo_dir, "CHCA-neg-ms2-10ev_MSMS306.0700.mgf")
)

out_dir <- file.path("results", "demo_ms2_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1) read_imzml
obj <- read_imzml(
  path_1 = imzml_file,
  path_2 = out_dir,
  peak_align_tolerance = 10,
  verbose = TRUE
)

# 2) mutate_ms2
obj <- mutate_ms2(
  object = obj,
  path = mgf_files,
  ms1.ms2.match.mz.tol = 10,
  precursor_tolerance = 0,
  remove_precursor_peak = TRUE
)

# 3) extract_ms2
a_ms2 <- extract_ms2(obj, mz = 191.02, tolerance_ppm = 10)
print(a_ms2$metadata)

# 4) infer_pseudo_ms2
obj <- infer_pseudo_ms2(
  object = obj,
  mz = 191.02,
  tolerance_ppm = 10,
  correlation_threshold = 0.8,
  correlation_method = "spearman",
  spatial_smoothing = TRUE,
  smoothing_k = 1,
  fragment_intensity_threshold = 0.1,
  replace_existing_pseudo = TRUE
)

# 5) plot_ms2
p_ms2 <- plot_ms2(obj, mz = 191.02, ms2_type = "pseudo", tolerance_ppm = 10)
ggplot2::ggsave(
  file.path(out_dir, "plot_ms2_pseudo_191.02.png"),
  p_ms2,
  width = 7,
  height = 4,
  dpi = 150
)

# 6) plot_spatial_ms2
p_spatial_ms2 <- plot_spatial_ms2(
  object = obj,
  mz = 191.02,
  ms2_type = "pseudo",
  section_id = NULL,
  roi_id = NULL,
  overlay_roi = FALSE,
  tolerance_ppm = 20,
  max_fragments = 12,
  scales = "fixed",
  section_scale = "section"
)
ggplot2::ggsave(
  file.path(out_dir, "plot_spatial_ms2_pseudo_191.02.png"),
  p_spatial_ms2,
  width = 10,
  height = 6,
  dpi = 150
)

# 7) build_ms2_network
g <- build_ms2_network(
  object = obj,
  mz = 191.02,
  ms2_type = "pseudo",
  edge_weight = "spatial_correlation",
  correlation_method = "spearman",
  tolerance_ppm = 20
)
print(g)

# 8) plot_ms2_network
if (requireNamespace("ggraph", quietly = TRUE)) {
  p_net <- plot_ms2_network(g)
  ggplot2::ggsave(
    file.path(out_dir, "plot_ms2_network_pseudo_191.02.png"),
    p_net,
    width = 7,
    height = 6,
    dpi = 150
  )
}

saveRDS(obj, file.path(out_dir, "demo_ms2_analysis_final.rds"))
message("MS2 analysis demo complete: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
