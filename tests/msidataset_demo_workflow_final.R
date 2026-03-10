# ============================================================================
# Complete msidataset demonstration workflow
# Data source: D:/Shen_lab/msidataset/data/demo data/
# ============================================================================

# 1. Load packages and set paths
# ============================================================================
library(msidataset)
library(methods)

cat("✓ Packages loaded successfully！\n\n")

# Data directory
demo_dir <- "D:/Shen_lab/msidataset/data/demo data"
output_dir <- "msidataset_results"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("Data directory: %s\n", demo_dir))
cat(sprintf("Output directory: %s\n\n", normalizePath(output_dir)))

render_step_progress <- function(step, label, percent) {
  percent <- max(0, min(100, as.integer(percent)))
  cli::cli_text(paste0("Step ", step, " ", label))
  pb <- cli::cli_progress_bar(
    total = 100,
    format = "{cli::pb_bar} {cli::pb_percent}",
    clear = TRUE
  )
  cli::cli_progress_update(pb, set = percent)
  cli::cli_progress_done(pb)
}

format_eta <- function(seconds) {
  seconds <- max(0, as.integer(seconds))
  if (seconds == 0) {
    return("0s")
  }
  minutes <- seconds %/% 60
  remainder <- seconds %% 60
  if (minutes > 0) {
    sprintf("%dm %02ds", minutes, remainder)
  } else {
    sprintf("%ds", remainder)
  }
}

report_alignment_eta <- function(done, total) {
  total <- max(1, as.integer(total))
  done <- min(max(0, as.integer(done)), total)
  remaining <- total - done
  eta_seconds <- remaining * 2
  eta_label <- format_eta(eta_seconds)
  cli::cli_alert_info(sprintf("Aligning spectra %d / %d (ETA: %s)", done, total, eta_label))
}

# ============================================================================
# 2. Read imzML files using read_imzml()
# ============================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Step 2: Reading imzML files...\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# 读取第一个imzML文件
cat("[1/2] Reading 03-total ion count.imzML...\n")
msi_obj_03 <- read_imzml(
  path_1 = file.path(demo_dir, "03-total ion count.imzML"),
  peak_align_tolerance = 50,
  path_2 = output_dir,
  resolution = 50,
  matrix_name = "CHCA"
)

if (methods::is(msi_obj_03, "msi_dataset")) {
  cat("✓ Successfully read!\n")
  cat(sprintf(" - Dimensions: %s\n", paste(dim(msi_obj_03), collapse = " × ")))
  cat(sprintf(" - Number of pixels: %d\n", n_pixels(msi_obj_03)))
  cat(sprintf(" - Number of m/z values: %d\n", n_variables(msi_obj_03)))
} else {
  cat("✗ Reading failed.\n")
  stop("Unable to read the first imzML file")
}

# Read the second imzML file
cat("\n[2/2] Reading 04-total ion count.imzML...\n")
msi_obj_04 <- read_imzml(
  path_1 = file.path(demo_dir, "04-total ion count.imzML"),
  peak_align_tolerance = 50,
  path_2 = output_dir,
  resolution = 50,
  matrix_name = "CHCA"
)

if (methods::is(msi_obj_04, "msi_dataset")) {
  cat("✓ Successfully read! \n")
  cat(sprintf("  - Dimensions: %s\n", paste(dim(msi_obj_04), collapse = " × ")))
  cat(sprintf("  - Number of pixels: %d\n", n_pixels(msi_obj_04)))
  cat(sprintf("  - Number of m/z values: %d\n", n_variables(msi_obj_04)))
} else {
  cat("✗ Reading failed\n")
  stop("Unable to read the second imzML file")
}


cat("\n")

# ============================================================================
# 3. Merge two datasets
# ============================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Step 3: Merging imzML data\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat("Merging 03 and 04 datasets...\n")
msi_obj <- merge_msi_datasets(msi_obj_03, msi_obj_04)

if (methods::is(msi_obj, "msi_dataset")) {
  cat("✓ Merge successful\n")
  cat(sprintf("  - Total dimensions:  %s\n", paste(dim(msi_obj), collapse = " × ")))
  cat(sprintf("  - Total pixels: %d\n", n_pixels(msi_obj)))
  cat(sprintf("  - Total m/z values: %d\n", n_variables(msi_obj)))
} else {
  stop("Merge failed")
}

cat("\n")

# ============================================================================
# 4. Basic data information
# ============================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Step 4: Basic data information\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat("Calling msidataset_logo():\n")
msidataset_logo()

cat("\nCalling msidataset_packages():\n")
print(msidataset_packages())

cat("\nCalling check_msi_dataset_class():\n")
print(check_msi_dataset_class(msi_obj))

cat("\nCalling show():\n")
show(msi_obj)

cat("\nBasic statistics:\n")
cat(sprintf("  - Dimensions: %s\n", paste(dim(msi_obj), collapse = " × ")))
cat(sprintf("  - Number of pixels: %d\n", n_pixels(msi_obj)))
cat(sprintf("  - Number of samples: %d\n", n_samples(msi_obj)))
cat(sprintf("  - Number of variables (m/z): %d\n", n_variables(msi_obj)))
cat(sprintf("  - m/z range: %.2f - %.2f\n", min(mz(msi_obj)), max(mz(msi_obj))))
cat(sprintf("  - Number of spatial coordinates: %d\n", nrow(get_spatial_coords(msi_obj))))

cat("\n")

# ============================================================================
# 5. Load MGF files (MS2 data) - mutate_ms2
# ============================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Step 5: Loading MGF files (MS2 data)\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

mgf_files <- c(
  file.path(demo_dir, "CHCA-neg-ms2-10ev_MSMS191.0200.mgf"),
  file.path(demo_dir, "NEDC-neg-ms2-10ev_MSMS130.0600.mgf"),
  file.path(demo_dir, "CHCA-neg-ms2-10ev_MSMS306.0700.mgf")
)

cat(sprintf("Checking MGF files:\n"))
for (i in seq_along(mgf_files)) {
  if (file.exists(mgf_files[i])) {
    cat(sprintf("  %d. ✓ %s\n", i, basename(mgf_files[i])))
  } else {
    cat(sprintf("  %d. ✗ %s (file does not exist)\n", i, basename(mgf_files[i])))
  }
}

cat("\nCalling mutate_ms2() to load MGF data...\n")
msi_obj <- mutate_ms2(
  object = msi_obj,
  path = mgf_files,
  ms1.ms2.match.mz.tol = 50
)

cat("✓ MGF data loaded successfully\n")
ms2_info <- extract_ms2(msi_obj, what = "ms2_info")
cat(sprintf("  - Number of MS2 spectra: %d\n", nrow(ms2_info)))

cat("\nFirst 5 MS2 entries:\n")
print(head(ms2_info, 5))

cat("\n")

# ============================================================================
# 6. Filter MS2 data
# ============================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Step 6: Filtering MS2 data\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat("Original MS2 data: \n")
print(utils::head(extract_ms2(msi_obj, what = "ms2_info"), 3))

cat("\nFiltering by precursor ion m/z range (100-400)...\n")
msi_obj <- filter_ms2(msi_obj, precursor_mz >= 100, precursor_mz <= 400)

cat("✓ Filtering complete\n")
ms2_info_filtered <- extract_ms2(msi_obj, what = "ms2_info")
report_alignment_eta(nrow(ms2_info_filtered), nrow(ms2_info))
cat(sprintf("  - MS2 count after filtering: %d\n", nrow(ms2_info_filtered)))

cat("\nFiltered precursor m/z range: %.2f - %.2f\n",
    min(ms2_info_filtered$precursor_mz, na.rm = TRUE),
    max(ms2_info_filtered$precursor_mz, na.rm = TRUE))

cat("\n")

# ============================================================================
# 7. Generate pseudo-MS2 - core functionality
# ============================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Step 7: Generating pseudo-MS2\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# example_mz <- c(151.02, 289.13)
example_mz <- c(151.02, 289.13)

cat(sprintf("Generating pseudo-MS2 for m/z = %s...\n", paste(example_mz, collapse = ", ")))
cat("Parameters: ppm=10, method='spearman', top_n=20\n\n")

msi_obj <- generate_pseudo_ms2(
  object = msi_obj,
  precursor_mz = example_mz,
  ppm = 10,
  top_n = 20,
  method = "spearman"
)

cat("✓ Pseudo-MS2 generation successful! \n\n")

# Display pseudo-MS2 results
pseudo_ms2_table <- msi_obj@ms2_data$pseudo_ms2_table

cat(sprintf("Pseudo-MS2 table: %d × %d\n", nrow(pseudo_ms2_table), ncol(pseudo_ms2_table)))
cat("\nPseudo-MS2 data preview:\n")
print(pseudo_ms2_table)

cat("\n")

# ============================================================================
# 8. Detect joint fragments
# ============================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Step 8: Detecting joint fragments\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat(sprintf("Detecting joint fragments for m/z = %s...\n", paste(example_mz, collapse = ", ")))
msi_obj <- detect_joint_fragments(
  object = msi_obj,
  precursor_mz = example_mz,
  ppm = 10,
  correlation_method = "spearman"
)

cat("✓ Detection complete\n\n")

if (!is.null(msi_obj@ms2_data$fragment_candidates)) {
  fragment_candidates <- msi_obj@ms2_data$fragment_candidates
  cat(sprintf("Fragment candidates found: %d × %d\n", nrow(fragment_candidates), ncol(fragment_candidates)))
  cat("\nFragment data preview:\n")
  print(utils::head(fragment_candidates, 10))
}

cat("\n")

# ============================================================================
# 9. Build fragment network
# ============================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Step 9: Building fragment network\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

if (requireNamespace("igraph", quietly = TRUE)) {
  cat("Calling build_fragment_network()...\n")
  msi_obj <- build_fragment_network(msi_obj, top_k = 10)
  cat("✓ Network built successfully!\n")
} else {
  cat("⚠ iigraph package not installed, skipping network constructionn")
}

cat("\n")

# ============================================================================
# 10. Visualize pseudo-MS2
# ============================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Step 10: Visualization\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Plot pseudo-MS2
for (target_mz in example_mz) {
  cat(sprintf("Plotting pseudo-MS2 for m/z=%.2f...\n", target_mz))
  
  p <- plot_pseudo_ms2(msi_obj, mz = target_mz, tolerance_ppm = 10)
  
  if (!is.null(p)) {
    # Save image
    filename <- file.path(output_dir, sprintf("pseudo_ms2_%.2f.png", target_mz))
    ggplot2::ggsave(
      filename = filename,
      plot = p,
      width = 10,
      height = 6,
      dpi = 150
    )
    cat(sprintf("  ✓ Saved: %s\n", basename(filename)))
  }
}

cat("\n")

# Plot spatial fragment distribution
if (requireNamespace("ggplot2", quietly = TRUE)) {
  for (target_mz in example_mz) {
    cat(sprintf("Plotting spatial fragment distribution for m/z=%.2f...\n", target_mz))
    
    p <- plot_spatial_fragment(msi_obj, precursor_mz = target_mz, ppm = 10)
    
    if (!is.null(p)) {
      filename <- file.path(output_dir, sprintf("spatial_fragment_%.2f.png", target_mz))
      ggplot2::ggsave(
        filename = filename,
        plot = p,
        width = 10,
        height = 8,
        dpi = 150
      )
      cat(sprintf("✓ Saved: %s\n", basename(filename)))
    }
  }
}

cat("\n")
if (requireNamespace("ggraph", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)) {
  net_plot <- tryCatch(
    plot_fragment_network(msi_obj),
    error = function(e) {
      cli::cli_alert_warning("plot_fragment_network() failed: {conditionMessage(e)}")
      NULL
    }
  )
  if (!is.null(net_plot)) {
    filename <- file.path(output_dir, "fragment_network_ggraph.png")
    ggplot2::ggsave(
      filename = filename,
      plot = net_plot,
      width = 10,
      height = 7,
      dpi = 150
    )
    cat(sprintf("✓ Saved: %s\n", basename(filename)))
  }
}

if (requireNamespace("igraph", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)) {
  for (layout in c("fr", "kk")) {
    net_plot <- tryCatch(
      plot_ms2_network_igraph(msi_obj, layout = layout),
      error = function(e) {
        cli::cli_alert_warning("plot_ms2_network_igraph(layout={layout}) failed: {conditionMessage(e)}")
        NULL
      }
    )
    if (!is.null(net_plot)) {
      filename <- file.path(output_dir, sprintf("fragment_network_igraph_%s.png", layout))
      ggplot2::ggsave(
        filename = filename,
        plot = net_plot,
        width = 10,
        height = 7,
        dpi = 150
      )
      cat(sprintf("✓ Saved: %s\n", basename(filename)))
    }
  }
}

cat("\n")

# ============================================================================
# 11. Pseudo-MS2 benchmark evaluation
# ============================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Step 11: Pseudo-MS2 benchmark evaluation\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat("Calling benchmark_pseudo_ms2()...\n")
msi_obj <- benchmark_pseudo_ms2(msi_obj, tolerance = 0.01)
bench_tbl <- msi_obj@other_files$pseudo_ms2_benchmark
if (!is.null(bench_tbl) && nrow(bench_tbl) > 0) {
  cat(sprintf("  - Number of compared spectra: %d\n", nrow(bench_tbl)))
  cat("  - First 3 benchmark results:\n")
  print(utils::head(bench_tbl, 3))
} else {
  cat("  - No pseudo-MS2 benchmark results found\n")
}

cat("\n")

# ============================================================================
# 12. Other visualizations and analyses
# ============================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Step 12: Other analyses\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# TIC map
cat("Plotting TIC map...\n")
p_tic <- tic_map(msi_obj, plot = TRUE)
if (!is.null(p_tic)) {
  ggplot2::ggsave(
    file.path(output_dir, "tic_map.png"),
    width = 8, height = 8, dpi = 150
  )
  cat("✓ Saved\n")
}

# Ion image
cat("\nPlotting ion image...\n")
target_ion_mz <- 191.02
p_ion <- ion_image(
  object = msi_obj,
  mz = target_ion_mz,
  tolerance_ppm = 20,
  plot = TRUE
)
if (!is.null(p_ion)) {
  ggplot2::ggsave(
    file.path(output_dir, sprintf("ion_image_%.2f_new.png", target_ion_mz)),
    width = 8, height = 8, dpi = 150
  )
  cat("✓ Saved\n")
}

# Mean spectrum
cat("\nPlotting mean spectrum...\n")
p_mean <- mean_spectrum(msi_obj, plot = TRUE)
if (!is.null(p_mean)) {
  ggplot2::ggsave(
    file.path(output_dir, "mean_spectrum_1.png"),
    width = 10, height = 6, dpi = 150
  )
  cat("✓ Saved\n")
}

# Coverage map
cat("\nPlotting coverage map...\n")
p_cov <- coverage_map(msi_obj, threshold = 0, plot = TRUE)
if (!is.null(p_cov)) {
  ggplot2::ggsave(
    file.path(output_dir, "coverage_map_1.png"),
    width = 8, height = 8, dpi = 150
  )
  cat("✓ Saved\n")
}

cat("\n")

# ============================================================================
# 13. Save results
# ============================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Step 13: Saving results\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# Save msi_dataset object
cat("Saving MSI dataset object...\n")
saveRDS(msi_obj, file.path(output_dir, "msi_dataset_final.rds"))
cat("✓ Saved to: msi_dataset_final.rds\n")

# Save pseudo-MS2 table
cat("\nSaving pseudo-MS2 table...\n")
if (!is.null(msi_obj@ms2_data$pseudo_ms2_table)) {
  write.csv(
    msi_obj@ms2_data$pseudo_ms2_table,
    file.path(output_dir, "pseudo_ms2_table.csv"),
    row.names = FALSE
  )
  cat("✓ Saved to: pseudo_ms2_table.csv\n")
}

# Save fragment candidates
cat("\nSaving fragment candidate data...\n")
if (!is.null(msi_obj@ms2_data$fragment_candidates)) {
  write.csv(
    msi_obj@ms2_data$fragment_candidates,
    file.path(output_dir, "fragment_candidates.csv"),
    row.names = FALSE
  )
  cat("✓ Saved to: fragment_candidates.csv\n")
}

# Save MS2 information
cat("\nSaving MS2 information...\n")
ms2_info_final <- extract_ms2(msi_obj, what = "ms2_info")
if (!is.null(ms2_info_final)) {
  write.csv(
    ms2_info_final,
    file.path(output_dir, "ms2_info.csv"),
    row.names = FALSE
  )
  cat("✓ Saved to: ms2_info.csv\n")
}

# Save pseudo-MS2 benchmark
cat("\nSaving pseudo-MS2 benchmark...\n")
bench_tbl <- msi_obj@other_files$pseudo_ms2_benchmark
if (!is.null(bench_tbl) && nrow(bench_tbl) > 0) {
  write.csv(
    bench_tbl,
    file.path(output_dir, "pseudo_ms2_benchmark.csv"),
    row.names = FALSE
  )
  cat("✓ Saved to: pseudo_ms2_benchmark.csv\n")
} else {
  cat(" - No pseudo-MS2 benchmark results to save\n")
}

# Save analysis summary
cat("\nSaving analysis summary...\n")
summary_file <- file.path(output_dir, "analysis_summary.txt")
sink(summary_file)

cat("msidataset Complete Demonstration Analysis Summary\n")
cat(strrep("=", 70), "\n\n")

cat("Data source information:\\n")
cat(strrep("─", 70), "\n")
cat(sprintf("Data directory: %s\n", demo_dir))
cat(sprintf("imzMLfiles:\n"))
cat("  - 03-total ion count.imzML\n")
cat("  - 04-total ion count.imzML\n")
cat(sprintf("MGF files (MS2 data):\n"))
cat("  - CHCA-neg-ms2-10ev_MSMS191.0200.mgf\n")
cat("  - NEDC-neg-ms2-10ev_MSMS130.0600.mgf\n")
cat("  - CHCA-neg-ms2-10ev_MSMS306.0700.mgf\n")

cat("\nMSI dataset statistics:\n")
cat(strrep("─", 70), "\n")
cat(sprintf("Merged dimensions: %s\n", paste(dim(msi_obj), collapse = " × ")))
cat(sprintf("Number of pixels: %d\n", n_pixels(msi_obj)))
cat(sprintf("Number of samples: %d\n", n_samples(msi_obj)))
cat(sprintf("Number of m/z variables: %d\n", n_variables(msi_obj)))
cat(sprintf("m/z range:  %.2f - %.2f\n", min(mz(msi_obj)), max(mz(msi_obj))))

cat("\nMS2 data statistics:\\n")
cat(strrep("─", 70), "\n")
cat(sprintf("Loaded MS2 spectra: %d\n", nrow(extract_ms2(msi_obj, what = "ms2_info"))))
cat(sprintf("Precursor m/z range (after filtering):  %.2f - %.2f\n",
            min(ms2_info_filtered$precursor_mz, na.rm = TRUE),
            max(ms2_info_filtered$precursor_mz, na.rm = TRUE)))

cat("\nPseudo-MS2 generation results:\n")
cat(strrep("─", 70), "\n")
cat(sprintf("Target m/z values: %s\n", paste(example_mz, collapse = ", ")))
cat(sprintf("Pseudo-MS2 table rows: %d\n", nrow(pseudo_ms2_table)))
cat(sprintf("Parameters: ppm=10, method='spearman', top_n=20\n"))

cat("\nFragment analysis:\n")
cat(strrep("─", 70), "\n")
if (!is.null(msi_obj@ms2_data$fragment_candidates)) {
  cat(sprintf("Joint fragment candidates: %d\n", nrow(msi_obj@ms2_data$fragment_candidates)))
}
cat(sprintf("Fragment network status: %s\n", 
            ifelse(!is.null(msi_obj@other_files$pseudo_fragment_network), "已构建", "未构建")))
cat(sprintf("Pseudo-MS2 benchmark entries: %d\n", if (!is.null(bench_tbl)) nrow(bench_tbl) else 0))

cat("\nOutput files:\n")
cat(strrep("─", 70), "\n")
cat("Generated files:\n")
cat("  - msi_dataset_final.rds: Complete MSI dataset object\n")
cat("  - pseudo_ms2_table.csv: Pseudo-MS2 spectra data\n")
cat("  - fragment_candidates.csv: Fragment candidate data\n")
cat("  - ms2_info.csv: MS2 information\n")
cat("  - pseudo_ms2_benchmark.csv: Pseudo-MS2 benchmark results\n")
cat("  - fragment_network_ggraph.png, fragment_network_igraph_fr.png, fragment_network_igraph_kk.png: Network visualizations\n")
cat("  - pseudo_ms2_*.png: Pseudo-MS2 mass spectra\n")
cat("  - spatial_fragment_*.png: Spatial fragment distribution maps\n")
cat("  - tic_map.png, ion_image_*.png, mean_spectrum.png, coverage_map.png: Other visualizations\n")

sink()

cat("✓ Saved to: analysis_summary.txt\n\n")

# ============================================================================
# 13. Demonstration complete
# ============================================================================
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("Demonstration complete!\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

cat("✓ Completed tasks:\n")
cat(" 1. ✓ Read two imzML files using read_imzml()\n")
cat(" 2. ✓ Merged data using merge_msi_datasets()\n")
cat(" 3. ✓ Loaded three MGF files (MS2 data) using mutate_ms2()\n")
cat(" 4. ✓ Filtered MS2 data using filter_ms2()\n")
cat(sprintf(" 5. ✓ Generated pseudo-MS2 for m/z=%.2f and %.2f\n", example_mz[1], example_mz[2]))
cat(" 6. ✓ Detected joint fragments using detect_joint_fragments()\n")
cat(" 7. ✓ Built fragment network using build_fragment_network()\n")
cat(" 8. ✓ Visualized network using plot_fragment_network() and plot_ms2_network_igraph()\n")
cat(" 9. ✓ Evaluated pseudo-MS2 spectra using benchmark_pseudo_ms2()\n")
cat(" 10. ✓ Saved all results and analysis data\n\n")

cat("📊 Output file location:\n")
cat(sprintf(" %s\n\n", normalizePath(output_dir)))

cat("📝 Available R objects:\n")
cat(" - msi_obj: Complete msi_dataset object (contains all data)\n")
cat(" - example_mz: Example m/z values [151.02, 289.13]\n\n")

cat("🎯 Suggestions for further analysis:\n")
cat(" - Use plot_spatial() for spatial analysis\n")
cat(" - Use spatial_clustering() for clustering analysis\n")
cat(" - Use ROI-related functions for region analysis\n")
cat(" - Use add_roi() and aggregate_by_roi() for ROI aggregation\n\n")

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
