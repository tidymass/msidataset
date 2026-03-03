###############################################################################
# MSI Demo Workflow (Comprehensive API smoke demo)
# Data path: D:/Shen_lab/msidataset/data/demo data/
###############################################################################
# devtools::document()
# devtools::install()

safe_run <- function(label, expr, default = NULL) {
  cat("\n---", label, "---\n")
  tryCatch(
    {
      value <- force(expr)
      cat("[OK]\n")
      value
    },
    error = function(e) {
      cat("[ERROR]", conditionMessage(e), "\n")
      default
    }
  )
}

extract_updated_object <- function(result, fallback_object) {
  if (is.list(result) && !is.null(result$object)) {
    return(result$object)
  }
  obj_attr <- attr(result, "object")
  if (!is.null(obj_attr)) {
    return(obj_attr)
  }
  fallback_object
}

safe_update_object <- function(label, result, object) {
  obj <- extract_updated_object(result, object)
  if (methods::is(obj, "msi_dataset")) {
    return(obj)
  }
  cat("[WARN]", label, "did not return an updated msi_dataset object.\n")
  object
}

safe_plot <- function(label, plot_obj) {
  if (inherits(plot_obj, "ggplot")) {
    safe_run(paste0("print(", label, ")"), print(plot_obj))
  }
}

safe_save_plot <- function(plot_obj, filename, width = 8, height = 6, dpi = 150) {
  if (inherits(plot_obj, "ggplot")) {
    safe_run(
      paste0("save plot: ", filename),
      ggplot2::ggsave(
        filename = file.path(output_dir, filename),
        plot = plot_obj,
        width = width,
        height = height,
        dpi = dpi
      )
    )
  }
}

# =========================
# Step 1. Build and install
# =========================
# safe_run("devtools::document()", devtools::document())
# safe_run("devtools::test()", devtools::test())
# safe_run("devtools::install()", devtools::install())

library(msidataset)
library(methods)

demo_dir <- "D:/Shen_lab/msidataset/data/demo data"
output_dir <- "D:/Shen_lab/msidataset/results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(output_dir, "test_demo_data_console_log.txt")
zz <- file(log_file, open = "wt")
sink(zz, split = TRUE)
on.exit({
  sink()
  close(zz)
}, add = TRUE)

# ====================
# Step 2. Data import
# ====================
msi_obj <- safe_run(
  "read_imzml()",
  read_imzml(
    path_1 = demo_dir,
    peak_align_tolerance = 50,
    path_2 = output_dir,
    resolution = 50,
    matrix_name = "CHCA"
  )
)

if (is.null(msi_obj) || !methods::is(msi_obj, "msi_dataset")) {
  stop("read_imzml failed; cannot continue demo.")
}

# ==========================================
# Step 3. Base checks + utility functions
# ==========================================
safe_run("msidataset_logo()", msidataset_logo())
safe_run("msidataset_packages()", print(msidataset_packages()))

conflicts_obj <- safe_run("msidataset_conflicts()", msidataset_conflicts())
if (!is.null(conflicts_obj)) {
  safe_run("print.msidataset_conflicts()", print(conflicts_obj))
}

safe_run(
  "check_column_name()",
  {
    tmp_df <- data.frame(x = 1:2, y = 2:3, y.1 = 3:4)
    print(check_column_name(tmp_df, "y"))
    print(check_column_name(tmp_df, "z"))
  }
)

safe_run("check_msi_dataset_class()", print(check_msi_dataset_class(msi_obj)))
safe_run("show(msi_obj)", show(msi_obj))
safe_run("dim(msi_obj)", print(dim(msi_obj)))
safe_run("dimnames(msi_obj)", print(lapply(dimnames(msi_obj), head)))

vinfo <- variable_info(msi_obj)
sinfo <- sample_info(msi_obj)
expr <- as.matrix(expression_data(msi_obj))

safe_run("n_variables()", print(n_variables(msi_obj)))
safe_run("n_samples()", print(n_samples(msi_obj)))
safe_run("n_pixels()", print(n_pixels(msi_obj)))
safe_run("mz()", print(utils::head(mz(msi_obj))))
safe_run("get_sample_id()", print(utils::head(get_sample_id(msi_obj))))
safe_run("get_variable_id()", print(utils::head(get_variable_id(msi_obj))))
safe_run("get_spatial_coords()", print(utils::head(get_spatial_coords(msi_obj))))
safe_run("get_section_ids()", print(get_section_ids(msi_obj)))

section_ids <- unique(as.character(sinfo$section_id))
if (length(section_ids) > 0) {
  safe_run(
    "get_section_pixels()",
    print(utils::head(get_section_pixels(msi_obj, section_id = section_ids[[1]])))
  )
}

safe_run("pixel_info()", print(utils::head(pixel_info(msi_obj))))
safe_run("pixel_info_note()", print(utils::head(pixel_info_note(msi_obj))))

safe_run("activate_msi_dataset(sample_info)", {
  msi_obj <<- activate_msi_dataset(msi_obj, what = "sample_info")
  print(msi_obj@activated)
})
safe_run("activate_msi_dataset(variable_info)", {
  msi_obj <<- activate_msi_dataset(msi_obj, what = "variable_info")
  print(msi_obj@activated)
})
safe_run("activate_msi_dataset(expression_data)", {
  msi_obj <<- activate_msi_dataset(msi_obj, what = "expression_data")
  print(msi_obj@activated)
})

# ==========================================
# Step 4. Accessor replacement methods
# ==========================================
safe_run("expression_data<-", {
  expression_data(msi_obj) <- expression_data(msi_obj)
  print(dim(expression_data(msi_obj)))
})
safe_run("sample_info<-", {
  sample_info(msi_obj) <- sample_info(msi_obj)
  print(dim(sample_info(msi_obj)))
})
safe_run("variable_info<-", {
  variable_info(msi_obj) <- variable_info(msi_obj)
  print(dim(variable_info(msi_obj)))
})
safe_run("sample_info_note<-", {
  sample_info_note(msi_obj) <- sample_info_note(msi_obj)
  print(dim(sample_info_note(msi_obj)))
})
safe_run("variable_info_note<-", {
  variable_info_note(msi_obj) <- variable_info_note(msi_obj)
  print(dim(variable_info_note(msi_obj)))
})
safe_run("ms2_data<-", {
  ms2_data(msi_obj) <- ms2_data(msi_obj)
  print(length(ms2_data(msi_obj)))
})

# ==========================================
# Step 5. Notes, IDs, constructor checks
# ==========================================
safe_run(
  "create_enhanced_sample_info_note()",
  print(utils::head(create_enhanced_sample_info_note(msi_obj)))
)
safe_run(
  "create_enhanced_variable_info_note()",
  print(utils::head(create_enhanced_variable_info_note(msi_obj)))
)
safe_run("update_msi_notes()", {
  msi_obj <<- update_msi_notes(msi_obj)
  print(dim(sample_info_note(msi_obj)))
  print(dim(variable_info_note(msi_obj)))
})

safe_run(
  "generate_variable_id()/generate_variable_id_extended()/validate_variable_id()",
  {
    mz_demo <- c(100.1234, 100.1235, 250.2233, 900.9988)
    var_id1 <- generate_variable_id(mz_demo)
    var_id2 <- generate_variable_id_extended(mz_demo, zero_pad = TRUE, pad_width = 4)
    print(var_id1)
    print(var_id2)
    print(validate_variable_id(var_id1, mz_demo))
  }
)

safe_run("create_msi_dataset() minimal", {
  n_row <- min(20, nrow(expression_data(msi_obj)))
  n_col <- min(20, ncol(expression_data(msi_obj)))
  mini_expr <- as.matrix(expression_data(msi_obj))[seq_len(n_row), seq_len(n_col), drop = FALSE]
  mini_sinfo <- sample_info(msi_obj)[seq_len(n_col), , drop = FALSE]
  mini_vinfo <- variable_info(msi_obj)[seq_len(n_row), , drop = FALSE]
  mini <- create_msi_dataset(
    expression_data = mini_expr,
    sample_info = mini_sinfo,
    variable_info = mini_vinfo
  )
  print(dim(mini))
})

safe_run("update_msi_dataset()", {
  msi_obj <<- update_msi_dataset(msi_obj)
  print(dim(msi_obj))
})

safe_run(
  "convert_cardinal2msi_dataset() smoke",
  {
    if (!requireNamespace("Cardinal", quietly = TRUE)) {
      cat("Skip: Cardinal not installed.\n")
    } else {
      imzml_candidates <- list.files(demo_dir, pattern = "\\.imzML$", ignore.case = TRUE, full.names = TRUE)
      if (length(imzml_candidates) == 0) {
        cat("Skip: no imzML file found in demo_dir.\n")
      } else {
        cardinal_object <- Cardinal::readMSIData(imzml_candidates[[1]])
        obj_from_cardinal <- convert_cardinal2msi_dataset(
          cardinal_object = cardinal_object,
          sample_name = "subject_demo",
          section_name = "section_demo",
          resolution = 50,
          matrix_name = "CHCA"
        )
        print(dim(obj_from_cardinal))
      }
    }
  }
)

# ==========================================
# Step 6. Annotation table check + setter
# ==========================================
safe_run(
  "check_annotation_table() and annotation_table<-",
  {
    ann <- data.frame(
      variable_id = head(variable_info(msi_obj)$variable_id, 2),
      ms2_files_id = c("file1", "file2"),
      ms2_spectrum_id = c("spec1", "spec2"),
      Compound.name = c("cmp1", "cmp2"),
      CAS.ID = c("cas1", "cas2"),
      HMDB.ID = c("hmdb1", "hmdb2"),
      KEGG.ID = c("kegg1", "kegg2"),
      Lab.ID = c("lab1", "lab2"),
      Adduct = c("[M+H]+", "[M+Na]+"),
      mz.error = c(0.001, 0.002),
      mz.match.score = c(0.95, 0.96),
      RT.error = c(0.1, 0.2),
      RT.match.score = c(0.9, 0.92),
      CE = c(20, 30),
      SS = c("ms2_1", "ms2_2"),
      Total.score = c(0.88, 0.77),
      Database = c("HMDB", "KEGG"),
      Level = c("Level 1", "Level 2"),
      stringsAsFactors = FALSE
    )
    check_annotation_table(ann)
    annotation_table(msi_obj) <- ann
    print(annotation_table(msi_obj))
  }
)

# ==========================================
# Step 7. Process log helpers
# ==========================================
safe_run("add_process_log()", {
  msi_obj <<- add_process_log(
    object = msi_obj,
    function_name = "manual_demo_step()",
    parameter = list(note = "demo")
  )
  print(length(process_info(msi_obj)))
})

safe_run("show_process_info()", show_process_info(msi_obj, detailed = FALSE, verbose = FALSE))
safe_run("show_process_timeline()", show_process_timeline(msi_obj))
proc_export <- safe_run("export_process_history()", export_process_history(msi_obj, include_details = TRUE))
if (!is.null(proc_export)) {
  safe_run("write process history csv", {
    utils::write.csv(proc_export, file.path(output_dir, "demo_process_history_export.csv"), row.names = FALSE)
  })
}

# ==========================================
# Step 8. QC and supplementary filter eval
# ==========================================
safe_run("show_filter_types()", show_filter_types())
safe_run("show_mz_statistics()", show_mz_statistics(msi_obj))
safe_run("recommend_filter_params()", recommend_filter_params(msi_obj))

safe_run("impute_mv()", {
  msi_obj <<- impute_mv(msi_obj, method = "median", by = "variable")
  print(dim(expression_data(msi_obj)))
})

safe_run("filter_out_zero_heavy_mz()", {
  obj_zero <- filter_out_zero_heavy_mz(msi_obj, zero_threshold = 0.95)
  print(dim(obj_zero))
})
safe_run("filter_high_frequency_mz_by_count()", {
  obj_hf <- filter_high_frequency_mz_by_count(msi_obj, min_pixels = 1)
  print(dim(obj_hf))
})
safe_run("filter_supplementary()", {
  obj_sup <- filter_supplementary(msi_obj, min_pixels = 1, zero_threshold = 0.95)
  print(dim(obj_sup))
})
safe_run("filter_msi_dataset_supplementary(high_frequency_count)", {
  obj_sup2 <- filter_msi_dataset_supplementary(
    msi_obj,
    target = "variable_info",
    filter_type = "high_frequency_count",
    min_pixels = 1
  )
  print(dim(obj_sup2))
})

safe_run("filter_msi_dataset(mz_range)", {
  mz_vals <- variable_info(msi_obj)$mz
  msi_filtered <- filter_msi_dataset(
    msi_obj,
    target = "variable_info",
    filter_type = "mz_range",
    mz_min = stats::quantile(mz_vals, 0.05, na.rm = TRUE),
    mz_max = stats::quantile(mz_vals, 0.95, na.rm = TRUE)
  )
  print(dim(msi_filtered))
})

# ==========================================
# Step 9. ROI lifecycle + ROI-aware analysis
# ==========================================
current_sinfo <- sample_info(msi_obj)
current_ids <- as.character(current_sinfo$sample_id)
roi_pixels_1 <- head(current_ids, min(30, length(current_ids)))
roi_pixels_2 <- tail(setdiff(current_ids, roi_pixels_1), min(30, length(current_ids)))
sid <- as.character(current_sinfo$section_id[[1]])

roi_obj_1 <- safe_run(
  "create_roi()",
  create_roi(
    roi_id = "roi_manual_1",
    roi_name = "Manual ROI 1",
    sample_id = roi_pixels_1,
    section_id = sid,
    metadata = list(source = "manual_head_pixels")
  )
)

safe_run("add_roi(roi object)", {
  msi_obj <<- add_roi(msi_obj, roi = roi_obj_1, overwrite = TRUE)
  print(list_roi(msi_obj))
})

safe_run("add_roi(direct args)", {
  msi_obj <<- add_roi(
    object = msi_obj,
    roi_id = "roi_manual_2",
    roi_name = "Manual ROI 2",
    sample_id = roi_pixels_2,
    section_id = sid,
    metadata = list(source = "manual_tail_pixels"),
    overwrite = TRUE
  )
  print(list_roi(msi_obj))
})

safe_run("add_roi_from_spatial()", {
  msi_obj <<- add_roi_from_spatial(
    object = msi_obj,
    roi_id = "roi_rect_1",
    roi_name = "Rectangle ROI",
    type = "rectangle",
    roi = list(
      xmin = min(current_sinfo$x, na.rm = TRUE),
      xmax = stats::median(current_sinfo$x, na.rm = TRUE),
      ymin = min(current_sinfo$y, na.rm = TRUE),
      ymax = stats::median(current_sinfo$y, na.rm = TRUE)
    ),
    section_id = sid,
    overwrite = TRUE
  )
  print(list_roi(msi_obj))
})

safe_run("rename_roi()", {
  msi_obj <<- rename_roi(msi_obj, roi_id = "roi_manual_2", new_name = "Renamed ROI 2")
  print(list_roi(msi_obj))
})

roi_summary <- safe_run("list_roi()", list_roi(msi_obj), default = data.frame())
roi_membership <- safe_run("roi_data()", roi_data(msi_obj), default = data.frame())
safe_run("validate_roi_integrity()", print(validate_roi_integrity(msi_obj)))

safe_run("aggregate_by_roi(variable_id)", {
  agg1 <- aggregate_by_roi(msi_obj, feature = NULL, feature_type = "variable_id", fun = "mean")
  print(utils::head(agg1))
})

target_mz <- variable_info(msi_obj)$mz[which.min(abs(variable_info(msi_obj)$mz - 191.0200))]
safe_run("aggregate_by_roi(mz)", {
  agg2 <- aggregate_by_roi(msi_obj, feature = target_mz, feature_type = "mz", fun = "median", tolerance_ppm = 20)
  print(agg2)
})

safe_run("subset_by_roi()", {
  subset_obj <- subset_by_roi(msi_obj, roi_id = c("roi_manual_1", "roi_rect_1"), keep_roi = FALSE)
  print(dim(subset_obj))
})

safe_run("filter_msi_dataset(roi)", {
  roi_ids_all <- list_roi(msi_obj)$roi_id
  obj_roi <- filter_msi_dataset(
    msi_obj,
    target = "sample_info",
    filter_type = "roi",
    roi_id = roi_ids_all
  )
  print(dim(obj_roi))
})

safe_run("roi_visualize(save_roi = TRUE)", {
  roi_vis <- roi_visualize(
    object = msi_obj,
    type = "rectangle",
    roi = list(
      xmin = min(current_sinfo$x, na.rm = TRUE),
      xmax = stats::quantile(current_sinfo$x, 0.25, na.rm = TRUE),
      ymin = min(current_sinfo$y, na.rm = TRUE),
      ymax = stats::quantile(current_sinfo$y, 0.25, na.rm = TRUE)
    ),
    plot = TRUE,
    save_roi = TRUE,
    roi_id = "roi_visual_1",
    roi_name = "ROI from roi_visualize",
    overwrite = TRUE
  )
  msi_obj <<- safe_update_object("roi_visualize()", roi_vis, msi_obj)
  print(utils::head(roi_vis))
})

safe_run("remove_roi()", {
  msi_obj <<- remove_roi(msi_obj, roi_id = "roi_manual_1")
  print(list_roi(msi_obj))
})

# ==========================================
# Step 10. Visualization and eval functions
# ==========================================
vinfo <- variable_info(msi_obj)
target_mz <- vinfo$mz[which.min(abs(vinfo$mz - 191.0200))]
multi_mz <- unique(round(stats::quantile(vinfo$mz, probs = c(0.2, 0.5, 0.8), na.rm = TRUE), 4))
sample_id_1 <- sample_info(msi_obj)$sample_id[[1]]
var_id_1 <- variable_info(msi_obj)$variable_id[[1]]
all_sections <- unique(as.character(sample_info(msi_obj)$section_id))

tic_res <- safe_run("tic_map()", tic_map(msi_obj, plot = TRUE))
msi_obj <- safe_update_object("tic_map()", tic_res, msi_obj)

ion_res <- safe_run(
  "ion_image()",
  ion_image(
    object = msi_obj,
    mz = target_mz,
    tolerance_ppm = 20,
    overlays = 1,
    color_scale = "YlOrRd",
    alpha = 1,
    xlab = "x",
    ylab = "y",
    plot = TRUE
  )
)
msi_obj <- safe_update_object("ion_image()", ion_res, msi_obj)

mean_res <- safe_run("mean_spectrum()", mean_spectrum(msi_obj, plot = TRUE))
msi_obj <- safe_update_object("mean_spectrum()", mean_res, msi_obj)

safe_run("pixel_spectrum()", {
  px <- pixel_spectrum(msi_obj, sample_id = sample_id_1, plot = TRUE)
  msi_obj <<- safe_update_object("pixel_spectrum()", px, msi_obj)
  print(utils::head(px))
})

overlay_res <- safe_run("overlay_image()", overlay_image(msi_obj, mz = multi_mz, alpha = 0.6, tolerance_ppm = 20))
msi_obj <- safe_update_object("overlay_image()", overlay_res, msi_obj)

safe_run("pixel_similarity_map()", {
  sim_res <- pixel_similarity_map(msi_obj, sample_id = sample_id_1, method = "pearson", plot = TRUE)
  msi_obj <<- safe_update_object("pixel_similarity_map()", sim_res, msi_obj)
  print(utils::head(sim_res))
})

safe_run("coverage_map()", {
  cov_res <- coverage_map(msi_obj, threshold = 0, plot = TRUE)
  msi_obj <<- safe_update_object("coverage_map()", cov_res, msi_obj)
  print(utils::head(cov_res))
})

safe_run("spatial_clustering()", {
  clu_res <- spatial_clustering(msi_obj, k = 4, n_pcs = 10, plot = TRUE)
  msi_obj <<- safe_update_object("spatial_clustering()", clu_res, msi_obj)
  print(utils::head(clu_res))
})

safe_run("embedding_visualize(PCA)", {
  emb_res <- embedding_visualize(msi_obj, method = "PCA", dims = 2, color_by = "section_id")
  msi_obj <<- safe_update_object("embedding_visualize()", emb_res, msi_obj)
  print(utils::head(emb_res))
})

safe_run("ms2_spatial_map(binary)", {
  ms2_res <- ms2_spatial_map(msi_obj, mode = "binary", plot = TRUE)
  msi_obj <<- safe_update_object("ms2_spatial_map()", ms2_res, msi_obj)
  print(utils::head(ms2_res))
})

safe_run("ms2_spatial_map(count)", {
  ms2_res2 <- ms2_spatial_map(msi_obj, mode = "count", plot = TRUE)
  msi_obj <<- safe_update_object("ms2_spatial_map()", ms2_res2, msi_obj)
  print(utils::head(ms2_res2))
})

p_single <- safe_run(
  "plot_msi_filtered()",
  plot_msi_filtered(msi_obj, variable_id = var_id_1, transform = "none", color_scale = "viridis")
)
safe_plot("plot_msi_filtered", p_single)
safe_save_plot(p_single, "plot_msi_filtered_single.png")

safe_run(
  "plot_msi_grid()",
  plot_msi_grid(msi_obj, variable_ids = head(variable_info(msi_obj)$variable_id, 6), ncol = 3)
)

corr_obj <- safe_run(
  "plot_spatial_correlation()",
  plot_spatial_correlation(msi_obj, n_vars = 20, method = "pearson", show_labels = FALSE)
)
if (!is.null(corr_obj)) {
  safe_run("print.spatial_correlation()", print(corr_obj))
  safe_run("summary.spatial_correlation()", summary(corr_obj, top_n = 5))
  safe_plot("spatial correlation plot_x", corr_obj$plot_x)
  safe_plot("spatial correlation plot_y", corr_obj$plot_y)
  safe_save_plot(corr_obj$plot_x, "plot_spatial_correlation_x.png")
  safe_save_plot(corr_obj$plot_y, "plot_spatial_correlation_y.png")
}

p_sp1 <- safe_run(
  "plot_spatial(single feature)",
  plot_spatial(
    object = msi_obj,
    feature = target_mz,
    feature_type = "mz",
    facet_by = c("feature"),
    scales = "fixed",
    aggregation = "none",
    section_scale = "global"
  )
)
safe_plot("plot_spatial(single)", p_sp1)
safe_save_plot(p_sp1, "plot_spatial_single.png")
msi_obj <- safe_update_object("plot_spatial(single)", p_sp1, msi_obj)

p_sp2 <- safe_run(
  "plot_spatial(multi feature)",
  plot_spatial(
    object = msi_obj,
    feature = multi_mz,
    feature_type = "mz",
    facet_by = c("feature"),
    scales = "free",
    aggregation = "none",
    section_scale = "global"
  )
)
safe_plot("plot_spatial(multi)", p_sp2)
safe_save_plot(p_sp2, "plot_spatial_multi.png")
msi_obj <- safe_update_object("plot_spatial(multi)", p_sp2, msi_obj)

p_sp3 <- safe_run(
  "plot_spatial(section facet)",
  plot_spatial(
    object = msi_obj,
    feature = target_mz,
    feature_type = "mz",
    section_id = all_sections,
    facet_by = c("section"),
    scales = "fixed",
    aggregation = "none",
    section_scale = "section"
  )
)
safe_plot("plot_spatial(section)", p_sp3)
safe_save_plot(p_sp3, "plot_spatial_section.png")
msi_obj <- safe_update_object("plot_spatial(section)", p_sp3, msi_obj)

p_sp4 <- safe_run(
  "plot_spatial(roi overlay)",
  {
    roi_ids_available <- list_roi(msi_obj)$roi_id
    if (length(roi_ids_available) == 0) {
      stop("No ROI available for roi overlay plot.")
    }
    plot_spatial(
      object = msi_obj,
      feature = target_mz,
      feature_type = "mz",
      roi_id = roi_ids_available[[1]],
      overlay_roi = TRUE,
      facet_by = c("roi"),
      scales = "fixed",
      aggregation = "mean",
      section_scale = "global"
    )
  }
)
safe_plot("plot_spatial(roi)", p_sp4)
safe_save_plot(p_sp4, "plot_spatial_roi_overlay.png")
msi_obj <- safe_update_object("plot_spatial(roi)", p_sp4, msi_obj)

# ==========================================
# Step 11. Save outputs
# ==========================================
safe_run("write outputs", {
  saveRDS(msi_obj, file.path(output_dir, "demo_workflow_msidataset_final.rds"))
  saveRDS(list_roi(msi_obj), file.path(output_dir, "demo_roi_summary_final.rds"))
  saveRDS(roi_data(msi_obj), file.path(output_dir, "demo_roi_membership_final.rds"))
  saveRDS(process_info(msi_obj), file.path(output_dir, "demo_process_info_final.rds"))
  saveRDS(
    list(
      p_single = p_single,
      corr_obj = corr_obj,
      p_sp1 = p_sp1,
      p_sp2 = p_sp2,
      p_sp3 = p_sp3,
      p_sp4 = p_sp4
    ),
    file.path(output_dir, "demo_plot_objects.rds")
  )
  if (is.data.frame(roi_summary) && nrow(roi_summary) > 0) {
    utils::write.csv(roi_summary, file.path(output_dir, "demo_roi_summary.csv"), row.names = FALSE)
  }
  if (is.data.frame(roi_membership) && nrow(roi_membership) > 0) {
    utils::write.csv(roi_membership, file.path(output_dir, "demo_roi_membership.csv"), row.names = FALSE)
  }
  proc_df <- try(export_process_history(msi_obj, include_details = TRUE), silent = TRUE)
  if (!inherits(proc_df, "try-error") && !is.null(proc_df)) {
    utils::write.csv(proc_df, file.path(output_dir, "demo_process_history_final.csv"), row.names = FALSE)
  }
  invisible(TRUE)
})

cat("\nWorkflow complete. Comprehensive API smoke demo finished.\n")
