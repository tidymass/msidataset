###############################################################################
# MSI Demo Workflow (msidataset)
# Data path (exact): D:/Shen_lab/msidataset/data/demo data/
###############################################################################

# =========================
# Step 1. Load environment
# =========================
devtools::document()
devtools::install()

library(msidataset)
library(methods)

# Optional plotting packages if you want additional custom plots
# library(ggplot2)

# Exact demo directory provided
demo_dir <- "D:/Shen_lab/msidataset/data/demo data"
output_dir <- file.path(demo_dir, "output1")

# ====================
# Utility: log chaining
# ====================
# Some functions return data + attach updated object as an attribute or list field.
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

# ====================
# Step 2. Data import
# ====================
msi_obj <- read_imzml(
  path_1 = demo_dir,
  peak_align_tolerance = 50,
  path_2 = output_dir
)

# Verify MS2 mapping was populated from MGF
vinfo <- variable_info(msi_obj)
cat("MS2 flag counts:\n")
print(table(vinfo$ms2, useNA = "ifany"))
stopifnot(all(vinfo$ms2 %in% c("Yes", "no")))

# ==========================================
# Step 3. Verify required sample_info columns
# ==========================================
required_cols <- c("sample_id", "x", "y", "z", "subject_id", "section_id", "resolution", "matrix")
sinfo <- sample_info(msi_obj)

missing_cols <- setdiff(required_cols, colnames(sinfo))
if (length(missing_cols) > 0) {
  stop("Missing required sample_info columns: ", paste(missing_cols, collapse = ", "))
}

cat("Required sample_info columns are present.\n")
print(head(sinfo[, required_cols, drop = FALSE]))

# ===========================================
# Step 4. QC + preprocessing (TIC + imputation)
# ===========================================
# TIC map
tic_df <- tic_map(msi_obj, plot = TRUE)
msi_obj <- extract_updated_object(tic_df, msi_obj)

# Missing value check + imputation
expr_mat <- as.matrix(expression_data(msi_obj))
has_missing <- anyNA(expr_mat) || any(!is.finite(expr_mat))

if (has_missing) {
  msi_obj <- impute_mv(msi_obj, method = "median", by = "variable")
  cat("Missing values detected and imputed using median (by variable).\n")
} else {
  cat("No missing values detected. Imputation skipped.\n")
}

# ==================================================
# Step 5. Visualization and MSI-specific analytics
# ==================================================

# 5.1 Ion image for one target m/z (use nearest to 191.0200)
target_mz <- vinfo$mz[which.min(abs(vinfo$mz - 191.0200))]
ion_res <- ion_image(
  msi_obj,
  mz = target_mz,
  tolerance_ppm = 20,
  overlays = 1,
  color_scale = "YlOrRd",
  alpha = 1,
  xlab = "x",
  ylab = "y",
  plot = TRUE
)
msi_obj <- extract_updated_object(ion_res, msi_obj)

# 5.2 Mean spectrum across all pixels
mean_spec <- mean_spectrum(msi_obj, plot = TRUE)
msi_obj <- extract_updated_object(mean_spec, msi_obj)

# 5.3 Multi-m/z overlay image
# Pick 3 representative m/z values from the current variable table
overlay_mz <- unique(round(quantile(vinfo$mz, probs = c(0.2, 0.5, 0.8), na.rm = TRUE), 4))
overlay_res <- overlay_image(msi_obj, mz = overlay_mz, alpha = 0.6, tolerance_ppm = 20)
msi_obj <- extract_updated_object(overlay_res, msi_obj)

# 5.4 Spatial clustering
cluster_df <- spatial_clustering(msi_obj, k = 4, n_pcs = 10, plot = TRUE)
msi_obj <- extract_updated_object(cluster_df, msi_obj)

# Optional: add a variable filtering step so log contains explicit filtering entry
msi_obj <- filter_msi_dataset(
  msi_obj,
  target = "variable_info",
  filter_type = "mz_range",
  mz_min = 100,
  mz_max = 1000
)

# =====================================
# Step 6. Process log / history checking
# =====================================
proc <- process_info(msi_obj)

cat("Process log keys:\n")
print(names(proc))

# Compact log summary
log_df <- do.call(
  rbind,
  lapply(seq_along(proc), function(i) {
    p <- proc[[i]]
    data.frame(
      step = i,
      function_name = p@function_name,
      package_name = p@pacakge_name,
      time = as.character(p@time),
      stringsAsFactors = FALSE
    )
  })
)
print(log_df)

# Optional detailed printout
# show_process_info(msi_obj, detailed = TRUE, verbose = FALSE)

# Save final object
saveRDS(msi_obj, file.path(output_dir, "demo_workflow_msidataset_final.rds"))
cat("Workflow complete. Final object saved.\n")