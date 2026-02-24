# ============================================================================
# msidataset — Complete Debug / Demo Script
# Author : Xiaoyu Xu  <XIAOYU021@e.ntu.edu.sg>
# Purpose: Full workflow test using real .imzML data via Cardinal.
#
# Workflow:
#   Step 1 : Read .imzML -> Cardinal -> msi_dataset
#   Step 2 : Validate & explore data structure
#   Step 3 : Accessor methods (getters / setters)
#   Step 4 : Subsetting with [
#   Step 5 : filter_pixels() / filter_variables()
#   Step 6 : activate() + process_info
#   Step 7 : update_msi_dataset()
#   Step 8 : msi_dataset <-> massdataset round-trip
#   Step 9 : msi_dataset -> Cardinal (back-conversion)
#   Step 10: Update metadata (tissue type, organism, date ...)
#   Step 11: Save & reload final check
# ============================================================================

# ----------------------------------------------------------------------------
# 0. LOAD / BUILD PACKAGE
# ----------------------------------------------------------------------------

setwd("D:/Shen_lab/msidataset")
devtools::document()
devtools::load_all()

# library(msidataset)   # use this after package is installed

library(S4Vectors)
library(Cardinal)
library(massdataset)
library(dplyr)
library(crayon)


# ============================================================================
# STEP 1 - Read real .imzML and convert to msi_dataset
# ============================================================================

## read data
cardinal_obj <- readMSIData(
  "data/demo data/03-total ion count.imzML",
  memory = FALSE
)
cat("Cardinal class:", class(cardinal_obj)[1], "\n")
# print(cardinal_obj)
# print(experimentData)

## aligned
cardinal_obj_aligned <- Cardinal::peakAlign(
  cardinal_obj,
  tolerance = 10,
  units = "ppm"
)

my_msi_full <- convert_cardinal2msi_dataset(
  cardinal_object = cardinal_obj_aligned,
  sample_name     = "renal_cancer",
  section_name    = "slice_01",
  resolution      = 50,
  matrix_name     = "NEDC"
)


## overview
print(my_msi_full)

## save
if (!dir.exists("data/output")) dir.create("data/output", recursive = TRUE)
saveRDS(my_msi_full, "data/output/my_msi_dataset_full.rds")
cat("Saved to data/output/my_msi_dataset_full.rds\n")

## reload (confirms file integrity; use this object from here on)
my_msi <- readRDS("data/output/my_msi_dataset_full.rds")
cat("Reloaded dim:", paste(dim(my_msi), collapse = " x "), "\n")


# ============================================================================
# STEP 2 - Validate & explore data structure
# ============================================================================

## -- 2.1 Validity check ------------------------------------------------------
check_result <- check_msi_dataset_class(my_msi)
if (isTRUE(check_result)) {
  cat(crayon::green("Object structure is valid.\n"))
} else {
  cat(crayon::red("Validation errors:\n"))
  print(check_result)
}

## -- 2.2 Expression data preview ---------------------------------------------
cat("\n--- Expression Data (first 5 variables x first 5 pixels) ---\n")
print(expression_data(my_msi)[1:5, 1:5])

## -- 2.3 Intensity statistics ------------------------------------------------
expr_mat <- as.matrix(expression_data(my_msi))
na_n     <- sum(is.na(expr_mat))
zero_n   <- sum(expr_mat == 0, na.rm = TRUE)
total_n  <- prod(dim(expr_mat))
cat("\nData quality summary:\n")
cat(sprintf("  Variables       : %d\n",   n_variables(my_msi)))
cat(sprintf("  Pixels          : %d\n",   n_pixels(my_msi)))
cat(sprintf("  Missing values  : %d\n",   na_n))
cat(sprintf("  Zero values     : %d  (%.1f%%)\n", zero_n, 100 * zero_n / total_n))
cat(sprintf("  Intensity min   : %.2f\n", min(expr_mat, na.rm = TRUE)))
cat(sprintf("  Intensity max   : %.2f\n", max(expr_mat, na.rm = TRUE)))
cat(sprintf("  Intensity mean  : %.2f\n", mean(expr_mat, na.rm = TRUE)))

## -- 2.4 Pixel spatial summary -----------------------------------------------
sinfo <- sample_info(my_msi)
cat("\nPixel spatial summary:\n")
cat("  x positions :", length(unique(sinfo$x)),
    "| range [", min(sinfo$x), "-", max(sinfo$x), "]\n")
cat("  y positions :", length(unique(sinfo$y)),
    "| range [", min(sinfo$y), "-", max(sinfo$y), "]\n")
if ("z" %in% colnames(sinfo))
  cat("  z slices    :", length(unique(sinfo$z)), "\n")

## -- 2.5 m/z summary ---------------------------------------------------------
vinfo <- variable_info(my_msi)
cat("\nm/z summary:\n")
cat(sprintf("  n features  : %d\n",     nrow(vinfo)))
cat(sprintf("  Min m/z     : %.4f\n",   min(vinfo$mz)))
cat(sprintf("  Max m/z     : %.4f\n",   max(vinfo$mz)))
cat(sprintf("  Median m/z  : %.4f\n",   median(vinfo$mz)))

## -- 2.6 sample_info / variable_info first rows ------------------------------
cat("\nsample_info (first 10 rows):\n")
print(head(sinfo, 10))

cat("\nvariable_info (first 20 rows):\n")
print(head(vinfo, 20))

cat("\nsample_info_note:\n")
print(sample_info_note(my_msi))

cat("\nvariable_info_note:\n")
print(variable_info_note(my_msi))


# ============================================================================
# STEP 3 - Accessor methods (getters & setters)
# ============================================================================

## -- 3.1 Count getters -------------------------------------------------------
cat("n_variables() :", n_variables(my_msi), "\n")
cat("n_pixels()    :", n_pixels(my_msi), "\n")
cat("n_samples()   :", n_samples(my_msi), "(= n_pixels)\n")
cat("dim()         :", paste(dim(my_msi), collapse = " x "), "\n")

## -- 3.2 Spatial / id helpers ------------------------------------------------
cat("\nget_spatial_coords() (first 6 rows):\n")
print(head(get_spatial_coords(my_msi), 6))

cat("\nget_section_ids()  :", paste(get_section_ids(my_msi), collapse = ", "), "\n")
cat("get_sample_id()   [1:5]:", paste(head(get_sample_id(my_msi), 5), collapse = " "), "\n")
cat("get_variable_id() [1:5]:", paste(head(get_variable_id(my_msi), 5), collapse = " "), "\n")

## -- 3.3 mz() ----------------------------------------------------------------
cat("\nmz() first 5:", paste(round(head(mz(my_msi), 5), 4), collapse = ", "), "\n")

## -- 3.4 Backward-compat pixel_info() ----------------------------------------
suppressMessages({
  cat("\npixel_info() [first 3 rows]:\n")
  print(head(pixel_info(my_msi), 3))
})

## -- 3.5 Setter: add metadata columns ----------------------------------------
my_msi2 <- my_msi

# Add annotation column to variable_info
new_vinfo              <- variable_info(my_msi2)
new_vinfo$annotation   <- NA_character_   # placeholder; fill after annotation step
variable_info(my_msi2) <- new_vinfo
cat("\nvariable_info after setter (cols):",
    paste(colnames(variable_info(my_msi2)), collapse = ", "), "\n")

# Add tissue metadata to sample_info
new_sinfo               <- sample_info(my_msi2)
new_sinfo$tissue_type   <- "renal_tumor"
new_sinfo$organism      <- "human"
new_sinfo$date          <- "2026-02-10"
sample_info(my_msi2)    <- new_sinfo
cat("sample_info after setter (cols):",
    paste(colnames(sample_info(my_msi2)), collapse = ", "), "\n")

rm(my_msi2)


# ============================================================================
# STEP 4 - Subsetting with [
# ============================================================================

## numeric index
sub_num <- my_msi[1:50, 1:10]
cat("my_msi[1:50, 1:10] dim:", paste(dim(sub_num), collapse = " x "), "\n")
cat("  valid:", if (isTRUE(check_msi_dataset_class(sub_num))) "yes" else "NO", "\n")

# ============================================================================
# STEP 5 - filter_pixels() and filter_variables()
# ============================================================================

x_mid  <- median(sample_info(my_msi)$x)
y_mid  <- median(sample_info(my_msi)$y)
mz_lo  <- quantile(mz(my_msi), 0.25)
mz_hi  <- quantile(mz(my_msi), 0.75)

cat(sprintf("Spatial midpoints: x=%.0f  y=%.0f\n", x_mid, y_mid))
cat(sprintf("m/z IQR window   : [%.2f, %.2f]\n", mz_lo, mz_hi))

## filter pixels (top-right spatial quadrant)
f_pix <- filter_pixels(my_msi, x > !!x_mid, y > !!y_mid)
cat(sprintf("filter_pixels(x > %.0f & y > %.0f): %d pixels remain\n",
            x_mid, y_mid, n_pixels(f_pix)))

## filter variables (m/z window)
f_var <- filter_variables(my_msi, mz >= !!mz_lo, mz <= !!mz_hi)
cat(sprintf("filter_variables(%.2f <= mz <= %.2f): %d variables remain\n",
            mz_lo, mz_hi, n_variables(f_var)))

## dplyr::filter() dispatched via filter.msi_dataset
f_dplyr <- dplyr::filter(my_msi, x <= !!x_mid, target = "sample_info")
cat(sprintf("dplyr::filter(x <= %.0f): %d pixels remain\n",
            x_mid, n_pixels(f_dplyr)))


# ============================================================================
# STEP 6 - activate() + process_info
# ============================================================================

my_msi_act <- activate_msi_dataset(my_msi, "variable_info")
cat("Activated slot:", my_msi_act@activated, "\n")

mz_cutoff  <- median(mz(my_msi_act))
my_msi_act <- dplyr::filter(my_msi_act, mz < !!mz_cutoff)
cat(sprintf("After filter(mz < %.4f): %d variables remain\n",
            mz_cutoff, n_variables(my_msi_act)))

cat("\nProcess info steps:", length(process_info(my_msi_act)), "\n")
cat("Step names:", paste(names(process_info(my_msi_act)), collapse = ", "), "\n")


# ============================================================================
# STEP 7 - update_msi_dataset()
# ============================================================================

my_msi_upd                   <- my_msi
my_msi_upd@sample_info_note  <- data.frame()   # simulate stale / missing note
my_msi_upd <- update_msi_dataset(my_msi_upd)
cat("sample_info_note rows after update:", nrow(my_msi_upd@sample_info_note), "\n")
cat("  column names:", paste(my_msi_upd@sample_info_note$name, collapse = ", "), "\n")
rm(my_msi_upd)


# ============================================================================
# STEP 8 - msi_dataset <-> massdataset round-trip
# ============================================================================

## -- 8.1 Convert TO massdataset ---------------------------------------------
my_mass <- convert_msi2mass_dataset(
  msi_object                 = my_msi,
  save_metadata              = TRUE,
  metadata_path              = "data/output/msi_metadata.rds",
  add_spatial_to_sample_info = FALSE
)

cat("mass_dataset class     :", class(my_mass), "\n")
cat("mass_dataset variables :", nrow(my_mass@variable_info), "\n")
cat("mass_dataset samples   :", nrow(my_mass@sample_info), "\n")
cat("sample_info columns    :", paste(colnames(my_mass@sample_info), collapse = ", "), "\n")
cat("\nmass_dataset sample_info (first 5 rows):\n")
print(head(my_mass@sample_info, 5))

## -- 8.2 Restore from .GlobalEnv (same session) -----------------------------
my_msi_restored_a <- convert_mass2msi_dataset(
  mass_object             = my_mass,
  metadata_source         = "auto",
  create_annotation_table = FALSE
)
cat("\nRestored (auto) dim:", paste(dim(my_msi_restored_a), collapse = " x "), "\n")
cat("  valid:", if (isTRUE(check_msi_dataset_class(my_msi_restored_a))) "yes" else "NO", "\n")

## -- 8.3 Restore from file (cross-session safe) -----------------------------
my_msi_restored_b <- convert_mass2msi_dataset(
  mass_object             = my_mass,
  metadata_source         = "file",
  metadata_path           = "data/output/msi_metadata.rds",
  create_annotation_table = FALSE
)
cat("Restored (file)  dim:", paste(dim(my_msi_restored_b), collapse = " x "), "\n")
cat("  valid:", if (isTRUE(check_msi_dataset_class(my_msi_restored_b))) "yes" else "NO", "\n")


# ============================================================================
# STEP 9 - msi_dataset -> Cardinal (back-conversion)
# ============================================================================
cat(crayon::cyan("\n====== STEP 9: msi_dataset -> Cardinal ======\n"))

card_back <- convert_msi_dataset2cardinal(my_msi)
cat("Cardinal class   :", class(card_back)[1], "\n")
cat("Cardinal features:", nrow(Cardinal::featureData(card_back)), "\n")
cat("Cardinal pixels  :", nrow(Cardinal::pixelData(card_back)), "\n")


# ============================================================================
# STEP 10 - Update metadata (resolution, matrix, tissue type, organism, date)
# ============================================================================
cat(crayon::cyan("\n====== STEP 10: Update metadata ======\n"))

cat("Before — sample_info columns:",
    paste(colnames(sample_info(my_msi)), collapse = ", "), "\n")

my_msi@sample_info$resolution  <- "50um"
my_msi@sample_info$matrix      <- "NEDC"
my_msi@sample_info$tissue_type <- "renal_tumor"
my_msi@sample_info$organism    <- "human"
my_msi@sample_info$date        <- "2026-02-10"

new_note_rows <- data.frame(
  name    = c("tissue_type", "organism", "date"),
  meaning = c("Type of tissue", "Source organism", "Collection date"),
  stringsAsFactors = FALSE
)
my_msi@sample_info_note <- rbind(my_msi@sample_info_note, new_note_rows)

my_msi <- update_msi_dataset(my_msi)

cat("After  — sample_info columns:",
    paste(colnames(sample_info(my_msi)), collapse = ", "), "\n")
print(my_msi)

saveRDS(my_msi, "data/output/my_msi_dataset_updated.rds")
cat("Saved updated object -> data/output/my_msi_dataset_updated.rds\n")


# ============================================================================
# STEP 11 - Save & reload final check
# ============================================================================
cat(crayon::cyan("\n====== STEP 11: Save & reload ======\n"))

tmp <- tempfile(fileext = ".rds")
saveRDS(my_msi, tmp)
my_msi_reload <- readRDS(tmp)
cat("Reloaded dim:", paste(dim(my_msi_reload), collapse = " x "), "\n")
check_reload <- check_msi_dataset_class(my_msi_reload)
cat("Valid       :",
    if (isTRUE(check_reload)) "yes" else paste(check_reload, collapse = "; "), "\n")
unlink(tmp)

cat(crayon::green("\n====== All steps completed successfully! ======\n"))