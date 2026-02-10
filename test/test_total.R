# ============================================================================
# Package: msidataset
# Author: Xiaoyu Xu
# Email: XIAOYU021@e.ntu.edu.sg
# Version: 1.2.0
# step 1: convert cardinal to msidataset
# step 2/3: convert msidataset to massdataset
# step 2/3: data preprocess & analysis & visualization 
# step 4: convert massdataset to msidataset
# ============================================================================

devtools::document()
setwd("D:/Shen_lab/msidataset")
devtools::load_all()

library(S4Vectors)
library(Cardinal)
library(massdataset)


# library(msidataset)

# step 1: convert cardinal to msidataset======================================
#=============================================================================

# read data
cardinal_obj <- readMSIData(
  "data/demo data/03-total ion count.imzML",
  memory = FALSE 
)


# mz style: m50_1
my_msi_full <- convert_cardinal2msi_dataset(
  cardinal_object = cardinal_obj,
  tolerance = 10,                      # ppm
  units = "ppm",                      
  sample_name = "renal_cancer", 
  section_name = "slice_01",   
  resolution = "50um",                 
  matrix_name = "NEDC"                  # MALDI matrix
)

# check data
print(my_msi_full)


cat("\n=== Example 3: Explore Data Structure ===\n")

# save data
saveRDS(my_msi_full, "data/output/my_msi_dataset_03.rds")

# reload data
my_msi <- readRDS("data/output/my_msi_dataset_03.rds")

# 1. Expression Data (20,5)
cat("\n--- Expression Data ---\n")
cat("Dimensions:", paste(dim(expression_data(my_msi_full)), collapse = " x "), "\n")
cat("First 20 variables, first 5 pixels:\n")
print(expression_data(my_msi_full)[1:20, 1:5])

# 2. Pixel Info
cat("\n--- Pixel Info ---\n")
cat("Columns:", paste(colnames(pixel_info(my_msi_full)), collapse = ", "), "\n")
cat("Number of pixels:", nrow(pixel_info(my_msi_full)), "\n")
print(head(pixel_info(my_msi_full), 10))

# spatial range
cat("\nSpatial ranges:\n")
cat("  x: [", min(pixel_info(my_msi_full)$x), ",", max(pixel_info(my_msi_full)$x), "]\n")
cat("  y: [", min(pixel_info(my_msi_full)$y), ",", max(pixel_info(my_msi_full)$y), "]\n")
cat("  z: [", min(pixel_info(my_msi_full)$z), ",", max(pixel_info(my_msi_full)$z), "]\n")

# 3. Sample Info
cat("\n--- Sample Info ---\n")
cat("Columns:", paste(colnames(sample_info(my_msi_full)), collapse = ", "), "\n")
print(sample_info(my_msi_full))

# 4. Variable Info
cat("\n--- Variable Info ---\n")
cat("Columns:", paste(colnames(variable_info(my_msi_full)), collapse = ", "), "\n")
cat("Number of variables:", nrow(variable_info(my_msi_full)), "\n")
print(head(variable_info(my_msi_full), 20))

# m/z range
cat("\nm/z statistics:\n")
cat("  Min m/z:", round(min(variable_info(my_msi_full)$mz), 4), "\n")
cat("  Max m/z:", round(max(variable_info(my_msi_full)$mz), 4), "\n")
cat("  Median m/z:", round(median(variable_info(my_msi_full)$mz), 4), "\n")


cat("\n=== Example 4: Quality Check ===\n")

# 1. Data Integrity
cat("\n--- Data Integrity Check ---\n")

# Missing values
na_count <- sum(is.na(expression_data(my_msi_full)))
cat("Missing values in expression_data:", na_count, "\n")

# Zero values
zero_count <- sum(expression_data(my_msi_full) == 0)
total_count <- nrow(expression_data(my_msi_full)) * ncol(expression_data(my_msi_full))
cat("Zero values:", zero_count, "(", 
    round(zero_count/total_count*100, 2), "%)\n")

# 2. data range
cat("\n--- Data Range Check ---\n")
cat("Expression data range:\n")
cat("  Min:", round(min(expression_data(my_msi_full)), 2), "\n")
cat("  Max:", round(max(expression_data(my_msi_full)), 2), "\n")
cat("  Mean:", round(mean(as.matrix(expression_data(my_msi_full))), 2), "\n")

# 3. Pixel Distribution
cat("\n--- Pixel Distribution ---\n")
pixel_coords <- pixel_info(my_msi_full)[, c("x", "y", "z")]
cat("Number of unique positions:\n")
cat("  x positions:", length(unique(pixel_coords$x)), "\n")
cat("  y positions:", length(unique(pixel_coords$y)), "\n")
cat("  z slices:", length(unique(pixel_coords$z)), "\n")

# 4. Object Structure Validation
cat("\n--- Object Structure Validation ---\n")
check_result <- check_msi_dataset_class(my_msi_full)

if (isTRUE(check_result)) {
  cat("✓ Object structure is valid!\n")
} else {
  cat("✗ Validation errors found:\n")
  print(check_result)
}

# Update Metadata examples
# ============================================================================

cat("\n=== Example 5: Update Metadata ===\n")

# 1. update sample_info
cat("\n--- Before Update ---\n")
print(sample_info(my_msi))

# resolution & matrix
my_msi@sample_info$resolution <- "10um"
my_msi@sample_info$matrix <- "CHCA"

cat("\n--- After Update ---\n")
print(sample_info(my_msi))

# 2. additional sample info
my_msi@sample_info$tissue_type <- "renal tumor"
my_msi@sample_info$organism <- "human"
my_msi@sample_info$date <- "2026-02-10"

# update sample_info_note
new_notes <- data.frame(
  name = c("tissue_type", "organism", "date"),
  meaning = c("Type of tissue", "Source organism", "Collection date"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

my_msi@sample_info_note <- rbind(
  my_msi@sample_info_note,
  new_notes
)

# check update
cat("\n--- Updated Object ---\n")
print(my_msi)

# save update data
saveRDS(my_msi, "data/output/my_msi_dataset_updated.rds")
cat("✓ Updated object saved\n")

# step 2: convert msidataset to massdataset===================================
#=============================================================================

# choice 1: basic
my_mass <- convert_msi2mass_dataset(
  msi_object = my_msi,
  save_metadata = TRUE,  # 保存元数据到环境
  add_spatial_to_sample_info = FALSE  # 不添加空间坐标到 sample_info
)

# choice 2: Persist metadata to a file
my_mass <- convert_msi2mass_dataset(
  msi_object = my_msi_full,
  save_metadata = TRUE,
  metadata_path = "data/output/msi_metadata.rds",  # 保存到文件
  add_spatial_to_sample_info = TRUE  # 可选：添加 x,y,z 到 sample_info
)

# check structure of mass_dataset
print(my_mass)

cat("\nMass Dataset summary:\n")
cat("  Variables:", nrow(my_mass@variable_info), "\n")
cat("  Samples:", nrow(my_mass@sample_info), "\n")

# 查看 sample_info
cat("\nSample info preview:\n")
print(head(my_mass@sample_info, 10))

# step 3: data preprocess & analysis & visualization==========================
#=============================================================================

# 3.1 ion spatial distribution
# 3.1.1 Plot first variable
plot_ion_spatial(
  my_msi, 
  variable_id = "m50_4",  
  color_scale = "viridis", 
  transform = "none"
)

# 3.1.2 Plot by m/z value
var_mz <- 50.3
plot_ion_spatial(
  my_msi,
  mz = var_mz,
  mz_tolerance = 0.05,
  color_scale = "magma",
  transform = "sqrt"
)

# 3.1.3 Select several interesting variables; ncol
real_ids <- rownames(expression_data(my_msi))[c(10, 15, 25, 35, 45, 55)]
print(real_ids)

plot_ion_grid(
  my_msi,
  variable_ids = real_ids,
  ncol = 2,
  color_scale = "plasma",
  transform = "log10"
)

# 3.1.4 Multi-ion spatial overlay visualization
real_ids <- rownames(expression_data(my_msi))[c(100, 200, 300, 400, 500, 600, 700)]
 
# 3.1.4-1 运行叠加绘图 (自动分配颜色)
plot_multi_ion_overlay(
  my_msi,
  variable_ids = real_ids,
  title = "Ion Spatial Co-localization"
)

# 3.1.4-2 自定义颜色
plot_multi_ion_overlay(
  my_msi,
  variable_ids = c("m50_5", "m50_15"),
  colors = c("cyan", "magenta"),
  title = "Duo-Ion Contrast"
)

# 3.2 Quality control - TIC spatial map
plot_tic_spatial(
  my_msi,
  transform = "log10",
  color_scale = "viridis"
)

# 3.3 Peak count spatial distribution
plot_peak_count_spatial(
  my_msi,
  threshold = 3000
)

# 3.4 Spatial clustering
# Perform k-means clustering
msi_clustered <- spatial_clustering(
  my_msi,
  n_clusters = 5,
  method = "kmeans",
  plot = TRUE
)

# Access cluster assignments
table(pixel_info(msi_clustered)$cluster)

# 3.5 Spatial PCA
pca_results <- plot_spatial_pca(
  msi_clustered,
  n_components = 2,
  color_by = "cluster"
)

# Access PCA scores
pc_scores <- pca_results$scores
variance_explained <- pca_results$variance_explained

cat("Variance explained by PC1:", variance_explained[1], "%\n")
cat("Variance explained by PC2:", variance_explained[2], "%\n")





