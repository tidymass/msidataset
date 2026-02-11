# ============================================================================
# Package: msidataset
# Author: Xiaoyu Xu
# Email: XIAOYU021@e.ntu.edu.sg
# Version: 1.2.0
# step 1: convert cardinal to msidataset & checkout
# step 2: data preprocess & analysis & visualization (msidataset)
# step 3: convert msidataset to massdataset
# step 4: ggplot2 & annotation and others)
# step 5: convert massdataset to msidataset
# ============================================================================

devtools::document()
setwd("D:/Shen_lab/msidataset")
devtools::load_all()

library(S4Vectors)
library(Cardinal)
library(massdataset)


# library(msidataset)

#=============================================================================
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

#=============================================================================
# step 2: data preprocess & analysis & visualization==========================
#=============================================================================

# 2.1 ion spatial distribution
# 2.1.1 Plot first variable
plot_ion_spatial(
  my_msi, 
  variable_id = "m50_4",  
  color_scale = "viridis", 
  transform = "none"
)

# 2.1.2 Plot by m/z value
var_mz <- 50.3
plot_ion_spatial(
  my_msi,
  mz = var_mz,
  mz_tolerance = 0.05,
  color_scale = "magma",
  transform = "sqrt"
)

# 2.1.3 Select several interesting variables; ncol
real_ids <- rownames(expression_data(my_msi))[c(10, 15, 25, 35, 45, 55)]
print(real_ids)

plot_ion_grid(
  my_msi,
  variable_ids = real_ids,
  ncol = 2,
  color_scale = "plasma",
  transform = "log10"
)

# 2.1.4 Multi-ion spatial overlay visualization
real_ids <- rownames(expression_data(my_msi))[c(100, 200, 300, 400, 500, 600, 700)]
 
# auto color
plot_multi_ion_overlay(
  my_msi,
  variable_ids = real_ids,
  title = "Ion Spatial Co-localization"
)

# user-defined color
plot_multi_ion_overlay(
  my_msi,
  variable_ids = c("m50_5", "m50_15"),
  colors = c("cyan", "magenta"),
  title = "Duo-Ion Contrast"
)

# 2.2 Quality control - TIC spatial map
plot_tic_spatial(
  my_msi,
  transform = "log10",
  color_scale = "viridis"
)

# 2.3 Peak count spatial distribution
plot_peak_count_spatial(
  my_msi,
  threshold = 3000
)

# 2.4 Spatial clustering
# Perform k-means clustering
msi_clustered <- spatial_clustering(
  my_msi,
  n_clusters = 5,
  method = "kmeans",
  plot = TRUE
)

# Access cluster assignments
table(pixel_info(msi_clustered)$cluster)

# 2.5 Spatial PCA
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


#=============================================================================
# step 3: convert msidataset to massdataset===================================
#=============================================================================

# choice 1: basic
my_mass <- convert_msi2mass_dataset(
  msi_object = my_msi_full,
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

#=============================================================================
# step 4: using massdataset func==============================================
#=============================================================================
# 4.1 intensity_plot
library(massdataset)

# choice 1: variable_index
intensity_plot(
  object = my_mass,
  variable_index = 1  # 绘制第1个变量
)

# choice 2: variable_id
intensity_plot(
  object = my_mass,
  variable_id = "m50_100"  # 使用具体的变量ID
)

# 4.2 ggplot_mass_dataset
library(ggplot2)
library(dplyr)
library(patchwork)

# 4.2.1 add group
my_mass@sample_info$x_region <- cut(
  my_mass@sample_info$x,
  breaks = 3,
  labels = c("Left", "Center", "Right")
)

my_mass@sample_info$y_region <- cut(
  my_mass@sample_info$y,
  breaks = 3,
  labels = c("Bottom", "Middle", "Top")
)

# 4.2.2 calculate
mean_int <- rowMeans(my_mass@expression_data)
cv <- apply(my_mass@expression_data, 1, sd) / mean_int

# high intensity
high_int_vars <- order(mean_int, decreasing = TRUE)[1:10]

# high cv
high_cv_vars <- order(cv, decreasing = TRUE)[1:10]


# 4.2.3 Highest Intensity Variable + Spatial Distribution
p1 <- ggplot_mass_dataset(
  my_mass,
  direction = "variable",
  variable_index = high_int_vars[1]
) +
  geom_boxplot(aes(x = x_region, fill = x_region)) +
  labs(title = "Highest Intensity Variable")

p2 <- ggplot_mass_dataset(
  my_mass,
  direction = "variable",
  variable_index = high_int_vars[1]
) +
  geom_point(aes(x = x, y = y, size = value, color = value)) +
  scale_color_viridis_c() +
  coord_equal() +
  labs(title = "Spatial Distribution")

p1 + p2

# 4.2.4 High Variability Variable Across Regions
ggplot_mass_dataset(
  my_mass,
  direction = "variable",
  variable_index = high_cv_vars[1]
) +
  geom_violin(aes(x = x_region, fill = x_region)) +
  geom_jitter(aes(x = x_region), width = 0.1, alpha = 0.3) +
  facet_wrap(~y_region) +
  labs(title = "High Variability Variable Across Regions")

# 4.2.5 
ggplot_mass_dataset(
  my_mass,
  direction = "sample",
  sample_index = 1
) +
  geom_line(aes(x = mz), color = "steelblue", size = 0.3) +
  geom_point(aes(x = mz), color = "red", size = 0.5, alpha = 0.3) +
  labs(
    title = "Mass Spectrum - Sample 1",
    x = "m/z",
    y = "Intensity"
  ) +
  theme_minimal()

# 4.2.6
samples_to_compare <- c(1, 10, 20, 30)

plot_list <- lapply(samples_to_compare, function(idx) {
  ggplot_mass_dataset(
    my_mass,
    direction = "sample",
    sample_index = idx
  ) +
    geom_line(aes(x = mz), color = "darkblue", size = 0.3) +
    labs(
      title = paste("Sample", idx),
      x = "m/z",
      y = "Intensity"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 8),
      plot.title = element_text(size = 10)
    )
})

wrap_plots(plot_list, ncol = 2)


# 4.4 Plot m/z vs Spatial Correlation for MSI Data

# 4.4.1 top 20 
result <- plot_mz_spatial_correlation(my_mass)
print(result)

summary(result)

correlation_data <- result$data
head(correlation_data)


# 4.4.2

# top 100
result_100 <- plot_mz_spatial_correlation(
  my_mass, 
  n_vars = 100
)

summary(result_100)


# 4.4.3 Spearman
result_spearman <- plot_mz_spatial_correlation(
  my_mass,
  n_vars = 50,
  method = "spearman"
)

print(result_spearman)

# 4.4.4 add m/z 
result_labeled <- plot_mz_spatial_correlation(
  my_mass,
  n_vars = 50,
  show_labels = TRUE,      # 显示标签
  label_top_n = 10         # 标记前10个相关性最强的
)

print(result_labeled)

# 4.4.5 Intensity
# more than 3,000
result_filtered <- plot_mz_spatial_correlation(
  my_mass,
  n_vars = 100,
  intensity_threshold = 3000
)


# ============================================
# 用例 6: 提取和分析强相关的变量
# ============================================

# 找出与 X 坐标强相关的变量 (|r| > 0.7)
strong_x_cor <- result$data %>%
  dplyr::filter(abs(cor_x) > 0.7 & sig_x == "Significant") %>%
  dplyr::arrange(desc(abs(cor_x)))

cat("Found", nrow(strong_x_cor), "variables strongly correlated with X\n")
print(strong_x_cor)

# 找出与 Y 坐标强相关的变量
strong_y_cor <- result$data %>%
  dplyr::filter(abs(cor_y) > 0.7 & sig_y == "Significant") %>%
  dplyr::arrange(desc(abs(cor_y)))

cat("Found", nrow(strong_y_cor), "variables strongly correlated with Y\n")
print(strong_y_cor)


# ============================================
# 用例 7: 可视化强相关变量的空间分布
# ============================================

# 结合之前的空间分析函数
source("msi_spatial_analysis.R")

# 如果有强相关的变量,绘制它们的空间分布
if (nrow(strong_x_cor) > 0) {
  
  # 获取最强相关的变量的索引
  top_var_id <- strong_x_cor$variable_id[1]
  top_var_idx <- which(my_mass@variable_info$variable_id == top_var_id)
  
  # 绘制空间分布
  plot_ion_spatial(
    convert_mass2msi_dataset(my_mass),  # 先转回 msi_dataset
    variable_index = top_var_idx,
    title = paste0("Spatial Distribution: m/z ", 
                   round(strong_x_cor$mz[1], 2),
                   " (r_x = ", round(strong_x_cor$cor_x[1], 3), ")")
  )
}


# ============================================
# 用例 8: 比较不同相关性方法
# ============================================

library(patchwork)

# Pearson 相关
res_pearson <- plot_mz_spatial_correlation(
  my_mass, 
  n_vars = 50, 
  method = "pearson"
)

# Spearman 相关
res_spearman <- plot_mz_spatial_correlation(
  my_mass, 
  n_vars = 50, 
  method = "spearman"
)

# 并排比较(只显示 X 相关性)
# 注意:需要从 result 对象中提取单个图
# (这里简化展示,实际可能需要重新绘制)


# ============================================
# 用例 9: 导出相关性数据
# ============================================

# 分析所有高强度变量
result_all <- plot_mz_spatial_correlation(
  my_mass,
  n_vars = 200,
  method = "pearson"
)

# 导出相关性数据到 CSV
write.csv(
  result_all$data,
  "mz_spatial_correlation_results.csv",
  row.names = FALSE
)

# 只导出显著相关的变量
sig_vars <- result_all$data %>%
  dplyr::filter(sig_x == "Significant" | sig_y == "Significant")

write.csv(
  sig_vars,
  "significant_spatial_correlations.csv",
  row.names = FALSE
)


# ============================================
# 用例 10: 创建相关性热图
# ============================================

library(ggplot2)

# 提取数据
cor_data <- result$data

# 创建散点图:X vs Y 相关性
ggplot(cor_data, aes(x = cor_x, y = cor_y)) +
  geom_point(aes(color = mean_intensity, size = total_cor), alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_viridis_c(name = "Mean\nIntensity") +
  scale_size_continuous(name = "Total\nCorrelation") +
  labs(
    title = "X vs Y Correlation Pattern",
    x = "Correlation with X",
    y = "Correlation with Y"
  ) +
  theme_bw() +
  coord_equal()


# ============================================
# 用例 11: 分区域分析
# ============================================

# 创建一个函数分析特定 m/z 范围
analyze_mz_range <- function(mass_obj, mz_min, mz_max, n_vars = 50) {
  
  # 筛选特定 m/z 范围的变量
  mz_values <- mass_obj@variable_info$mz
  in_range <- which(mz_values >= mz_min & mz_values <= mz_max)
  
  if (length(in_range) == 0) {
    stop("No variables in the specified m/z range")
  }
  
  # 创建临时对象,只包含该范围的变量
  temp_obj <- mass_obj
  temp_obj@expression_data <- mass_obj@expression_data[in_range, ]
  temp_obj@variable_info <- mass_obj@variable_info[in_range, ]
  
  # 分析
  result <- plot_mz_spatial_correlation(
    temp_obj,
    n_vars = min(n_vars, length(in_range))
  )
  
  return(result)
}

# 分析脂质范围 (m/z 700-900)
lipid_result <- analyze_mz_range(my_mass, 700, 900, n_vars = 50)
print(lipid_result)

# 分析代谢物范围 (m/z 100-500)
metabolite_result <- analyze_mz_range(my_mass, 100, 500, n_vars = 50)
print(metabolite_result)


# ============================================
# 用例 12: 完整分析流程
# ============================================

# 完整的空间相关性分析工作流

# 1. 全局分析
cat("Step 1: Global correlation analysis\n")
global_result <- plot_mz_spatial_correlation(
  my_mass,
  n_vars = 100,
  method = "pearson",
  show_labels = FALSE
)

print(global_result)
summary(global_result)

# 2. 提取显著变量
cat("\nStep 2: Extract significant variables\n")
sig_data <- global_result$data %>%
  dplyr::filter(sig_x == "Significant" | sig_y == "Significant")

cat("Found", nrow(sig_data), "significant variables\n")

# 3. 分类显著变量
cat("\nStep 3: Classify variables by correlation pattern\n")

# X方向正相关
x_positive <- sig_data %>% 
  dplyr::filter(cor_x > 0.5 & sig_x == "Significant")
cat("X-positive correlations:", nrow(x_positive), "\n")

# X方向负相关
x_negative <- sig_data %>% 
  dplyr::filter(cor_x < -0.5 & sig_x == "Significant")
cat("X-negative correlations:", nrow(x_negative), "\n")

# Y方向正相关
y_positive <- sig_data %>% 
  dplyr::filter(cor_y > 0.5 & sig_y == "Significant")
cat("Y-positive correlations:", nrow(y_positive), "\n")

# Y方向负相关
y_negative <- sig_data %>% 
  dplyr::filter(cor_y < -0.5 & sig_y == "Significant")
cat("Y-negative correlations:", nrow(y_negative), "\n")

# 4. 保存结果
cat("\nStep 4: Save results\n")

# 创建输出目录
if (!dir.exists("correlation_analysis")) {
  dir.create("correlation_analysis")
}

# 保存图形
ggsave(
  "correlation_analysis/mz_spatial_correlation.pdf",
  plot = global_result$plot,
  width = 12,
  height = 10
)

# 保存数据
write.csv(
  global_result$data,
  "correlation_analysis/correlation_data.csv",
  row.names = FALSE
)

# 保存分类结果
write.csv(x_positive, "correlation_analysis/x_positive.csv", row.names = FALSE)
write.csv(x_negative, "correlation_analysis/x_negative.csv", row.names = FALSE)
write.csv(y_positive, "correlation_analysis/y_positive.csv", row.names = FALSE)
write.csv(y_negative, "correlation_analysis/y_negative.csv", row.names = FALSE)

cat("\nAnalysis complete! Results saved to 'correlation_analysis/'\n")















# 4.5 annotation?







#=============================================================================
# step 5: convert massdataset to msidataset===================================
#=============================================================================
# 方法 A: 自动检测元数据（如果在同一 R session）
my_msi_restored <- convert_mass2msi_dataset(
  mass_object = my_mass_annotated,
  metadata_source = "auto",  # 自动检测
  create_annotation_table = TRUE  # 创建 annotation_table
)

# 或方法 B: 从文件加载元数据（如果跨 session）
my_msi_restored <- convert_mass2msi_dataset(
  mass_object = my_mass_annotated,
  metadata_source = "file",
  metadata_path = "2_data/output/msi_metadata.rds",
  create_annotation_table = TRUE
)

# 查看恢复后的 msi_dataset
print(my_msi_restored)

cat("\nRestored MSI Dataset summary:\n")
cat("  Variables:", n_variables(my_msi_restored), "\n")
cat("  Pixels:", n_pixels(my_msi_restored), "\n")








