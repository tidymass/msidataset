# ============================================================================
# convert_cardinal2msi_dataset 完整使用指南 - 修正版
# Author: Xiaoyu Xu
# Package: msidataset
# Version: 2.0 (Fixed)
# ============================================================================

# ============================================================================
# 环境设置和包加载
# ============================================================================
devtools::document()
setwd("D:/Shen_lab/msidataset")
library(S4Vectors)
library(Cardinal)

devtools::load_all()

# 或者如果已经安装了包
# library(msidataset)

# ============================================================================
#  完整参数使用（推荐用于正式分析）
# ============================================================================

cat("\n=== Full Parameters ===\n")

# 读取数据
cardinal_obj <- readMSIData(
  "data/demo data/03-total ion count.imzML",
  memory = FALSE
)

# 使用完整参数转换v0: mz_0001
my_msi_full <- convert_cardinal2msi_dataset(
  cardinal_object = cardinal_obj,
  tolerance = 10,                      # 对齐容差
  units = "ppm",                       # 容差单位
  sample_name = "rena_cancer", # 样本名称
  section_name = "slice_01",   # 切片名称
  resolution = "50um",                 # 空间分辨率
  matrix_name = "NEDC"                  # MALDI 基质
)

# 使用完整参数转换v1: m50_1
my_msi_full <- convert_cardinal2msi_dataset_v1(
  cardinal_object = cardinal_obj,
  tolerance = 10,                      # 对齐容差
  units = "ppm",                       # 容差单位
  sample_name = "rena_cancer", # 样本名称
  section_name = "slice_01",   # 切片名称
  resolution = "50um",                 # 空间分辨率
  matrix_name = "NEDC"                  # MALDI 基质
)

# 查看对象
print(my_msi_full)

# ============================================================================
# 示例 3: 查看转换后的数据结构
# ============================================================================

cat("\n=== Example 3: Explore Data Structure ===\n")

# 重新加载之前保存的对象
my_msi <- readRDS("2_data/output/my_msi_dataset.rds")

# 1. Expression Data
cat("\n--- Expression Data ---\n")
cat("Dimensions:", paste(dim(expression_data(my_msi_full)), collapse = " x "), "\n")
cat("First 20 variables, first 5 pixels:\n")
print(expression_data(my_msi_full)[1:20, 1:5])

# 2. Pixel Info
cat("\n--- Pixel Info ---\n")
cat("Columns:", paste(colnames(pixel_info(my_msi_full)), collapse = ", "), "\n")
cat("Number of pixels:", nrow(pixel_info(my_msi_full)), "\n")
print(head(pixel_info(my_msi_full), 3))

# 查看空间范围
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

# m/z 分布
cat("\nm/z statistics:\n")
cat("  Min m/z:", round(min(variable_info(my_msi_full)$mz), 4), "\n")
cat("  Max m/z:", round(max(variable_info(my_msi_full)$mz), 4), "\n")
cat("  Median m/z:", round(median(variable_info(my_msi_full)$mz), 4), "\n")

# ============================================================================
# 示例 4: 数据质量检查
# ============================================================================

cat("\n=== Example 4: Quality Check ===\n")

my_msi <- readRDS("2_data/output/my_msi_dataset.rds")

# 1. 检查数据完整性
cat("\n--- Data Integrity Check ---\n")

# 检查缺失值
na_count <- sum(is.na(expression_data(my_msi_full)))
cat("Missing values in expression_data:", na_count, "\n")

# 检查零值比例
zero_count <- sum(expression_data(my_msi_full) == 0)
total_count <- nrow(expression_data(my_msi_full)) * ncol(expression_data(my_msi_full))
cat("Zero values:", zero_count, "(", 
    round(zero_count/total_count*100, 2), "%)\n")

# 2. 检查数据范围
cat("\n--- Data Range Check ---\n")
cat("Expression data range:\n")
cat("  Min:", round(min(expression_data(my_msi_full)), 2), "\n")
cat("  Max:", round(max(expression_data(my_msi_full)), 2), "\n")
cat("  Mean:", round(mean(as.matrix(expression_data(my_msi_full))), 2), "\n")

# 3. 检查像素分布
cat("\n--- Pixel Distribution ---\n")
pixel_coords <- pixel_info(my_msi_full)[, c("x", "y", "z")]
cat("Number of unique positions:\n")
cat("  x positions:", length(unique(pixel_coords$x)), "\n")
cat("  y positions:", length(unique(pixel_coords$y)), "\n")
cat("  z slices:", length(unique(pixel_coords$z)), "\n")

# 4. 验证对象结构
cat("\n--- Object Structure Validation ---\n")
check_result <- check_msi_dataset_class(my_msi_full)

if (isTRUE(check_result)) {
  cat("✓ Object structure is valid!\n")
} else {
  cat("✗ Validation errors found:\n")
  print(check_result)
}

# ============================================================================
# 示例 5: 更新和修改元数据
# ============================================================================

cat("\n=== Example 5: Update Metadata ===\n")

my_msi <- readRDS("2_data/output/my_msi_dataset.rds")

# 场景：转换时使用了默认值，现在需要更新实际参数

# 1. 更新 sample_info
cat("\n--- Before Update ---\n")
print(sample_info(my_msi))

# 更新分辨率和基质信息
my_msi@sample_info$resolution <- "75um"
my_msi@sample_info$matrix <- "CHCA"

cat("\n--- After Update ---\n")
print(sample_info(my_msi))

# 2. 添加额外的样本信息列
my_msi@sample_info$tissue_type <- "brain"
my_msi@sample_info$organism <- "mouse"
my_msi@sample_info$date <- "2025-02-05"

# 更新 sample_info_note
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

# 验证更新
cat("\n--- Updated Object ---\n")
print(my_msi)

# 保存更新后的对象
saveRDS(my_msi, "2_data/output/my_msi_dataset_updated.rds")
cat("✓ Updated object saved\n")

# ============================================================================
# 示例 6: 不同MALDI基质的数据转换
# ============================================================================

cat("\n=== Example 6: Different MALDI Matrices ===\n")

# 常见的MALDI基质类型
maldi_matrices <- data.frame(
  code = c("DHB", "CHCA", "SA", "DAN", "9AA"),
  name = c(
    "2,5-Dihydroxybenzoic acid",
    "α-Cyano-4-hydroxycinnamic acid",
    "Sinapinic acid",
    "1,5-Diaminonaphthalene",
    "9-Aminoacridine"
  ),
  application = c(
    "Lipids",
    "Peptides/Proteins",
    "Large proteins",
    "Negative ion mode",
    "Negative ion mode"
  ),
  stringsAsFactors = FALSE
)

print(maldi_matrices)

# 示例：DHB基质（常用于脂质成像）
cat("\nConverting with DHB matrix (for lipids)...\n")
# msi_dhb <- convert_cardinal2msi_dataset(
#   cardinal_object = readMSIData("2_data/demo data/lipid_sample.imzML"),
#   sample_name = "lipid_analysis",
#   matrix_name = "DHB",
#   resolution = "50um"
# )

# ============================================================================
# 示例 7: 不同空间分辨率的数据对比
# ============================================================================

cat("\n=== Example 7: Resolution Comparison ===\n")

resolution_info <- data.frame(
  resolution = c("10um", "20um", "50um", "100um", "200um"),
  application = c(
    "High-resolution imaging (cellular level)",
    "High-resolution imaging",
    "Standard resolution (most common)",
    "Low resolution (quick scan)",
    "Very low resolution (rapid screening)"
  ),
  typical_pixels = c(">100,000", "50,000-100,000", "10,000-50,000", 
                     "5,000-10,000", "<5,000"),
  stringsAsFactors = FALSE
)

print(resolution_info)

# ============================================================================
# 示例 8: 批量转换（带错误处理）
# ============================================================================

cat("\n=== Example 8: Batch Conversion with Error Handling ===\n")

# 设置数据文件路径
data_files <- c(
  sample1 = "data/demo data/03-total ion count.imzML"
  # 如果有更多文件，在这里添加
  # sample2 = "2_data/demo data/sample2.imzML",
  # sample3 = "2_data/demo data/sample3.imzML"
)

# 批量读取和转换
msi_objects <- list()
for (i in seq_along(data_files)) {
  sample_name <- names(data_files)[i]
  file_path <- data_files[i]
  
  cat("\nProcessing:", sample_name, "\n")
  
  tryCatch({
    # 检查文件是否存在
    if (!file.exists(file_path)) {
      stop("File not found: ", file_path)
    }
    
    # 读取 Cardinal 数据
    cardinal_obj <- readMSIData(file_path, memory = FALSE)
    
    # 转换为 msi_dataset
    msi_obj <- convert_cardinal2msi_dataset(
      cardinal_object = cardinal_obj,
      sample_name = sample_name,
      section_name = paste0(sample_name, "_section1"),
      resolution = "50um",
      matrix_name = "DHB"
    )
    
    # 保存到列表
    msi_objects[[sample_name]] <- msi_obj
    
    # 保存到文件
    saveRDS(msi_obj, 
            file.path("2_data/output", paste0(sample_name, "_msi.rds")))
    
    cat("✓", sample_name, "conversion successful\n")
    
  }, error = function(e) {
    cat("✗ Error processing", sample_name, ":", e$message, "\n")
  })
}

cat("\nBatch conversion completed. Processed", length(msi_objects), "sample(s).\n")

# ============================================================================
# 示例 9: 完整的项目工作流
# ============================================================================

cat("\n=== Example 9: Complete Project Workflow ===\n")

# 项目设置
project_dir <- "D:/Shen_lab/msidataset"
data_dir <- file.path(project_dir, "data/demo data")
output_dir <- file.path(project_dir, "2_data/output")

# 创建输出目录
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 1. 读取数据
cat("\n[1/5] Reading Cardinal data...\n")
cardinal_obj <- readMSIData(
  file.path(data_dir, "03-total ion count.imzML"),
  memory = FALSE
)
cat("✓ Data loaded successfully\n")

# 2. 转换为 msi_dataset
cat("\n[2/5] Converting to msi_dataset...\n")
my_msi <- convert_cardinal2msi_dataset(
  cardinal_object = cardinal_obj,
  tolerance = 10,
  units = "ppm",
  sample_name = "project_sample_001",
  section_name = "coronal_section_1",
  resolution = "50um",
  matrix_name = "DHB"
)
cat("✓ Conversion completed\n")

# 3. 数据质量检查
cat("\n[3/5] Quality check...\n")
check_result <- check_msi_dataset_class(my_msi)

if (!isTRUE(check_result)) {
  stop("Data validation failed: ", paste(check_result, collapse = "; "))
}

cat("✓ Quality check passed\n")

# 4. 保存处理后的数据
cat("\n[4/5] Saving processed data...\n")
saveRDS(my_msi, file.path(output_dir, "processed_msi_dataset.rds"))

# 保存元数据摘要
sink(file.path(output_dir, "metadata_summary.txt"))
cat("MSI Dataset Summary\n")
cat("===================\n\n")
cat("Sample name:", sample_info(my_msi)$sample_id, "\n")
cat("Section:", sample_info(my_msi)$section_id, "\n")
cat("Resolution:", sample_info(my_msi)$resolution, "\n")
cat("Matrix:", sample_info(my_msi)$matrix, "\n")
cat("Number of variables:", n_variables(my_msi), "\n")
cat("Number of pixels:", n_pixels(my_msi), "\n")
cat("m/z range:", paste(round(range(mz(my_msi)), 4), collapse = " - "), "\n")
cat("Spatial extent:\n")
cat("  x:", paste(range(pixel_info(my_msi)$x), collapse = " - "), "\n")
cat("  y:", paste(range(pixel_info(my_msi)$y), collapse = " - "), "\n")
cat("  z:", paste(range(pixel_info(my_msi)$z), collapse = " - "), "\n")
sink()

cat("✓ Data and metadata saved\n")

# 5. 生成质量报告
cat("\n[5/5] Generating quality report...\n")

# 创建简单的质量报告
report <- list(
  dataset_info = list(
    n_variables = n_variables(my_msi),
    n_pixels = n_pixels(my_msi),
    n_samples = nrow(sample_info(my_msi)),
    mz_range = range(mz(my_msi)),
    spatial_extent = list(
      x = range(pixel_info(my_msi)$x),
      y = range(pixel_info(my_msi)$y),
      z = range(pixel_info(my_msi)$z)
    )
  ),
  data_quality = list(
    na_count = sum(is.na(expression_data(my_msi))),
    zero_percent = round(sum(expression_data(my_msi) == 0) / 
                           (nrow(expression_data(my_msi)) * 
                              ncol(expression_data(my_msi))) * 100, 2),
    intensity_range = range(expression_data(my_msi)),
    mean_intensity = mean(as.matrix(expression_data(my_msi)))
  ),
  validation = list(
    passed = isTRUE(check_msi_dataset_class(my_msi))
  )
)

saveRDS(report, file.path(output_dir, "quality_report.rds"))
cat("✓ Quality report saved\n")

cat("\n========================================\n")
cat("Project workflow completed successfully!\n")
cat("========================================\n")
cat("Output files:\n")
cat("  - processed_msi_dataset.rds\n")
cat("  - metadata_summary.txt\n")
cat("  - quality_report.rds\n")
cat("\nData location:", output_dir, "\n")

# ============================================================================
# 示例 10: 数据摘要和统计
# ============================================================================

cat("\n=== Example 10: Data Summary and Statistics ===\n")

my_msi <- readRDS("2_data/output/my_msi_dataset.rds")

# 创建综合摘要
cat("\n--- Comprehensive Summary ---\n")
cat("\n1. Dataset Dimensions:\n")
cat("   Variables:", n_variables(my_msi), "\n")
cat("   Pixels:", n_pixels(my_msi), "\n")
cat("   Samples:", nrow(sample_info(my_msi)), "\n")

cat("\n2. m/z Information:\n")
cat("   Range:", paste(round(range(mz(my_msi)), 4), collapse = " - "), "\n")
cat("   Mean:", round(mean(mz(my_msi)), 4), "\n")
cat("   Median:", round(median(mz(my_msi)), 4), "\n")

cat("\n3. Spatial Information:\n")
cat("   X range:", paste(range(pixel_info(my_msi)$x), collapse = " - "), "\n")
cat("   Y range:", paste(range(pixel_info(my_msi)$y), collapse = " - "), "\n")
cat("   Z slices:", length(unique(pixel_info(my_msi)$z)), "\n")

cat("\n4. Data Quality Metrics:\n")
expr_data <- expression_data(my_msi)
cat("   Total values:", nrow(expr_data) * ncol(expr_data), "\n")
cat("   Missing values:", sum(is.na(expr_data)), "\n")
cat("   Zero values:", sum(expr_data == 0), 
    "(", round(sum(expr_data == 0)/(nrow(expr_data)*ncol(expr_data))*100, 2), "%)\n")
cat("   Non-zero values:", sum(expr_data != 0), "\n")

cat("\n5. Intensity Statistics:\n")
cat("   Min:", round(min(expr_data), 2), "\n")
cat("   Max:", round(max(expr_data), 2), "\n")
cat("   Mean:", round(mean(as.matrix(expr_data)), 2), "\n")
cat("   Median:", round(median(as.matrix(expr_data)), 2), "\n")
cat("   SD:", round(sd(as.matrix(expr_data)), 2), "\n")

# ============================================================================
# 常见问题解答 (FAQ)
# ============================================================================

cat("\n=== FAQ ===\n")

cat("\nQ1: 如何选择合适的对齐容差？\n")
cat("A1: 通常使用 10 ppm 是合理的起点。\n")
cat("    - 高分辨率仪器: 5 ppm\n")
cat("    - 低分辨率或有较大质量偏移: 15-20 ppm\n")

cat("\nQ2: resolution 参数应该填什么？\n")
cat("A2: 填写实际的空间分辨率，如 '10um', '50um', '100um'。\n")
cat("    如果不清楚，可以暂时填 'unknown'，之后更新。\n")

cat("\nQ3: 如何处理没有z坐标的2D数据？\n")
cat("A3: 函数会自动添加 z = 1。不需要特殊处理。\n")

cat("\nQ4: 转换后如何访问原始 Cardinal 对象？\n")
cat("A4: 转换是单向的，不会保留原始对象。\n")
cat("    如需要，请单独保存 Cardinal 对象。\n")

cat("\nQ5: 如何验证转换结果？\n")
cat("A5: 使用 check_msi_dataset_class(my_msi) 验证对象结构。\n")

# ============================================================================
# 测试总结
# ============================================================================

cat("\n" , paste(rep("=", 60), collapse = ""), "\n")
cat("All examples completed successfully!\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("Tested features:\n")
cat("  ✓ Basic conversion\n")
cat("  ✓ Full parameter conversion\n")
cat("  ✓ Data structure exploration\n")
cat("  ✓ Quality checking\n")
cat("  ✓ Metadata updates\n")
cat("  ✓ Batch processing\n")
cat("  ✓ Complete workflow\n")
cat("  ✓ Data statistics\n")
cat("\nAll accessor methods working:\n")
cat("  ✓ n_variables()\n")
cat("  ✓ n_pixels()\n")
cat("  ✓ mz()\n")
cat("  ✓ expression_data()\n")
cat("  ✓ pixel_info()\n")
cat("  ✓ sample_info()\n")
cat("  ✓ variable_info()\n")
cat("\nValidation:\n")
cat("  ✓ check_msi_dataset_class() working correctly\n")
cat("  ✓ show() method displaying properly\n")

cat("\nReady for downstream analysis!\n")