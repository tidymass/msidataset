# ============================================================================
# convert_cardinal2msi_dataset 完整使用指南
# Author: Xiaoyu Xu
# Package: msidataset
# ============================================================================

# ============================================================================
# 环境设置和包加载
# ============================================================================
setwd("D:/Shen_lab/msidataset")
library(S4Vectors)
library(Cardinal)
devtools::load_all()

# 或者如果已经安装了包
# library(msidataset)

# ============================================================================
# 示例 1: 基础工作流（推荐用于入门）
# ============================================================================

# Step 1: 读取 imzML 文件（使用 Cardinal）
cardinal_obj <- readMSIData(
  "data/demo data/03-total ion count.imzML",
  memory = FALSE
)

# 查看 Cardinal 对象
print(cardinal_obj)

# Step 2: 转换为 msi_dataset 对象
my_msi <- convert_cardinal2msi_dataset(
  cardinal_object = cardinal_obj,
  sample_name = "demo_sample",
  section_name = "section_1"
)

# 查看转换后的对象
print(my_msi)

# Step 3: 验证数据
n_variables(my_msi)  # 变量数
n_pixels(my_msi)     # 像素数
range(mz(my_msi))    # m/z 范围

# Step 4: 保存对象
saveRDS(my_msi, "2_data/output/my_msi_dataset.rds")

# Step 5: 后续分析
# - 数据过滤
# - 归一化
# - 统计分析
# - 可视化等

# ============================================================================
# 示例 2: 完整参数使用（推荐用于正式分析）
# ============================================================================

# 读取数据
cardinal_obj <- readMSIData(
  "2_data/demo data/03-total ion count.imzML",
  memory = FALSE
)

# 使用完整参数转换
my_msi_full <- convert_cardinal2msi_dataset(
  cardinal_object = cardinal_obj,
  tolerance = 10,                      # 对齐容差
  units = "ppm",                       # 容差单位
  sample_name = "mouse_brain_section", # 样本名称
  section_name = "coronal_slice_01",   # 切片名称
  resolution = "50um",                 # 空间分辨率
  matrix_name = "DHB"                  # MALDI 基质
)

# 查看对象
print(my_msi_full)

# ============================================================================
# 示例 3: 查看转换后的数据结构
# ============================================================================

# 重新加载之前保存的对象
my_msi <- readRDS("2_data/output/my_msi_dataset.rds")

# 1. Expression Data
cat("\n=== Expression Data ===\n")
cat("Dimensions:", dim(expression_data(my_msi)), "\n")
cat("First 5 variables, first 3 pixels:\n")
print(expression_data(my_msi)[1:5, 1:3])

# 2. Pixel Info
cat("\n=== Pixel Info ===\n")
cat("Columns:", colnames(pixel_info(my_msi)), "\n")
cat("Number of pixels:", nrow(pixel_info(my_msi)), "\n")
print(head(pixel_info(my_msi)))

# 查看空间范围
cat("\nSpatial ranges:\n")
cat("  x: [", min(pixel_info(my_msi)$x), ",", max(pixel_info(my_msi)$x), "]\n")
cat("  y: [", min(pixel_info(my_msi)$y), ",", max(pixel_info(my_msi)$y), "]\n")
cat("  z: [", min(pixel_info(my_msi)$z), ",", max(pixel_info(my_msi)$z), "]\n")

# 3. Sample Info
cat("\n=== Sample Info ===\n")
cat("Columns:", colnames(sample_info(my_msi)), "\n")
print(sample_info(my_msi))

# 4. Variable Info
cat("\n=== Variable Info ===\n")
cat("Columns:", colnames(variable_info(my_msi)), "\n")
cat("Number of variables:", nrow(variable_info(my_msi)), "\n")
print(head(variable_info(my_msi)))

# m/z 分布
cat("\nm/z statistics:\n")
cat("  Min m/z:", min(variable_info(my_msi)$mz), "\n")
cat("  Max m/z:", max(variable_info(my_msi)$mz), "\n")
cat("  Median m/z:", median(variable_info(my_msi)$mz), "\n")

# ============================================================================
# 示例 4: 批量转换多个文件
# ============================================================================

# 设置数据文件路径
data_files <- list(
  sample1 = "2_data/demo data/sample1.imzML",
  sample2 = "2_data/demo data/sample2.imzML",
  sample3 = "2_data/demo data/sample3.imzML"
)

# 批量读取和转换
msi_objects <- lapply(names(data_files), function(sample_name) {
  
  cat("\nProcessing:", sample_name, "\n")
  
  # 读取 Cardinal 数据
  cardinal_obj <- readMSIData(data_files[[sample_name]], memory = FALSE)
  
  # 转换为 msi_dataset
  convert_cardinal2msi_dataset(
    cardinal_object = cardinal_obj,
    sample_name = sample_name,
    section_name = paste0(sample_name, "_section1"),
    resolution = "50um",
    matrix_name = "DHB"
  )
})

names(msi_objects) <- names(data_files)

# 查看所有对象
lapply(msi_objects, print)

# 批量保存
for (i in seq_along(msi_objects)) {
  saveRDS(
    msi_objects[[i]], 
    file.path("2_data/output", paste0(names(msi_objects)[i], "_msi.rds"))
  )
}

# ============================================================================
# 示例 5: 处理 MSImagingArrays（需要自动对齐）
# ============================================================================

# 如果读取的是未对齐的 Processed 数据
cardinal_arrays <- readMSIData(
  "2_data/demo data/processed_data.imzML",
  memory = FALSE
)

# 检查对象类型
cat("Object class:", class(cardinal_arrays), "\n")

# 转换（会自动进行峰对齐）
my_msi_aligned <- convert_cardinal2msi_dataset(
  cardinal_object = cardinal_arrays,
  tolerance = 10,              # 可以调整对齐容差
  units = "ppm",
  sample_name = "aligned_sample",
  section_name = "section_A"
)

# 查看对齐后的结果
print(my_msi_aligned)
cat("Aligned m/z values:", length(mz(my_msi_aligned)), "\n")

# ============================================================================
# 示例 6: 数据质量检查
# ============================================================================

my_msi <- readRDS("2_data/output/my_msi_dataset.rds")

# 1. 检查数据完整性
cat("\n=== Data Integrity Check ===\n")

# 检查缺失值
na_count <- sum(is.na(expression_data(my_msi)))
cat("Missing values in expression_data:", na_count, "\n")

# 检查零值比例
zero_count <- sum(expression_data(my_msi) == 0)
total_count <- nrow(expression_data(my_msi)) * ncol(expression_data(my_msi))
cat("Zero values:", zero_count, "(", 
    round(zero_count/total_count*100, 2), "%)\n")

# 2. 检查数据范围
cat("\n=== Data Range Check ===\n")
cat("Expression data range:\n")
cat("  Min:", min(expression_data(my_msi)), "\n")
cat("  Max:", max(expression_data(my_msi)), "\n")
cat("  Mean:", mean(as.matrix(expression_data(my_msi))), "\n")

# 3. 检查像素分布
cat("\n=== Pixel Distribution ===\n")
pixel_coords <- pixel_info(my_msi)[, c("x", "y", "z")]
cat("Number of unique positions:\n")
cat("  x positions:", length(unique(pixel_coords$x)), "\n")
cat("  y positions:", length(unique(pixel_coords$y)), "\n")
cat("  z slices:", length(unique(pixel_coords$z)), "\n")

# 4. 验证对象结构
cat("\n=== Object Structure Validation ===\n")
check_result <- check_msi_dataset(
  expression_data = my_msi@expression_data,
  pixel_info = my_msi@pixel_info,
  sample_info = my_msi@sample_info,
  variable_info = my_msi@variable_info
)

if (check_result == TRUE) {
  cat("✓ Object structure is valid!\n")
} else {
  cat("✗ Validation errors found:\n")
  print(check_result)
}

# ============================================================================
# 示例 7: 处理3D数据（多z切片）
# ============================================================================

# 读取包含多个z切片的3D数据
cardinal_3d <- readMSIData(
  "2_data/demo data/3d_volume.imzML",
  memory = FALSE
)

# 转换3D数据
my_msi_3d <- convert_cardinal2msi_dataset(
  cardinal_object = cardinal_3d,
  sample_name = "3d_tissue_volume",
  section_name = "3d_scan",
  resolution = "25um",
  matrix_name = "CHCA"
)

# 查看3D数据结构
cat("\n=== 3D Data Structure ===\n")

# 查看z切片分布
z_table <- table(pixel_info(my_msi_3d)$z)
cat("Z-slice distribution:\n")
print(z_table)

# 查看sample_info（每个z切片可能是独立的样本）
cat("\nSample info for 3D data:\n")
print(sample_info(my_msi_3d))

# 提取单个z切片
z_slice_1 <- pixel_info(my_msi_3d)[pixel_info(my_msi_3d)$z == 1, ]
cat("\nPixels in z-slice 1:", nrow(z_slice_1), "\n")

# ============================================================================
# 示例 8: 更新和修改元数据
# ============================================================================

my_msi <- readRDS("2_data/output/my_msi_dataset.rds")

# 场景：转换时使用了默认值，现在需要更新实际参数

# 1. 更新 sample_info
cat("\n=== Before Update ===\n")
print(sample_info(my_msi))

# 更新分辨率和基质信息
my_msi@sample_info$resolution <- "75um"
my_msi@sample_info$matrix <- "CHCA"

cat("\n=== After Update ===\n")
print(sample_info(my_msi))

# 2. 添加额外的样本信息列
my_msi@sample_info$tissue_type <- "brain"
my_msi@sample_info$organism <- "mouse"
my_msi@sample_info$date <- "2025-02-05"

# 更新 sample_info_note
new_notes <- data.frame(
  name = c("tissue_type", "organism", "date"),
  meaning = c("Type of tissue", "Source organism", "Collection date"),
  stringsAsFactors = FALSE
)

my_msi@sample_info_note <- rbind(
  my_msi@sample_info_note,
  new_notes
)

# 验证更新
print(my_msi)

# 保存更新后的对象
saveRDS(my_msi, "2_data/output/my_msi_dataset_updated.rds")

# ============================================================================
# 示例 9: 不同MALDI基质的数据转换
# ============================================================================

# 常见的MALDI基质类型
maldi_matrices <- c(
  "DHB",      # 2,5-Dihydroxybenzoic acid
  "CHCA",     # α-Cyano-4-hydroxycinnamic acid
  "SA",       # Sinapinic acid
  "DAN",      # 1,5-Diaminonaphthalene
  "9AA"       # 9-Aminoacridine
)

# 示例：DHB基质（常用于脂质成像）
msi_dhb <- convert_cardinal2msi_dataset(
  cardinal_object = readMSIData("2_data/demo data/lipid_sample.imzML"),
  sample_name = "lipid_analysis",
  matrix_name = "DHB",
  resolution = "50um"
)

# 示例：CHCA基质（常用于肽段和蛋白质成像）
msi_chca <- convert_cardinal2msi_dataset(
  cardinal_object = readMSIData("2_data/demo data/protein_sample.imzML"),
  sample_name = "protein_analysis",
  matrix_name = "CHCA",
  resolution = "100um"
)

# ============================================================================
# 示例 10: 不同空间分辨率的数据
# ============================================================================

# 高分辨率数据（10um）
msi_high_res <- convert_cardinal2msi_dataset(
  cardinal_object = readMSIData("2_data/demo data/high_res.imzML"),
  sample_name = "high_resolution_scan",
  resolution = "10um",
  matrix_name = "DHB"
)

# 中等分辨率数据（50um）
msi_mid_res <- convert_cardinal2msi_dataset(
  cardinal_object = readMSIData("2_data/demo data/mid_res.imzML"),
  sample_name = "mid_resolution_scan",
  resolution = "50um",
  matrix_name = "DHB"
)

# 低分辨率数据（100um）
msi_low_res <- convert_cardinal2msi_dataset(
  cardinal_object = readMSIData("2_data/demo data/low_res.imzML"),
  sample_name = "low_resolution_scan",
  resolution = "100um",
  matrix_name = "DHB"
)

# 比较不同分辨率的数据
cat("\n=== Resolution Comparison ===\n")
cat("High res pixels:", n_pixels(msi_high_res), "\n")
cat("Mid res pixels:", n_pixels(msi_mid_res), "\n")
cat("Low res pixels:", n_pixels(msi_low_res), "\n")

# ============================================================================
# 示例 11: 错误处理和故障排除
# ============================================================================

# 1. 处理读取失败的情况
tryCatch({
  cardinal_obj <- readMSIData("2_data/demo data/nonexistent.imzML")
  my_msi <- convert_cardinal2msi_dataset(cardinal_obj)
}, error = function(e) {
  cat("Error reading file:", e$message, "\n")
  cat("Please check the file path and format.\n")
})

# 2. 处理内存不足的情况
# 对于大文件，使用 memory = FALSE
cardinal_obj <- readMSIData(
  "2_data/demo data/large_file.imzML",
  memory = FALSE  # 不将数据全部加载到内存
)

# 3. 检查对齐结果
# 如果对齐后变量数异常，可以调整容差
my_msi_strict <- convert_cardinal2msi_dataset(
  cardinal_object = cardinal_obj,
  tolerance = 5,   # 更严格的对齐
  units = "ppm"
)

my_msi_loose <- convert_cardinal2msi_dataset(
  cardinal_object = cardinal_obj,
  tolerance = 20,  # 更宽松的对齐
  units = "ppm"
)

cat("Strict alignment variables:", n_variables(my_msi_strict), "\n")
cat("Loose alignment variables:", n_variables(my_msi_loose), "\n")

# ============================================================================
# 示例 12: 完整的项目工作流
# ============================================================================

# 项目设置
project_dir <- "D:/Shen_lab/msidataset"
data_dir <- file.path(project_dir, "2_data/demo data")
output_dir <- file.path(project_dir, "2_data/output")

# 创建输出目录
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 1. 读取数据
cat("\n[1/5] Reading Cardinal data...\n")
cardinal_obj <- readMSIData(
  file.path(data_dir, "03-total ion count.imzML"),
  memory = FALSE
)

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

# 3. 数据质量检查
cat("\n[3/5] Quality check...\n")
check_result <- check_msi_dataset(
  expression_data = my_msi@expression_data,
  pixel_info = my_msi@pixel_info,
  sample_info = my_msi@sample_info,
  variable_info = my_msi@variable_info
)

if (check_result != TRUE) {
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

# 5. 后续分析准备完成
cat("\n[5/5] Setup complete!\n")
cat("Data saved to:", file.path(output_dir, "processed_msi_dataset.rds"), "\n")
cat("Ready for downstream analysis.\n\n")

# ============================================================================
# 常见问题解答 (FAQ)
# ============================================================================

# Q1: 如何选择合适的对齐容差？
# A1: 通常使用 10 ppm 是合理的起点。对于高分辨率仪器可以用 5 ppm，
#     对于低分辨率或有较大质量偏移的数据可以用 15-20 ppm。

# Q2: resolution 参数应该填什么？
# A2: 填写实际的空间分辨率，如 "10um", "50um", "100um"。
#     如果不清楚，可以暂时填 "unknown"，之后更新。

# Q3: 如何处理没有z坐标的2D数据？
# A3: 函数会自动添加 z = 1。不需要特殊处理。

# Q4: 转换后如何访问原始 Cardinal 对象？
# A4: 转换是单向的，不会保留原始对象。如需要，请单独保存 Cardinal 对象。

# Q5: 如何合并多个 msi_dataset 对象？
# A5: 这需要专门的合并函数（未来功能）。目前建议在 Cardinal 层面先合并。

cat("\n=== All examples completed! ===\n")