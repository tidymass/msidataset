# convert_cardinal2msi_dataset 使用示例
# Author: Xiaoyu Xu
# Date: 2025

# ============================================================
# 示例 1: 基本使用（最小参数）
# ============================================================

library(Cardinal)
library(msidataset)  # 假设你的包名

# 加载 Cardinal 数据
data("example_msi_data")  # 这是一个假设的示例数据

# 基本转换
msi_obj1 <- convert_cardinal2msi_dataset(
  cardinal_object = example_msi_data
)

# 查看对象
msi_obj1

# ============================================================
# 示例 2: 完整参数使用（推荐）
# ============================================================

msi_obj2 <- convert_cardinal2msi_dataset(
  cardinal_object = example_msi_data,
  tolerance = 10,                    # 对齐容差
  units = "ppm",                     # 容差单位
  sample_name = "mouse_brain_01",    # 样本名称
  section_name = "coronal_slice_1",  # 切片名称
  resolution = "50um",               # 空间分辨率
  matrix_name = "DHB"                # MALDI 基质
)

# 查看对象
msi_obj2

# ============================================================
# 示例 3: 查看转换后的数据结构
# ============================================================

# 查看 expression_data
cat("\n=== Expression Data ===\n")
cat("Dimensions:", dim(msi_obj2@expression_data), "\n")
cat("First 5 variables, first 3 pixels:\n")
print(msi_obj2@expression_data[1:5, 1:3])

# 查看 pixel_info
cat("\n=== Pixel Info ===\n")
cat("Columns:", colnames(msi_obj2@pixel_info), "\n")
print(head(msi_obj2@pixel_info))

# 查看 sample_info
cat("\n=== Sample Info ===\n")
cat("Columns:", colnames(msi_obj2@sample_info), "\n")
print(msi_obj2@sample_info)

# 查看 variable_info
cat("\n=== Variable Info ===\n")
cat("Columns:", colnames(msi_obj2@variable_info), "\n")
print(head(msi_obj2@variable_info))

# ============================================================
# 示例 4: 3D 数据转换
# ============================================================

# 假设有一个包含多个z切片的3D Cardinal数据
# 3d_cardinal_data <- ...

msi_obj_3d <- convert_cardinal2msi_dataset(
  cardinal_object = 3d_cardinal_data,
  sample_name = "3d_tissue_scan",
  section_name = "3d_volume",
  resolution = "25um",
  matrix_name = "CHCA"
)

# 检查多个z切片
cat("\n=== 3D Data - Sample Info ===\n")
print(msi_obj_3d@sample_info)

cat("\n=== 3D Data - Unique Z values in pixel_info ===\n")
print(table(msi_obj_3d@pixel_info$z))

# ============================================================
# 示例 5: 使用 accessor 函数访问数据
# ============================================================

# 如果定义了 accessor 函数
pixel_info_df <- pixel_info(msi_obj2)
sample_info_df <- sample_info(msi_obj2)
variable_info_df <- variable_info(msi_obj2)
expr_data <- expression_data(msi_obj2)

# ============================================================
# 示例 6: 转换 MSImagingArrays（需要自动对齐）
# ============================================================

# 假设有一个未对齐的 MSImagingArrays 对象
# arrays_data <- ...

msi_obj_aligned <- convert_cardinal2msi_dataset(
  cardinal_object = arrays_data,
  tolerance = 15,        # 使用更大的容差
  units = "ppm",
  sample_name = "unaligned_sample",
  section_name = "section_A",
  resolution = "100um",
  matrix_name = "DHB"
)

# 自动对齐会在第一步执行
# 输出会显示: "MSImagingArrays detected, performing automatic alignment..."

# ============================================================
# 示例 7: 数据验证
# ============================================================

# 验证对象结构
check_result <- check_msi_dataset(
  expression_data = msi_obj2@expression_data,
  pixel_info = msi_obj2@pixel_info,
  sample_info = msi_obj2@sample_info,
  variable_info = msi_obj2@variable_info
)

if (check_result == TRUE) {
  cat("\n✓ Object structure is valid!\n")
} else {
  cat("\n✗ Validation errors:\n")
  print(check_result)
}

# ============================================================
# 示例 8: 批量转换多个文件
# ============================================================

# 假设有多个 Cardinal 对象
cardinal_files <- list(
  sample1 = cardinal_data1,
  sample2 = cardinal_data2,
  sample3 = cardinal_data3
)

# 批量转换
msi_objects <- lapply(names(cardinal_files), function(sample_name) {
  convert_cardinal2msi_dataset(
    cardinal_object = cardinal_files[[sample_name]],
    sample_name = sample_name,
    section_name = paste0(sample_name, "_section1"),
    resolution = "50um",
    matrix_name = "DHB"
  )
})

names(msi_objects) <- names(cardinal_files)

# ============================================================
# 示例 9: 保存和加载转换后的对象
# ============================================================

# 保存对象
saveRDS(msi_obj2, "msi_dataset_sample1.rds")

# 加载对象
msi_obj_loaded <- readRDS("msi_dataset_sample1.rds")

# 验证加载的对象
msi_obj_loaded

# ============================================================
# 示例 10: 常见问题处理
# ============================================================

# 问题1: 如果不知道原始数据的分辨率和基质
msi_obj_unknown <- convert_cardinal2msi_dataset(
  cardinal_object = example_msi_data,
  sample_name = "unknown_sample"
  # resolution 和 matrix_name 会自动设置为 "unknown"
)

# 问题2: 手动更新元数据
msi_obj_unknown@sample_info$resolution <- "75um"
msi_obj_unknown@sample_info$matrix <- "CHCA"

# 问题3: 检查对齐前后的m/z范围
cat("\n=== m/z Range ===\n")
cat("Min m/z:", min(msi_obj2@variable_info$mz), "\n")
cat("Max m/z:", max(msi_obj2@variable_info$mz), "\n")
cat("Number of m/z values:", nrow(msi_obj2@variable_info), "\n")

# ============================================================
# 完整的工作流示例
# ============================================================

# 1. 加载 Cardinal 数据
library(Cardinal)
cardinal_data <- readMSIData("path/to/imzml_file.imzML")

# 2. 转换为 msi_dataset
msi_obj <- convert_cardinal2msi_dataset(
  cardinal_object = cardinal_data,
  sample_name = "tissue_sample_A",
  section_name = "section_001",
  resolution = "50um",
  matrix_name = "DHB"
)

# 3. 检查对象
msi_obj

# 4. 验证数据
check_msi_dataset(
  expression_data = msi_obj@expression_data,
  pixel_info = msi_obj@pixel_info,
  sample_info = msi_obj@sample_info,
  variable_info = msi_obj@variable_info
)

# 5. 进行下游分析
# ... 你的分析代码 ...

# 6. 保存结果
saveRDS(msi_obj, "processed_msi_dataset.rds")

cat("\n=== Conversion complete! ===\n")