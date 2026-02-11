# Quick Start Guide - convert_cardinal2msi_dataset
# ================================================

# ============================================
# 最小化示例（适合快速测试）
# ============================================

# Step 0: 环境设置
setwd("D:/Shen_lab/msidataset")
library(S4Vectors)
library(Cardinal)
devtools::load_all()

# Step 1: 读取 imzML 文件
cardinal_obj <- readMSIData(
  "data/demo data/03-total ion count.imzML",
  memory = FALSE
)
print(cardinal_obj)

# Step 2: 转换为 msi_dataset（最简单调用）
my_msi <- convert_cardinal2msi_dataset(
  cardinal_object = cardinal_obj,
  sample_name = "demo_sample",
  section_name = "section_1"
)
print(my_msi)

# Step 3: 快速验证
n_variables(my_msi)  # 变量数
n_pixels(my_msi)     # 像素数
range(mz(my_msi))    # m/z 范围

# Step 4: 保存
saveRDS(my_msi, "2_data/output/my_msi_dataset.rds")

# ============================================
# 推荐示例（用于正式分析）
# ============================================

# Step 0: 环境设置
setwd("D:/Shen_lab/msidataset")
library(S4Vectors)
library(Cardinal)
devtools::load_all()

# Step 1: 读取数据
cardinal_obj <- readMSIData(
  "2_data/demo data/03-total ion count.imzML",
  memory = FALSE
)

# Step 2: 转换（使用完整参数）
my_msi <- convert_cardinal2msi_dataset(
  cardinal_object = cardinal_obj,
  tolerance = 10,                    # ppm对齐容差
  units = "ppm",                     # 容差单位
  sample_name = "mouse_brain_01",    # 样本名称
  section_name = "coronal_slice_1",  # 切片名称
  resolution = "50um",               # 空间分辨率
  matrix_name = "DHB"                # MALDI基质
)

# Step 3: 查看数据结构
print(my_msi)

# 查看各个组件
head(pixel_info(my_msi))      # 像素空间信息
sample_info(my_msi)            # 样本元数据
head(variable_info(my_msi))   # 变量（m/z）信息

# Step 4: 数据验证
check_result <- check_msi_dataset(
  expression_data = my_msi@expression_data,
  pixel_info = my_msi@pixel_info,
  sample_info = my_msi@sample_info,
  variable_info = my_msi@variable_info
)

if (check_result == TRUE) {
  cat("✓ 数据验证通过！\n")
} else {
  cat("✗ 数据验证失败:\n", check_result, "\n")
}

# Step 5: 保存
saveRDS(my_msi, "2_data/output/my_msi_dataset_full.rds")

# Step 6: 后续分析
# 现在可以进行：
# - 数据过滤
# - 归一化
# - 统计分析
# - 空间可视化
# - 等等...

# ============================================
# 参数说明
# ============================================

# 必需参数：
# - cardinal_object: Cardinal的MSImagingExperiment或MSImagingArrays对象

# 推荐参数：
# - sample_name: 给样本一个有意义的名称（如 "mouse_brain_01"）
# - section_name: 切片标识符（如 "coronal_slice_1"）
# - resolution: 空间分辨率（如 "50um", "10um", "100um"）
# - matrix_name: MALDI基质（如 "DHB", "CHCA", "SA"）

# 可选参数：
# - tolerance: 对齐容差，默认10（一般不需要改）
# - units: 容差单位，默认"ppm"（一般不需要改）

# ============================================
# 常见MALDI基质类型
# ============================================

# matrix_name 参数的常见值：
# - "DHB"   : 2,5-Dihydroxybenzoic acid (常用于脂质)
# - "CHCA"  : α-Cyano-4-hydroxycinnamic acid (常用于肽段/蛋白)
# - "SA"    : Sinapinic acid (常用于大蛋白)
# - "DAN"   : 1,5-Diaminonaphthalene (用于负离子模式)
# - "9AA"   : 9-Aminoacridine (用于负离子模式)
# - "unknown": 如果不确定或没有使用基质

# ============================================
# 常见空间分辨率
# ============================================

# resolution 参数的常见值：
# - "10um"  : 高分辨率成像
# - "20um"  : 高分辨率成像
# - "50um"  : 标准分辨率成像（最常见）
# - "75um"  : 中等分辨率成像
# - "100um" : 低分辨率成像
# - "200um" : 快速扫描
# - "unknown": 如果不确定实际分辨率

# ============================================
# 故障排除
# ============================================

# 问题1: "Error in readMSIData..."
# 解决: 检查文件路径是否正确，文件是否存在

# 问题2: 转换后变量数太少或太多
# 解决: 调整 tolerance 参数
#   - tolerance = 5   (更严格，变量更多)
#   - tolerance = 20  (更宽松，变量更少)

# 问题3: 内存不足
# 解决: 在readMSIData时使用 memory = FALSE

# 问题4: 不知道分辨率或基质
# 解决: 暂时使用 "unknown"，之后可以手动更新
#   my_msi@sample_info$resolution <- "50um"
#   my_msi@sample_info$matrix <- "DHB"

# ============================================
# 数据检查清单
# ============================================

# 转换后应该检查：
# 1. 变量数是否合理
n_variables(my_msi)

# 2. 像素数是否正确
n_pixels(my_msi)

# 3. m/z范围是否符合预期
range(mz(my_msi))

# 4. 空间范围是否正确
range(pixel_info(my_msi)$x)
range(pixel_info(my_msi)$y)

# 5. 样本信息是否完整
sample_info(my_msi)

# 6. 数据结构是否有效
check_msi_dataset(
  expression_data = my_msi@expression_data,
  pixel_info = my_msi@pixel_info,
  sample_info = my_msi@sample_info,
  variable_info = my_msi@variable_info
)

cat("\n=== Quick Start Guide Complete! ===\n")