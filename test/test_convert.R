# ============================================================================
# MSI Dataset ↔ Mass Dataset
# ============================================================================

# ============================================================================
# environment
# ============================================================================

setwd("D:/Shen_lab/msidataset")

# 加载必需的包
library(S4Vectors)
library(Cardinal)
library(massdataset)

# 加载 msidataset
devtools::document()
devtools::load_all()

# ============================================================================
# 步骤 1: 创建或加载 msi_dataset 对象
# ============================================================================

cat("\n=== Step 1: Prepare MSI Dataset ===\n")

# 方法 A: 从 Cardinal 对象转换
cardinal_obj <- readMSIData(
  "data/demo data/03-total ion count.imzML",
  memory = FALSE
)

my_msi <- convert_cardinal2msi_dataset_v1(
  cardinal_object = cardinal_obj,
  sample_name = "demo_sample",
  section_name = "section_1",
  resolution = "50um",
  matrix_name = "DHB",
  variable_id_mode = "mz_based"  # 使用基于 m/z 的命名
)

# 或方法 B: 加载已有的 msi_dataset
# my_msi <- readRDS("2_data/output/my_msi_dataset.rds")

# 查看 msi_dataset 结构
print(my_msi)

cat("\nMSI Dataset summary:\n")
cat("  Variables:", n_variables(my_msi), "\n")
cat("  Pixels:", n_pixels(my_msi), "\n")
cat("  m/z range:", paste(round(range(mz(my_msi)), 4), collapse = " - "), "\n")
cat("  Spatial extent:\n")
cat("    x:", paste(range(pixel_info(my_msi)$x), collapse = " - "), "\n")
cat("    y:", paste(range(pixel_info(my_msi)$y), collapse = " - "), "\n")
cat("    z:", paste(range(pixel_info(my_msi)$z), collapse = " - "), "\n")

# ============================================================================
# 步骤 2: 转换 msi_dataset → mass_dataset
# ============================================================================

cat("\n=== Step 2: Convert MSI to Mass Dataset ===\n")

# 基础转换（推荐）
my_mass <- convert_msi2mass_dataset(
  msi_object = my_msi,
  save_metadata = TRUE,  # 保存元数据到环境
  add_spatial_to_sample_info = FALSE  # 不添加空间坐标到 sample_info
)

# 或者：保存元数据到文件（用于长期保存或跨 session）
my_mass <- convert_msi2mass_dataset(
  msi_object = my_msi,
  save_metadata = TRUE,
  metadata_path = "2_data/output/msi_metadata.rds",  # 保存到文件
  add_spatial_to_sample_info = TRUE  # 可选：添加 x,y,z 到 sample_info
)

# 查看 mass_dataset 结构
print(my_mass)

cat("\nMass Dataset summary:\n")
cat("  Variables:", nrow(my_mass@variable_info), "\n")
cat("  Samples:", nrow(my_mass@sample_info), "\n")

# 查看 sample_info
cat("\nSample info preview:\n")
print(head(my_mass@sample_info, 5))

# ============================================================================
# 步骤 3: 使用 massdataset 进行分析
# ============================================================================

cat("\n=== Step 3: Perform Analysis with massdataset ===\n")

# 3.1 数据预处理
cat("\n--- Data Preprocessing ---\n")

# 变量过滤（示例：保留强度较高的变量）
cat("Filtering variables...\n")
my_mass_filtered <- filter_variables(
  object = my_mass,
  min_intensity = 1000,  # 最小强度阈值
  min_freq = 0.5  # 至少在 50% 的样本中检出
)

cat("  Variables before:", nrow(my_mass@variable_info), "\n")
cat("  Variables after:", nrow(my_mass_filtered@variable_info), "\n")

# 数据归一化（示例：TIC 归一化）
cat("\nNormalizing data...\n")
my_mass_normalized <- normalize_data(
  object = my_mass_filtered,
  method = "total"  # Total Ion Current (TIC) normalization
)

# 缺失值填补（如果需要）
if (sum(is.na(my_mass_normalized@expression_data)) > 0) {
  cat("\nImputing missing values...\n")
  my_mass_imputed <- impute_mvs(
    object = my_mass_normalized,
    method = "knn"  # K-nearest neighbors imputation
  )
  my_mass_processed <- my_mass_imputed
} else {
  cat("\nNo missing values detected\n")
  my_mass_processed <- my_mass_normalized
}

# 3.2 统计分析（示例）
cat("\n--- Statistical Analysis ---\n")

# 数据缩放
cat("Scaling data...\n")
my_mass_scaled <- scale_variables(
  object = my_mass_processed,
  method = "auto"  # Auto-scaling (mean-centering and unit variance)
)

# 如果有分组信息，可以进行差异分析
# 这里我们创建一个简单的分组（基于空间位置）
if ("x" %in% colnames(my_mass@sample_info)) {
  cat("\nCreating spatial groups...\n")
  
  # 根据 x 坐标分为左右两组
  median_x <- median(my_mass@sample_info$x)
  my_mass_scaled@sample_info$class <- ifelse(
    my_mass@sample_info$x < median_x,
    "left",
    "right"
  )
  
  cat("  Left group:", sum(my_mass_scaled@sample_info$class == "left"), "samples\n")
  cat("  Right group:", sum(my_mass_scaled@sample_info$class == "right"), "samples\n")
}

# 3.3 代谢物注释（示例）
cat("\n--- Metabolite Annotation (Example) ---\n")

# 注意：实际注释需要数据库，这里仅作演示
# 假设我们添加一些注释信息到 variable_info

# 示例：基于 m/z 值添加模拟注释
my_mass_annotated <- my_mass_scaled

# 创建一些示例注释
set.seed(123)
n_vars <- nrow(my_mass_annotated@variable_info)

my_mass_annotated@variable_info$compound_name <- NA
my_mass_annotated@variable_info$formula <- NA
my_mass_annotated@variable_info$adduct <- "[M+H]+"
my_mass_annotated@variable_info$confidence_level <- NA

# 模拟：为一些变量添加注释
annotated_idx <- sample(1:n_vars, min(50, floor(n_vars * 0.1)))
my_mass_annotated@variable_info$compound_name[annotated_idx] <- 
  paste0("Compound_", 1:length(annotated_idx))
my_mass_annotated@variable_info$confidence_level[annotated_idx] <- 
  sample(1:3, length(annotated_idx), replace = TRUE)

cat("  Total variables:", n_vars, "\n")
cat("  Annotated:", sum(!is.na(my_mass_annotated@variable_info$compound_name)), "\n")

# 3.4 主成分分析 (PCA)
cat("\n--- PCA Analysis (Example) ---\n")

# 提取表达矩阵用于 PCA
expr_matrix <- t(as.matrix(my_mass_scaled@expression_data))

# 执行 PCA
pca_result <- prcomp(expr_matrix, center = TRUE, scale. = FALSE)

cat("  PCA completed\n")
cat("  PC1 variance explained:", 
    round(summary(pca_result)$importance[2,1] * 100, 2), "%\n")
cat("  PC2 variance explained:", 
    round(summary(pca_result)$importance[2,2] * 100, 2), "%\n")

# 保存 PCA 结果到 process_info
pca_parameter <- new(
  Class = "tidymass_parameter",
  pacakge_name = "stats",
  function_name = "prcomp()",
  parameter = list(
    pc1_var = summary(pca_result)$importance[2,1],
    pc2_var = summary(pca_result)$importance[2,2]
  ),
  time = Sys.time()
)

my_mass_annotated@process_info$pca_analysis <- pca_parameter

# ============================================================================
# 步骤 4: 转换回 msi_dataset
# ============================================================================

cat("\n=== Step 4: Convert Mass Dataset back to MSI Dataset ===\n")

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

# ============================================================================
# 步骤 5: 检查恢复的数据
# ============================================================================

cat("\n=== Step 5: Verify Restored Data ===\n")

# 5.1 检查空间信息是否恢复
cat("\n--- Spatial Information ---\n")
cat("pixel_info dimensions:", dim(pixel_info(my_msi_restored)), "\n")
print(head(pixel_info(my_msi_restored), 3))

cat("\nsample_info dimensions:", dim(sample_info(my_msi_restored)), "\n")
print(sample_info(my_msi_restored))

# 5.2 检查 annotation_table
cat("\n--- Annotation Table ---\n")
cat("annotation_table dimensions:", dim(annotation_table(my_msi_restored)), "\n")
if (nrow(annotation_table(my_msi_restored)) > 0) {
  print(head(annotation_table(my_msi_restored), 5))
}

# 5.3 检查 variable_info
cat("\n--- Variable Info ---\n")
cat("Annotated variables:", 
    sum(!is.na(variable_info(my_msi_restored)$mz)), "\n")
print(head(variable_info(my_msi_restored), 5))

# 5.4 检查 process_info
cat("\n--- Processing History ---\n")
cat("Number of processing steps:", length(my_msi_restored@process_info), "\n")
cat("Processing steps:\n")
print(names(my_msi_restored@process_info))

# ============================================================================
# 步骤 6: 保存结果
# ============================================================================

cat("\n=== Step 6: Save Results ===\n")

# 保存更新后的 msi_dataset
dir.create("2_data/output", showWarnings = FALSE, recursive = TRUE)

saveRDS(my_msi_restored, "2_data/output/my_msi_annotated.rds")
cat("✓ MSI dataset saved to: 2_data/output/my_msi_annotated.rds\n")

# 保存 annotation_table 为 CSV
if (nrow(annotation_table(my_msi_restored)) > 0) {
  write.csv(
    annotation_table(my_msi_restored),
    "2_data/output/annotation_table.csv",
    row.names = FALSE
  )
  cat("✓ Annotation table saved to: 2_data/output/annotation_table.csv\n")
}

# ============================================================================
# 步骤 7: 查看转换历史
# ============================================================================

cat("\n=== Step 7: View Conversion History ===\n")

# 查看转换历史
conversion_summary <- get_conversion_summary(my_msi_restored)

# 列出缓存中的元数据
cat("\n--- Metadata Cache ---\n")
list_msi_metadata()

# ============================================================================
# 辅助函数：比较原始和恢复的数据
# ============================================================================

compare_msi_objects <- function(original, restored) {
  cat("\n========================================\n")
  cat("Comparison: Original vs Restored\n")
  cat("========================================\n\n")
  
  # 维度比较
  cat("Dimensions:\n")
  cat("  Variables - Original:", n_variables(original), 
      "| Restored:", n_variables(restored), "\n")
  cat("  Pixels    - Original:", n_pixels(original), 
      "| Restored:", n_pixels(restored), "\n")
  
  # 数据范围比较
  cat("\nData ranges:\n")
  cat("  Original - min:", round(min(expression_data(original)), 2),
      "max:", round(max(expression_data(original)), 2), "\n")
  cat("  Restored - min:", round(min(expression_data(restored)), 2),
      "max:", round(max(expression_data(restored)), 2), "\n")
  
  # 空间信息比较
  cat("\nSpatial extent:\n")
  cat("  Original - x:", paste(range(pixel_info(original)$x), collapse = "-"),
      "y:", paste(range(pixel_info(original)$y), collapse = "-"), "\n")
  cat("  Restored - x:", paste(range(pixel_info(restored)$x), collapse = "-"),
      "y:", paste(range(pixel_info(restored)$y), collapse = "-"), "\n")
  
  # 新增内容
  if (nrow(annotation_table(restored)) > 0) {
    cat("\nNew content in restored object:\n")
    cat("  Annotation table:", nrow(annotation_table(restored)), "rows\n")
  }
  
  cat("\n========================================\n")
}

# 执行比较
compare_msi_objects(my_msi, my_msi_restored)

# ============================================================================
# 总结
# ============================================================================

cat("\n========================================\n")
cat("Workflow Completed Successfully!\n")
cat("========================================\n\n")

cat("Summary of workflow:\n")
cat("1. ✓ Created/loaded msi_dataset\n")
cat("2. ✓ Converted to mass_dataset\n")
cat("3. ✓ Performed analysis with massdataset\n")
cat("4. ✓ Converted back to msi_dataset\n")
cat("5. ✓ Verified data integrity\n")
cat("6. ✓ Saved results\n\n")

cat("Output files:\n")
cat("  - 2_data/output/my_msi_annotated.rds (updated msi_dataset)\n")
cat("  - 2_data/output/annotation_table.csv (annotation results)\n")
cat("  - 2_data/output/msi_metadata.rds (conversion metadata)\n\n")

cat("Next steps:\n")
cat("  - Use updated msi_dataset for MSI-specific visualization\n")
cat("  - Generate spatial distribution maps\n")
cat("  - Perform region-of-interest analysis\n")
cat("  - Integrate with imaging data\n\n")

cat("========================================\n")