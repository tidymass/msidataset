# ============================================================================
# Massdataset 完整使用指南（基于 GitHub 官方文档）
# 适用于从 convert_msi2mass_dataset() 转换后的对象
# ============================================================================

# ============================================================================
# 环境设置
# ============================================================================

library(massdataset)
library(tidyverse)
library(ggplot2)

# 假设你已经有了转换后的 my_mass 对象
# my_mass <- convert_msi2mass_dataset(my_msi, save_metadata = TRUE)

# ============================================================================
# 第一部分：数据过滤（使用正确的函数）
# ============================================================================

cat("\n=== Part 1: Data Filtering ===\n")

# ----------------------------------------------------------------------------
# 1.1 mutate_variable_na_freq() - 计算缺失值频率
# ----------------------------------------------------------------------------

cat("\n--- 1.1 Calculate NA frequency ---\n")

# 为所有样本计算缺失值频率
my_mass <- my_mass %>%
  mutate_variable_na_freq()

# 查看结果
head(extract_variable_info(my_mass))

# 查看新增的 na_freq 列
# variable_id        mz   rt   na_freq
# m50_1        50.00125    0  0.000000
# m50_2        50.00250    0  0.428571
# m50_3        50.00375    0  0.000000

# ----------------------------------------------------------------------------
# 1.2 filter_variables() - 基于缺失值过滤变量（正确用法）
# 函数‘n_variables’标签‘object = "mass_dataset"’找不到继承方法
# ----------------------------------------------------------------------------

cat("\n--- 1.2 Filter variables by NA frequency ---\n")

# 方法 A: 使用 filter_variables() 函数
my_mass_filtered <- filter_variables(
  object = my_mass,
  flist = function(x) {
    sum(is.na(x)) / length(x) < 0.5  # 缺失值 < 50%
  },
  prune = TRUE
)

# 方法 B: 使用 tidyverse 风格（推荐）
my_mass_filtered <- my_mass %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(na_freq < 0.5)  # 保留缺失值 < 50% 的变量

cat("Variables before:", n_variables(my_mass), "\n")
cat("Variables after:", n_variables(my_mass_filtered), "\n")

# ----------------------------------------------------------------------------
# 1.3 filter_samples() - 过滤样本
# ----------------------------------------------------------------------------

cat("\n--- 1.3 Filter samples ---\n")

# 基于缺失值过滤样本
my_mass_filtered2 <- filter_samples(
  object = my_mass,
  flist = function(x) {
    sum(is.na(x)) / length(x) < 0.4  # 缺失值 < 40%
  },
  prune = TRUE
)

# ----------------------------------------------------------------------------
# 1.4 基于空间分组过滤
# ----------------------------------------------------------------------------

cat("\n--- 1.4 Filter by spatial groups ---\n")

# 如果你已经创建了空间分组（spatial_group）
if ("spatial_group" %in% colnames(my_mass@sample_info)) {
  # 只保留 "left" 区域的样本
  my_mass_left <- my_mass %>%
    activate_mass_dataset(what = "sample_info") %>%
    filter(spatial_group == "left")
  
  cat("Samples in left region:", ncol(my_mass_left), "\n")
}

# ============================================================================
# 第二部分：数据预处理
# ============================================================================

cat("\n=== Part 2: Data Preprocessing ===\n")

# ----------------------------------------------------------------------------
# 2.1 mutate_rsd() - 计算相对标准偏差
# ----------------------------------------------------------------------------

cat("\n--- 2.1 Calculate RSD ---\n")

# 对所有样本计算 RSD
my_mass_filtered <- my_mass_filtered %>%
  mutate_rsd()

# 对特定样本计算 RSD（例如只计算某个空间区域）
if ("spatial_group" %in% colnames(my_mass_filtered@sample_info)) {
  left_sample_ids <- my_mass_filtered %>%
    activate_mass_dataset(what = "sample_info") %>%
    filter(spatial_group == "left") %>%
    pull(sample_id)
  
  my_mass_filtered <- my_mass_filtered %>%
    mutate_rsd(according_to_samples = left_sample_ids)
}

# 查看 RSD 结果
head(extract_variable_info(my_mass_filtered))

# ----------------------------------------------------------------------------
# 2.2 数据归一化（没有专门的函数，需要手动实现）
# ----------------------------------------------------------------------------

cat("\n--- 2.2 Normalize data ---\n")

# TIC (Total Ion Current) 归一化
normalize_tic <- function(object) {
  expr_data <- object@expression_data
  
  # 计算每个样本的总离子流
  tic <- colSums(expr_data, na.rm = TRUE)
  median_tic <- median(tic, na.rm = TRUE)
  
  # 归一化
  expr_data_normalized <- sweep(expr_data, 2, tic, "/") * median_tic
  
  # 更新对象
  object@expression_data <- expr_data_normalized
  
  return(object)
}

# 应用 TIC 归一化
my_mass_normalized <- normalize_tic(my_mass_filtered)

# ----------------------------------------------------------------------------
# 2.3 缺失值填补
# ----------------------------------------------------------------------------

cat("\n--- 2.3 Impute missing values ---\n")

# 检查缺失值
na_count <- sum(is.na(my_mass_normalized@expression_data))
cat("Missing values:", na_count, "\n")

if (na_count > 0) {
  # 方法 A: 中位数填补
  impute_median <- function(object) {
    expr_data <- object@expression_data
    
    for (i in 1:nrow(expr_data)) {
      if (any(is.na(expr_data[i, ]))) {
        median_val <- median(as.numeric(expr_data[i, ]), na.rm = TRUE)
        expr_data[i, is.na(expr_data[i, ])] <- median_val
      }
    }
    
    object@expression_data <- expr_data
    return(object)
  }
  
  # 方法 B: 最小值填补
  impute_min <- function(object) {
    expr_data <- object@expression_data
    
    for (i in 1:nrow(expr_data)) {
      if (any(is.na(expr_data[i, ]))) {
        min_val <- min(as.numeric(expr_data[i, ]), na.rm = TRUE) / 2
        expr_data[i, is.na(expr_data[i, ])] <- min_val
      }
    }
    
    object@expression_data <- expr_data
    return(object)
  }
  
  # 应用中位数填补
  my_mass_imputed <- impute_median(my_mass_normalized)
  
  cat("Missing values after imputation:", 
      sum(is.na(my_mass_imputed@expression_data)), "\n")
}

# ----------------------------------------------------------------------------
# 2.4 数据缩放
# 错误: 找不到对象'my_mass_imputed'
# ----------------------------------------------------------------------------

cat("\n--- 2.4 Scale data ---\n")

# Auto-scaling (mean-centering + unit variance)
scale_auto <- function(object) {
  expr_data <- object@expression_data
  
  # Mean centering
  row_means <- rowMeans(expr_data, na.rm = TRUE)
  expr_data_centered <- expr_data - row_means
  
  # Unit variance
  row_sds <- apply(expr_data, 1, sd, na.rm = TRUE)
  expr_data_scaled <- expr_data_centered / row_sds
  
  object@expression_data <- expr_data_scaled
  return(object)
}

my_mass_scaled <- scale_auto(my_mass_imputed)

# ============================================================================
# 第三部分：可视化
# ============================================================================

cat("\n=== Part 3: Visualization ===\n")

# ----------------------------------------------------------------------------
# 3.1 intensity_plot() - 强度图
# ----------------------------------------------------------------------------

cat("\n--- 3.1 Intensity plots ---\n")

# 基础用法
intensity_plot(object = my_mass, variable_index = 1)

# 按分组着色（需要先创建分组）
if ("spatial_group" %in% colnames(my_mass@sample_info)) {
  intensity_plot(
    object = my_mass, 
    variable_index = 1,
    color_by = "spatial_group"
  )
}

# 绘制多个变量
for (i in 1:5) {
  p <- intensity_plot(
    object = my_mass,
    variable_index = i,
    color_by = "spatial_group"
  )
  print(p)
}

# ----------------------------------------------------------------------------
# 3.2 show_sample_missing_values() - 缺失值热图
# ----------------------------------------------------------------------------

cat("\n--- 3.2 Missing values heatmap ---\n")

# 基础用法
show_sample_missing_values(
  object = my_mass,
  percentage = TRUE
) +
  theme(axis.text.x = element_text(size = 8))

# 按注射顺序排列（如果有 injection.order 列）
if ("injection.order" %in% colnames(my_mass@sample_info)) {
  show_sample_missing_values(
    object = my_mass,
    order_by = "injection.order",
    percentage = TRUE
  )
}

# ----------------------------------------------------------------------------
# 3.3 show_missing_values() - 变量缺失值图
# ----------------------------------------------------------------------------

cat("\n--- 3.3 Variable missing values ---\n")

# 显示所有缺失值
show_missing_values(object = my_mass)

# 只显示特定 m/z 范围
my_mass %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(mz >= 100 & mz <= 300) %>%
  show_missing_values(
    cell_color = "white",
    show_row_names = TRUE
  )

# ============================================================================
# 第四部分：统计分析
# ============================================================================

cat("\n=== Part 4: Statistical Analysis ===\n")

# ----------------------------------------------------------------------------
# 4.1 mutate_fc() - 计算 Fold Change
# 错误于h(simpleError(msg, call)): 
# 在为函数“%in%”选择方法时计算参数“table”时出错：在为函数“colnames”选择方法时计算参数“x”时出错：找不到对象'my_mass_scaled'
# ----------------------------------------------------------------------------

cat("\n--- 4.1 Calculate Fold Change ---\n")

# 首先需要定义对照组和实验组
# 如果有空间分组
if ("spatial_group" %in% colnames(my_mass_scaled@sample_info)) {
  # 左边作为对照组，右边作为实验组
  control_sample_id <- my_mass_scaled %>%
    activate_mass_dataset(what = "sample_info") %>%
    filter(spatial_group == "left") %>%
    pull(sample_id)
  
  case_sample_id <- my_mass_scaled %>%
    activate_mass_dataset(what = "sample_info") %>%
    filter(spatial_group == "right") %>%
    pull(sample_id)
  
  # 计算 Fold Change
  my_mass_scaled <- mutate_fc(
    object = my_mass_scaled,
    control_sample_id = control_sample_id,
    case_sample_id = case_sample_id,
    mean_median = "mean"
  )
  
  # 查看结果
  head(extract_variable_info(my_mass_scaled))
}

# ----------------------------------------------------------------------------
# 4.2 mutate_p_value() - 计算 p 值
# ----------------------------------------------------------------------------

cat("\n--- 4.2 Calculate p-values ---\n")

if ("spatial_group" %in% colnames(my_mass_scaled@sample_info)) {
  my_mass_scaled <- mutate_p_value(
    object = my_mass_scaled,
    control_sample_id = control_sample_id,
    case_sample_id = case_sample_id,
    method = "t.test",  # 或 "wilcox.test"
    p_adjust_methods = "BH"  # Benjamini-Hochberg 校正
  )
  
  # 查看结果
  var_info <- extract_variable_info(my_mass_scaled)
  head(var_info[, c("variable_id", "mz", "fc", "p_value", "p_value_adjust")])
}

# ----------------------------------------------------------------------------
# 4.3 volcano_plot() - 火山图
# ----------------------------------------------------------------------------

cat("\n--- 4.3 Volcano plot ---\n")

if ("fc" %in% colnames(extract_variable_info(my_mass_scaled)) &&
    "p_value" %in% colnames(extract_variable_info(my_mass_scaled))) {
  
  volcano_plot(
    object = my_mass_scaled,
    add_text = TRUE,
    text_from = "variable_id",
    point_size_scale = "p_value"
  ) +
    scale_size_continuous(range = c(0.5, 5))
}

# ============================================================================
# 第五部分：使用 tidyverse 操作
# ============================================================================

cat("\n=== Part 5: Tidyverse Operations ===\n")

# ----------------------------------------------------------------------------
# 5.1 activate_mass_dataset() - 激活不同数据层
# ----------------------------------------------------------------------------

cat("\n--- 5.1 Activate different data layers ---\n")

# 激活 sample_info
my_mass %>%
  activate_mass_dataset(what = "sample_info") %>%
  glimpse()

# 激活 variable_info
my_mass %>%
  activate_mass_dataset(what = "variable_info") %>%
  glimpse()

# 激活 expression_data
my_mass %>%
  activate_mass_dataset(what = "expression_data") %>%
  glimpse()

# ----------------------------------------------------------------------------
# 5.2 filter() - 过滤
# ----------------------------------------------------------------------------

cat("\n--- 5.2 Filter operations ---\n")

# 过滤特定 class 的样本
my_mass %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class == "pixel")

# 过滤特定 m/z 范围的变量
my_mass %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(mz >= 100 & mz <= 200)

# 过滤 RSD > 30 的变量
my_mass %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(rsd > 30)

# ----------------------------------------------------------------------------
# 5.3 mutate() - 添加新列
# ----------------------------------------------------------------------------

cat("\n--- 5.3 Add new columns ---\n")

# 在 sample_info 中添加新列
my_mass2 <- my_mass %>%
  activate_mass_dataset(what = "sample_info") %>%
  mutate(sample_id2 = paste0(sample_id, "_copy"))

# 在 variable_info 中添加新列
my_mass2 <- my_mass %>%
  activate_mass_dataset(what = "variable_info") %>%
  mutate(mz_rounded = round(mz, 2))

# ----------------------------------------------------------------------------
# 5.4 select() - 选择特定样本
# ----------------------------------------------------------------------------

cat("\n--- 5.4 Select samples ---\n")

# 选择前 10 个样本
my_mass_subset <- my_mass %>%
  activate_mass_dataset(what = "expression_data") %>%
  select(1:10)

# 选择名称包含特定字符的样本
my_mass_subset <- my_mass %>%
  activate_mass_dataset(what = "expression_data") %>%
  select(starts_with("pixel_000"))

# ============================================================================
# 第六部分：数据导出
# ============================================================================

cat("\n=== Part 6: Data Export ===\n")

# ----------------------------------------------------------------------------
# 6.1 extract_*() - 提取数据
# ----------------------------------------------------------------------------

cat("\n--- 6.1 Extract data ---\n")

# 提取 expression_data
expr_data <- extract_expression_data(my_mass)
dim(expr_data)

# 提取 sample_info
sample_info <- extract_sample_info(my_mass)
head(sample_info)

# 提取 variable_info
variable_info <- extract_variable_info(my_mass)
head(variable_info)

# ----------------------------------------------------------------------------
# 6.2 export_mass_dataset() - 导出为文件
# ----------------------------------------------------------------------------

cat("\n--- 6.2 Export to files ---\n")

# 导出整个 mass_dataset
dir.create("2_data/output/massdataset_export", showWarnings = FALSE, recursive = TRUE)

export_mass_dataset(
  object = my_mass_scaled,
  path = "2_data/output/massdataset_export"
)

# 会生成以下文件：
# - expression_data.csv
# - sample_info.csv
# - variable_info.csv
# - sample_info_note.csv
# - variable_info_note.csv

# ============================================================================
# 第七部分：完整工作流示例
# ============================================================================

cat("\n=== Part 7: Complete Workflow Example ===\n")

complete_workflow <- function(mass_object) {
  
  cat("\n========================================\n")
  cat("Complete Mass Dataset Analysis Workflow\n")
  cat("========================================\n\n")
  
  # Step 1: 计算缺失值频率
  cat("[1/7] Calculating NA frequencies...\n")
  mass_object <- mass_object %>%
    mutate_variable_na_freq()
  
  # Step 2: 过滤变量
  cat("[2/7] Filtering variables...\n")
  mass_object <- mass_object %>%
    activate_mass_dataset(what = "variable_info") %>%
    filter(na_freq < 0.5)
  
  # Step 3: 归一化
  cat("[3/7] Normalizing data...\n")
  mass_object <- normalize_tic(mass_object)
  
  # Step 4: 缺失值填补
  cat("[4/7] Imputing missing values...\n")
  if (sum(is.na(mass_object@expression_data)) > 0) {
    mass_object <- impute_median(mass_object)
  }
  
  # Step 5: 计算 RSD
  cat("[5/7] Calculating RSD...\n")
  mass_object <- mass_object %>%
    mutate_rsd()
  
  # Step 6: 数据缩放
  cat("[6/7] Scaling data...\n")
  mass_object <- scale_auto(mass_object)
  
  # Step 7: 统计分析（如果有分组）
  cat("[7/7] Statistical analysis...\n")
  if ("spatial_group" %in% colnames(mass_object@sample_info)) {
    control_ids <- mass_object %>%
      activate_mass_dataset(what = "sample_info") %>%
      filter(spatial_group == "left") %>%
      pull(sample_id)
    
    case_ids <- mass_object %>%
      activate_mass_dataset(what = "sample_info") %>%
      filter(spatial_group == "right") %>%
      pull(sample_id)
    
    if (length(control_ids) > 0 && length(case_ids) > 0) {
      mass_object <- mutate_fc(
        mass_object,
        control_sample_id = control_ids,
        case_sample_id = case_ids,
        mean_median = "mean"
      )
      
      mass_object <- mutate_p_value(
        mass_object,
        control_sample_id = control_ids,
        case_sample_id = case_ids,
        method = "t.test",
        p_adjust_methods = "BH"
      )
    }
  }
  
  cat("\n========================================\n")
  cat("Workflow Complete!\n")
  cat("========================================\n")
  
  return(mass_object)
}

# 应用完整工作流
my_mass_processed <- complete_workflow(my_mass)

# 保存结果
saveRDS(my_mass_processed, "2_data/output/mass_dataset_processed.rds")

cat("\n=== Analysis Complete ===\n")
cat("Processed data saved to: 2_data/output/mass_dataset_processed.rds\n")