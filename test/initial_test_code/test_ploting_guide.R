# ============================================================================
# Massdataset 画图和分析函数使用指南
# 基于转换后的 mass_dataset 对象
# ============================================================================

# ============================================================================
# 环境设置
# ============================================================================

setwd("D:/Shen_lab/msidataset")
library(S4Vectors)
library(Cardinal)
library(massdataset)
library(tidyverse)  # 用于数据处理

devtools::load_all()

# ============================================================================
# 前提：已经转换好的 mass_dataset 对象
# ============================================================================

# 假设你已经有了转换好的 my_mass 对象
# 如果没有，先进行转换：

# my_msi <- readRDS("2_data/output/my_msi_dataset.rds")
# my_mass <- convert_msi2mass_dataset(
#   msi_object = my_msi,
#   save_metadata = TRUE,
#   metadata_path = "2_data/output/msi_metadata.rds"
# )

# 查看当前对象
cat("Current mass_dataset object:\n")
print(my_mass)

# ============================================================================
# 1. intensity_plot() - 变量强度可视化
# ============================================================================

cat("\n=== Example 1: Intensity Plot ===\n")

# 1.1 基础用法：绘制第一个变量的强度
cat("\n--- Basic intensity plot ---\n")
intensity_plot(object = my_mass, variable_index = 1)

# 1.2 绘制特定变量（按索引）
cat("\n--- Plot specific variable by index ---\n")
intensity_plot(object = my_mass, variable_index = 10)

# 1.3 按类别着色（如果有 class 列）
cat("\n--- Plot with color by class ---\n")

# 首先检查 sample_info 是否有 class 列
if ("class" %in% colnames(my_mass@sample_info)) {
  intensity_plot(
    object = my_mass, 
    variable_index = 1, 
    color_by = "class"
  )
} else {
  cat("No 'class' column found in sample_info\n")
}

# 1.4 创建有意义的分组后绘图
cat("\n--- Create spatial groups and plot ---\n")

# 从元数据恢复空间坐标（如果之前保存了）
if (exists(".msi_conversion_cache", envir = .GlobalEnv)) {
  cache_env <- get(".msi_conversion_cache", envir = .GlobalEnv)
  cache_ids <- ls(cache_env)
  
  if (length(cache_ids) > 0) {
    metadata <- cache_env[[cache_ids[1]]]
    pixel_coords <- metadata$pixel_info
    
    # 创建空间分组
    # 例如：根据 x 坐标分为左中右三组
    x_coords <- pixel_coords$x[match(my_mass@sample_info$sample_id, 
                                     pixel_coords$pixel_id)]
    
    # 计算三分位数
    quantiles <- quantile(x_coords, probs = c(1/3, 2/3), na.rm = TRUE)
    
    my_mass@sample_info$spatial_group <- cut(
      x_coords,
      breaks = c(-Inf, quantiles[1], quantiles[2], Inf),
      labels = c("left", "center", "right")
    )
    
    cat("Created spatial groups:\n")
    print(table(my_mass@sample_info$spatial_group))
    
    # 绘图
    intensity_plot(
      object = my_mass, 
      variable_index = 1, 
      color_by = "spatial_group"
    )
    
    # 多个变量比较
    cat("\nPlotting multiple variables:\n")
    for (i in 1:5) {
      p <- intensity_plot(
        object = my_mass, 
        variable_index = i, 
        color_by = "spatial_group"
      )
      print(p)
    }
  }
}

# 1.5 按 y 坐标分组
if (exists("pixel_coords")) {
  y_coords <- pixel_coords$y[match(my_mass@sample_info$sample_id, 
                                   pixel_coords$pixel_id)]
  
  y_quantiles <- quantile(y_coords, probs = c(1/3, 2/3), na.rm = TRUE)
  
  my_mass@sample_info$y_group <- cut(
    y_coords,
    breaks = c(-Inf, y_quantiles[1], y_quantiles[2], Inf),
    labels = c("bottom", "middle", "top")
  )
  
  intensity_plot(
    object = my_mass, 
    variable_index = 1, 
    color_by = "y_group"
  )
}

# 1.6 按变量名绘图（如果知道变量 ID）
cat("\n--- Plot by variable ID ---\n")

# 找到特定 m/z 的变量
target_mz <- 100.0  # 例如：寻找 m/z 100 附近的变量
mz_values <- my_mass@variable_info$mz
closest_idx <- which.min(abs(mz_values - target_mz))

cat("Closest variable to m/z", target_mz, ":\n")
cat("  Variable ID:", my_mass@variable_info$variable_id[closest_idx], "\n")
cat("  Actual m/z:", mz_values[closest_idx], "\n")

intensity_plot(
  object = my_mass, 
  variable_index = closest_idx
)

# ============================================================================
# 2. show_missing_values() - 缺失值可视化
# ============================================================================

cat("\n=== Example 2: Missing Values Visualization ===\n")

# 2.1 基础用法：显示所有缺失值
cat("\n--- Show all missing values ---\n")
show_missing_values(object = my_mass)

# 2.2 自定义单元格颜色
cat("\n--- Custom cell color ---\n")
show_missing_values(
  object = my_mass,
  cell_color = "white"
)

# 2.3 只显示部分变量（前 20 个）
cat("\n--- Show first 20 variables ---\n")
show_missing_values(
  object = my_mass[1:20, ],
  cell_color = "lightblue"
)

# 2.4 过滤特定 class 的样本
cat("\n--- Filter by class ---\n")

# 如果有 class 列
if ("class" %in% colnames(my_mass@sample_info)) {
  # 显示特定 class 的缺失值
  my_mass %>%
    activate_mass_dataset(what = "sample_info") %>%
    filter(class == "pixel") %>%
    show_missing_values()
}

# 2.5 过滤特定空间区域的样本
cat("\n--- Filter by spatial group ---\n")

if ("spatial_group" %in% colnames(my_mass@sample_info)) {
  # 只显示 "left" 区域的缺失值
  my_mass %>%
    activate_mass_dataset(what = "sample_info") %>%
    filter(spatial_group == "left") %>%
    show_missing_values(cell_color = "white")
}

# 2.6 选择特定样本
cat("\n--- Select specific samples ---\n")

# 选择前 10 个样本
my_mass %>%
  activate_mass_dataset(what = "expression_data") %>%
  dplyr::select(1:10) %>%
  show_missing_values()

# 2.7 过滤特定 m/z 范围的变量
cat("\n--- Filter by m/z range ---\n")

# 只显示 m/z < 200 的变量
my_mass %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(mz < 200) %>%
  show_missing_values(
    cell_color = "white",
    show_row_names = TRUE,
    row_names_side = "left"
  )

# 2.8 显示百分比并设置缺失值阈值
cat("\n--- Show percentage with cutoffs ---\n")

my_mass %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(mz >= 100 & mz <= 300) %>%
  show_missing_values(
    cell_color = "white",
    show_row_names = TRUE,
    row_names_side = "left",
    percentage = TRUE,           # 显示百分比
    sample_na_cutoff = 50,       # 样本缺失值阈值 50%
    variable_na_cutoff = 20      # 变量缺失值阈值 20%
  )

# 2.9 检查不同区域的缺失值模式
cat("\n--- Compare missing values across spatial groups ---\n")

if ("spatial_group" %in% colnames(my_mass@sample_info)) {
  for (group in c("left", "center", "right")) {
    cat("\nMissing values for", group, "region:\n")
    
    my_mass %>%
      activate_mass_dataset(what = "sample_info") %>%
      filter(spatial_group == group) %>%
      show_missing_values(
        cell_color = "lightblue",
        percentage = TRUE
      )
  }
}

# ============================================================================
# 3. 准备 annotation_table 用于后续分析
# ============================================================================

cat("\n=== Example 3: Create Annotation Table ===\n")

# 3.1 创建基础 annotation_table 模板
cat("\n--- Creating annotation table template ---\n")

# 获取所有变量
all_variable_ids <- my_mass@variable_info$variable_id

# 创建空的 annotation_table
annotation_table <- data.frame(
  variable_id = all_variable_ids,
  ms2_files_id = NA_character_,
  ms2_spectrum_id = NA_character_,
  Compound.name = NA_character_,
  CAS.ID = NA_character_,
  HMDB.ID = NA_character_,
  KEGG.ID = NA_character_,
  Lab.ID = NA_character_,
  Adduct = "[M+H]+",  # 默认加合离子
  mz.error = NA_real_,
  mz.match.score = NA_real_,
  RT.error = NA_real_,
  RT.match.score = NA_real_,
  CE = NA_character_,
  SS = NA_real_,
  Total.score = NA_real_,
  Database = NA_character_,
  Level = NA_integer_,
  stringsAsFactors = FALSE
)

cat("Annotation table created with", nrow(annotation_table), "rows\n")
print(head(annotation_table, 3))

# 3.2 模拟添加一些注释（实际应该从数据库匹配）
cat("\n--- Adding mock annotations ---\n")

# 随机选择一些变量进行"注释"
set.seed(123)
n_annotated <- min(50, nrow(annotation_table))
annotated_idx <- sample(1:nrow(annotation_table), n_annotated)

# 添加模拟注释
annotation_table$Compound.name[annotated_idx] <- 
  paste0("Metabolite_", sprintf("%03d", 1:n_annotated))

annotation_table$HMDB.ID[annotated_idx] <- 
  paste0("HMDB", sprintf("%07d", sample(1:9999999, n_annotated)))

annotation_table$Level[annotated_idx] <- 
  sample(1:3, n_annotated, replace = TRUE, prob = c(0.2, 0.5, 0.3))

annotation_table$Total.score[annotated_idx] <- 
  round(runif(n_annotated, 0.5, 1.0), 3)

annotation_table$mz.match.score[annotated_idx] <- 
  round(runif(n_annotated, 0.8, 1.0), 3)

annotation_table$Database[annotated_idx] <- "HMDB"

cat("Added annotations for", n_annotated, "variables\n")
cat("\nAnnotation levels:\n")
print(table(annotation_table$Level, useNA = "always"))

# 3.3 将 annotation_table 添加到 mass_dataset
cat("\n--- Adding annotation to mass_dataset ---\n")

# 方法 A: 使用 mutate_annotation_table()（如果该函数可用）
if ("mutate_annotation_table" %in% getNamespaceExports("massdataset")) {
  my_mass_annotated <- mutate_annotation_table(
    object = my_mass,
    annotation_table = annotation_table
  )
  cat("✓ Annotation table added using mutate_annotation_table()\n")
} else {
  # 方法 B: 直接添加注释到 variable_info
  cat("Adding annotations directly to variable_info...\n")
  
  # 合并注释信息到 variable_info
  my_mass@variable_info <- left_join(
    my_mass@variable_info,
    annotation_table,
    by = "variable_id"
  )
  
  my_mass_annotated <- my_mass
  cat("✓ Annotation info merged to variable_info\n")
}

# 3.4 查看注释结果
cat("\n--- Annotation Summary ---\n")
cat("Total variables:", nrow(annotation_table), "\n")
cat("Annotated variables:", sum(!is.na(annotation_table$Compound.name)), "\n")
cat("Level 1 annotations:", sum(annotation_table$Level == 1, na.rm = TRUE), "\n")
cat("Level 2 annotations:", sum(annotation_table$Level == 2, na.rm = TRUE), "\n")
cat("Level 3 annotations:", sum(annotation_table$Level == 3, na.rm = TRUE), "\n")

# 查看注释的变量
annotated_vars <- annotation_table[!is.na(annotation_table$Compound.name), ]
cat("\nExample annotated variables:\n")
print(head(annotated_vars[, c("variable_id", "Compound.name", 
                              "HMDB.ID", "Level", "Total.score")], 5))

# ============================================================================
# 4. 综合示例：过滤 + 可视化
# ============================================================================

cat("\n=== Example 4: Combined Filtering and Visualization ===\n")

# 4.1 只显示高置信度注释的变量
cat("\n--- High confidence annotations (Level 1) ---\n")

if (sum(annotation_table$Level == 1, na.rm = TRUE) > 0) {
  my_mass_annotated %>%
    activate_mass_dataset(what = "variable_info") %>%
    filter(Level == 1) %>%
    show_missing_values(
      cell_color = "white",
      show_row_names = TRUE,
      row_names_side = "left"
    )
}

# 4.2 特定 m/z 范围 + 特定空间区域
cat("\n--- m/z 100-200 in left region ---\n")

if ("spatial_group" %in% colnames(my_mass@sample_info)) {
  my_mass_annotated %>%
    activate_mass_dataset(what = "variable_info") %>%
    filter(mz >= 100 & mz <= 200) %>%
    activate_mass_dataset(what = "sample_info") %>%
    filter(spatial_group == "left") %>%
    show_missing_values(cell_color = "lightgreen")
}

# 4.3 绘制已注释变量的强度
cat("\n--- Intensity plots for annotated variables ---\n")

annotated_indices <- which(!is.na(annotation_table$Compound.name))

if (length(annotated_indices) > 0) {
  # 绘制前 5 个已注释的变量
  for (i in head(annotated_indices, 5)) {
    var_name <- annotation_table$Compound.name[i]
    var_id <- annotation_table$variable_id[i]
    
    cat("\nPlotting:", var_name, "(", var_id, ")\n")
    
    if ("spatial_group" %in% colnames(my_mass_annotated@sample_info)) {
      p <- intensity_plot(
        object = my_mass_annotated,
        variable_index = i,
        color_by = "spatial_group"
      ) + 
        labs(title = paste(var_name, "-", var_id))
      
      print(p)
    }
  }
}

# ============================================================================
# 5. 保存结果
# ============================================================================

cat("\n=== Saving Results ===\n")

# 保存带注释的 mass_dataset
dir.create("2_data/output", showWarnings = FALSE, recursive = TRUE)

saveRDS(my_mass_annotated, "2_data/output/mass_dataset_annotated.rds")
cat("✓ Annotated mass_dataset saved\n")

# 保存 annotation_table
write.csv(
  annotation_table,
  "2_data/output/annotation_table_mass.csv",
  row.names = FALSE
)
cat("✓ Annotation table saved as CSV\n")

# 保存只包含已注释变量的表格
annotated_only <- annotation_table[!is.na(annotation_table$Compound.name), ]
write.csv(
  annotated_only,
  "2_data/output/annotated_variables_only.csv",
  row.names = FALSE
)
cat("✓ Annotated variables list saved\n")

# ============================================================================
# 6. 创建可视化报告
# ============================================================================

cat("\n=== Creating Visualization Report ===\n")

# 保存多个图形到 PDF
pdf("2_data/output/massdataset_plots.pdf", width = 10, height = 8)

# 缺失值总览
cat("Generating missing values plot...\n")
show_missing_values(my_mass_annotated)

# 不同 m/z 范围的缺失值
for (mz_range in list(c(50, 200), c(200, 400), c(400, 600))) {
  cat("Generating plot for m/z", mz_range[1], "-", mz_range[2], "\n")
  
  my_mass_annotated %>%
    activate_mass_dataset(what = "variable_info") %>%
    filter(mz >= mz_range[1] & mz <= mz_range[2]) %>%
    show_missing_values(cell_color = "white")
}

# 已注释变量的强度图
if (length(annotated_indices) > 0) {
  for (i in head(annotated_indices, 10)) {
    if ("spatial_group" %in% colnames(my_mass_annotated@sample_info)) {
      p <- intensity_plot(
        object = my_mass_annotated,
        variable_index = i,
        color_by = "spatial_group"
      )
      print(p)
    }
  }
}

dev.off()
cat("✓ Plots saved to massdataset_plots.pdf\n")

# ============================================================================
# 总结
# ============================================================================

cat("\n========================================\n")
cat("Analysis Complete!\n")
cat("========================================\n\n")

cat("Generated outputs:\n")
cat("  1. mass_dataset_annotated.rds - Annotated mass_dataset object\n")
cat("  2. annotation_table_mass.csv - Complete annotation table\n")
cat("  3. annotated_variables_only.csv - Only annotated variables\n")
cat("  4. massdataset_plots.pdf - Visualization report\n\n")

cat("Summary:\n")
cat("  Total variables:", nrow(my_mass_annotated@variable_info), "\n")
cat("  Total samples:", nrow(my_mass_annotated@sample_info), "\n")
cat("  Annotated variables:", sum(!is.na(annotation_table$Compound.name)), "\n")

if ("spatial_group" %in% colnames(my_mass_annotated@sample_info)) {
  cat("\nSpatial groups:\n")
  print(table(my_mass_annotated@sample_info$spatial_group))
}

cat("\n========================================\n")