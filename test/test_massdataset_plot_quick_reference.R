# ============================================================================
# Massdataset 函数快速参考指南
# ============================================================================

# ============================================================================
# 前提：已有 my_mass 对象
# ============================================================================

# 假设你已经转换好了 mass_dataset
# my_mass <- convert_msi2mass_dataset(my_msi, save_metadata = TRUE)

# ============================================================================
# 1. intensity_plot() - 强度可视化
# ============================================================================

# --- 基础用法 ---
# 绘制第1个变量
intensity_plot(object = my_mass, variable_index = 1)

# 绘制第10个变量
intensity_plot(object = my_mass, variable_index = 10)

# --- 按分组着色 ---
# 前提：sample_info 中需要有分组列（如 class, spatial_group 等）
# 从转换缓存获取空间坐标
cache_env <- get(".msi_conversion_cache", envir = .GlobalEnv)
metadata <- cache_env[[ls(cache_env)[1]]]
pixel_coords <- metadata$pixel_info

# 创建左中右分组
x_coords <- pixel_coords$x[match(
  my_mass@sample_info$sample_id, 
  pixel_coords$pixel_id
)]

quantiles <- quantile(x_coords, probs = c(1/3, 2/3))

my_mass@sample_info$spatial_group <- cut(
  x_coords,
  breaks = c(-Inf, quantiles[1], quantiles[2], Inf),
  labels = c("left", "center", "right")
)

# 现在可以使用 color_by
intensity_plot(
  object = my_mass, 
  variable_index = 1, 
  color_by = "spatial_group"
)


# 按 class 着色
intensity_plot(
  object = my_mass, 
  variable_index = 1, 
  color_by = "class"
)

# 按空间分组着色
intensity_plot(
  object = my_mass, 
  variable_index = 1, 
  color_by = "spatial_group"
)

# --- 找特定 m/z 的变量并绘图 ---
target_mz <- 200.5
mz_values <- my_mass@variable_info$mz
closest_idx <- which.min(abs(mz_values - target_mz))

intensity_plot(object = my_mass, variable_index = closest_idx)

# --- 批量绘制多个变量 ---
for (i in 1:10) {
  p <- intensity_plot(
    object = my_mass, 
    variable_index = i, 
    color_by = "spatial_group"
  )
  print(p)
}

# ============================================================================
# 2. show_missing_values() - 缺失值可视化
# ============================================================================

# --- 基础用法 ---
# 显示所有缺失值
show_missing_values(object = my_mass)

# 自定义颜色
show_missing_values(object = my_mass, cell_color = "white")

# --- 只显示部分变量 ---
# 显示前20个变量
show_missing_values(object = my_mass[1:20, ])

# --- 按条件过滤样本 ---
# 只显示特定 class 的样本
my_mass %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class == "pixel") %>%
  show_missing_values()

# 只显示特定空间区域
my_mass %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(spatial_group == "left") %>%
  show_missing_values()

# --- 选择特定样本列 ---
# 选择前10个样本
my_mass %>%
  activate_mass_dataset(what = "expression_data") %>%
  dplyr::select(1:10) %>%
  show_missing_values()

# --- 按 m/z 范围过滤变量 ---
# 只显示 m/z < 200 的变量
my_mass %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(mz < 200) %>%
  show_missing_values(
    cell_color = "white",
    show_row_names = TRUE,
    row_names_side = "left"
  )

# 显示特定 m/z 范围
my_mass %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(mz >= 100 & mz <= 300) %>%
  show_missing_values()

# --- 显示百分比和设置阈值 ---
my_mass %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(mz >= 100 & mz <= 300) %>%
  show_missing_values(
    cell_color = "white",
    show_row_names = TRUE,
    row_names_side = "left",
    percentage = TRUE,           # 显示百分比
    sample_na_cutoff = 50,       # 样本缺失≥50%的显示
    variable_na_cutoff = 20      # 变量缺失≥20%的显示
  )

# ============================================================================
# 3. 创建空间分组（用于着色）
# ============================================================================

# --- 从转换缓存中获取空间坐标 ---
if (exists(".msi_conversion_cache", envir = .GlobalEnv)) {
  cache_env <- get(".msi_conversion_cache", envir = .GlobalEnv)
  cache_ids <- ls(cache_env)
  
  if (length(cache_ids) > 0) {
    metadata <- cache_env[[cache_ids[1]]]
    pixel_coords <- metadata$pixel_info
    
    # 获取 x 坐标
    x_coords <- pixel_coords$x[match(
      my_mass@sample_info$sample_id, 
      pixel_coords$pixel_id
    )]
    
    # 创建左中右分组
    quantiles <- quantile(x_coords, probs = c(1/3, 2/3), na.rm = TRUE)
    
    my_mass@sample_info$spatial_group <- cut(
      x_coords,
      breaks = c(-Inf, quantiles[1], quantiles[2], Inf),
      labels = c("left", "center", "right")
    )
    
    # 查看分组
    table(my_mass@sample_info$spatial_group)
  }
}

# --- 创建 y 坐标分组 ---
if (exists("pixel_coords")) {
  y_coords <- pixel_coords$y[match(
    my_mass@sample_info$sample_id, 
    pixel_coords$pixel_id
  )]
  
  y_quantiles <- quantile(y_coords, probs = c(1/3, 2/3), na.rm = TRUE)
  
  my_mass@sample_info$y_group <- cut(
    y_coords,
    breaks = c(-Inf, y_quantiles[1], y_quantiles[2], Inf),
    labels = c("bottom", "middle", "top")
  )
}

# ============================================================================
# 4. 创建 annotation_table
# ============================================================================

# --- 创建空模板 ---
annotation_table <- data.frame(
  variable_id = my_mass@variable_info$variable_id,
  ms2_files_id = NA_character_,
  ms2_spectrum_id = NA_character_,
  Compound.name = NA_character_,
  CAS.ID = NA_character_,
  HMDB.ID = NA_character_,
  KEGG.ID = NA_character_,
  Lab.ID = NA_character_,
  Adduct = "[M+H]+",
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

# --- 添加示例注释 ---
# 随机选择50个变量添加注释
set.seed(123)
n_annotated <- min(50, nrow(annotation_table))
annotated_idx <- sample(1:nrow(annotation_table), n_annotated)

annotation_table$Compound.name[annotated_idx] <- 
  paste0("Metabolite_", sprintf("%03d", 1:n_annotated))

annotation_table$HMDB.ID[annotated_idx] <- 
  paste0("HMDB", sprintf("%07d", sample(1:9999999, n_annotated)))

annotation_table$Level[annotated_idx] <- 
  sample(1:3, n_annotated, replace = TRUE)

annotation_table$Total.score[annotated_idx] <- 
  round(runif(n_annotated, 0.5, 1.0), 3)

# --- 将注释添加到 variable_info ---
library(dplyr)
my_mass@variable_info <- left_join(
  my_mass@variable_info,
  annotation_table,
  by = "variable_id"
)

# ============================================================================
# 5. 组合使用示例
# ============================================================================

# --- 只显示已注释且高置信度的变量 ---
my_mass %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(!is.na(Compound.name) & Level == 1) %>%
  show_missing_values(cell_color = "white")

# --- 特定 m/z 范围 + 特定区域 ---
my_mass %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(mz >= 100 & mz <= 200) %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(spatial_group == "left") %>%
  show_missing_values()

# --- 绘制已注释变量的强度 ---
annotated_vars <- which(!is.na(annotation_table$Compound.name))

for (i in head(annotated_vars, 5)) {
  p <- intensity_plot(
    object = my_mass,
    variable_index = i,
    color_by = "spatial_group"
  )
  print(p)
}

# ============================================================================
# 6. 保存结果
# ============================================================================

# 保存带注释的对象
saveRDS(my_mass, "2_data/output/mass_dataset_annotated.rds")

# 保存注释表
write.csv(
  annotation_table,
  "2_data/output/annotation_table.csv",
  row.names = FALSE
)

# 保存图形到 PDF
pdf("2_data/output/plots.pdf", width = 10, height = 8)

# 缺失值图
show_missing_values(my_mass)

# 强度图
for (i in 1:10) {
  p <- intensity_plot(my_mass, variable_index = i, color_by = "spatial_group")
  print(p)
}

dev.off()

# ============================================================================
# 常见问题
# ============================================================================

# Q1: 如何找到特定 m/z 的变量索引？
# A1:
target_mz <- 100.0
idx <- which.min(abs(my_mass@variable_info$mz - target_mz))

# Q2: 如何查看哪些样本有 spatial_group？
# A2:
if ("spatial_group" %in% colnames(my_mass@sample_info)) {
  table(my_mass@sample_info$spatial_group)
}

# Q3: 如何统计注释情况？
# A3:
cat("Total variables:", nrow(annotation_table), "\n")
cat("Annotated:", sum(!is.na(annotation_table$Compound.name)), "\n")
cat("Level 1:", sum(annotation_table$Level == 1, na.rm = TRUE), "\n")
cat("Level 2:", sum(annotation_table$Level == 2, na.rm = TRUE), "\n")
cat("Level 3:", sum(annotation_table$Level == 3, na.rm = TRUE), "\n")

# Q4: 如何只提取已注释的变量？
# A4:
annotated_only <- annotation_table[!is.na(annotation_table$Compound.name), ]

# Q5: 如何按注释level过滤？
# A5:
my_mass %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(Level == 1) %>%
  show_missing_values()