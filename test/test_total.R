# ============================================================================
# msidataset — Complete Debug / Demo Script
# Author : Xiaoyu Xu  <XIAOYU021@e.ntu.edu.sg>
# Purpose: Full workflow test using real .imzML data via Cardinal.
#
# Workflow:
#   Step 1 : create
#   Step 2 : convert [Read .imzML -> Cardinal -> msi_dataset]
#   Step 3 : filter
#   Step 4 : show
#   Step 5 : 
#   Step 6 : 
#   Step 7 : 
#   Step 8 : 
#   Step 9 : 
#   Step 10: 
#   Step 11: Save & reload final check
# ============================================================================

# ----------------------------------------------------------------------------
# 0. LOAD / BUILD PACKAGE
# ----------------------------------------------------------------------------

setwd("D:/Shen_lab/msidataset")
devtools::document()
devtools::load_all()

# library(msidataset)   # use this after package is installed

library(S4Vectors)
library(Cardinal)
library(massdataset)
library(dplyr)
library(crayon)


# ============================================================================
# STEP 1 - Read real .imzML and convert to msi_dataset
# ============================================================================

## read data
cardinal_obj <- readMSIData(
  "data/demo data/03-total ion count.imzML",
  memory = FALSE
)

# print(cardinal_obj)
# print(experimentData)

## aligned
cardinal_obj_aligned <- Cardinal::peakAlign(
  cardinal_obj,
  tolerance = 50,
  units = "ppm"
)

my_msi_full <- convert_cardinal2msi_dataset(
  cardinal_object = cardinal_obj_aligned,
  sample_name     = "rcc_1",
  section_name    = "tumor",
  resolution      = 50,
  matrix_name     = "NEDC"
)

## add notes
my_msi_full <- update_msi_notes(my_msi_full)
View(my_msi_full@sample_info_note)  
View(my_msi_full@variable_info_note)

## filter (sample & variable)

## sample filter
# 1. filter_spatial()
msi_spatial <- filter_msi_dataset(
  my_msi_full, target = "sample_info", filter_type = "spatial",
  x_min = 2, x_max = 8, y_min = 1, y_max = 5
)

# 2. filter_by_section()
msi_section <- filter_msi_dataset(
  my_msi_full, target = "sample_info", filter_type = "section",
  section_ids = c("slice_01")
)

# 3. filter_by_run()
msi_run <- filter_msi_dataset(
  my_msi_full, target = "sample_info", filter_type = "run",
  run_ids = c("run1", "run2")
)

# 4. filter_by_intensity()
msi_intensity <- filter_msi_dataset(
  my_msi_full, target = "sample_info", filter_type = "intensity",
  intensity_type = "median", intensity_threshold = 10
)

# 5. filter_samples_combined()
msi_samples_combined <- filter_msi_dataset(
  my_msi_full, target = "sample_info", filter_type = "combined",
  spatial = list(x_min = 2, x_max = 8, y_min = 1, y_max = 5),
  section = list(section_ids = c("slice_01")),
  intensity = list(intensity_type = "median", intensity_threshold = 10)
)

## variable filter
# 1. filter_mz_range()
msi_mz_range <- filter_msi_dataset(
  my_msi_full, target = "variable_info", filter_type = "mz_range",
  mz_min = 300, mz_max = 500
)

# 2. filter_by_frequency()
msi_freq <- filter_msi_dataset(
  my_msi_full, target = "variable_info", filter_type = "frequency",
  freq_min = 0.4,   
  freq_max = 1.0    
)

# 3. filter_by_count()
msi_count <- filter_msi_dataset(
  my_msi_full, target = "variable_info", filter_type = "count",
  count_min = 5,    
  count_max = 54    
)

# 4. filter_variable_intensity()
msi <- filter_msi_dataset(
  my_msi_full, target = "variable_info", filter_type = "intensity",
  intensity_type = "median", intensity_threshold = 50
)

## **filter supplementary**
show_mz_statistics(my_msi_full)

params <- recommend_filter_params(my_msi_full)


# 1. high freq mz (by non-zero count)
msi_high_freq <- filter_msi_dataset_supplementary(
  my_msi_full, filter_type = 'high_frequency_count', min_pixels = 17)

print(msi_high_freq)

# 2. mz with excessive zeros
msi_no_zeros <- filter_msi_dataset_supplementary(
  msi_high_freq, filter_type = 'zero_heavy', zero_threshold = 0.5)

print(msi_no_zeros)

# 3. intensity (> 3%)?
msi_intensity <- filter_msi_dataset(
  msi_no_zeros, target = "sample_info", filter_type = "intensity",
  intensity_type = "median", intensity_threshold = 90
)

print(msi_intensity)

# 4. combined
msi_filtered <- filter_msi_dataset_supplementary(
  my_msi_full, filter_type = 'supplementary',
  min_pixels = 17, zero_threshold = 0.5)

print(msi_filtered)

## plot

# EXAMPLE 1
plot_msi_filtered(
  msi_filtered,           # 筛选后的msi_dataset对象
  variable_id = NULL,     # 变量ID，如"var_1"
  mz = 885.55,              # m/z值，如885.55
  mz_tolerance = 0.1,    # m/z匹配的容差
  title = NULL,           # 自定义标题
  color_scale = "viridis",# 颜色方案：viridis、magma、plasma、inferno、cividis
  transform = "none",     # 数据转换：none、log10、sqrt
  add_grid = TRUE         # 是否添加网格线
)

# BASIC
p1 <- plot_msi_filtered(msi_filtered, variable_id = "var_1")
print(p1)

# CUSTOM
p2 <- plot_msi_filtered(
  msi_filtered, 
  mz = 600.5, 
  color_scale = "magma",
  transform = "log10",
  title = "Ion m/z 600.5 Distribution"
)
print(p2)

# SAVE PDF
pdf("ion_distribution.pdf", width = 8, height = 6)
print(p2)
dev.off()

# EXAMPLE 2 plot_msi_grid()
plot_msi_grid(
  msi_filtered,              # 筛选后的msi_dataset对象
  variable_ids,              # 要绘制的变量ID向量
  ncol = 3,                  # 网格列数
  color_scale = "viridis",   # 所有图的颜色方案
  transform = "none"         # 所有图的数据转换
)

# 获取前6个变量ID
var_info <- variable_info(msi_filtered)
top_6_vars <- head(var_info$variable_id, 6)

# 创建2x3网格
plot_msi_grid(msi_filtered, top_6_vars, ncol = 3)

# 自定义参数
plot_msi_grid(
  msi_filtered,
  top_6_vars,
  ncol = 2,
  color_scale = "plasma",
  transform = "log10"
)

# 保存为PDF
pdf("ion_grid.pdf", width = 12, height = 8)
plot_msi_grid(msi_filtered, top_6_vars, ncol = 3)
dev.off()

# EXAMPLE 3 
plot_spatial_correlation(
  msi_filtered,         # 筛选后的msi_dataset对象
  n_vars = 20,          # 分析的顶部变量数
  method = "pearson",   # 相关性方法：pearson、spearman、kendall
  show_labels = FALSE,  # 是否显示标签
  label_top_n = 5       # 标签的顶部变量数
)

# 基本用法
result <- plot_spatial_correlation(msi_filtered)

# 查看X相关性图
print(result$plot_x)

# 查看Y相关性图
print(result$plot_y)

# 查看相关性数据
head(result$data)

# 查看摘要统计
summary(result)

# 高级用法：Spearman相关性，带标签
result2 <- plot_spatial_correlation(
  msi_filtered,
  n_vars = 50,
  method = "spearman",
  show_labels = TRUE,
  label_top_n = 10
)

# 保存两个图
pdf("spatial_correlation.pdf", width = 14, height = 6)
gridExtra::grid.arrange(result2$plot_x, result2$plot_y, ncol = 2)
dev.off()



## save
if (!dir.exists("data/output")) dir.create("data/output", recursive = TRUE)
saveRDS(my_msi_full, "data/output/my_msi_dataset_full.rds")
cat("Saved to data/output/my_msi_dataset_full.rds\n")

## reload (confirms file integrity; use this object from here on)
my_msi <- readRDS("data/output/my_msi_dataset_full.rds")
cat("Reloaded dim:", paste(dim(my_msi), collapse = " x "), "\n")

