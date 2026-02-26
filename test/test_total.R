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

# 3. intensity (> 3%)
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


# # 1. 查看处理信息摘要（简洁版）
show_process_info(msi_filtered, detailed = FALSE, verbose = FALSE)

# # 2. 查看完整的处理信息（推荐）
show_process_info(msi_filtered, detailed = TRUE, verbose = TRUE)

# # 3. 查看处理时间线
show_process_timeline(msi_filtered)

# # 4. 导出处理历史为表格
process_history <- export_process_history(msi_filtered, include_details = TRUE)
print(process_history)

# # 5. 保存处理历史到CSV
write.csv(process_history, "msi_processing_history.csv", row.names = FALSE)

print()

## save
if (!dir.exists("data/output")) dir.create("data/output", recursive = TRUE)
saveRDS(my_msi_full, "data/output/my_msi_dataset_full.rds")
cat("Saved to data/output/my_msi_dataset_full.rds\n")

## reload (confirms file integrity; use this object from here on)
my_msi <- readRDS("data/output/my_msi_dataset_full.rds")
cat("Reloaded dim:", paste(dim(my_msi), collapse = " x "), "\n")

