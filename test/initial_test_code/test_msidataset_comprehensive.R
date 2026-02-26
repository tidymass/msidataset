# ============================================================================
# Comprehensive Test Suite for msidataset Package
# Author: Xiaoyu Xu
# Purpose: Test conversion functions and massdataset integration
# ============================================================================

# ============================================================================
# Comprehensive Test Suite for msidataset Package
# ============================================================================

# 辅助函数：生成模拟数据 (保留在脚本内供测试使用)
generate_test_msi_data <- function(n_variables = 50, n_pixels_x = 10, n_pixels_y = 10) {
  n_pixels_total <- n_pixels_x * n_pixels_y
  
  variable_info <- data.frame(
    variable_id = paste0("var_", 1:n_variables),
    mz = seq(100, 1000, length.out = n_variables),
    stringsAsFactors = FALSE
  )
  
  pixel_info <- data.frame(
    pixel_id = paste0("pixel_", 1:n_pixels_total),
    x = rep(1:n_pixels_x, times = n_pixels_y),
    y = rep(1:n_pixels_y, each = n_pixels_x),
    z = 1,
    stringsAsFactors = FALSE
  )
  
  sample_info <- data.frame(
    sample_id = "sample_1",
    stringsAsFactors = FALSE
  )
  
  expression_data <- matrix(runif(n_variables * n_pixels_total, 0, 1000), 
                            nrow = n_variables)
  rownames(expression_data) <- variable_info$variable_id
  colnames(expression_data) <- pixel_info$pixel_id
  
  # 调用你包里的构造函数
  create_msi_dataset(
    expression_data = as.data.frame(expression_data),
    pixel_info = pixel_info,
    sample_info = sample_info,
    variable_info = variable_info
  )
}

# ============================================================================
# 实际测试用例 (testthat 会自动寻找这些 test_that 块)
# ============================================================================

test_that("Conversion: MSI to mass_dataset roundtrip works", {
  # 1. 准备数据
  msi_original <- generate_test_msi_data(n_variables = 30, n_pixels_x = 5, n_pixels_y = 5)
  
  # 2. 测试 MSI -> Mass
  mass_obj <- convert_msi2mass_dataset(msi_original, save_metadata = TRUE)
  
  expect_s4_class(mass_obj, "mass_dataset")
  expect_equal(nrow(mass_obj@expression_data), 30)
  expect_equal(ncol(mass_obj@expression_data), 25)
  expect_true("x" %in% colnames(mass_obj@sample_info))
  
  # 3. 测试 Mass -> MSI
  msi_restored <- convert_mass2msi_dataset(mass_obj, metadata_source = "auto")
  
  expect_s4_class(msi_restored, "msi_dataset")
  expect_equal(n_variables(msi_restored), n_variables(msi_original))
  expect_equal(mz(msi_restored), mz(msi_original), tolerance = 1e-6)
})

test_that("Integration: massdataset processing works", {
  # 只有安装了 massdataset 才跑这个测试
  skip_if_not_installed("massdataset")
  
  msi_obj <- generate_test_msi_data(n_variables = 100, n_pixels_x = 10, n_pixels_y = 10)
  mass_obj <- convert_msi2mass_dataset(msi_obj)
  
  # 使用 massdataset 的功能进行过滤
  # 注意：明确指定命名空间以防 filter 冲突
  mass_filtered <- mass_obj %>%
    massdataset::activate_mass_dataset(what = "variable_info") %>%
    massdataset::filter(mz > 500)
  
  expect_lt(nrow(mass_filtered@variable_info), nrow(mass_obj@variable_info))
  expect_true(all(mass_filtered@variable_info$mz > 500))
})

test_that("Edge Cases: Minimal and specific datasets", {
  # 极小数据集测试
  msi_tiny <- generate_test_msi_data(n_variables = 1, n_pixels_x = 2, n_pixels_y = 2)
  expect_silent(mass_tiny <- convert_msi2mass_dataset(msi_tiny))
  expect_equal(ncol(mass_tiny@expression_data), 4)
  
  # 缓存缺失报错测试
  # 清理可能存在的缓存环境
  if (exists(".msi_conversion_cache", envir = .GlobalEnv)) {
    rm(".msi_conversion_cache", envir = .GlobalEnv)
  }
  
  # 故意不提供缓存去转换，看是否报预期的错误
  # 注意：这里假设你的函数在找不到缓存时会 stop()
  expect_error(convert_mass2msi_dataset(mass_tiny, metadata_source = "env"))
})