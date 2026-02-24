#' @title Convert Cardinal Object to msi_dataset (Fixed Version)
#' @description 
#' 将 Cardinal 对象转换为 msi_dataset。
#' 修复了 iData 弃用警告和 generate_variable_id 缺失问题。
#'
#' @param cardinal_object Cardinal 镜像对象
#' @param sample_name 样本名称前缀
#' @param section_name 切片名称
#' @param resolution 空间分辨率 (um)
#' @param matrix_name 基质名称
#' @export

convert_cardinal2msi_dataset <- function(cardinal_object,
                                         sample_name = "sample_1",
                                         section_name = "section_1",
                                         resolution = NA,
                                         matrix_name = NA) {
  
  # 1. 检查依赖
  if (!requireNamespace("Cardinal", quietly = TRUE)) {
    stop("请安装 Cardinal 包: BiocManager::install('Cardinal')")
  }
  
  # 2. 提取强度矩阵 (Expression Data)
  # 解决警告：使用 spectra() 替代 iData()
  # as.matrix 确保数据从磁盘加载到内存中，符合 massdataset 的矩阵要求
  cat(crayon::blue("正在提取强度矩阵...\n"))
  expression_data <- as.matrix(Cardinal::spectra(cardinal_object))
  
  # 3. 提取像素信息 (Sample Info)
  # 使用 pData 避开 colData 的 S4 分发问题
  sample_info <- as.data.frame(Cardinal::pData(cardinal_object))
  
  # 统一 ID 命名，确保与 massdataset 风格一致
  sample_info$sample_id <- paste0(sample_name, "_pixel_", seq_len(nrow(sample_info)))
  sample_info$section_id <- section_name
  if (!is.na(resolution)) sample_info$resolution <- resolution
  if (!is.na(matrix_name)) sample_info$matrix <- matrix_name
  
  # 4. 提取变量信息 (Variable Info)
  variable_info <- as.data.frame(Cardinal::fData(cardinal_object))
  
  # 确保有 mz 列
  if (!"mz" %in% colnames(variable_info)) {
    variable_info$mz <- Cardinal::mz(cardinal_object)
  }
  
  # 解决报错：手动定义或调用 generate_variable_id
  # 这里假设你已经加载了 generate_variable_id.R 里的函数
  if (!exists("generate_variable_id")) {
    stop("未找到 generate_variable_id 函数。请确保已 source('generate_variable_id.R')")
  }
  
  variable_info$variable_id <- generate_variable_id(variable_info$mz)
  
  # 5. 同步行列名
  rownames(expression_data) <- as.character(variable_info$variable_id)
  colnames(expression_data) <- as.character(sample_info$sample_id)
  
  # 6. 创建 msi_dataset 对象
  # 这里调用你在 0-msi_dataset_class.R 中定义的 create_msi_dataset
  object <- create_msi_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )
  
  # 记录处理日志
  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "msidataset",
    function_name = "convert_cardinal2msi_dataset()",
    parameter = list(
      sample_name = sample_name,
      section_name = section_name,
      cardinal_class = class(cardinal_object)[1]
    ),
    time = Sys.time()
  )
  object@process_info$convert_cardinal2msi_dataset <- parameter
  
  cat(crayon::green("转换成功！\n"))
  return(object)
}