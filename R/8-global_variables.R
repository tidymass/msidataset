#' @title Global Variables for msidataset
#' @description This file defines global variables to avoid "no visible binding for global variable" notes during R CMD check.

# Define global variables
msidataset_version <- NULL
update_date <- NULL

globalVariables(names = c(
  # Version information
  "msidataset_version",
  "update_date",
  
  # Base Core Variables
  ".", 
  "value", 
  "name", 
  "index",
  
  # Mass Spectrometry Core
  "mz", 
  "intensity", 
  "variable_id", 
  "sample_id",
  
  # MSI Specific Spatial Variables
  "x", 
  "y", 
  "z", 
  "pixel_id", 
  "subject_id",
  "section_id",
  "resolution",
  "matrix",
  "coordinate_x", 
  "coordinate_y",
  "sum_intensity",
  "mean_intensity",
  "median_intensity",
  "tic", # Total Ion Count
  "bpi", # Base Peak Intensity
  
  # Data Cleaning & Metadata
  "na_number",
  "na_freq",
  "p_value",
  "p_value_adjust",
  "vip",
  "fdr",
  
  # Activation & Slots
  "sample_info",
  "variable_info",
  "expression_data",
  "annotation_table",
  "activated"
))
