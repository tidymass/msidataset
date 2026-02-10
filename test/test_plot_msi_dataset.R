# Plot first variable
plot_ion_spatial(
  msi_test,
  variable_id = "var_1",
  color_scale = "viridis",
  transform = "none"
)

# Plot by m/z value
var_mz <- mz(msi_test)[10]
plot_ion_spatial(
  msi_test,
  mz = var_mz,
  mz_tolerance = 0.1,
  color_scale = "magma",
  transform = "sqrt"
)


## Example 4: Multiple ions grid view
## --------------------------------------------------------

# Select 6 interesting variables
selected_vars <- paste0("var_", c(5, 15, 25, 35, 45, 55))

plot_ion_grid(
  msi_test,
  variable_ids = selected_vars,
  ncol = 3,
  color_scale = "plasma",
  transform = "log10"
)


## Example 5: RGB overlay
## --------------------------------------------------------

plot_rgb_overlay(
  msi_test,
  red_var = "var_10",
  green_var = "var_50",
  blue_var = "var_90",
  title = "Three Ion Overlay"
)


## Example 6: Quality control - TIC spatial map
## --------------------------------------------------------

plot_tic_spatial(
  msi_test,
  transform = "log10",
  color_scale = "viridis"
)


## Example 7: Peak count spatial distribution
## --------------------------------------------------------

plot_peak_count_spatial(
  msi_test,
  threshold = 100
)


# ============================================================================
# PART 4: Spatial Statistical Analysis
# ============================================================================

## Example 8: Spatial autocorrelation (Moran's I)
## --------------------------------------------------------

# Test spatial autocorrelation for a specific variable
morans_result <- calculate_morans_i(
  msi_test,
  variable_id = "var_15",
  distance_threshold = 3
)

# Print results
print(morans_result)

# Test multiple variables
autocorr_results <- data.frame(
  variable_id = character(),
  morans_i = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:20) {
  var_id <- paste0("var_", i * 10)
  result <- calculate_morans_i(msi_test, var_id)
  
  autocorr_results <- rbind(
    autocorr_results,
    data.frame(
      variable_id = var_id,
      morans_i = result$morans_i,
      p_value = result$p_value,
      stringsAsFactors = FALSE
    )
  )
}

# Visualize
ggplot(autocorr_results, aes(x = variable_id, y = morans_i, fill = p_value < 0.05)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Spatial Autocorrelation Across Variables",
    x = "Variable ID",
    y = "Moran's I",
    fill = "Significant"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## Example 9: Spatial clustering
## --------------------------------------------------------

# Perform k-means clustering
msi_clustered <- spatial_clustering(
  msi_test,
  n_clusters = 5,
  method = "kmeans",
  plot = TRUE
)

# Access cluster assignments
table(pixel_info(msi_clustered)$cluster)


## Example 10: ROI comparison
## --------------------------------------------------------

# Define two regions of interest
roi_results <- compare_roi(
  msi_test,
  roi1_condition = "x <= 15",
  roi2_condition = "x > 15",
  roi1_name = "Left Half",
  roi2_name = "Right Half",
  test_method = "wilcox.test",
  top_n = 20
)

# View significant variables
significant_vars <- roi_results[roi_results$p_adjusted < 0.05, ]
cat("Number of significant variables:", nrow(significant_vars), "\n")

# Plot top differential variable
if (nrow(significant_vars) > 0) {
  top_var <- significant_vars$variable_id[1]
  plot_ion_spatial(msi_test, variable_id = top_var)
}


## Example 11: Spatial PCA
## --------------------------------------------------------

pca_results <- plot_spatial_pca(
  msi_clustered,
  n_components = 2,
  color_by = "cluster"
)

# Access PCA scores
pc_scores <- pca_results$scores
variance_explained <- pca_results$variance_explained

cat("Variance explained by PC1:", variance_explained[1], "%\n")
cat("Variance explained by PC2:", variance_explained[2], "%\n")


# ============================================================================
# PART 5: Conversion to massdataset and Analysis
# ============================================================================

## Example 12: Convert to mass_dataset
## --------------------------------------------------------

source("conversion_functions_v2.R")

# Enhanced conversion with metadata preservation
mass_obj <- convert_msi2mass_dataset_v2(
  msi_test,
  save_metadata = TRUE,
  add_spatial_to_sample_info = TRUE,
  preserve_pixel_metadata = TRUE,
  preserve_variable_metadata = TRUE
)

# Check conversion
print(mass_obj)
dim(mass_obj@expression_data)
head(mass_obj@sample_info)
head(mass_obj@variable_info)


## Example 13: Data filtering with massdataset
## --------------------------------------------------------

# Filter by m/z range
mass_filtered <- mass_obj %>%
  activate_mass_dataset("variable_info") %>%
  filter(mz > 300 & mz < 700)

cat("Original variables:", nrow(mass_obj@variable_info), "\n")
cat("Filtered variables:", nrow(mass_filtered@variable_info), "\n")

# Spatial filtering
mass_spatial_filtered <- mass_obj %>%
  activate_mass_dataset("sample_info") %>%
  filter(x >= 10 & x <= 20 & y >= 10 & y <= 20)

cat("Original pixels:", nrow(mass_obj@sample_info), "\n")
cat("Filtered pixels:", nrow(mass_spatial_filtered@sample_info), "\n")


## Example 14: Normalization
## --------------------------------------------------------

# TIC normalization
mass_normalized <- normalize_data(
  mass_obj,
  method = "sum"
)

# Compare TIC before and after
tic_before <- colSums(mass_obj@expression_data)
tic_after <- colSums(mass_normalized@expression_data)

par(mfrow = c(1, 2))
hist(tic_before, main = "TIC Before Normalization", xlab = "TIC")
hist(tic_after, main = "TIC After Normalization", xlab = "TIC")
par(mfrow = c(1, 1))


## Example 15: Missing value imputation
## --------------------------------------------------------

# Introduce missing values for demonstration
mass_with_na <- mass_obj
na_count <- 500
na_indices <- sample(1:length(mass_obj@expression_data), na_count)
mass_with_na@expression_data[na_indices] <- NA

cat("Missing values introduced:", sum(is.na(mass_with_na@expression_data)), "\n")

# Impute using KNN
mass_imputed <- impute_mv(
  mass_with_na,
  method = "knn"
)

cat("Missing values after imputation:", 
    sum(is.na(mass_imputed@expression_data)), "\n")


## Example 16: Statistical testing
## --------------------------------------------------------

# Add group information based on spatial location
mass_obj@sample_info$spatial_group <- ifelse(
  mass_obj@sample_info$x <= 15,
  "left",
  "right"
)

# Perform t-test
mass_tested <- perform_t_test(
  mass_obj,
  group = "spatial_group",
  group1 = "left",
  group2 = "right"
)

# View results
sig_vars <- mass_tested@variable_info %>%
  filter(p_value < 0.05) %>%
  arrange(p_value)

cat("Significant variables:", nrow(sig_vars), "\n")
print(head(sig_vars))

# Create volcano plot
ggplot(
  mass_tested@variable_info,
  aes(x = fold_change, y = -log10(p_value))
) +
  geom_point(
    aes(color = p_value < 0.05 & abs(fold_change) > 1.5),
    alpha = 0.6
  ) +
  scale_color_manual(values = c("gray", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(1/1.5, 1.5), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Left vs Right",
    x = "Fold Change",
    y = "-log10(P-value)",
    color = "Significant"
  )


## Example 17: PCA analysis
## --------------------------------------------------------

# Perform PCA
mass_pca <- perform_pca(mass_obj)

# Plot PCA
plot_pca_score(
  mass_pca,
  color_by = "spatial_group"
)

plot_pca_loading(
  mass_pca,
  pc_x = 1,
  pc_y = 2
)


# ============================================================================
# PART 6: Convert back to msi_dataset
# ============================================================================

## Example 18: Add annotations and convert back
## --------------------------------------------------------

# Add some annotations (e.g., from database matching)
mass_obj@variable_info$lipid_class <- sample(
  c("PC", "PE", "SM", "TAG", "Unknown"),
  nrow(mass_obj@variable_info),
  replace = TRUE
)

mass_obj@variable_info$annotation_confidence <- sample(
  c("High", "Medium", "Low"),
  nrow(mass_obj@variable_info),
  replace = TRUE
)

# Convert back to msi_dataset
msi_restored <- convert_mass2msi_dataset_v2(
  mass_obj,
  metadata_source = "auto",
  create_annotation_table = TRUE
)

# Check annotations
print(msi_restored@annotation_table)
head(msi_restored@annotation_table)


## Example 19: Verify roundtrip conversion
## --------------------------------------------------------

# Compare original and restored
cat("Original msi_dataset:\n")
cat("  Variables:", n_variables(msi_test), "\n")
cat("  Pixels:", n_pixels(msi_test), "\n")

cat("\nRestored msi_dataset:\n")
cat("  Variables:", n_variables(msi_restored), "\n")
cat("  Pixels:", n_pixels(msi_restored), "\n")

# Check data integrity
max_diff <- max(abs(
  as.matrix(expression_data(msi_test)) - 
    as.matrix(expression_data(msi_restored))
))
cat("\nMaximum difference in expression data:", max_diff, "\n")

# Check spatial coordinates
coord_match <- all(
  pixel_info(msi_test)$x == pixel_info(msi_restored)$x &
    pixel_info(msi_test)$y == pixel_info(msi_restored)$y
)
cat("Spatial coordinates match:", coord_match, "\n")


# ============================================================================
# PART 7: Advanced Workflows
# ============================================================================

## Example 20: Complete analysis pipeline
## --------------------------------------------------------

complete_msi_analysis <- function(msi_object, output_dir = "msi_analysis_results") {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  cat("\n=== Starting Complete MSI Analysis Pipeline ===\n\n")
  
  # 1. Quality Control
  cat("Step 1: Quality Control\n")
  
  pdf(file.path(output_dir, "01_qc_plots.pdf"), width = 10, height = 8)
  
  print(plot_tic_spatial(msi_object, transform = "log10"))
  print(plot_peak_count_spatial(msi_object, threshold = 100))
  
  dev.off()
  
  # 2. Spatial Distribution of Top Variables
  cat("Step 2: Spatial Distribution Analysis\n")
  
  pdf(file.path(output_dir, "02_spatial_distributions.pdf"), width = 12, height = 10)
  
  top_vars <- paste0("var_", seq(1, 100, by = 10))
  plot_ion_grid(msi_object, top_vars, ncol = 3, transform = "log10")
  
  dev.off()
  
  # 3. Spatial Clustering
  cat("Step 3: Spatial Clustering\n")
  
  pdf(file.path(output_dir, "03_spatial_clustering.pdf"), width = 10, height = 8)
  
  msi_clustered <- spatial_clustering(
    msi_object,
    n_clusters = 5,
    method = "kmeans",
    plot = TRUE
  )
  
  dev.off()
  
  # 4. Convert to mass_dataset and perform statistical analysis
  cat("Step 4: Statistical Analysis\n")
  
  mass_obj <- convert_msi2mass_dataset_v2(
    msi_clustered,
    save_metadata = TRUE,
    add_spatial_to_sample_info = TRUE
  )
  
  # Normalize
  mass_norm <- normalize_data(mass_obj, method = "sum")
  
  # Add group based on cluster
  mass_norm@sample_info$group <- pixel_info(msi_clustered)$cluster
  
  # Save results
  saveRDS(mass_norm, file.path(output_dir, "mass_dataset_normalized.rds"))
  
  # 5. Export summary
  cat("Step 5: Generating Summary Report\n")
  
  summary_text <- sprintf(
    "MSI Analysis Summary
====================

Input Data:
- Variables: %d
- Pixels: %d (%d × %d)
- m/z range: %.2f - %.2f

Analysis Steps:
1. Quality control plots generated
2. Spatial distribution visualized for %d variables
3. Spatial clustering performed (k=%d)
4. Data normalized using TIC normalization
5. Mass dataset saved for further analysis

Output Files:
- 01_qc_plots.pdf
- 02_spatial_distributions.pdf
- 03_spatial_clustering.pdf
- mass_dataset_normalized.rds

Completed: %s
",
    n_variables(msi_object),
    n_pixels(msi_object),
    max(pixel_info(msi_object)$x),
    max(pixel_info(msi_object)$y),
    min(mz(msi_object)),
    max(mz(msi_object)),
    length(top_vars),
    5,
    Sys.time()
  )
  
  writeLines(summary_text, file.path(output_dir, "ANALYSIS_SUMMARY.txt"))
  
  cat("\n=== Analysis Pipeline Complete ===\n")
  cat("Results saved to:", output_dir, "\n\n")
  
  return(list(
    msi_clustered = msi_clustered,
    mass_normalized = mass_norm
  ))
}

# Run complete analysis
# results <- complete_msi_analysis(msi_test, output_dir = "my_msi_analysis")


# ============================================================================
# PART 8: Best Practices and Tips
# ============================================================================

cat("
================================================================================
BEST PRACTICES FOR USING msidataset
================================================================================

1. DATA IMPORT
   - Always validate your data after import
   - Check for missing values and outliers
   - Verify spatial coordinates are correct

2. QUALITY CONTROL
   - Always inspect TIC spatial distribution
   - Check for regional signal dropouts
   - Monitor peak count variations

3. NORMALIZATION
   - TIC normalization is recommended for MSI data
   - Consider tissue-specific normalization if analyzing heterogeneous samples
   - Always save normalized data separately

4. SPATIAL ANALYSIS
   - Use spatial autocorrelation to identify meaningful patterns
   - Perform clustering to identify tissue regions
   - Compare ROIs rather than individual pixels when possible

5. CONVERSION TO mass_dataset
   - Always save metadata when converting
   - Use preserve_pixel_metadata = TRUE to keep all information
   - Check conversion summary for any warnings

6. STATISTICAL ANALYSIS
   - Account for spatial correlation in statistical tests
   - Use appropriate multiple testing correction
   - Validate findings with spatial visualization

7. DATA EXPORT
   - Save intermediate results at each analysis step
   - Export high-quality figures for publications
   - Keep metadata files for reproducibility

8. PERFORMANCE
   - For large datasets, consider filtering before intensive operations
   - Use parallel processing when available
   - Save processed data to avoid re-computation

================================================================================
")


# ============================================================================
# PART 9: Troubleshooting Common Issues
# ============================================================================

cat("
================================================================================
TROUBLESHOOTING GUIDE
================================================================================

ISSUE 1: Conversion fails with dimension mismatch
SOLUTION: Check that expression_data columns match pixel_info rows
  nrow(pixel_info) == ncol(expression_data)
  
ISSUE 2: Missing values after conversion
SOLUTION: Check if variables were filtered in mass_dataset
  Compare n_variables before and after conversion

ISSUE 3: Spatial plots look strange
SOLUTION: Verify x, y coordinates are correct
  Check for duplicated coordinates
  Ensure coordinate system is consistent

ISSUE 4: Statistical tests not significant
SOLUTION: Check if spatial correlation is accounted for
  Verify group assignments are correct
  Consider using non-parametric tests

ISSUE 5: Out of memory errors
SOLUTION: Process data in chunks
  Filter variables/pixels before analysis
  Increase system memory allocation

ISSUE 6: Metadata lost after conversion
SOLUTION: Ensure save_metadata = TRUE
  Check .msi_conversion_cache exists
  Save metadata to file as backup

================================================================================
")


# ============================================================================
# End of Usage Guide
# ============================================================================

cat("
================================================================================
For more information:
- Package documentation: ?msidataset
- Bug reports: [GitHub issues page]
- Contact: XIAOYU021@e.ntu.edu.sg
================================================================================
")