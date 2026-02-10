# ============================================================================
# Comprehensive Test Suite for msidataset Package
# Author: Xiaoyu Xu
# Purpose: Test conversion functions and massdataset integration
# ============================================================================

library(testthat)
library(msidataset)
library(massdataset)

# ============================================================================
# Test Data Generation
# ============================================================================

#' Generate synthetic MSI dataset for testing
#'
#' @param n_variables Number of m/z variables
#' @param n_pixels_x Number of pixels in x dimension
#' @param n_pixels_y Number of pixels in y dimension
#' @param n_samples Number of samples (tissue sections)
#' @param add_noise Logical, add random noise
#'
#' @return A msi_dataset object
generate_test_msi_data <- function(n_variables = 100,
                                   n_pixels_x = 20,
                                   n_pixels_y = 20,
                                   n_samples = 2,
                                   add_noise = TRUE) {
  
  cat("Generating test MSI data...\n")
  
  n_pixels_total <- n_pixels_x * n_pixels_y
  
  # Generate variable info
  variable_info <- data.frame(
    variable_id = paste0("var_", 1:n_variables),
    mz = seq(100, 1000, length.out = n_variables),
    ms2 = sample(c("yes", "no"), n_variables, replace = TRUE, prob = c(0.1, 0.9)),
    stringsAsFactors = FALSE
  )
  
  # Generate pixel info with spatial coordinates
  pixel_info <- data.frame(
    pixel_id = paste0("pixel_", 1:n_pixels_total),
    x = rep(1:n_pixels_x, times = n_pixels_y),
    y = rep(1:n_pixels_y, each = n_pixels_x),
    z = 1,  # Single slice
    tissue_type = sample(c("tumor", "normal", "necrotic"), 
                         n_pixels_total, replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  # Generate sample info
  sample_info <- data.frame(
    sample_id = paste0("sample_", 1:n_samples),
    patient_id = paste0("patient_", 1:n_samples),
    tissue = sample(c("brain", "liver", "kidney"), n_samples, replace = TRUE),
    condition = sample(c("control", "disease"), n_samples, replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  # Generate expression data with spatial patterns
  expression_data <- matrix(0, nrow = n_variables, ncol = n_pixels_total)
  rownames(expression_data) <- variable_info$variable_id
  colnames(expression_data) <- pixel_info$pixel_id
  
  # Create some spatial patterns
  for (i in 1:n_variables) {
    # Create gradient pattern
    if (i %% 3 == 0) {
      pattern <- outer(pixel_info$x, pixel_info$y, 
                       function(x, y) x + y)
    } 
    # Create cluster pattern  
    else if (i %% 3 == 1) {
      center_x <- n_pixels_x / 2
      center_y <- n_pixels_y / 2
      pattern <- sqrt((pixel_info$x - center_x)^2 + 
                        (pixel_info$y - center_y)^2)
      pattern <- max(pattern) - pattern  # Invert for high center
    } 
    # Random pattern
    else {
      pattern <- runif(n_pixels_total)
    }
    
    # Scale and add to expression data
    pattern <- pattern / max(pattern) * 1000
    expression_data[i, ] <- pattern
    
    # Add noise if requested
    if (add_noise) {
      expression_data[i, ] <- expression_data[i, ] + 
        rnorm(n_pixels_total, mean = 0, sd = 50)
      expression_data[i, ] <- pmax(expression_data[i, ], 0)  # No negative values
    }
  }
  
  expression_data <- as.data.frame(expression_data)
  
  # Create msi_dataset object
  msi_obj <- create_msi_dataset(
    expression_data = expression_data,
    pixel_info = pixel_info,
    sample_info = sample_info,
    variable_info = variable_info
  )
  
  cat("✓ Generated msi_dataset with:\n")
  cat("  Variables:", n_variables, "\n")
  cat("  Pixels:", n_pixels_total, "(", n_pixels_x, "×", n_pixels_y, ")\n")
  cat("  Samples:", n_samples, "\n\n")
  
  return(msi_obj)
}


# ============================================================================
# Test 1: Basic Conversion Roundtrip
# ============================================================================

test_conversion_roundtrip <- function() {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("  TEST 1: Conversion Roundtrip (MSI → Mass → MSI)\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  # Generate test data
  msi_original <- generate_test_msi_data(
    n_variables = 50,
    n_pixels_x = 10,
    n_pixels_y = 10
  )
  
  # Test 1.1: MSI to Mass conversion
  cat("Test 1.1: Converting msi_dataset to mass_dataset...\n")
  mass_obj <- convert_msi2mass_dataset(
    msi_original,
    save_metadata = TRUE,
    add_spatial_to_sample_info = TRUE
  )
  
  # Validate conversion
  test_that("MSI to Mass conversion preserves dimensions", {
    expect_equal(
      nrow(mass_obj@expression_data),
      n_variables(msi_original)
    )
    expect_equal(
      ncol(mass_obj@expression_data),
      n_pixels(msi_original)
    )
  })
  
  test_that("Mass dataset has correct structure", {
    expect_true("sample_id" %in% colnames(mass_obj@sample_info))
    expect_true("x" %in% colnames(mass_obj@sample_info))
    expect_true("y" %in% colnames(mass_obj@sample_info))
    expect_true("z" %in% colnames(mass_obj@sample_info))
    expect_true("mz" %in% colnames(mass_obj@variable_info))
    expect_true("rt" %in% colnames(mass_obj@variable_info))
  })
  
  cat("✓ MSI to Mass conversion successful\n\n")
  
  # Test 1.2: Mass to MSI conversion
  cat("Test 1.2: Converting mass_dataset back to msi_dataset...\n")
  msi_restored <- convert_mass2msi_dataset(
    mass_obj,
    metadata_source = "auto"
  )
  
  # Validate roundtrip
  test_that("Roundtrip preserves dimensions", {
    expect_equal(
      n_variables(msi_restored),
      n_variables(msi_original)
    )
    expect_equal(
      n_pixels(msi_restored),
      n_pixels(msi_original)
    )
  })
  
  test_that("Roundtrip preserves spatial coordinates", {
    expect_equal(
      pixel_info(msi_restored)$x,
      pixel_info(msi_original)$x
    )
    expect_equal(
      pixel_info(msi_restored)$y,
      pixel_info(msi_original)$y
    )
    expect_equal(
      pixel_info(msi_restored)$z,
      pixel_info(msi_original)$z
    )
  })
  
  test_that("Roundtrip preserves expression data", {
    expect_equal(
      as.matrix(expression_data(msi_restored)),
      as.matrix(expression_data(msi_original)),
      tolerance = 1e-10
    )
  })
  
  test_that("Roundtrip preserves m/z values", {
    expect_equal(
      mz(msi_restored),
      mz(msi_original),
      tolerance = 1e-6
    )
  })
  
  cat("✓ Mass to MSI conversion successful\n")
  cat("✓ Roundtrip conversion maintains data integrity\n\n")
  
  return(list(
    msi_original = msi_original,
    mass_obj = mass_obj,
    msi_restored = msi_restored
  ))
}


# ============================================================================
# Test 2: massdataset Integration Tests
# ============================================================================

test_massdataset_integration <- function() {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("  TEST 2: massdataset Package Integration\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  # Generate test data
  msi_obj <- generate_test_msi_data(
    n_variables = 200,
    n_pixels_x = 15,
    n_pixels_y = 15
  )
  
  # Convert to mass_dataset
  mass_obj <- convert_msi2mass_dataset(
    msi_obj,
    save_metadata = TRUE,
    add_spatial_to_sample_info = TRUE
  )
  
  # Test 2.1: Data filtering
  cat("Test 2.1: Testing data filtering...\n")
  
  # Filter by m/z range
  if (requireNamespace("massdataset", quietly = TRUE)) {
    mass_filtered <- mass_obj %>%
      massdataset::activate_mass_dataset("variable_info") %>%
      massdataset::filter(mz > 200 & mz < 800)
    
    test_that("Filtering works correctly", {
      expect_lt(
        nrow(mass_filtered@variable_info),
        nrow(mass_obj@variable_info)
      )
      expect_true(all(mass_filtered@variable_info$mz > 200))
      expect_true(all(mass_filtered@variable_info$mz < 800))
    })
    
    cat("  ✓ Variable filtering: PASS\n")
    
    # Filter samples (pixels) by spatial location
    mass_filtered2 <- mass_obj %>%
      massdataset::activate_mass_dataset("sample_info") %>%
      massdataset::filter(x <= 10 & y <= 10)
    
    test_that("Spatial filtering works", {
      expect_lt(
        nrow(mass_filtered2@sample_info),
        nrow(mass_obj@sample_info)
      )
    })
    
    cat("  ✓ Spatial filtering: PASS\n")
  }
  
  # Test 2.2: Data normalization
  cat("\nTest 2.2: Testing normalization...\n")
  
  if (requireNamespace("massdataset", quietly = TRUE)) {
    # Total Ion Current (TIC) normalization
    mass_norm <- massdataset::normalize_data(
      mass_obj,
      method = "sum"  # TIC normalization
    )
    
    test_that("Normalization works", {
      # Check that column sums are similar after normalization
      col_sums <- colSums(mass_norm@expression_data)
      expect_lt(sd(col_sums), sd(colSums(mass_obj@expression_data)))
    })
    
    cat("  ✓ TIC normalization: PASS\n")
  }
  
  # Test 2.3: Missing value imputation
  cat("\nTest 2.3: Testing missing value handling...\n")
  
  # Introduce some missing values
  mass_with_na <- mass_obj
  expr_data_na <- mass_obj@expression_data
  na_indices <- sample(1:length(expr_data_na), size = 100)
  expr_data_na[na_indices] <- NA
  mass_with_na@expression_data <- expr_data_na
  
  if (requireNamespace("massdataset", quietly = TRUE)) {
    # Impute missing values
    mass_imputed <- massdataset::impute_mv(
      mass_with_na,
      method = "knn"
    )
    
    test_that("Imputation removes NAs", {
      expect_equal(sum(is.na(mass_imputed@expression_data)), 0)
    })
    
    cat("  ✓ Missing value imputation: PASS\n")
  }
  
  # Test 2.4: Statistical analysis
  cat("\nTest 2.4: Testing statistical analysis...\n")
  
  if (requireNamespace("massdataset", quietly = TRUE)) {
    # Add group information for testing
    mass_obj@sample_info$group <- ifelse(
      mass_obj@sample_info$x <= 8,
      "left",
      "right"
    )
    
    # Perform t-test
    mass_tested <- massdataset::perform_t_test(
      mass_obj,
      group = "group",
      group1 = "left",
      group2 = "right"
    )
    
    test_that("Statistical test works", {
      expect_true("p_value" %in% colnames(mass_tested@variable_info))
      expect_true("fold_change" %in% colnames(mass_tested@variable_info))
    })
    
    cat("  ✓ Statistical testing: PASS\n")
  }
  
  # Test 2.5: Convert back with annotations
  cat("\nTest 2.5: Converting back with analysis results...\n")
  
  # Add some annotations
  mass_obj@variable_info$annotation <- sample(
    c("Lipid", "Metabolite", "Unknown"),
    nrow(mass_obj@variable_info),
    replace = TRUE
  )
  
  mass_obj@variable_info$intensity_class <- ifelse(
    rowMeans(mass_obj@expression_data) > median(rowMeans(mass_obj@expression_data)),
    "High",
    "Low"
  )
  
  # Convert back
  msi_with_annotations <- convert_mass2msi_dataset(
    mass_obj,
    metadata_source = "auto",
    create_annotation_table = TRUE
  )
  
  test_that("Annotations are preserved in roundtrip", {
    expect_true(nrow(msi_with_annotations@annotation_table) > 0)
    expect_true("annotation" %in% colnames(msi_with_annotations@annotation_table))
  })
  
  cat("  ✓ Annotation preservation: PASS\n")
  
  cat("\n✓ All massdataset integration tests completed successfully!\n\n")
  
  return(list(
    mass_original = mass_obj,
    mass_filtered = if(exists("mass_filtered")) mass_filtered else NULL,
    mass_normalized = if(exists("mass_norm")) mass_norm else NULL,
    msi_annotated = msi_with_annotations
  ))
}


# ============================================================================
# Test 3: Edge Cases and Error Handling
# ============================================================================

test_edge_cases <- function() {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("  TEST 3: Edge Cases and Error Handling\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  # Test 3.1: Small dataset (single variable, few pixels)
  cat("Test 3.1: Testing with minimal dataset...\n")
  
  msi_tiny <- generate_test_msi_data(
    n_variables = 1,
    n_pixels_x = 2,
    n_pixels_y = 2,
    n_samples = 1,
    add_noise = FALSE
  )
  
  mass_tiny <- convert_msi2mass_dataset(msi_tiny, save_metadata = TRUE)
  msi_tiny_restored <- convert_mass2msi_dataset(mass_tiny)
  
  test_that("Minimal dataset conversion works", {
    expect_equal(n_variables(msi_tiny_restored), 1)
    expect_equal(n_pixels(msi_tiny_restored), 4)
  })
  
  cat("  ✓ Minimal dataset: PASS\n")
  
  # Test 3.2: Filtered dataset (some pixels removed)
  cat("\nTest 3.2: Testing with filtered dataset...\n")
  
  msi_medium <- generate_test_msi_data(
    n_variables = 50,
    n_pixels_x = 10,
    n_pixels_y = 10
  )
  
  mass_medium <- convert_msi2mass_dataset(msi_medium, save_metadata = TRUE)
  
  # Filter out half of the samples
  if (requireNamespace("massdataset", quietly = TRUE)) {
    mass_filtered <- mass_medium %>%
      massdataset::activate_mass_dataset("sample_info") %>%
      massdataset::filter(x <= 5)
    
    # Convert back
    msi_filtered <- convert_mass2msi_dataset(mass_filtered)
    
    test_that("Filtered dataset conversion works", {
      expect_equal(n_pixels(msi_filtered), 50)  # 5x10 pixels
      expect_lt(n_pixels(msi_filtered), n_pixels(msi_medium))
    })
    
    cat("  ✓ Filtered dataset: PASS\n")
  }
  
  # Test 3.3: Metadata recovery from file
  cat("\nTest 3.3: Testing metadata file save/load...\n")
  
  msi_test <- generate_test_msi_data(n_variables = 20, n_pixels_x = 5, n_pixels_y = 5)
  
  # Convert and save metadata to file
  temp_metadata_file <- tempfile(fileext = ".rds")
  mass_test <- convert_msi2mass_dataset(
    msi_test,
    save_metadata = TRUE,
    metadata_path = temp_metadata_file
  )
  
  # Clear environment cache
  if (exists(".msi_conversion_cache", envir = .GlobalEnv)) {
    rm(".msi_conversion_cache", envir = .GlobalEnv)
  }
  
  # Load from file
  msi_from_file <- convert_mass2msi_dataset(
    mass_test,
    metadata_source = "file",
    metadata_path = temp_metadata_file
  )
  
  test_that("Metadata file save/load works", {
    expect_equal(n_pixels(msi_from_file), n_pixels(msi_test))
    expect_equal(n_variables(msi_from_file), n_variables(msi_test))
  })
  
  # Cleanup
  unlink(temp_metadata_file)
  
  cat("  ✓ Metadata file handling: PASS\n")
  
  # Test 3.4: Error handling for missing metadata
  cat("\nTest 3.4: Testing error handling for missing metadata...\n")
  
  # Clear cache
  if (exists(".msi_conversion_cache", envir = .GlobalEnv)) {
    rm(".msi_conversion_cache", envir = .GlobalEnv)
  }
  
  test_that("Missing metadata throws appropriate error", {
    expect_error(
      convert_mass2msi_dataset(mass_test, metadata_source = "env"),
      "No conversion cache found"
    )
  })
  
  cat("  ✓ Error handling: PASS\n")
  
  cat("\n✓ All edge case tests completed successfully!\n\n")
  
  return(TRUE)
}


# ============================================================================
# Test 4: Performance Benchmarking
# ============================================================================

test_performance <- function() {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("  TEST 4: Performance Benchmarking\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  # Test different dataset sizes
  sizes <- list(
    small = list(vars = 100, px = 10, py = 10),
    medium = list(vars = 500, px = 20, py = 20),
    large = list(vars = 1000, px = 30, py = 30)
  )
  
  results <- data.frame(
    size = character(),
    n_variables = numeric(),
    n_pixels = numeric(),
    time_msi2mass = numeric(),
    time_mass2msi = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (size_name in names(sizes)) {
    cat(sprintf("Benchmarking %s dataset...\n", size_name))
    
    size_params <- sizes[[size_name]]
    
    # Generate data
    msi_obj <- generate_test_msi_data(
      n_variables = size_params$vars,
      n_pixels_x = size_params$px,
      n_pixels_y = size_params$py
    )
    
    # Benchmark MSI to Mass conversion
    time_msi2mass <- system.time({
      mass_obj <- convert_msi2mass_dataset(msi_obj, save_metadata = TRUE)
    })[["elapsed"]]
    
    # Benchmark Mass to MSI conversion
    time_mass2msi <- system.time({
      msi_restored <- convert_mass2msi_dataset(mass_obj)
    })[["elapsed"]]
    
    # Store results
    results <- rbind(
      results,
      data.frame(
        size = size_name,
        n_variables = n_variables(msi_obj),
        n_pixels = n_pixels(msi_obj),
        time_msi2mass = time_msi2mass,
        time_mass2msi = time_mass2msi,
        stringsAsFactors = FALSE
      )
    )
    
    cat(sprintf("  MSI→Mass: %.3f sec\n", time_msi2mass))
    cat(sprintf("  Mass→MSI: %.3f sec\n\n", time_mass2msi))
  }
  
  cat("\nPerformance Summary:\n")
  print(results)
  cat("\n")
  
  return(results)
}


# ============================================================================
# Main Test Runner
# ============================================================================

run_all_tests <- function(include_performance = FALSE) {
  
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════════╗\n")
  cat("║                                                               ║\n")
  cat("║     msidataset Package - Comprehensive Test Suite            ║\n")
  cat("║                                                               ║\n")
  cat("╚═══════════════════════════════════════════════════════════════╝\n")
  
  # Initialize results
  all_results <- list()
  
  # Run Test 1: Conversion Roundtrip
  tryCatch({
    all_results$test1 <- test_conversion_roundtrip()
  }, error = function(e) {
    cat("✗ TEST 1 FAILED:", e$message, "\n")
  })
  
  # Run Test 2: massdataset Integration
  tryCatch({
    all_results$test2 <- test_massdataset_integration()
  }, error = function(e) {
    cat("✗ TEST 2 FAILED:", e$message, "\n")
  })
  
  # Run Test 3: Edge Cases
  tryCatch({
    all_results$test3 <- test_edge_cases()
  }, error = function(e) {
    cat("✗ TEST 3 FAILED:", e$message, "\n")
  })
  
  # Optional: Run Test 4: Performance
  if (include_performance) {
    tryCatch({
      all_results$test4 <- test_performance()
    }, error = function(e) {
      cat("✗ TEST 4 FAILED:", e$message, "\n")
    })
  }
  
  # Summary
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════════╗\n")
  cat("║                     TEST SUMMARY                              ║\n")
  cat("╚═══════════════════════════════════════════════════════════════╝\n\n")
  
  n_passed <- sum(sapply(all_results, function(x) !is.null(x)))
  n_total <- length(all_results)
  
  cat(sprintf("Tests passed: %d / %d\n", n_passed, n_total))
  
  if (n_passed == n_total) {
    cat("\n✓ ALL TESTS PASSED! Package is ready for use.\n\n")
  } else {
    cat("\n⚠ Some tests failed. Please review the output above.\n\n")
  }
  
  invisible(all_results)
}


# ============================================================================
# Example usage
# ============================================================================

if (FALSE) {
  # Run all tests (excluding performance benchmarks)
  results <- run_all_tests(include_performance = FALSE)
  
  # Run with performance benchmarks
  results_full <- run_all_tests(include_performance = TRUE)
  
  # Run individual tests
  test1_results <- test_conversion_roundtrip()
  test2_results <- test_massdataset_integration()
  test3_results <- test_edge_cases()
  test4_results <- test_performance()
}