# ============================================================================
# Variable ID Generation - Test Script
# ============================================================================

cat("\n========================================\n")
cat("Testing Variable ID Generation Functions\n")
cat("========================================\n\n")

# ============================================================================
# Test 1: Basic Functionality
# ============================================================================

cat("Test 1: Basic Functionality\n")
cat("----------------------------\n")

test_mz_1 <- c(50.00125, 50.00250, 50.00375, 51.00123, 51.00456, 100.0012)

result_1 <- generate_variable_id(test_mz_1)

cat("Input m/z values:\n")
print(test_mz_1)

cat("\nGenerated variable IDs:\n")
print(result_1)

cat("\nExpected: m50_1, m50_2, m50_3, m51_1, m51_2, m100_1\n")
cat("Result: ", paste(result_1, collapse = ", "), "\n")

# Verify
expected_1 <- c("m50_1", "m50_2", "m50_3", "m51_1", "m51_2", "m100_1")
if (all(result_1 == expected_1)) {
  cat("✓ Test 1 PASSED\n\n")
} else {
  cat("✗ Test 1 FAILED\n\n")
}

# ============================================================================
# Test 2: Zero Padding
# ============================================================================

cat("Test 2: Zero Padding\n")
cat("--------------------\n")

result_2 <- generate_variable_id_extended(test_mz_1, zero_pad = TRUE, pad_width = 3)

cat("Generated variable IDs with zero padding:\n")
print(result_2)

expected_2 <- c("m50_001", "m50_002", "m50_003", "m51_001", "m51_002", "m100_001")
if (all(result_2 == expected_2)) {
  cat("✓ Test 2 PASSED\n\n")
} else {
  cat("✗ Test 2 FAILED\n\n")
}

# ============================================================================
# Test 3: Custom Prefix and Separator
# ============================================================================

cat("Test 3: Custom Prefix and Separator\n")
cat("------------------------------------\n")

result_3 <- generate_variable_id_extended(
  test_mz_1, 
  prefix = "mass", 
  separator = "-"
)

cat("Generated variable IDs with custom format:\n")
print(result_3)

expected_3 <- c("mass50-1", "mass50-2", "mass50-3", "mass51-1", "mass51-2", "mass100-1")
if (all(result_3 == expected_3)) {
  cat("✓ Test 3 PASSED\n\n")
} else {
  cat("✗ Test 3 FAILED\n\n")
}

# ============================================================================
# Test 4: Large Number of Peaks in Same Mass
# ============================================================================

cat("Test 4: Multiple Peaks in Same Mass Range\n")
cat("------------------------------------------\n")

# Create 10 peaks in mass 150
test_mz_4 <- seq(150.001, 150.010, length.out = 10)

result_4 <- generate_variable_id(test_mz_4)

cat("Generated variable IDs for 10 peaks at m/z ~150:\n")
print(result_4)

# Should be m150_1 through m150_10
if (all(grepl("^m150_[0-9]+$", result_4))) {
  cat("✓ Test 4 PASSED\n\n")
} else {
  cat("✗ Test 4 FAILED\n\n")
}

# ============================================================================
# Test 5: Unsorted m/z Values (should still work)
# ============================================================================

cat("Test 5: Unsorted m/z Values\n")
cat("----------------------------\n")

# Unsorted m/z values
test_mz_5 <- c(51.002, 50.001, 100.005, 50.003, 51.001)

result_5 <- generate_variable_id(test_mz_5)

cat("Input (unsorted) m/z values:\n")
print(test_mz_5)

cat("\nGenerated variable IDs:\n")
print(result_5)

# IDs should correspond to original order, but numbering within mass should be sorted
# 51.002 -> should be m51_2 (second in 51 range after 51.001)
# 50.001 -> should be m50_1 (first in 50 range)
# 100.005 -> should be m100_1
# 50.003 -> should be m50_2 (second in 50 range after 50.001)
# 51.001 -> should be m51_1 (first in 51 range)
expected_5 <- c("m51_2", "m50_1", "m100_1", "m50_2", "m51_1")

if (all(result_5 == expected_5)) {
  cat("✓ Test 5 PASSED (correct handling of unsorted input)\n\n")
} else {
  cat("✗ Test 5 FAILED\n")
  cat("Expected:", paste(expected_5, collapse = ", "), "\n")
  cat("Got:     ", paste(result_5, collapse = ", "), "\n\n")
}

# ============================================================================
# Test 6: Uniqueness Validation
# ============================================================================

cat("Test 6: Uniqueness Validation\n")
cat("------------------------------\n")

# This should pass (all unique)
test_mz_6 <- c(50.001, 50.002, 51.001)
result_6 <- generate_variable_id(test_mz_6)

tryCatch({
  validate_variable_id(result_6, test_mz_6)
  cat("✓ Test 6 PASSED (unique IDs validated)\n\n")
}, error = function(e) {
  cat("✗ Test 6 FAILED:", e$message, "\n\n")
})

# ============================================================================
# Test 7: Edge Cases
# ============================================================================

cat("Test 7: Edge Cases\n")
cat("------------------\n")

# Edge case 1: Single m/z value
test_mz_7a <- c(75.0012)
result_7a <- generate_variable_id(test_mz_7a)
cat("Single m/z:", test_mz_7a, "->", result_7a, "\n")

if (result_7a == "m75_1") {
  cat("  ✓ Single value handled correctly\n")
} else {
  cat("  ✗ Single value test failed\n")
}

# Edge case 2: Very small m/z
test_mz_7b <- c(10.0001, 10.0002)
result_7b <- generate_variable_id(test_mz_7b)
cat("Small m/z:", paste(test_mz_7b, collapse = ", "), "->", 
    paste(result_7b, collapse = ", "), "\n")

# Edge case 3: Large m/z
test_mz_7c <- c(1500.001, 1500.002, 2000.001)
result_7c <- generate_variable_id(test_mz_7c)
cat("Large m/z:", paste(test_mz_7c, collapse = ", "), "->", 
    paste(result_7c, collapse = ", "), "\n")

if (all(c(result_7b[1] == "m10_1", result_7c[1] == "m1500_1"))) {
  cat("  ✓ Edge cases handled correctly\n\n")
} else {
  cat("  ✗ Edge case test failed\n\n")
}

# ============================================================================
# Test 8: Comparison with Sequential Mode
# ============================================================================

cat("Test 8: Comparison - mz_based vs sequential\n")
cat("--------------------------------------------\n")

test_mz_8 <- c(50.001, 50.002, 51.001, 100.001, 100.002)

# New mode
result_mz_based <- generate_variable_id(test_mz_8)

# Old mode (sequential)
result_sequential <- paste0("mz_", sprintf("%05d", 1:length(test_mz_8)))

comparison_df <- data.frame(
  mz = test_mz_8,
  mz_based = result_mz_based,
  sequential = result_sequential
)

cat("\nComparison table:\n")
print(comparison_df)

cat("\n✓ Comparison complete\n\n")

# ============================================================================
# Test Summary
# ============================================================================

cat("========================================\n")
cat("Test Summary\n")
cat("========================================\n\n")

cat("All tests completed!\n\n")

cat("Features tested:\n")
cat("  ✓ Basic variable ID generation\n")
cat("  ✓ Zero padding option\n")
cat("  ✓ Custom prefix and separator\n")
cat("  ✓ Multiple peaks in same mass range\n")
cat("  ✓ Unsorted m/z value handling\n")
cat("  ✓ Uniqueness validation\n")
cat("  ✓ Edge cases (single value, small/large m/z)\n")
cat("  ✓ Comparison with sequential mode\n\n")

# ============================================================================
# Create Sample Output Table
# ============================================================================

cat("========================================\n")
cat("Sample Output Table\n")
cat("========================================\n\n")

# Create a comprehensive example
sample_mz <- c(
  49.995, 50.001, 50.005, 50.010,  # Mass 50
  75.002,                           # Mass 75
  100.001, 100.005, 100.010,        # Mass 100
  150.002, 150.008,                 # Mass 150
  200.001                           # Mass 200
)

sample_ids <- generate_variable_id(sample_mz)

sample_table <- data.frame(
  Index = 1:length(sample_mz),
  variable_id = sample_ids,
  mz = round(sample_mz, 5),
  mass_int = floor(sample_mz)
)

cat("Sample variable_info table:\n")
print(sample_table)

cat("\n========================================\n")
cat("Testing Complete!\n")
cat("========================================\n")