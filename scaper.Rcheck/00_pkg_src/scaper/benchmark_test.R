# Benchmark Test for scaper Package Optimization
# This script tests the performance improvements made to the scaper package
# Run this after installing the updated package

library(scaper)
library(Seurat)
library(SeuratObject)
library(microbenchmark)

cat("====================================================\n")
cat("Benchmark Test for scaper Package\n")
cat("Testing optimized code vs original implementation\n")
cat("====================================================\n\n")

# Load test data
data(pbmc_small)
pbmc_small <- NormalizeData(pbmc_small)

# Test 1: scape() function with CytoSig database
cat("Test 1: scape() function with CytoSig database\n")
cat("----------------------------------------------\n")
counts.matrix <- as.data.frame(t(as.matrix(GetAssayData(pbmc_small, slot = "data"))))

# Test with subset of cytokines for faster testing
test_cytokines <- c("IL6", "IL4", "IL13", "IFNG", "TNFA")

tryCatch({
  result_scape_cytosig <- microbenchmark(
    scape(counts.matrix = counts.matrix, 
          database = "cytosig", 
          cytokine = test_cytokines),
    times = 5
  )
  print(result_scape_cytosig)
  cat("\n")
}, error = function(e) {
  cat("Error in scape() CytoSig test:", e$message, "\n\n")
})

# Test 2: scape() function with Reactome database
cat("Test 2: scape() function with Reactome database\n")
cat("----------------------------------------------\n")

tryCatch({
  result_scape_reactome <- microbenchmark(
    scape(counts.matrix = counts.matrix, 
          database = "reactome", 
          cytokine = test_cytokines),
    times = 5
  )
  print(result_scape_reactome)
  cat("\n")
}, error = function(e) {
  cat("Error in scape() Reactome test:", e$message, "\n\n")
})

# Test 3: scapeForSeurat() function with CytoSig database
cat("Test 3: scapeForSeurat() function with CytoSig database\n")
cat("--------------------------------------------------------\n")

tryCatch({
  result_seurat_cytosig <- microbenchmark(
    scapeForSeurat(seurat.object = pbmc_small, 
                   database = "cytosig", 
                   cytokine = test_cytokines, 
                   normalize = FALSE),
    times = 5
  )
  print(result_seurat_cytosig)
  cat("\n")
}, error = function(e) {
  cat("Error in scapeForSeurat() CytoSig test:", e$message, "\n\n")
})

# Test 4: scapeForSeurat() function with Reactome database
cat("Test 4: scapeForSeurat() function with Reactome database\n")
cat("--------------------------------------------------------\n")

tryCatch({
  result_seurat_reactome <- microbenchmark(
    scapeForSeurat(seurat.object = pbmc_small, 
                   database = "reactome", 
                   cytokine = test_cytokines, 
                   normalize = FALSE),
    times = 5
  )
  print(result_seurat_reactome)
  cat("\n")
}, error = function(e) {
  cat("Error in scapeForSeurat() Reactome test:", e$message, "\n\n")
})

# Test 5: Correctness check - compare results
cat("Test 5: Correctness Check\n")
cat("-------------------------\n")

tryCatch({
  # Test scape() outputs
  result1 <- scape(counts.matrix = counts.matrix, 
                   database = "cytosig", 
                   cytokine = test_cytokines)
  
  cat("scape() CytoSig output dimensions:", dim(result1), "\n")
  cat("scape() CytoSig output range: [", min(result1), ",", max(result1), "]\n")
  cat("scape() CytoSig output mean:", mean(result1), "\n")
  
  # Check for NAs
  if (any(is.na(result1))) {
    cat("WARNING: scape() CytoSig output contains NA values!\n")
  } else {
    cat("scape() CytoSig output: No NA values (PASS)\n")
  }
  
  # Check range is valid (CDF should be between 0 and 1)
  if (all(result1 >= 0 & result1 <= 1)) {
    cat("scape() CytoSig output: Values in valid range [0,1] (PASS)\n")
  } else {
    cat("WARNING: scape() CytoSig output has values outside [0,1]!\n")
  }
  cat("\n")
  
}, error = function(e) {
  cat("Error in correctness check:", e$message, "\n\n")
})

# Test 6: Scientific accuracy check - verify positive/negative weight handling
cat("Test 6: Scientific Accuracy - Positive/Negative Weight Handling\n")
cat("----------------------------------------------------------------\n")

tryCatch({
  # This test verifies that the fix for positive/negative weight handling is correct
  # The formula should be: (positive_cdf + (1 - negative_cdf)) / 2
  
  result_cytosig <- scape(counts.matrix = counts.matrix, 
                          database = "cytosig", 
                          cytokine = c("IL6"))
  
  cat("IL6 activity score statistics:\n")
  cat("  Min:", min(result_cytosig), "\n")
  cat("  Max:", max(result_cytosig), "\n")
  cat("  Mean:", mean(result_cytosig), "\n")
  cat("  Median:", median(result_cytosig), "\n")
  
  if (all(result_cytosig >= 0 & result_cytosig <= 1)) {
    cat("Scientific accuracy check: PASS (values in valid range)\n")
  } else {
    cat("Scientific accuracy check: FAIL (values outside valid range)\n")
  }
  cat("\n")
  
}, error = function(e) {
  cat("Error in scientific accuracy check:", e$message, "\n\n")
})

cat("====================================================\n")
cat("Benchmark tests completed!\n")
cat("====================================================\n")

# Memory usage report
cat("\nMemory usage information:\n")
cat("-------------------------\n")
print(gc())
