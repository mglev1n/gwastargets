# Integration tests using inst/testdata/ parquet files
# These tests require the actual parquet files to be present.

testdata_dir <- system.file("testdata", package = "gwastargets")

skip_if_no_testdata <- function() {
  if (!file.exists(file.path(testdata_dir, "binary_study1.parquet"))) {
    skip("inst/testdata parquet files not found — run inst/testdata/create_testdata.R")
  }
}

test_that("integration: summarize_sumstats binary has correct schema", {
  skip_if_no_testdata()

  f1 <- file.path(testdata_dir, "binary_study1.parquet")
  f2 <- file.path(testdata_dir, "binary_study2.parquet")

  result <- summarize_sumstats(c(f1, f2), trait_type = "binary")
  expect_equal(nrow(result), 2)
  expect_true(all(c("n_cases", "n_controls", "n_effective", "n_snps") %in% names(result)))
  expect_true(all(result$n_snps > 0))
})

test_that("integration: meta_analyze_ivw binary has expected columns and valid p-values", {
  skip_if_no_testdata()

  f1 <- file.path(testdata_dir, "binary_study1.parquet")
  f2 <- file.path(testdata_dir, "binary_study2.parquet")

  result <- meta_analyze_ivw(c(f1, f2), trait_type = "binary",
                             chromosomes = 1:2, min_mac = 10)
  expect_true(nrow(result) > 0)
  expect_true(all(c("CHR", "POS_38", "RSID", "EffectAllele", "OtherAllele",
                    "B", "SE", "p_value", "CaseN", "ControlN", "EffectiveN",
                    "Q", "I2") %in% names(result)))
  expect_true(all(result$p_value >= 0 & result$p_value <= 1))
})

test_that("integration: extract_common_variants returns variants in all files", {
  skip_if_no_testdata()

  f1 <- file.path(testdata_dir, "binary_study1.parquet")
  f2 <- file.path(testdata_dir, "binary_study2.parquet")

  common <- extract_common_variants(c(f1, f2), trait_type = "binary",
                                    maf_threshold = 0.05)
  expect_true(nrow(common) > 0)

  # Every variant must appear in both files
  pos_by_file <- split(common$POS_38, common$file)
  if (length(pos_by_file) == 2) {
    expect_identical(sort(unique(pos_by_file[[1]])),
                     sort(unique(pos_by_file[[2]])))
  }
})

test_that("integration: calculate_mr_mega_mds returns valid MDS output", {
  skip_if_no_testdata()

  f1 <- file.path(testdata_dir, "binary_study1.parquet")
  f2 <- file.path(testdata_dir, "binary_study2.parquet")

  common <- extract_common_variants(c(f1, f2), trait_type = "binary",
                                    maf_threshold = 0.05)
  skip_if(nrow(common) < 5, "too few common variants for MDS")

  mds <- calculate_mr_mega_mds(common, pcCount = 2)
  expect_equal(nrow(mds$pc_coordinates), 2)
  expect_equal(ncol(mds$pc_coordinates), 2)
  expect_true(mds$summary$n_variants_used > 0)
})
