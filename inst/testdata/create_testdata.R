# Script to generate inst/testdata parquet files
# Run once from the package root: source("inst/testdata/create_testdata.R")

library(arrow)

set.seed(2024)

make_binary_parquet <- function(path, n = 100, chrs = c(1L, 2L), seed_offset = 0) {
  set.seed(2024 + seed_offset)
  n_per_chr <- n %/% length(chrs)
  df <- do.call(rbind, lapply(chrs, function(chr) {
    data.frame(
      CHR          = rep(chr, n_per_chr),
      POS_38       = seq(1000000L, by = 50000L, length.out = n_per_chr),
      RSID         = paste0("rs", chr, "_", seq_len(n_per_chr)),
      EffectAllele = sample(c("A", "C"), n_per_chr, replace = TRUE),
      OtherAllele  = sample(c("G", "T"), n_per_chr, replace = TRUE),
      B            = rnorm(n_per_chr, 0, 0.1),
      SE           = runif(n_per_chr, 0.02, 0.06),
      EAF          = runif(n_per_chr, 0.05, 0.95),
      Z            = rnorm(n_per_chr),
      CaseN        = rep(5000L, n_per_chr),
      ControlN     = rep(10000L, n_per_chr),
      EffectiveN   = rep(6667L, n_per_chr),
      N            = rep(15000L, n_per_chr),
      P            = runif(n_per_chr, 0, 1),
      SE_ldsc_adjusted = runif(n_per_chr, 0.02, 0.06),
      Z_ldsc_adjusted  = rnorm(n_per_chr),
      P_ldsc_adjusted  = runif(n_per_chr, 0, 1),
      ldsc_adjustment  = FALSE,
      intercept        = 0.98,
      stringsAsFactors = FALSE
    )
  }))
  arrow::write_parquet(df, path)
  cat("Written:", path, "\n")
  invisible(path)
}

make_quant_parquet <- function(path, n = 100, chrs = c(1L, 2L), seed_offset = 0) {
  set.seed(2024 + seed_offset)
  n_per_chr <- n %/% length(chrs)
  df <- do.call(rbind, lapply(chrs, function(chr) {
    data.frame(
      CHR          = rep(chr, n_per_chr),
      POS_38       = seq(1000000L, by = 50000L, length.out = n_per_chr),
      RSID         = paste0("rs", chr, "_", seq_len(n_per_chr)),
      EffectAllele = sample(c("A", "C"), n_per_chr, replace = TRUE),
      OtherAllele  = sample(c("G", "T"), n_per_chr, replace = TRUE),
      B            = rnorm(n_per_chr, 0, 0.1),
      SE           = runif(n_per_chr, 0.02, 0.06),
      EAF          = runif(n_per_chr, 0.05, 0.95),
      Z            = rnorm(n_per_chr),
      N            = rep(20000L, n_per_chr),
      EffectiveN   = rep(20000L, n_per_chr),
      P            = runif(n_per_chr, 0, 1),
      SE_ldsc_adjusted = runif(n_per_chr, 0.02, 0.06),
      Z_ldsc_adjusted  = rnorm(n_per_chr),
      P_ldsc_adjusted  = runif(n_per_chr, 0, 1),
      ldsc_adjustment  = FALSE,
      intercept        = 0.98,
      stringsAsFactors = FALSE
    )
  }))
  arrow::write_parquet(df, path)
  cat("Written:", path, "\n")
  invisible(path)
}

pkg_root <- here::here()
testdata_dir <- file.path(pkg_root, "inst", "testdata")
dir.create(testdata_dir, showWarnings = FALSE, recursive = TRUE)

make_binary_parquet(file.path(testdata_dir, "binary_study1.parquet"), seed_offset = 0)
make_binary_parquet(file.path(testdata_dir, "binary_study2.parquet"), seed_offset = 10)
make_quant_parquet(file.path(testdata_dir, "quant_study1.parquet"),   seed_offset = 20)
make_quant_parquet(file.path(testdata_dir, "quant_study2.parquet"),   seed_offset = 30)

cat("All testdata files created.\n")
