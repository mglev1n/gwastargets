make_common_df <- function(n = 40, chr = 1L, seed = 1L) {
  set.seed(seed)
  data.frame(
    CHR          = rep(chr, n),
    POS_38       = seq(1e6L, by = 2e4L, length.out = n),
    RSID         = paste0("rs", seq_len(n)),
    EffectAllele = rep("A", n),
    OtherAllele  = rep("G", n),
    EAF          = runif(n, 0.1, 0.9),
    N            = rep(10000L, n),
    EffectiveN   = rep(10000L, n),
    B            = rnorm(n, 0, 0.05),
    SE           = runif(n, 0.02, 0.08),
    CaseN        = rep(5000L, n),
    ControlN     = rep(5000L, n),
    stringsAsFactors = FALSE
  )
}

test_that("filters variants below maf_threshold", {
  df <- make_common_df(40)
  # set 5 variants to very low MAF
  df$EAF[1:5] <- 0.01

  tmp1 <- withr::local_tempfile(fileext = ".parquet")
  tmp2 <- withr::local_tempfile(fileext = ".parquet")
  arrow::write_parquet(df, tmp1)
  arrow::write_parquet(df, tmp2)

  result <- extract_common_variants(c(tmp1, tmp2), trait_type = "binary",
                                    maf_threshold = 0.05)
  # low-MAF variants should not appear
  expect_false(any(result$minor_af < 0.05))
})

test_that("returns only variants present in ALL files", {
  df1 <- make_common_df(30)
  df2 <- make_common_df(20)  # only first 20 of df1's variants

  tmp1 <- withr::local_tempfile(fileext = ".parquet")
  tmp2 <- withr::local_tempfile(fileext = ".parquet")
  arrow::write_parquet(df1, tmp1)
  arrow::write_parquet(df2, tmp2)

  result <- extract_common_variants(c(tmp1, tmp2), trait_type = "binary",
                                    maf_threshold = 0.0)
  # Every position in result must appear in both files
  pos_in_result <- unique(result$POS_38)
  pos_in_df1    <- df1$POS_38
  pos_in_df2    <- df2$POS_38
  expect_true(all(pos_in_result %in% pos_in_df1))
  expect_true(all(pos_in_result %in% pos_in_df2))
})

test_that("select_per_mb_window: at most 1 variant per window per chromosome", {
  df <- make_common_df(40)
  tmp1 <- withr::local_tempfile(fileext = ".parquet")
  tmp2 <- withr::local_tempfile(fileext = ".parquet")
  arrow::write_parquet(df, tmp1)
  arrow::write_parquet(df, tmp2)

  result <- extract_common_variants(
    c(tmp1, tmp2), trait_type = "binary",
    maf_threshold = 0.0, select_per_mb_window = TRUE, window_size_mb = 1,
    seed = 42L
  )
  # Per file, per chromosome: count distinct variants per 1Mb window
  result_one_file <- result[result$file == basename(tmp1), ]
  result_one_file$window <- floor(result_one_file$POS_38 / 1e6)
  counts <- table(result_one_file$CHR, result_one_file$window)
  expect_true(all(counts <= 1))
})

test_that("seed is respected for reproducibility", {
  df <- make_common_df(40)
  tmp1 <- withr::local_tempfile(fileext = ".parquet")
  tmp2 <- withr::local_tempfile(fileext = ".parquet")
  arrow::write_parquet(df, tmp1)
  arrow::write_parquet(df, tmp2)

  r1 <- extract_common_variants(c(tmp1, tmp2), trait_type = "binary",
                                 maf_threshold = 0.0,
                                 select_per_mb_window = TRUE, seed = 99L)
  r2 <- extract_common_variants(c(tmp1, tmp2), trait_type = "binary",
                                 maf_threshold = 0.0,
                                 select_per_mb_window = TRUE, seed = 99L)
  expect_identical(r1$POS_38, r2$POS_38)
})
