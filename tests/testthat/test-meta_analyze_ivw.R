make_meta_binary_df <- function(n = 60, chr = 1L) {
  set.seed(123)
  data.frame(
    CHR          = rep(chr, n),
    POS_38       = seq(1e6L, by = 1e4L, length.out = n),
    RSID         = paste0("rs", seq_len(n)),
    EffectAllele = rep("A", n),
    OtherAllele  = rep("G", n),
    B            = rnorm(n, 0, 0.05),
    SE           = runif(n, 0.02, 0.08),
    EAF          = runif(n, 0.1, 0.9),
    CaseN        = rep(5000L, n),
    ControlN     = rep(10000L, n),
    EffectiveN   = rep(6667L, n),
    N            = rep(15000L, n),
    stringsAsFactors = FALSE
  )
}

make_meta_quant_df <- function(n = 60, chr = 1L) {
  set.seed(123)
  data.frame(
    CHR          = rep(chr, n),
    POS_38       = seq(1e6L, by = 1e4L, length.out = n),
    RSID         = paste0("rs", seq_len(n)),
    EffectAllele = rep("A", n),
    OtherAllele  = rep("G", n),
    B            = rnorm(n, 0, 0.05),
    SE           = runif(n, 0.02, 0.08),
    EAF          = runif(n, 0.1, 0.9),
    N            = rep(20000L, n),
    EffectiveN   = rep(20000L, n),
    stringsAsFactors = FALSE
  )
}

test_that("binary: output has CaseN, ControlN, EffectiveN columns", {
  tmp1 <- withr::local_tempfile(fileext = ".parquet")
  tmp2 <- withr::local_tempfile(fileext = ".parquet")
  arrow::write_parquet(make_meta_binary_df(60, 1L), tmp1)
  arrow::write_parquet(make_meta_binary_df(60, 1L), tmp2)

  result <- meta_analyze_ivw(c(tmp1, tmp2), trait_type = "binary",
                             chromosomes = 1L, min_mac = 10)
  expect_true("CaseN" %in% names(result))
  expect_true("ControlN" %in% names(result))
  expect_true("EffectiveN" %in% names(result))
})

test_that("quantitative: output does NOT have CaseN, ControlN but has EffectiveN", {
  tmp1 <- withr::local_tempfile(fileext = ".parquet")
  tmp2 <- withr::local_tempfile(fileext = ".parquet")
  arrow::write_parquet(make_meta_quant_df(60, 1L), tmp1)
  arrow::write_parquet(make_meta_quant_df(60, 1L), tmp2)

  result <- meta_analyze_ivw(c(tmp1, tmp2), trait_type = "quantitative",
                             chromosomes = 1L, min_mac = 10)
  expect_false("CaseN" %in% names(result))
  expect_false("ControlN" %in% names(result))
  expect_true("EffectiveN" %in% names(result))
})

test_that("all p-values in [0, 1]", {
  tmp1 <- withr::local_tempfile(fileext = ".parquet")
  tmp2 <- withr::local_tempfile(fileext = ".parquet")
  arrow::write_parquet(make_meta_binary_df(60, 1L), tmp1)
  arrow::write_parquet(make_meta_binary_df(60, 1L), tmp2)

  result <- meta_analyze_ivw(c(tmp1, tmp2), trait_type = "binary",
                             chromosomes = 1L, min_mac = 10)
  expect_true(all(result$p_value >= 0 & result$p_value <= 1))
})

test_that("weighted beta/SE mathematically correct for 2-study case", {
  # Construct two studies with known B and SE for a single variant
  # IVW: B_meta = (B1/SE1^2 + B2/SE2^2) / (1/SE1^2 + 1/SE2^2)
  b1 <- 0.2; se1 <- 0.1
  b2 <- 0.4; se2 <- 0.2
  w1 <- 1 / se1^2; w2 <- 1 / se2^2
  expected_b  <- (b1 * w1 + b2 * w2) / (w1 + w2)
  expected_se <- 1 / sqrt(w1 + w2)

  make_single_variant <- function(b, se, n = 20000L, chr = 1L) {
    data.frame(
      CHR = chr, POS_38 = 1000000L, RSID = "rs1",
      EffectAllele = "A", OtherAllele = "G",
      B = b, SE = se, EAF = 0.3,
      N = n, EffectiveN = n,
      stringsAsFactors = FALSE
    )
  }

  tmp1 <- withr::local_tempfile(fileext = ".parquet")
  tmp2 <- withr::local_tempfile(fileext = ".parquet")
  arrow::write_parquet(make_single_variant(b1, se1), tmp1)
  arrow::write_parquet(make_single_variant(b2, se2), tmp2)

  result <- meta_analyze_ivw(c(tmp1, tmp2), trait_type = "quantitative",
                             chromosomes = 1L, min_mac = 0)

  expect_equal(result$B[1], expected_b, tolerance = 1e-6)
  expect_equal(result$SE[1], expected_se, tolerance = 1e-6)
})

test_that("I2 values in [0, 1]", {
  tmp1 <- withr::local_tempfile(fileext = ".parquet")
  tmp2 <- withr::local_tempfile(fileext = ".parquet")
  arrow::write_parquet(make_meta_binary_df(60, 1L), tmp1)
  arrow::write_parquet(make_meta_binary_df(60, 1L), tmp2)

  result <- meta_analyze_ivw(c(tmp1, tmp2), trait_type = "binary",
                             chromosomes = 1L, min_mac = 10)
  expect_true(all(result$I2 >= 0 & result$I2 <= 1))
})

test_that("MAC filter removes low-count variants", {
  tmp1 <- withr::local_tempfile(fileext = ".parquet")
  df <- make_meta_binary_df(10, 1L)
  df$N <- 10L  # very small N → low MAC
  arrow::write_parquet(df, tmp1)

  result <- meta_analyze_ivw(tmp1, trait_type = "binary",
                             chromosomes = 1L, min_mac = 10000)
  expect_equal(nrow(result), 0)
})

test_that("warning when no variants remain for a chromosome", {
  tmp1 <- withr::local_tempfile(fileext = ".parquet")
  df <- make_meta_binary_df(10, 1L)
  df$N <- 10L
  arrow::write_parquet(df, tmp1)

  expect_warning(
    meta_analyze_ivw(tmp1, trait_type = "binary",
                     chromosomes = 1L, min_mac = 10000),
    "No variants"
  )
})
