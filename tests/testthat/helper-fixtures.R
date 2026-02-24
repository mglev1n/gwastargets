# Helper fixtures for gwastargets tests
# Provides synthetic data factories and temp-file helpers

#' Make a synthetic binary-trait summary statistics data frame
#' Columns: CHR, POS_38, RSID, EffectAllele, OtherAllele, B, SE, EAF, Z,
#'          CaseN, ControlN, EffectiveN, N, P,
#'          SE_ldsc_adjusted, Z_ldsc_adjusted, P_ldsc_adjusted,
#'          ldsc_adjustment, intercept
make_binary_sumstats_df <- function(n = 50, chr = 1L) {
  set.seed(42)
  data.frame(
    CHR          = rep(chr, n),
    POS_38       = seq(1000000L, by = 10000L, length.out = n),
    RSID         = paste0("rs", seq_len(n)),
    EffectAllele = sample(c("A", "C", "G", "T"), n, replace = TRUE),
    OtherAllele  = sample(c("A", "C", "G", "T"), n, replace = TRUE),
    B            = rnorm(n, 0, 0.1),
    SE           = runif(n, 0.01, 0.05),
    EAF          = runif(n, 0.05, 0.95),
    Z            = rnorm(n),
    CaseN        = rep(5000L, n),
    ControlN     = rep(10000L, n),
    EffectiveN   = rep(6667L, n),
    N            = rep(15000L, n),
    P            = runif(n, 0, 1),
    SE_ldsc_adjusted = runif(n, 0.01, 0.05),
    Z_ldsc_adjusted  = rnorm(n),
    P_ldsc_adjusted  = runif(n, 0, 1),
    ldsc_adjustment  = FALSE,
    intercept        = 0.98,
    stringsAsFactors = FALSE
  )
}

#' Make a synthetic quantitative-trait summary statistics data frame
#' Same as binary minus CaseN, ControlN, EffectiveN
make_quant_sumstats_df <- function(n = 50, chr = 1L) {
  set.seed(42)
  data.frame(
    CHR          = rep(chr, n),
    POS_38       = seq(1000000L, by = 10000L, length.out = n),
    RSID         = paste0("rs", seq_len(n)),
    EffectAllele = sample(c("A", "C", "G", "T"), n, replace = TRUE),
    OtherAllele  = sample(c("A", "C", "G", "T"), n, replace = TRUE),
    B            = rnorm(n, 0, 0.1),
    SE           = runif(n, 0.01, 0.05),
    EAF          = runif(n, 0.05, 0.95),
    Z            = rnorm(n),
    N            = rep(20000L, n),
    P            = runif(n, 0, 1),
    SE_ldsc_adjusted = runif(n, 0.01, 0.05),
    Z_ldsc_adjusted  = rnorm(n),
    P_ldsc_adjusted  = runif(n, 0, 1),
    ldsc_adjustment  = FALSE,
    intercept        = 0.98,
    stringsAsFactors = FALSE
  )
}

#' Write a binary parquet file to a path
write_binary_parquet <- function(path, n = 50, chr = 1L) {
  df <- make_binary_sumstats_df(n = n, chr = chr)
  arrow::write_parquet(df, path)
  invisible(path)
}

#' Write a quantitative parquet file to a path
write_quant_parquet <- function(path, n = 50, chr = 1L) {
  df <- make_quant_sumstats_df(n = n, chr = chr)
  arrow::write_parquet(df, path)
  invisible(path)
}
