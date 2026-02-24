make_loci_df <- function(n_signif = 0) {
  set.seed(42)
  df <- data.frame(
    CHR     = rep(1L, 100),
    POS_38  = seq(1e6L, by = 1e4L, length.out = 100),
    RSID    = paste0("rs", seq_len(100)),
    EAF     = runif(100, 0.05, 0.95),
    B       = rnorm(100, 0, 0.05),
    SE      = runif(100, 0.02, 0.08),
    p_value = runif(100, 0.01, 0.99),
    stringsAsFactors = FALSE
  )
  if (n_signif > 0) {
    df$p_value[seq_len(n_signif)] <- runif(n_signif, 1e-10, 1e-8)
  }
  df
}

test_that("returns empty tibble when no variants pass p-threshold", {
  df <- make_loci_df(n_signif = 0)
  result <- extract_loci(df)
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
})

test_that("warns on missing chr/pos/p values", {
  df <- make_loci_df()
  df$CHR[1:5]     <- NA
  df$POS_38[6:10] <- NA
  # No genome-wide significant so it'll return empty, but warnings fire first
  expect_warning(
    extract_loci(df),
    "missing"
  )
})

test_that("errors when required columns missing", {
  df <- make_loci_df()
  df$B <- NULL
  expect_error(extract_loci(df), "Missing required columns")
})

test_that("calls gwasRtools functions and returns annotated tibble on success", {
  skip_if_not_installed("mockery")

  df <- make_loci_df(n_signif = 3)
  mock_get_loci <- mockery::mock(
    df[df$p_value < 5e-8, ] |> dplyr::mutate(lead = TRUE)
  )
  mock_get_nearest_gene <- mockery::mock(
    df[df$p_value < 5e-8, ] |>
      dplyr::mutate(lead = TRUE, nearest_gene = "GENE1")
  )

  mockery::stub(extract_loci, "gwasRtools::get_loci", mock_get_loci)
  mockery::stub(extract_loci, "gwasRtools::get_nearest_gene", mock_get_nearest_gene)

  result <- extract_loci(df)
  expect_s3_class(result, "tbl_df")
  expect_true("nearest_gene" %in% names(result))
})

