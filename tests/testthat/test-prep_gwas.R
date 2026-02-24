test_that("errors when n_col missing", {
  expect_error(
    prep_gwas("file.txt.gz", hm3 = NULL, ancestry = "EUR",
              output_path = tempdir(), trait_type = "binary"),
    "n_col"
  )
})

test_that("errors when trait_type='binary' and n_col=N", {
  expect_error(
    prep_gwas("file.txt.gz", hm3 = NULL, ancestry = "EUR",
              output_path = tempdir(), trait_type = "binary", n_col = N),
    "EffectiveN"
  )
})

test_that("warns when trait_type='quantitative' and n_col=EffectiveN", {
  skip_if_not_installed("mockery")

  fake_df <- make_quant_sumstats_df(n = 10)
  mock_clean  <- mockery::mock(fake_df)
  mock_munge  <- mockery::mock(
    data.frame(SNP = paste0("rs", 1:5), N = 1000, Z = rnorm(5),
               A1 = "A", A2 = "G", EAF = 0.3, stringsAsFactors = FALSE)
  )
  mock_ldsc <- mockery::mock(list(intercept = 0.98))

  mockery::stub(prep_gwas, "clean_gwas", mock_clean)
  mockery::stub(prep_gwas, "snp_match_munge", mock_munge)
  mockery::stub(prep_gwas, "ldscr::ldsc_h2", mock_ldsc)
  mockery::stub(prep_gwas, "arrow::write_parquet", function(...) invisible(NULL))

  expect_warning(
    prep_gwas("file.txt.gz", hm3 = data.frame(SNP = "rs1", A1 = "A", A2 = "G"),
              ancestry = "EUR", output_path = tempdir(),
              trait_type = "quantitative", n_col = EffectiveN),
    "unusual"
  )
})

test_that("writes parquet with ldsc_adjustment=TRUE when intercept>1", {
  skip_if_not_installed("mockery")

  fake_df <- make_binary_sumstats_df(n = 10)
  mock_clean <- mockery::mock(fake_df)
  mock_munge <- mockery::mock(
    data.frame(SNP = paste0("rs", 1:5), N = 1000, Z = rnorm(5),
               A1 = "A", A2 = "G", EAF = 0.3, stringsAsFactors = FALSE)
  )
  mock_ldsc <- mockery::mock(list(intercept = 1.05))

  written_data <- NULL
  mock_write <- function(x, ...) { written_data <<- x; invisible(NULL) }

  mockery::stub(prep_gwas, "clean_gwas", mock_clean)
  mockery::stub(prep_gwas, "snp_match_munge", mock_munge)
  mockery::stub(prep_gwas, "ldscr::ldsc_h2", mock_ldsc)
  mockery::stub(prep_gwas, "arrow::write_parquet", mock_write)

  withr::with_tempdir({
    prep_gwas("cohort.txt.gz",
              hm3 = data.frame(SNP = "rs1", A1 = "A", A2 = "G"),
              ancestry = "EUR", output_path = ".",
              trait_type = "binary", n_col = EffectiveN)
    expect_true(all(written_data$ldsc_adjustment))
  })
})

test_that("writes parquet with ldsc_adjustment=FALSE when intercept<=1", {
  skip_if_not_installed("mockery")

  fake_df <- make_binary_sumstats_df(n = 10)
  mock_clean <- mockery::mock(fake_df)
  mock_munge <- mockery::mock(
    data.frame(SNP = paste0("rs", 1:5), N = 1000, Z = rnorm(5),
               A1 = "A", A2 = "G", EAF = 0.3, stringsAsFactors = FALSE)
  )
  mock_ldsc <- mockery::mock(list(intercept = 0.95))

  written_data <- NULL
  mock_write <- function(x, ...) { written_data <<- x; invisible(NULL) }

  mockery::stub(prep_gwas, "clean_gwas", mock_clean)
  mockery::stub(prep_gwas, "snp_match_munge", mock_munge)
  mockery::stub(prep_gwas, "ldscr::ldsc_h2", mock_ldsc)
  mockery::stub(prep_gwas, "arrow::write_parquet", mock_write)

  withr::with_tempdir({
    prep_gwas("cohort.txt.gz",
              hm3 = data.frame(SNP = "rs1", A1 = "A", A2 = "G"),
              ancestry = "EUR", output_path = ".",
              trait_type = "binary", n_col = EffectiveN)
    expect_false(any(written_data$ldsc_adjustment))
  })
})
