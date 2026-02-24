test_that("clean_gwas calls harmonize_sumstats_headers then tidyGWAS", {
  skip_if_not_installed("mockery")

  fake_cleaned <- make_binary_sumstats_df(n = 10)

  mock_harmonize <- mockery::mock(
    data.frame(EFF_ALL_FREQ = rep(0.3, 10), stringsAsFactors = FALSE)
  )
  mock_tidygwas <- mockery::mock(fake_cleaned)

  withr::local_tempdir() |> (\(d) {
    tmp_sumstats <- file.path(d, "test.txt.gz")
    write.csv(data.frame(x = 1), tmp_sumstats)
    log_dir <- file.path(d, "logs")

    mockery::stub(clean_gwas, "harmonize_sumstats_headers", mock_harmonize)
    mockery::stub(clean_gwas, "tidyGWAS::tidyGWAS", mock_tidygwas)
    mockery::stub(clean_gwas, "fs::file_copy", function(...) invisible(NULL))

    result <- clean_gwas(tmp_sumstats, log_dir)
    expect_equal(result, fake_cleaned)
    mockery::expect_called(mock_harmonize, 1)
    mockery::expect_called(mock_tidygwas, 1)
  })()
})
