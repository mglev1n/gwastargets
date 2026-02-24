test_that("binary: output has n_cases, n_controls, n_effective columns", {
  tmp1 <- withr::local_tempfile(fileext = ".parquet")
  tmp2 <- withr::local_tempfile(fileext = ".parquet")
  write_binary_parquet(tmp1, n = 20)
  write_binary_parquet(tmp2, n = 20)

  result <- summarize_sumstats(c(tmp1, tmp2), trait_type = "binary")
  expect_true("n_cases" %in% names(result))
  expect_true("n_controls" %in% names(result))
  expect_true("n_effective" %in% names(result))
})

test_that("quantitative: output does NOT have n_cases, n_controls", {
  tmp1 <- withr::local_tempfile(fileext = ".parquet")
  write_quant_parquet(tmp1, n = 20)

  result <- summarize_sumstats(tmp1, trait_type = "quantitative")
  expect_false("n_cases" %in% names(result))
  expect_false("n_controls" %in% names(result))
})

test_that("errors when binary columns missing for trait_type='binary'", {
  tmp <- withr::local_tempfile(fileext = ".parquet")
  write_quant_parquet(tmp, n = 20)  # no CaseN/ControlN/EffectiveN
  expect_error(
    summarize_sumstats(tmp, trait_type = "binary"),
    "CaseN|required"
  )
})

test_that("correct row count: one row per file", {
  tmp1 <- withr::local_tempfile(fileext = ".parquet")
  tmp2 <- withr::local_tempfile(fileext = ".parquet")
  write_quant_parquet(tmp1, n = 20)
  write_quant_parquet(tmp2, n = 20)

  result <- summarize_sumstats(c(tmp1, tmp2), trait_type = "quantitative")
  expect_equal(nrow(result), 2)
})

test_that("min_eaf <= max_eaf for each file", {
  tmp <- withr::local_tempfile(fileext = ".parquet")
  write_binary_parquet(tmp, n = 50)

  result <- summarize_sumstats(tmp, trait_type = "binary")
  expect_true(all(result$min_eaf <= result$max_eaf))
})
