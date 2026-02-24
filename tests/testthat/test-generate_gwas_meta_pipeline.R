test_that("returns a character string", {
  result <- generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                        n_col = "EffectiveN")
  expect_type(result, "character")
  expect_length(result, 1)
})

test_that("contains trait name in uppercase and lowercase", {
  result <- generate_gwas_meta_pipeline("BMI", trait_type = "quantitative",
                                        n_col = "N")
  expect_true(grepl("BMI", result))
  expect_true(grepl("bmi", result))
})

test_that("binary: output contains EffectiveN, n_cases, n_controls", {
  result <- generate_gwas_meta_pipeline("T2D", trait_type = "binary",
                                        n_col = "EffectiveN")
  expect_true(grepl("EffectiveN", result))
  expect_true(grepl("n_cases", result))
  expect_true(grepl("n_controls", result))
})

test_that("quantitative: output contains n_total but NOT n_cases", {
  result <- generate_gwas_meta_pipeline("HDL", trait_type = "quantitative",
                                        n_col = "N")
  expect_true(grepl("n_total", result))
  expect_false(grepl("n_cases", result))
})

test_that("errors when trait_type missing", {
  expect_error(
    generate_gwas_meta_pipeline("CAD", n_col = "EffectiveN"),
    "trait_type"
  )
})

test_that("errors when trait_type invalid", {
  expect_error(
    generate_gwas_meta_pipeline("CAD", trait_type = "continuous", n_col = "N"),
    "trait_type"
  )
})

test_that("errors when n_col missing", {
  expect_error(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary"),
    "n_col"
  )
})

test_that("errors when n_col invalid", {
  expect_error(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary", n_col = "BadCol"),
    "n_col"
  )
})

test_that("errors when binary + N (should use EffectiveN)", {
  expect_error(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary", n_col = "N"),
    "EffectiveN"
  )
})

test_that("warns when quantitative + EffectiveN", {
  expect_warning(
    generate_gwas_meta_pipeline("BMI", trait_type = "quantitative",
                                n_col = "EffectiveN"),
    "unusual"
  )
})
