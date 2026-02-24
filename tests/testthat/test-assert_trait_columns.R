test_that("binary: errors when CaseN/ControlN/EffectiveN missing", {
  expect_error(
    assert_trait_columns(c("CHR", "POS_38", "B", "SE"), "binary"),
    "CaseN"
  )
})

test_that("binary: errors when only some binary columns missing", {
  expect_error(
    assert_trait_columns(c("CHR", "CaseN", "ControlN"), "binary"),
    "EffectiveN"
  )
})

test_that("binary: passes silently when all three binary columns present", {
  expect_silent(
    assert_trait_columns(c("CHR", "POS_38", "B", "SE", "CaseN", "ControlN", "EffectiveN"), "binary")
  )
})

test_that("quantitative: warns when binary columns detected", {
  expect_warning(
    assert_trait_columns(c("CHR", "B", "SE", "N", "CaseN", "ControlN"), "quantitative"),
    "binary columns"
  )
})

test_that("quantitative: silent when no binary columns present", {
  expect_silent(
    assert_trait_columns(c("CHR", "POS_38", "B", "SE", "N"), "quantitative")
  )
})
