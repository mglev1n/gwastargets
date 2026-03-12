test_that("returns NULL when all args are NA", {
  expect_null(build_column_map())
})

test_that("returns named vector with correct source->target mapping", {
  result <- build_column_map(col_eaf = "MY_FREQ", col_beta = "BETA_VAL")
  expect_equal(result[["MY_FREQ"]], "EFF_ALL_FREQ")
  expect_equal(result[["BETA_VAL"]], "B")
  expect_length(result, 2)
})

test_that("handles mix of NA and non-NA args", {
  result <- build_column_map(col_chr = "CHROM", col_p = NA_character_, col_n = "SAMPLE_N")
  expect_equal(result, c("CHROM" = "CHR", "SAMPLE_N" = "N"))
  expect_false("P" %in% result)
})

test_that("single column override works", {
  result <- build_column_map(col_eaf = "AF")
  expect_equal(result, c("AF" = "EFF_ALL_FREQ"))
  expect_length(result, 1)
})

test_that("all 12 parameters produce correct targets", {
  result <- build_column_map(
    col_chr = "c1", col_pos = "c2", col_rsid = "c3",
    col_effect_allele = "c4", col_other_allele = "c5",
    col_beta = "c6", col_se = "c7", col_p = "c8",
    col_eaf = "c9", col_n = "c10",
    col_n_cases = "c11", col_n_controls = "c12"
  )
  expect_length(result, 12)
  expect_equal(unname(result),
    c("CHR", "POS", "RSID", "EffectAllele", "OtherAllele",
      "B", "SE", "P", "EFF_ALL_FREQ", "N", "N_CASES", "N_CONTROLS"))
  expect_equal(names(result), paste0("c", 1:12))
})
