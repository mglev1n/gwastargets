test_that("maps BETA to B", {
  tmp <- withr::local_tempfile(fileext = ".txt")
  write.table(
    data.frame(BETA = 0.1, P = 0.05, SNP = "rs1"),
    tmp, row.names = FALSE, quote = FALSE, sep = "\t"
  )
  result <- harmonize_sumstats_headers(tmp)
  expect_true("B" %in% names(result))
  expect_false("BETA" %in% names(result))
})

test_that("maps PVAL to P", {
  tmp <- withr::local_tempfile(fileext = ".txt")
  write.table(
    data.frame(PVAL = 0.05, B = 0.1, SNP = "rs1"),
    tmp, row.names = FALSE, quote = FALSE, sep = "\t"
  )
  result <- harmonize_sumstats_headers(tmp)
  expect_true("P" %in% names(result))
  expect_false("PVAL" %in% names(result))
})

test_that("maps N_CASE to N_CASES and N_CONTROL to N_CONTROLS", {
  tmp <- withr::local_tempfile(fileext = ".txt")
  write.table(
    data.frame(N_CASE = 1000, N_CONTROL = 2000, SNP = "rs1"),
    tmp, row.names = FALSE, quote = FALSE, sep = "\t"
  )
  result <- harmonize_sumstats_headers(tmp)
  expect_true("N_CASES" %in% names(result))
  expect_true("N_CONTROLS" %in% names(result))
})

test_that("maps FREQ to EFF_ALL_FREQ", {
  tmp <- withr::local_tempfile(fileext = ".txt")
  write.table(
    data.frame(FREQ = 0.3, SNP = "rs1"),
    tmp, row.names = FALSE, quote = FALSE, sep = "\t"
  )
  result <- harmonize_sumstats_headers(tmp)
  expect_true("EFF_ALL_FREQ" %in% names(result))
})

test_that("leaves unmapped columns unchanged", {
  tmp <- withr::local_tempfile(fileext = ".txt")
  write.table(
    data.frame(CHR = 1, CUSTOM_COL = "foo", SNP = "rs1"),
    tmp, row.names = FALSE, quote = FALSE, sep = "\t"
  )
  result <- harmonize_sumstats_headers(tmp)
  expect_true("CUSTOM_COL" %in% names(result))
  expect_true("CHR" %in% names(result))
})

test_that("handles file with no mappable columns", {
  tmp <- withr::local_tempfile(fileext = ".txt")
  write.table(
    data.frame(FOO = 1, BAR = 2),
    tmp, row.names = FALSE, quote = FALSE, sep = "\t"
  )
  result <- harmonize_sumstats_headers(tmp)
  expect_equal(sort(names(result)), c("BAR", "FOO"))
})
