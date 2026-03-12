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

# --- column_map tests ---

test_that("column_map = NULL behaves like no argument (backward compat)", {
  tmp <- withr::local_tempfile(fileext = ".txt")
  write.table(
    data.frame(BETA = 0.1, P = 0.05, SNP = "rs1"),
    tmp, row.names = FALSE, quote = FALSE, sep = "\t"
  )
  result <- harmonize_sumstats_headers(tmp, column_map = NULL)
  expect_true("B" %in% names(result))
})

test_that("custom column_map renames applied before hardcoded dictionary", {
  tmp <- withr::local_tempfile(fileext = ".txt")
  write.table(
    data.frame(MY_FREQ = 0.3, SNP = "rs1"),
    tmp, row.names = FALSE, quote = FALSE, sep = "\t"
  )
  result <- harmonize_sumstats_headers(tmp, column_map = c("MY_FREQ" = "EFF_ALL_FREQ"))
  expect_true("EFF_ALL_FREQ" %in% names(result))
  expect_false("MY_FREQ" %in% names(result))
})

test_that("custom map + hardcoded dictionary work together", {
  tmp <- withr::local_tempfile(fileext = ".txt")
  write.table(
    data.frame(MY_FREQ = 0.3, BETA = 0.1, SNP = "rs1"),
    tmp, row.names = FALSE, quote = FALSE, sep = "\t"
  )
  # Custom map renames MY_FREQ, hardcoded dict renames BETA
  result <- harmonize_sumstats_headers(tmp, column_map = c("MY_FREQ" = "EFF_ALL_FREQ"))
  expect_true("EFF_ALL_FREQ" %in% names(result))
  expect_true("B" %in% names(result))
})

test_that("custom map overrides hardcoded dictionary for same source column", {
  tmp <- withr::local_tempfile(fileext = ".txt")
  # FREQ would be mapped to EFF_ALL_FREQ by hardcoded dict,

  # but custom map renames it first to something else
  write.table(
    data.frame(FREQ = 0.3, SNP = "rs1"),
    tmp, row.names = FALSE, quote = FALSE, sep = "\t"
  )
  result <- harmonize_sumstats_headers(tmp, column_map = c("FREQ" = "MY_CUSTOM"))
  expect_true("MY_CUSTOM" %in% names(result))
  expect_false("EFF_ALL_FREQ" %in% names(result))
  expect_false("FREQ" %in% names(result))
})

test_that("column_map ignores source names not present in file", {
  tmp <- withr::local_tempfile(fileext = ".txt")
  write.table(
    data.frame(BETA = 0.1, SNP = "rs1"),
    tmp, row.names = FALSE, quote = FALSE, sep = "\t"
  )
  result <- harmonize_sumstats_headers(tmp, column_map = c("NOT_HERE" = "EFF_ALL_FREQ"))
  # BETA still gets mapped by hardcoded dict, NOT_HERE is silently ignored
  expect_true("B" %in% names(result))
  expect_false("EFF_ALL_FREQ" %in% names(result))
})
