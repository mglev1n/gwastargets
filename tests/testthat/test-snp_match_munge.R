make_sumstats <- function() {
  data.frame(
    chr  = c(1L, 1L, 2L, 2L),
    rsid = c("rs1", "rs2", "rs3", "rs4"),
    a0   = c("A", "C", "G", "T"),
    a1   = c("G", "T", "A", "C"),
    beta = c(0.1, -0.2, 0.3, -0.1),
    N    = c(1000L, 1000L, 1000L, 1000L),
    EAF  = c(0.3, 0.4, 0.5, 0.2),
    stringsAsFactors = FALSE
  )
}

make_info_snp <- function() {
  data.frame(
    rsid = c("rs1", "rs2", "rs3", "rs4"),
    a0   = c("A", "C", "G", "T"),
    a1   = c("G", "T", "A", "C"),
    stringsAsFactors = FALSE
  )
}

test_that("matches variants on rsid+a0+a1", {
  result <- snp_match_munge(make_sumstats(), make_info_snp(),
                            strand_flip = FALSE, match.min.prop = 0)
  expect_true(nrow(result) >= 4)
  expect_true(all(c("rsid", "beta") %in% names(result)))
})

test_that("removes ambiguous SNPs when strand_flip = TRUE", {
  ss <- data.frame(
    chr  = c(1L, 1L),
    rsid = c("rs_amb", "rs_ok"),
    a0   = c("A", "C"),    # A/T is ambiguous
    a1   = c("T", "G"),
    beta = c(0.1, 0.2),
    N    = c(1000L, 1000L),
    EAF  = c(0.3, 0.4),
    stringsAsFactors = FALSE
  )
  info <- data.frame(
    rsid = c("rs_amb", "rs_ok"),
    a0   = c("A", "C"),
    a1   = c("T", "G"),
    stringsAsFactors = FALSE
  )
  result <- snp_match_munge(ss, info, strand_flip = TRUE, match.min.prop = 0)
  expect_false("rs_amb" %in% result$rsid)
  expect_true("rs_ok" %in% result$rsid)
})

test_that("reverses alleles and negates beta when needed", {
  ss <- data.frame(
    chr  = 1L,
    rsid = "rs1",
    a0   = "G",   # swapped vs info_snp
    a1   = "A",
    beta = 0.5,
    N    = 1000L,
    EAF  = 0.3,
    stringsAsFactors = FALSE
  )
  info <- data.frame(
    rsid = "rs1",
    a0   = "A",
    a1   = "G",
    stringsAsFactors = FALSE
  )
  result <- snp_match_munge(ss, info, strand_flip = FALSE,
                            match.min.prop = 0, return_flip_and_rev = TRUE)
  expect_equal(result$beta[result$`_REV_`], -0.5)
})

test_that("errors when required columns missing from sumstats", {
  bad_ss <- data.frame(chr = 1, rsid = "rs1", a0 = "A", a1 = "G")  # no beta
  expect_error(
    snp_match_munge(bad_ss, make_info_snp(), match.min.prop = 0),
    "proper names"
  )
})

test_that("errors when fewer than min_match variants matched", {
  expect_error(
    snp_match_munge(make_sumstats(), make_info_snp(),
                    strand_flip = FALSE, match.min.prop = 0.99),
    "Not enough variants"
  )
})
