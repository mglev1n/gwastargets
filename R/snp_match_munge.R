#' Match and munge SNPs against a reference panel
#'
#' @title Match GWAS SNPs against a reference SNP list
#'
#' @description
#' Matches variants in `sumstats` against `info_snp` by rsid, effect allele
#' (`a1`), and other allele (`a0`). Optionally removes ambiguous SNPs
#' (A/T and C/G), performs strand flipping, and reverses allele coding.
#'
#' @param sumstats A data frame of GWAS summary statistics. Must contain
#'   columns `rsid`, `a0`, `a1`, and `beta`.
#' @param info_snp A data frame of reference SNPs. Must contain columns
#'   `rsid`, `a0`, and `a1`.
#' @param strand_flip Logical; if `TRUE` (default), remove ambiguous SNPs
#'   (A/T, C/G) and attempt strand flipping.
#' @param join_by_pos Logical; if `TRUE`, also join by position. Default
#'   `TRUE`.
#' @param remove_dups Logical; if `TRUE` (default), remove duplicate
#'   rsid entries after matching.
#' @param match.min.prop Minimum proportion of variants that must match.
#'   Default `0.2`. An error is raised if fewer variants match.
#' @param return_flip_and_rev Logical; if `TRUE`, retain the `_FLIP_` and
#'   `_REV_` indicator columns in the output. Default `FALSE`.
#'
#' @return A data frame of matched variants ordered by chromosome.
#'
#' @examples
#' \dontrun{
#' matched <- snp_match_munge(
#'   sumstats = my_sumstats,
#'   info_snp = hm3_snps
#' )
#' }
#'
#' @importFrom bigassertr stop2 message2
#' @importFrom bigparallelr rows_along
#' @importFrom data.table as.data.table
#' @importFrom vctrs vec_in vec_duplicate_detect
#' @export
snp_match_munge <- function(sumstats, info_snp,
                            strand_flip         = TRUE,
                            join_by_pos         = TRUE,
                            remove_dups         = TRUE,
                            match.min.prop      = 0.2,
                            return_flip_and_rev = FALSE) {
  sumstats <- as.data.frame(sumstats)
  info_snp <- as.data.frame(info_snp)

  sumstats$`_NUM_ID_` <- bigparallelr::rows_along(sumstats)
  info_snp$`_NUM_ID_` <- bigparallelr::rows_along(info_snp)

  min_match <- match.min.prop * min(nrow(sumstats), nrow(info_snp))

  join_by <- c("rsid", "a0", "a1")

  if (!all(c(join_by, "beta") %in% names(sumstats))) {
    bigassertr::stop2(
      "Please use proper names for variables in 'sumstats'. Expected '%s'.",
      paste(c(join_by, "beta"), collapse = ", ")
    )
  }
  if (!all(c(join_by) %in% names(info_snp))) {
    bigassertr::stop2(
      "Please use proper names for variables in 'info_snp'. Expected '%s'.",
      paste(unique(c(join_by)), collapse = ", ")
    )
  }

  bigassertr::message2("%s variants to be matched.", format(nrow(sumstats), big.mark = ","))

  sumstats <- sumstats[vctrs::vec_in(
    sumstats[, join_by[1:2]],
    info_snp[, join_by[1:2]]
  ), ]
  if (nrow(sumstats) == 0) bigassertr::stop2("No variant has been matched.")

  if (strand_flip) {
    is_ambiguous <- with(sumstats, paste(a0, a1) %in% c("A T", "T A", "C G", "G C"))
    bigassertr::message2(
      "%s ambiguous SNPs have been removed.",
      format(sum(is_ambiguous), big.mark = ",")
    )
    sumstats2 <- sumstats[!is_ambiguous, ]
    sumstats3 <- sumstats2
    sumstats2$`_FLIP_` <- FALSE
    sumstats3$`_FLIP_` <- TRUE
    # bigsnpr:::flip_strand is an internal function used intentionally here;
    # this generates an R CMD CHECK NOTE which is acceptable.
    sumstats3$a0 <- bigsnpr:::flip_strand(sumstats2$a0)
    sumstats3$a1 <- bigsnpr:::flip_strand(sumstats2$a1)
    sumstats3 <- rbind(sumstats2, sumstats3)
  } else {
    sumstats3 <- sumstats
    sumstats3$`_FLIP_` <- FALSE
  }

  sumstats4          <- sumstats3
  sumstats3$`_REV_` <- FALSE
  sumstats4$`_REV_` <- TRUE
  sumstats4$a0      <- sumstats3$a1
  sumstats4$a1      <- sumstats3$a0
  sumstats4$beta    <- -sumstats3$beta
  sumstats4          <- rbind(sumstats3, sumstats4)

  matched <- merge(
    data.table::as.data.table(sumstats4), data.table::as.data.table(info_snp),
    by  = join_by, all = FALSE, suffixes = c(".ss", "")
  )

  if (remove_dups) {
    dups <- vctrs::vec_duplicate_detect(matched[, c("rsid")])
    if (any(dups)) {
      matched <- matched[!dups, ]
      bigassertr::message2("Some duplicates were removed.")
    }
  }

  bigassertr::message2(
    "%s variants have been matched; %s were flipped and %s were reversed.",
    format(nrow(matched), big.mark = ","),
    format(sum(matched$`_FLIP_`), big.mark = ","),
    format(sum(matched$`_REV_`), big.mark = ",")
  )

  if (nrow(matched) < min_match) bigassertr::stop2("Not enough variants have been matched.")

  if (!return_flip_and_rev) {
    matched[, c("_FLIP_", "_REV_") := NULL]
  }

  as.data.frame(matched[order(chr)])
}
