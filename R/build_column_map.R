#' Build a column rename map from per-cohort manifest values
#'
#' @description
#' Constructs a named character vector mapping non-standard column names
#' (as they appear in the raw summary statistics file) to harmonised names.
#' Each `col_*` argument corresponds to a single harmonised column. Pass the
#' cohort-specific raw column name (a string) to override the default mapping,
#' or leave as `NA` to skip.
#'
#' The returned vector is suitable for the `column_map` argument of
#' [harmonize_sumstats_headers()].
#'
#' @param col_chr Column name mapping to `CHR`.
#' @param col_pos Column name mapping to `POS`.
#' @param col_rsid Column name mapping to `RSID`.
#' @param col_effect_allele Column name mapping to `EffectAllele`.
#' @param col_other_allele Column name mapping to `OtherAllele`.
#' @param col_beta Column name mapping to `B`.
#' @param col_se Column name mapping to `SE`.
#' @param col_p Column name mapping to `P`.
#' @param col_eaf Column name mapping to `EFF_ALL_FREQ`.
#' @param col_n Column name mapping to `N`.
#' @param col_n_cases Column name mapping to `N_CASES`.
#' @param col_n_controls Column name mapping to `N_CONTROLS`.
#'
#' @return A named character vector where names are source column names and
#'   values are target harmonised names, or `NULL` if all arguments are `NA`.
#'
#' @examples
#' build_column_map(col_eaf = "MY_FREQ", col_beta = "BETA_VAL")
#' # c("MY_FREQ" = "EFF_ALL_FREQ", "BETA_VAL" = "B")
#'
#' build_column_map()
#' # NULL
#'
#' @export
build_column_map <- function(col_chr = NA_character_,
                             col_pos = NA_character_,
                             col_rsid = NA_character_,
                             col_effect_allele = NA_character_,
                             col_other_allele = NA_character_,
                             col_beta = NA_character_,
                             col_se = NA_character_,
                             col_p = NA_character_,
                             col_eaf = NA_character_,
                             col_n = NA_character_,
                             col_n_cases = NA_character_,
                             col_n_controls = NA_character_) {
  targets <- c(
    col_chr            = "CHR",
    col_pos            = "POS",
    col_rsid           = "RSID",
    col_effect_allele  = "EffectAllele",
    col_other_allele   = "OtherAllele",
    col_beta           = "B",
    col_se             = "SE",
    col_p              = "P",
    col_eaf            = "EFF_ALL_FREQ",
    col_n              = "N",
    col_n_cases        = "N_CASES",
    col_n_controls     = "N_CONTROLS"
  )

  sources <- c(
    col_chr            = col_chr,
    col_pos            = col_pos,
    col_rsid           = col_rsid,
    col_effect_allele  = col_effect_allele,
    col_other_allele   = col_other_allele,
    col_beta           = col_beta,
    col_se             = col_se,
    col_p              = col_p,
    col_eaf            = col_eaf,
    col_n              = col_n,
    col_n_cases        = col_n_cases,
    col_n_controls     = col_n_controls
  )

  keep <- !is.na(sources)
  if (!any(keep)) return(NULL)

  setNames(targets[names(sources[keep])], sources[keep])
}
