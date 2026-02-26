#' Prepare GWAS summary statistics for meta-analysis
#'
#' @title Prepare and QC GWAS summary statistics
#'
#' @description
#' A high-level wrapper that runs the full per-cohort preparation pipeline:
#' (1) cleans the raw file via [clean_gwas()], (2) validates trait-type columns
#' via `assert_trait_columns()`, (3) matches to HapMap3 SNPs via
#' [snp_match_munge()], (4) estimates the LDSC intercept via
#' [ldscr::ldsc_h2()], and (5) writes an `.parquet` file with LDSC-adjusted
#' (or unadjusted) SE, Z, and P columns.
#'
#' @param sumstats_file Path to a raw GWAS summary statistics file.
#' @param hm3 Data frame of HapMap3 SNPs with columns `SNP`, `A1`, `A2`.
#' @param ancestry Ancestry label (e.g. `"EUR"`). Passed to
#'   [ldscr::ldsc_h2()].
#' @param output_path Directory where the output parquet file is written.
#' @param trait_type `"binary"` or `"quantitative"` (required, no default).
#' @param n_col Bare column name for sample size. Use `EffectiveN` for binary
#'   traits and `N` for quantitative traits (required, no default).
#' @param dbsnp_path Path to the dbSNP155 reference. Default
#'   `"Data/dbSNP155"`.
#' @param logging_path Directory for tidyGWAS log files. Default
#'   `"Data/logs"`.
#' @param ... Additional arguments forwarded to [clean_gwas()].
#'
#' @return The file path of the written parquet file (invisibly).
#'
#' @examples
#' \dontrun{
#' out <- prep_gwas(
#'   sumstats_file = "cohort1.txt.gz",
#'   hm3           = hm3_snps,
#'   ancestry      = "EUR",
#'   output_path   = "Data/prepped",
#'   trait_type    = "binary",
#'   n_col         = EffectiveN
#' )
#' }
#'
#' @importFrom cli cli_abort cli_warn cli_alert_info
#' @importFrom dplyr mutate select
#' @importFrom fs dir_create path
#' @importFrom rlang enquo as_name missing_arg
#' @importFrom stringr str_to_upper str_replace
#' @importFrom arrow write_parquet
#' @importFrom stats pnorm
#' @export
prep_gwas <- function(sumstats_file, hm3, ancestry, output_path,
                      trait_type,
                      n_col,
                      dbsnp_path   = "Data/dbSNP155",
                      logging_path = "Data/logs",
                      ...) {

  trait_type <- match.arg(trait_type, choices = c("binary", "quantitative"))

  # n_col must be explicitly provided
  if (missing(n_col)) {
    cli::cli_abort(c(
      "{.arg n_col} must be specified explicitly.",
      "i" = "Use {.code n_col = EffectiveN} for binary traits",
      "i" = "Use {.code n_col = N} for quantitative traits"
    ))
  }

  # Cross-validate trait_type and n_col
  n_col_name <- rlang::as_name(rlang::enquo(n_col))
  if (trait_type == "binary" && n_col_name == "N") {
    cli::cli_abort(c(
      "{.arg trait_type} = {.val binary} expects {.arg n_col} = {.val EffectiveN}, not {.val N}",
      "i" = "For quantitative traits, set {.arg trait_type} = {.val quantitative}"
    ))
  }
  if (trait_type == "quantitative" && n_col_name == "EffectiveN") {
    cli::cli_warn(c(
      "{.arg trait_type} = {.val quantitative} with {.arg n_col} = {.val EffectiveN}",
      "i" = "This is unusual -- did you mean {.arg trait_type} = {.val binary}?"
    ))
  }

  ancestry <- stringr::str_to_upper(ancestry)
  fs::dir_create(output_path, recurse = TRUE)

  cli::cli_alert_info("Preparing {sumstats_file}")
  gwas_cleaned <- clean_gwas(
    sumstats_file = sumstats_file,
    dbsnp_path    = dbsnp_path,
    logging_path  = logging_path,
    ...
  )

  # Assert that the cleaned data has the columns expected for this trait_type
  assert_trait_columns(names(gwas_cleaned), trait_type)

  gwas_munged <- snp_match_munge(
    sumstats = gwas_cleaned |>
      dplyr::select(chr = CHR, rsid = RSID, a0 = OtherAllele, a1 = EffectAllele,
                    beta = Z, N = {{ n_col }}, EAF),
    info_snp    = hm3 |> dplyr::select(rsid = SNP, a1 = A1, a0 = A2),
    join_by_pos = FALSE,
    match.min.prop = 0
  ) |>
    dplyr::select(SNP = rsid, N, Z = beta, A1 = a1, A2 = a0, EAF)

  # LDSC intercept correction -- identical for binary and quantitative.
  # sample_prev / pop_prev default to NA in ldscr::ldsc_h2(), so
  # observed-scale h2 is computed for both trait types. The primary
  # purpose here is the intercept, not the h2 estimate itself.
  gwas_ldsc_res <- ldscr::ldsc_h2(
    munged_sumstats = gwas_munged,
    ancestry        = ancestry
  )

  out_path <- fs::path(
    output_path,
    stringr::str_replace(basename(sumstats_file), "\\.txt\\.gz$", ""),
    ext = "parquet"
  )

  if (trait_type == "quantitative") {
    gwas_cleaned <- gwas_cleaned |> dplyr::mutate(EffectiveN = N)
  }

  if (gwas_ldsc_res$intercept > 1) {
    cli::cli_alert_info("LDSC intercept > 1: Adjusting SE, Z, and P-values")
    gwas_cleaned |>
      dplyr::mutate(intercept = gwas_ldsc_res$intercept) |>
      dplyr::mutate(
        SE_ldsc_adjusted = SE * sqrt(gwas_ldsc_res$intercept),
        Z_ldsc_adjusted  = B / SE_ldsc_adjusted,
        P_ldsc_adjusted  = 2 * stats::pnorm(-abs(Z_ldsc_adjusted)),
        ldsc_adjustment  = TRUE
      ) |>
      arrow::write_parquet(out_path)
  } else {
    cli::cli_alert_info("LDSC intercept <= 1: No adjustment necessary")
    gwas_cleaned |>
      dplyr::mutate(intercept = gwas_ldsc_res$intercept) |>
      dplyr::mutate(
        SE_ldsc_adjusted = SE,
        Z_ldsc_adjusted  = Z,
        P_ldsc_adjusted  = P,
        ldsc_adjustment  = FALSE
      ) |>
      arrow::write_parquet(out_path)
  }

  out_path
}
