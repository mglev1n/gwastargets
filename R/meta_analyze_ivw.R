#' IVW meta-analysis across cohorts
#'
#' @title Inverse-variance weighted GWAS meta-analysis
#'
#' @description
#' Performs a fixed-effects inverse-variance weighted (IVW) meta-analysis
#' across multiple cohort parquet files. Processes each chromosome via
#' `meta_analyze_chromosome()`, applies a minor allele count (MAC) filter,
#' and returns heterogeneity statistics (Q, I2) alongside the weighted effect
#' estimates.
#'
#' @param parquet_files Character vector of paths to prepared parquet files
#'   (output of [prep_gwas()]).
#' @param trait_type `"binary"` or `"quantitative"` (required).
#' @param chromosomes Integer vector of chromosomes to process. Default `1:22`.
#' @param min_mac Minimum minor allele count. Variants below this threshold are
#'   excluded. Default `100`.
#' @param se_col Bare column name for SE. Default `SE`.
#' @param beta_col Bare column name for effect size. Default `B`.
#' @param eaf_col Bare column name for EAF. Default `EAF`.
#' @param chr_col Bare column name for chromosome. Default `CHR`.
#' @param pos_col Bare column name for position. Default `POS_38`.
#' @param rsid_col Bare column name for rsid. Default `RSID`.
#' @param ea_col Bare column name for effect allele. Default `EffectAllele`.
#' @param oa_col Bare column name for other allele. Default `OtherAllele`.
#' @param n_col Bare column name for total N. Default `N`.
#' @param effective_n_col Bare column name for effective N. Default
#'   `EffectiveN`.
#' @param case_n_col Bare column name for case N. Default `CaseN`.
#' @param control_n_col Bare column name for control N. Default `ControlN`.
#'
#' @return A [tibble::tibble()] with columns: `CHR`, `POS_38`, `RSID`,
#'   `EffectAllele`, `OtherAllele`, `B`, `SE`, `p_value`, `z_score`, `EAF`,
#'   `n_contributions`, `N`, `EffectiveN`, `Q`, `Q_df`, `Q_pval`, `I2`.
#'   Binary traits additionally include `CaseN` and `ControlN`.
#'
#' @examples
#' \dontrun{
#' meta <- meta_analyze_ivw(
#'   parquet_files = c("study1.parquet", "study2.parquet"),
#'   trait_type    = "binary",
#'   chromosomes   = 1:22
#' )
#' }
#'
#' @importFrom arrow open_dataset schema
#' @importFrom cli cli_alert_info cli_abort
#' @importFrom dplyr filter select mutate
#' @importFrom purrr map_dfr
#' @importFrom rlang enquo as_name
#' @export
meta_analyze_ivw <- function(parquet_files,
                             trait_type,
                             chromosomes     = 1:22,
                             min_mac         = 100,
                             se_col          = SE,
                             beta_col        = B,
                             eaf_col         = EAF,
                             chr_col         = CHR,
                             pos_col         = POS_38,
                             rsid_col        = RSID,
                             ea_col          = EffectAllele,
                             oa_col          = OtherAllele,
                             n_col           = N,
                             effective_n_col = EffectiveN,
                             case_n_col      = CaseN,
                             control_n_col   = ControlN) {

  trait_type <- match.arg(trait_type, choices = c("binary", "quantitative"))

  se_col_name     <- rlang::as_name(rlang::enquo(se_col))
  beta_col_name   <- rlang::as_name(rlang::enquo(beta_col))
  eaf_col_name    <- rlang::as_name(rlang::enquo(eaf_col))
  chr_col_name    <- rlang::as_name(rlang::enquo(chr_col))
  pos_col_name    <- rlang::as_name(rlang::enquo(pos_col))
  rsid_col_name   <- rlang::as_name(rlang::enquo(rsid_col))
  ea_col_name     <- rlang::as_name(rlang::enquo(ea_col))
  oa_col_name     <- rlang::as_name(rlang::enquo(oa_col))
  n_col_name      <- rlang::as_name(rlang::enquo(n_col))
  eff_n_col_name  <- rlang::as_name(rlang::enquo(effective_n_col))
  case_n_col_name <- rlang::as_name(rlang::enquo(case_n_col))
  ctrl_n_col_name <- rlang::as_name(rlang::enquo(control_n_col))

  ds_raw <- arrow::open_dataset(parquet_files)
  schema_names <- names(arrow::schema(ds_raw))

  # Assert columns match declared trait_type before doing any work
  assert_trait_columns(schema_names, trait_type)

  # Check all required columns are present
  required_cols <- c(
    se_col_name, beta_col_name, eaf_col_name, chr_col_name,
    pos_col_name, rsid_col_name, ea_col_name, oa_col_name, n_col_name
  )
  if (trait_type == "binary") {
    required_cols <- c(required_cols, eff_n_col_name, case_n_col_name, ctrl_n_col_name)
  }
  missing_cols <- setdiff(required_cols, schema_names)
  if (length(missing_cols) > 0) {
    cli::cli_abort("Missing required columns: {.field {missing_cols}}")
  }

  # Standardise column names for meta_analyze_chromosome
  if (trait_type == "binary") {
    ds <- ds_raw |>
      dplyr::select(
        CHR           = {{ chr_col }},
        POS           = {{ pos_col }},
        RSID          = {{ rsid_col }},
        EFFECT_ALLELE = {{ ea_col }},
        OTHER_ALLELE  = {{ oa_col }},
        BETA          = {{ beta_col }},
        SE            = {{ se_col }},
        EAF           = {{ eaf_col }},
        CASE_N        = {{ case_n_col }},
        CONTROL_N     = {{ control_n_col }},
        N             = {{ n_col }},
        EFFECTIVE_N   = {{ effective_n_col }}
      )
  } else {
    ds <- ds_raw |>
      dplyr::select(
        CHR           = {{ chr_col }},
        POS           = {{ pos_col }},
        RSID          = {{ rsid_col }},
        EFFECT_ALLELE = {{ ea_col }},
        OTHER_ALLELE  = {{ oa_col }},
        BETA          = {{ beta_col }},
        SE            = {{ se_col }},
        EAF           = {{ eaf_col }},
        N             = {{ n_col }},
        EFFECTIVE_N   = {{ effective_n_col }}
      )
  }

  cli::cli_alert_info("Processing {length(chromosomes)} chromosome{?s}...")
  cli::cli_alert_info("Using SE column: {se_col_name}")
  cli::cli_alert_info("Trait type: {trait_type}")

  results <- purrr::map_dfr(chromosomes, function(chr) {
    cli::cli_alert_info("Processing chromosome {chr}")
    meta_analyze_chromosome(ds, chr, min_mac, trait_type)
  })

  return(results)
}
