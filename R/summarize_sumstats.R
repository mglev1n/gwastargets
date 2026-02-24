#' Summarize GWAS summary statistics files
#'
#' @title Summarize GWAS summary statistics
#'
#' @description
#' Opens one or more parquet files as an Arrow dataset, computes per-file
#' summary statistics (SNP count, EAF range, median SE, etc.), and appends
#' binary-trait columns (`n_cases`, `n_controls`, `n_effective`) when
#' `trait_type = "binary"`.
#'
#' @param parquet_files Character vector of paths to parquet files.
#' @param trait_type `"binary"` or `"quantitative"` (required).
#' @param min_mac Minimum minor allele count filter. Default `100`.
#' @param se_col Bare column name for SE. Default `SE`.
#' @param beta_col Bare column name for effect size. Default `B`.
#' @param eaf_col Bare column name for effect allele frequency. Default `EAF`.
#' @param chr_col Bare column name for chromosome. Default `CHR`.
#' @param pos_col Bare column name for position. Default `POS_38`.
#' @param rsid_col Bare column name for rsid. Default `RSID`.
#' @param n_col Bare column name for total N. Default `N`.
#' @param effective_n_col Bare column name for effective N. Default
#'   `EffectiveN`.
#' @param case_n_col Bare column name for case N. Default `CaseN`.
#' @param control_n_col Bare column name for control N. Default `ControlN`.
#'
#' @return A [tibble::tibble()] with one row per file and summary columns.
#'
#' @examples
#' \dontrun{
#' summary <- summarize_sumstats(
#'   parquet_files = c("study1.parquet", "study2.parquet"),
#'   trait_type    = "binary"
#' )
#' }
#'
#' @importFrom arrow open_dataset schema add_filename
#' @importFrom dplyr mutate group_by summarize n n_distinct collect left_join
#' @export
summarize_sumstats <- function(parquet_files,
                               trait_type,
                               min_mac         = 100,
                               se_col          = SE,
                               beta_col        = B,
                               eaf_col         = EAF,
                               chr_col         = CHR,
                               pos_col         = POS_38,
                               rsid_col        = RSID,
                               n_col           = N,
                               effective_n_col = EffectiveN,
                               case_n_col      = CaseN,
                               control_n_col   = ControlN) {

  trait_type <- match.arg(trait_type, choices = c("binary", "quantitative"))

  ds <- arrow::open_dataset(parquet_files)
  assert_trait_columns(names(arrow::schema(ds)), trait_type)

  ds <- ds |>
    dplyr::mutate(file = arrow::add_filename()) |>
    dplyr::mutate(minor_af = pmin({{ eaf_col }}, 1 - {{ eaf_col }}))

  # Columns summarised for all trait types
  base_summary <- ds |>
    dplyr::group_by(file) |>
    dplyr::summarize(
      min_eaf      = min({{ eaf_col }},   na.rm = TRUE),
      max_eaf      = max({{ eaf_col }},   na.rm = TRUE),
      n_chromosome = dplyr::n_distinct({{ chr_col }}),
      n_snps       = dplyr::n(),
      n_total      = max({{ n_col }},     na.rm = TRUE),
      median_se    = median({{ se_col }}, na.rm = TRUE),
      min_beta     = min({{ beta_col }},  na.rm = TRUE),
      max_beta     = max({{ beta_col }},  na.rm = TRUE),
      median_beta  = median({{ beta_col }}, na.rm = TRUE)
    ) |>
    dplyr::collect()

  if (trait_type == "binary") {
    binary_summary <- ds |>
      dplyr::group_by(file) |>
      dplyr::summarize(
        n_cases     = max({{ case_n_col }},    na.rm = TRUE),
        n_controls  = max({{ control_n_col }}, na.rm = TRUE),
        n_effective = max({{ effective_n_col }}, na.rm = TRUE)
      ) |>
      dplyr::collect()

    base_summary <- dplyr::left_join(base_summary, binary_summary, by = "file")
  }

  base_summary
}
