#' Extract common variants present in all studies
#'
#' @title Extract variants common to all GWAS studies
#'
#' @description
#' Opens multiple parquet files, filters to variants with minor allele
#' frequency above `maf_threshold`, and returns only variants present in
#' every file. Optionally selects one variant per megabase window per
#' chromosome for LD-thinned analyses.
#'
#' @param parquet_files Character vector of paths to parquet files.
#' @param trait_type `"binary"` or `"quantitative"` (required).
#' @param maf_threshold Minimum minor allele frequency. Default `0.05`.
#' @param select_per_mb_window Logical; if `TRUE`, select one variant per
#'   `window_size_mb` Mb window per chromosome. Default `FALSE`.
#' @param window_size_mb Window size in megabases for thinning. Default `1`.
#' @param seed Random seed for reproducible window selection. Default `NULL`.
#' @param eaf_col Bare column name for EAF. Default `EAF`.
#' @param chr_col Bare column name for chromosome. Default `CHR`.
#' @param pos_col Bare column name for position. Default `POS_38`.
#' @param rsid_col Bare column name for rsid. Default `RSID`.
#' @param ea_col Bare column name for effect allele. Default `EffectAllele`.
#' @param oa_col Bare column name for other allele. Default `OtherAllele`.
#' @param n_col Bare column name for N. Default `N`.
#'
#' @return A [tibble::tibble()] with columns `file`, chromosome, position,
#'   rsid, effect allele, other allele, and `minor_af`. One row per
#'   variant × study combination.
#'
#' @examples
#' \dontrun{
#' common <- extract_common_variants(
#'   parquet_files = c("study1.parquet", "study2.parquet"),
#'   trait_type    = "binary"
#' )
#' }
#'
#' @importFrom arrow open_dataset schema
#' @importFrom cli cli_alert_info
#' @importFrom dplyr mutate filter group_by count select collect inner_join
#' @importFrom rlang enquo as_name
#' @export
extract_common_variants <- function(parquet_files,
                                    trait_type,
                                    maf_threshold        = 0.05,
                                    select_per_mb_window = FALSE,
                                    window_size_mb       = 1,
                                    seed                 = NULL,
                                    eaf_col              = EAF,
                                    chr_col              = CHR,
                                    pos_col              = POS_38,
                                    rsid_col             = RSID,
                                    ea_col               = EffectAllele,
                                    oa_col               = OtherAllele,
                                    n_col                = N) {

  trait_type <- match.arg(trait_type, choices = c("binary", "quantitative"))

  ds_raw <- arrow::open_dataset(parquet_files)
  assert_trait_columns(names(arrow::schema(ds_raw)), trait_type)

  n_files <- length(parquet_files)
  cli::cli_alert_info("Extracting variants with MAF > {maf_threshold}")

  ds <- ds_raw |>
    dplyr::mutate(
      file     = arrow:::add_filename(),
      minor_af = pmin({{ eaf_col }}, 1 - {{ eaf_col }})
    ) |>
    dplyr::filter(minor_af > maf_threshold)

  cli::cli_alert_info("Finding variants present in all {n_files} files")
  variants_in_all_files <- ds |>
    dplyr::group_by({{ chr_col }}, {{ pos_col }}, {{ ea_col }}, {{ oa_col }}) |>
    dplyr::count(name = "file_count") |>
    dplyr::filter(file_count == n_files) |>
    dplyr::select({{ chr_col }}, {{ pos_col }}, {{ ea_col }}, {{ oa_col }}) |>
    dplyr::collect()

  cli::cli_alert_info("Filtering to variants present in all studies")
  common_variants <- ds |>
    dplyr::inner_join(
      variants_in_all_files,
      by = c(
        rlang::as_name(rlang::enquo(chr_col)),
        rlang::as_name(rlang::enquo(pos_col)),
        rlang::as_name(rlang::enquo(ea_col)),
        rlang::as_name(rlang::enquo(oa_col))
      )
    ) |>
    dplyr::select(file, {{ chr_col }}, {{ pos_col }}, {{ rsid_col }},
                  {{ ea_col }}, {{ oa_col }}, minor_af) |>
    dplyr::collect() |>
    dplyr::mutate(file = basename(file))

  if (select_per_mb_window) {
    if (!is.null(seed)) set.seed(seed)
    cli::cli_alert_info("Selecting 1 variant per {window_size_mb}MB window per chromosome")
    window_size_bp <- window_size_mb * 1e6

    selected_variants <- common_variants |>
      tidytable::mutate(window = floor({{ pos_col }} / window_size_bp)) |>
      tidytable::distinct({{ chr_col }}, {{ pos_col }}, {{ ea_col }}, {{ oa_col }}, window) |>
      tidytable::group_by({{ chr_col }}, window) |>
      tidytable::slice_sample(n = 1) |>
      tidytable::ungroup() |>
      tidytable::select(-window)

    common_variants <- common_variants |>
      tidytable::inner_join(
        selected_variants,
        by = c(
          rlang::as_name(rlang::enquo(chr_col)),
          rlang::as_name(rlang::enquo(pos_col)),
          rlang::as_name(rlang::enquo(ea_col)),
          rlang::as_name(rlang::enquo(oa_col))
        )
      )
  }

  return(common_variants)
}
