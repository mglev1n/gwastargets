#' Extract genome-wide significant loci
#'
#' @title Extract genome-wide significant independent loci
#'
#' @description
#' Filters a meta-analysis data frame to genome-wide significant variants,
#' identifies independent loci via [gwasRtools::get_loci()], and annotates
#' each lead variant with the nearest gene via [gwasRtools::get_nearest_gene()].
#' Returns an empty [tibble::tibble()] when no variants pass the p-value
#' threshold.
#'
#' @param df A data frame of meta-analysis results (e.g. output of
#'   [meta_analyze_ivw()]).
#' @param snp_col Bare column name for rsid. Default `RSID`.
#' @param chr_col Bare column name for chromosome. Default `CHR`.
#' @param pos_col Bare column name for position. Default `POS_38`.
#' @param maf_col Bare column name for minor allele frequency. Default `EAF`.
#' @param beta_col Bare column name for effect size. Default `B`.
#' @param se_col Bare column name for SE. Default `SE`.
#' @param p_col Bare column name for p-value. Default `p_value`.
#' @param p_threshold Genome-wide significance threshold. Default `5e-8`.
#' @param build Genome build for gene annotation (`37` or `38`). Default `38`.
#' @param ... Additional arguments forwarded to [gwasRtools::get_loci()].
#'
#' @return A [tibble::tibble()] of independent lead variants with gene
#'   annotations, or an empty tibble if no variants pass the threshold.
#'
#' @examples
#' \dontrun{
#' loci <- extract_loci(meta_results, p_threshold = 5e-8)
#' }
#'
#' @importFrom cli cli_abort cli_alert_info cli_alert_warning
#'   cli_alert_success cli_progress_step
#' @importFrom dplyr filter mutate
#' @importFrom rlang enquo as_name
#' @importFrom scales comma
#' @importFrom tibble as_tibble tibble
#' @export
extract_loci <- function(df,
                         snp_col     = RSID,
                         chr_col     = CHR,
                         pos_col     = POS_38,
                         maf_col     = EAF,
                         beta_col    = B,
                         se_col      = SE,
                         p_col       = p_value,
                         p_threshold = 5e-8,
                         build       = 38,
                         ...) {

  snp_col  <- rlang::enquo(snp_col)
  chr_col  <- rlang::enquo(chr_col)
  pos_col  <- rlang::enquo(pos_col)
  maf_col  <- rlang::enquo(maf_col)
  beta_col <- rlang::enquo(beta_col)
  se_col   <- rlang::enquo(se_col)
  p_col    <- rlang::enquo(p_col)

  snp_col_str  <- rlang::as_name(snp_col)
  chr_col_str  <- rlang::as_name(chr_col)
  pos_col_str  <- rlang::as_name(pos_col)
  maf_col_str  <- rlang::as_name(maf_col)
  beta_col_str <- rlang::as_name(beta_col)
  se_col_str   <- rlang::as_name(se_col)
  p_col_str    <- rlang::as_name(p_col)

  required_cols <- c(snp_col_str, chr_col_str, pos_col_str, maf_col_str,
                     beta_col_str, se_col_str, p_col_str)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    cli::cli_abort("Missing required columns: {.field {missing_cols}}")
  }

  initial_rows <- nrow(df)
  cli::cli_alert_info("Starting with {.val {scales::comma(initial_rows)}} rows")

  na_chr <- sum(is.na(df[[chr_col_str]]))
  na_pos <- sum(is.na(df[[pos_col_str]]))
  na_p   <- sum(is.na(df[[p_col_str]]))
  if (na_chr > 0) cli::cli_warn("Found {.val {scales::comma(na_chr)}} rows with missing {.field {chr_col_str}}")
  if (na_pos > 0) cli::cli_warn("Found {.val {scales::comma(na_pos)}} rows with missing {.field {pos_col_str}}")
  if (na_p   > 0) cli::cli_warn("Found {.val {scales::comma(na_p)}} rows with missing {.field {p_col_str}}")

  valid_p <- df[[p_col_str]][!is.na(df[[p_col_str]])]
  if (any(valid_p < 0 | valid_p > 1, na.rm = TRUE)) {
    cli::cli_alert_warning("Some p-values are outside the expected range [0, 1]")
  }

  cli::cli_progress_step("Filtering by p-value threshold ({.val {p_threshold}})")
  df_filtered    <- df[df[[p_col_str]] < p_threshold & !is.na(df[[p_col_str]]), ]
  after_p_filter <- nrow(df_filtered)
  removed_p      <- initial_rows - after_p_filter

  if (removed_p > 0) {
    cli::cli_alert_info("Removed {.val {scales::comma(removed_p)}} rows; {.val {scales::comma(after_p_filter)}} remaining")
  }
  if (after_p_filter == 0) {
    cli::cli_alert_warning("No genome-wide significant variants found (p < {.val {p_threshold}})")
    return(tibble::tibble())
  }

  cli::cli_progress_step("Removing rows with missing chromosome or position values")
  df_filtered      <- df_filtered[!is.na(df_filtered[[chr_col_str]]) & !is.na(df_filtered[[pos_col_str]]), ]
  after_pos_filter <- nrow(df_filtered)

  if (after_pos_filter == 0) {
    cli::cli_abort("No rows remaining after filtering. Check your p-value threshold and data quality.")
  }

  cli::cli_progress_step("Identifying independent loci")
  df_loci <- df_filtered |>
    gwasRtools::get_loci(
      snp_col  = snp_col_str,
      chr_col  = chr_col_str,
      pos_col  = pos_col_str,
      maf_col  = maf_col_str,
      beta_col = beta_col_str,
      se_col   = se_col_str,
      p_col    = p_col_str,
      ...
    ) |>
    dplyr::filter(lead)

  n_loci <- nrow(df_loci)
  if (n_loci > 0) {
    cli::cli_alert_success("Found {.val {n_loci}} independent loci")
  } else {
    cli::cli_alert_warning("No independent loci identified.")
    return(df_loci)
  }

  cli::cli_progress_step("Annotating with nearest genes (build {.val {build}})")
  result <- df_loci |>
    gwasRtools::get_nearest_gene(
      snp_col = snp_col_str,
      chr_col = chr_col_str,
      pos_col = pos_col_str,
      build   = build
    )

  cli::cli_alert_success("Analysis complete: {.val {nrow(result)}} loci with gene annotations")
  return(tibble::as_tibble(result))
}
