#' Clean GWAS summary statistics
#'
#' @title Clean GWAS summary statistics with tidyGWAS
#'
#' @description
#' Harmonises column headers via [harmonize_sumstats_headers()], applies a
#' minor-allele frequency pre-filter (EFF_ALL_FREQ between 0.001 and 0.999),
#' and passes the result through [tidyGWAS::tidyGWAS()] for full QC. The
#' tidyGWAS log file is copied to `logging_path`.
#'
#' @param sumstats_file Path to a raw GWAS summary statistics file (plain text
#'   or `.gz`).
#' @param logging_path Directory where the tidyGWAS log file will be saved.
#' @param ... Additional arguments forwarded to [tidyGWAS::tidyGWAS()].
#'
#' @return The cleaned summary statistics object returned by
#'   [tidyGWAS::tidyGWAS()].
#'
#' @examples
#' \dontrun{
#' cleaned <- clean_gwas(
#'   sumstats_file = "my_gwas.txt.gz",
#'   logging_path  = "logs/"
#' )
#' }
#'
#' @importFrom cli cli_alert_info
#' @importFrom dplyr filter
#' @importFrom fs file_temp dir_create path file_copy
#' @export
clean_gwas <- function(sumstats_file, logging_path, ...) {
  cli::cli_alert_info("Cleaning GWAS summary statistics: {sumstats_file}")

  out_dir <- fs::file_temp()

  sumstats_cleaned <- tidyGWAS::tidyGWAS(
    tbl = harmonize_sumstats_headers(sumstats_file) |>
      dplyr::filter(EFF_ALL_FREQ > 0.001 & 1 - EFF_ALL_FREQ > 0.001),
    output_dir = out_dir,
    logfile    = TRUE,
    ...
  )

  fs::dir_create(logging_path, recurse = TRUE)
  log_file <- fs::path(logging_path, basename(sumstats_file), ext = "log")
  fs::file_copy(fs::path(out_dir, "tidyGWAS_logfile.txt"), log_file, overwrite = TRUE)

  return(sumstats_cleaned)
}
