#' Harmonize summary statistics column headers
#'
#' @title Harmonize GWAS summary statistics headers
#'
#' @description
#' Reads a GWAS summary statistics text file (plain or `.gz`) and renames
#' columns to a standardised vocabulary. Unmapped columns are left unchanged.
#'
#' @param input_file Path to a plain-text or gzip-compressed summary statistics
#'   file readable by [data.table::fread()].
#'
#' @return A [data.table::data.table()] with harmonised column names.
#'
#' @examples
#' \dontrun{
#' dt <- harmonize_sumstats_headers("my_gwas.txt.gz")
#' }
#'
#' @importFrom cli cli_alert_info cli_alert_success
#' @importFrom data.table fread
#' @export
harmonize_sumstats_headers <- function(input_file) {
  header_map <- c(
    # Case / control N
    "N_CASE"    = "N_CASES",
    "N_CAS"     = "N_CASES",
    "N_CONTROL" = "N_CONTROLS",
    "N_CTRL"    = "N_CONTROLS",
    "N_CON"     = "N_CONTROLS",

    # Allele frequency
    "AF_CASE"          = "AF_CASES",
    "AF_CTRL"          = "AF_CONTROLS",
    "AC_EFFECT_ALLELE" = "AC_CASES",
    "AC_ALLELE2"       = "AC_CONTROLS",

    # Imputation quality
    "RSQ" = "IMP_QUALITY",

    # Total N (quantitative)
    "TOTAL_N"   = "N",
    "N_SAMPLES" = "N",

    # Effect size
    "BETA"      = "B",
    "LOG_ODDS"  = "B",

    # P-value
    "PVAL"    = "P",
    "P_VALUE" = "P",

    # Effect allele frequency
    "FREQ"      = "EFF_ALL_FREQ",
    "A1FREQ"    = "EFF_ALL_FREQ",
    "EFFECT_AF" = "EFF_ALL_FREQ"
  )

  cli::cli_alert_info("Reading file: {input_file}")
  data <- data.table::fread(input_file)
  original_cols <- colnames(data)
  cli::cli_alert_info("Original columns: {paste(original_cols, collapse = ', ')}")

  new_cols <- original_cols
  for (old_name in names(header_map)) {
    if (old_name %in% new_cols) {
      new_cols[new_cols == old_name] <- header_map[old_name]
      cli::cli_alert_success("Mapped: {old_name} -> {header_map[old_name]}")
    }
  }

  colnames(data) <- new_cols
  data
}
