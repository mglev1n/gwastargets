#' Generate Quarto report chunks for a GWAS meta-analysis trait
#'
#' @description
#' Generates quarto-compatible R code chunks as a character string for
#' inclusion in multi-trait reports. Manhattan plots use
#' `knitr::include_graphics()` with PDF file targets for fast rendering.
#' Precision plots and LDSC heatmaps use native `tar_read()` calls.
#'
#' @param trait Trait name (e.g. `"CAD"`).
#' @param manifest_df A data frame with at least `ancestry` and `cohort`
#'   columns. Same manifest passed to [generate_gwas_meta_pipeline()].
#'   Used to determine which ancestries have per-ancestry meta-analysis
#'   targets (only ancestries with 2+ cohorts).
#' @param include_loci Logical; whether to include genome-wide significant
#'   loci table chunks. Default `TRUE`.
#'
#' @return A length-1 character string containing quarto markdown with
#'   embedded R code chunks.
#'
#' @examples
#' \dontrun{
#' manifest <- data.frame(
#'   path     = c("/data/UKBB_EUR.txt.gz", "/data/MVP_AFR.txt.gz",
#'                "/data/BioVU_EUR.txt.gz"),
#'   file     = c("UKBB_EUR.txt.gz", "MVP_AFR.txt.gz", "BioVU_EUR.txt.gz"),
#'   cohort   = c("UKBB", "MVP", "BioVU"),
#'   ancestry = c("EUR", "AFR", "EUR"),
#'   study    = c("UKBB", "MVP", "BioVU"),
#'   stringsAsFactors = FALSE
#' )
#' cat(generate_report_chunks("CAD", manifest))
#' }
#'
#' @importFrom cli cli_abort
#' @importFrom glue glue
#' @importFrom dplyr add_count filter distinct pull
#' @export
generate_report_chunks <- function(trait, manifest_df, include_loci = TRUE) {

  # Validate trait

  if (missing(trait) || !is.character(trait) || length(trait) != 1 || !nzchar(trait)) {
    cli::cli_abort(c(
      "{.arg trait} must be a single non-empty character string"
    ))
  }

  # Validate manifest_df

  if (missing(manifest_df)) {
    cli::cli_abort(c(
      "{.arg manifest_df} must be provided",
      "i" = "Supply a data frame with at least {.val ancestry} and {.val cohort} columns"
    ))
  }

  required_cols <- c("ancestry", "cohort")
  missing_cols  <- setdiff(required_cols, names(manifest_df))
  if (length(missing_cols) > 0) {
    cli::cli_abort(c(
      "{.arg manifest_df} is missing required column{?s}: {.val {missing_cols}}",
      "i" = "Required columns: {.val {required_cols}}"
    ))
  }

  trait_lower <- tolower(trait)
  trait_upper <- toupper(trait)

  # Ancestries with 2+ cohorts (same logic as pipeline generator)
  ancs <- manifest_df |>
    dplyr::add_count(ancestry) |>
    dplyr::filter(n > 1) |>
    dplyr::distinct(ancestry) |>
    dplyr::pull(ancestry)

  # Build chunks
  chunks <- character()

  # Header + cohort overview + precision plot
  chunks <- c(chunks, glue::glue(
    '## <<trait_upper>> Results\n\n### Cohort Overview\n```{r}\ntargets::tar_read(<<trait_lower>>_cohort_overview)\n```\n\n### Precision Plot\n```{r}\ntargets::tar_read(<<trait_lower>>_precision_plot)\n```',
    .open = "<<", .close = ">>"
  ))

  # LDSC rg heatmaps per ancestry
  for (anc in ancs) {
    chunks <- c(chunks, glue::glue(
      '\n\n### LDSC Genetic Correlation QC - <<anc>>\n```{r}\ntargets::tar_read(<<trait_lower>>_ldsc_rg_qc_heatmap_<<anc>>)\n```',
      .open = "<<", .close = ">>"
    ))
  }

  # ALL population manhattan
  chunks <- c(chunks, glue::glue(
    '\n\n### Manhattan Plot - All Populations\n```{r}\nknitr::include_graphics(targets::tar_read(<<trait_lower>>_meta_manhattan_pdf_ALL))\n```',
    .open = "<<", .close = ">>"
  ))

  # ALL population loci
  if (include_loci) {
    chunks <- c(chunks, glue::glue(
      '\n\n### Genome-wide Significant Loci - All Populations\n```{r}\ntargets::tar_read(<<trait_lower>>_meta_loci_ALL)\n```',
      .open = "<<", .close = ">>"
    ))
  }

  # Per-ancestry manhattan + loci
  for (anc in ancs) {
    chunks <- c(chunks, glue::glue(
      '\n\n### Manhattan Plot - <<anc>>\n```{r}\nknitr::include_graphics(targets::tar_read(<<trait_lower>>_meta_manhattan_pdf_<<anc>>))\n```',
      .open = "<<", .close = ">>"
    ))
    if (include_loci) {
      chunks <- c(chunks, glue::glue(
        '\n\n### Genome-wide Significant Loci - <<anc>>\n```{r}\ntargets::tar_read(<<trait_lower>>_meta_loci_<<anc>>)\n```',
        .open = "<<", .close = ">>"
      ))
    }
  }

  paste(chunks, collapse = "")
}
