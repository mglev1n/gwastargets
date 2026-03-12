#' Generate Quarto report chunks for a GWAS meta-analysis trait
#'
#' @description
#' Generates quarto-compatible R code chunks as a character string for
#' inclusion in multi-trait reports. Manhattan plots use
#' `knitr::include_graphics()` with PDF file targets for fast rendering.
#' Precision plots and LDSC heatmaps use native `tar_read()` calls.
#' Sections with per-ancestry results (LDSC heatmaps, Manhattan plots)
#' are organized using Quarto tabset panels.
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
#' @importFrom glue glue glue_collapse
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

  # Build LDSC tabset (only if there are per-ancestry sections)
  ldsc_section <- ""
  if (length(ancs) > 0) {
    ldsc_panels <- vapply(ancs, function(anc) {
      glue::glue(
'#### <<anc>>

```{r}
targets::tar_read(<<trait_lower>>_ldsc_rg_qc_heatmap_<<anc>>)
```
', .open = "<<", .close = ">>")
    }, character(1))

    ldsc_section <- glue::glue(
'
### LDSC Genetic Correlation QC

::: {.panel-tabset}

<<glue::glue_collapse(ldsc_panels, sep = "\n")>>
:::
', .open = "<<", .close = ">>")
  }

  # Build Manhattan + Loci tabset
  # ALL populations tab
  all_loci <- ""
  if (include_loci) {
    all_loci <- glue::glue(
'
```{r}
targets::tar_read(<<trait_lower>>_meta_loci_ALL)
```
', .open = "<<", .close = ">>")
  }

  all_tab <- glue::glue(
'#### All Populations

```{r}
knitr::include_graphics(targets::tar_read(<<trait_lower>>_meta_manhattan_pdf_ALL))
```
<<all_loci>>', .open = "<<", .close = ">>")

  # Per-ancestry tabs
  ancestry_tabs <- ""
  if (length(ancs) > 0) {
    ancestry_panels <- vapply(ancs, function(anc) {
      anc_loci <- ""
      if (include_loci) {
        anc_loci <- glue::glue(
'
```{r}
targets::tar_read(<<trait_lower>>_meta_loci_<<anc>>)
```
', .open = "<<", .close = ">>")
      }
      glue::glue(
'#### <<anc>>

```{r}
knitr::include_graphics(targets::tar_read(<<trait_lower>>_meta_manhattan_pdf_<<anc>>))
```
<<anc_loci>>', .open = "<<", .close = ">>")
    }, character(1))

    ancestry_tabs <- glue::glue(
'
<<glue::glue_collapse(ancestry_panels, sep = "\n")>>', .open = "<<", .close = ">>")
  }

  glue::glue(
'## <<trait_upper>> Results

### Cohort Overview

```{r}
targets::tar_read(<<trait_lower>>_cohort_overview)
```

### Precision Plot

```{r}
targets::tar_read(<<trait_lower>>_precision_plot)
```
<<ldsc_section>>
### Manhattan Plots & Genome-wide Significant Loci

::: {.panel-tabset}

<<all_tab>>
<<ancestry_tabs>>
:::
', .open = "<<", .close = ">>")
}
