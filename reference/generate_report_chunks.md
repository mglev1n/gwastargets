# Generate Quarto report chunks for a GWAS meta-analysis trait

Generates quarto-compatible R code chunks as a character string for
inclusion in multi-trait reports. Manhattan plots use
[`knitr::include_graphics()`](https://rdrr.io/pkg/knitr/man/include_graphics.html)
with high-DPI PNG file targets for native HTML rendering at full width.
Loci tables are rendered with
[`DT::datatable()`](https://rdrr.io/pkg/DT/man/datatable.html) with a
CSV export button. Precision plots and LDSC heatmaps use native
`tar_read()` calls. Sections with per-ancestry results (LDSC heatmaps,
Manhattan plots) are organized using Quarto tabset panels.

## Usage

``` r
generate_report_chunks(
  trait,
  manifest_df,
  include_loci = TRUE,
  output_file = NULL
)
```

## Arguments

- trait:

  Trait name (e.g. `"CAD"`).

- manifest_df:

  A data frame with at least `ancestry` and `cohort` columns. Same
  manifest passed to
  [`generate_gwas_meta_pipeline()`](https://mglev1n.github.io/gwastargets/reference/generate_gwas_meta_pipeline.md).
  Used to determine which ancestries have per-ancestry meta-analysis
  targets (only ancestries with 2+ cohorts).

- include_loci:

  Logical; whether to include genome-wide significant loci table chunks.
  Default `TRUE`.

- output_file:

  Optional file path to write the generated chunks to. When provided,
  parent directories are created automatically and the string is
  returned invisibly. When `NULL` (default), the string is returned
  visibly.

## Value

A length-1 character string containing quarto markdown with embedded R
code chunks.

## Examples

``` r
if (FALSE) { # \dontrun{
manifest <- data.frame(
  path     = c("/data/UKBB_EUR.txt.gz", "/data/MVP_AFR.txt.gz",
               "/data/BioVU_EUR.txt.gz"),
  file     = c("UKBB_EUR.txt.gz", "MVP_AFR.txt.gz", "BioVU_EUR.txt.gz"),
  cohort   = c("UKBB", "MVP", "BioVU"),
  ancestry = c("EUR", "AFR", "EUR"),
  study    = c("UKBB", "MVP", "BioVU"),
  stringsAsFactors = FALSE
)
cat(generate_report_chunks("CAD", manifest))
} # }
```
