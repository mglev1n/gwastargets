# Generate targets pipeline code for GWAS meta-analysis

Generates a complete `targets` / `tarchetypes` pipeline as a character
string, ready to be pasted into a `_targets.R` file. The generated code
covers per-cohort preparation, LDSC heritability estimation,
per-ancestry and all-populations IVW meta-analysis, loci extraction,
Manhattan plots, and MR-MEGA MDS. Validates `trait_type`, `n_col`, and
`manifest_df` before generating any code. The manifest is serialized as
a
[`tibble::tribble()`](https://tibble.tidyverse.org/reference/tribble.html)
call in the generated output so the pipeline is self-contained.

## Usage

``` r
generate_gwas_meta_pipeline(
  trait,
  trait_type,
  n_col,
  manifest_df,
  output_base_dir = "Data"
)
```

## Arguments

- trait:

  Trait name used in target names (e.g. `"CAD"`).

- trait_type:

  `"binary"` or `"quantitative"` (required, no default).

- n_col:

  Sample size column: `"EffectiveN"` for binary traits, `"N"` for
  quantitative traits (required, no default).

- manifest_df:

  A data frame describing the cohort files to include. Required columns:
  `path` (full path to raw summary statistics file), `file` (basename of
  the file), `cohort` (cohort label), `ancestry` (ancestry label, e.g.
  `"EUR"`), `study` (unique identifier used as the `tar_map` target name
  – must contain no spaces or special characters). Additional columns
  are preserved in the serialized output.

- output_base_dir:

  Base output directory for intermediate files. Default `"Data"`.

## Value

A length-1 character string containing valid R code defining a `targets`
pipeline.

## Details

Generate a targets GWAS meta-analysis pipeline

## Examples

``` r
if (FALSE) { # \dontrun{
manifest <- data.frame(
  path     = c("/data/UKBB_EUR.txt.gz", "/data/MVP_AFR.txt.gz"),
  file     = c("UKBB_EUR.txt.gz",       "MVP_AFR.txt.gz"),
  cohort   = c("UKBB",                  "MVP"),
  ancestry = c("EUR",                   "AFR"),
  study    = c("UKBB_EUR",              "MVP_AFR"),
  stringsAsFactors = FALSE
)
cat(generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = manifest))
} # }
```
