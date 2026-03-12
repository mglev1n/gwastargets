# Generate targets pipeline code for GWAS meta-analysis

Generates a complete `targets` / `tarchetypes` pipeline as a character
string, ready to be pasted into a `_targets.R` file. The generated code
covers per-cohort preparation, LDSC heritability estimation,
per-ancestry and all-populations IVW meta-analysis, loci extraction,
Manhattan plots (with PDF and high-DPI PNG export targets for fast
report rendering), and MR-MEGA MDS. Validates `trait_type`, `n_col`, and
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
  hm3_path,
  dbsnp_path,
  crew_controller = NULL,
  output_base_dir = "Data",
  manhattan_width = 16,
  manhattan_height = 6,
  manhattan_dpi = 300,
  output_file = NULL
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
  `"EUR"`), `study` (unique human-readable identifier; must contain no
  spaces or special characters; does not need to include the ancestry).
  A `tar_name` column (`{study}_{ancestry}`) is automatically added and
  used as the `tar_map` target name so that ancestry-based target
  selectors work reliably. Additional columns are preserved in the
  serialized output.

- hm3_path:

  Path to the LDSC HapMap3 SNP list file (e.g. `w_hm3.snplist`).
  Required.

- dbsnp_path:

  Path to the dbSNP155 reference directory passed to
  [`prep_gwas()`](https://mglev1n.github.io/gwastargets/reference/prep_gwas.md).
  Required.

- crew_controller:

  Name of the `crew` controller to use for resource-intensive targets
  (`meta_ALL` and `meta_common_variants_ALL`). If `NULL` (default), no
  `resources` block is added to those targets.

- output_base_dir:

  Base output directory for intermediate files. Default `"Data"`.

- manhattan_width:

  Width in inches for Manhattan plot export. Default `16`.

- manhattan_height:

  Height in inches for Manhattan plot export. Default `6`.

- manhattan_dpi:

  DPI for high-resolution Manhattan plot PNG export. Default `300`.

- output_file:

  Optional file path to write the generated code to. When provided,
  parent directories are created automatically, the code is wrapped in a
  ```` ```{targets <trait>-pipeline} ```` chunk, and the raw pipeline
  code string is returned invisibly. When `NULL` (default), the string
  is returned visibly.

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
                                n_col = "EffectiveN", manifest_df = manifest,
                                hm3_path   = "/path/to/w_hm3.snplist",
                                dbsnp_path = "/path/to/dbSNP155"))
} # }
```
