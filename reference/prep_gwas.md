# Prepare and QC GWAS summary statistics

A high-level wrapper that runs the full per-cohort preparation pipeline:
(1) cleans the raw file via
[`clean_gwas()`](https://mglev1n.github.io/gwastargets/reference/clean_gwas.md),
(2) validates trait-type columns via `assert_trait_columns()`, (3)
matches to HapMap3 SNPs via
[`snp_match_munge()`](https://mglev1n.github.io/gwastargets/reference/snp_match_munge.md),
(4) estimates the LDSC intercept via
[`ldscr::ldsc_h2()`](https://mglev1n.github.io/ldscr/reference/ldsc_h2.html),
and (5) writes an `.parquet` file with LDSC-adjusted (or unadjusted) SE,
Z, and P columns.

## Usage

``` r
prep_gwas(
  sumstats_file,
  hm3,
  ancestry,
  output_path,
  trait_type,
  n_col,
  dbsnp_path = "Data/dbSNP155",
  logging_path = "Data/logs",
  ...
)
```

## Arguments

- sumstats_file:

  Path to a raw GWAS summary statistics file.

- hm3:

  Data frame of HapMap3 SNPs with columns `SNP`, `A1`, `A2`.

- ancestry:

  Ancestry label (e.g. `"EUR"`). Passed to
  [`ldscr::ldsc_h2()`](https://mglev1n.github.io/ldscr/reference/ldsc_h2.html).

- output_path:

  Directory where the output parquet file is written.

- trait_type:

  `"binary"` or `"quantitative"` (required, no default).

- n_col:

  Bare column name for sample size. Use `EffectiveN` for binary traits
  and `N` for quantitative traits (required, no default).

- dbsnp_path:

  Path to the dbSNP155 reference. Default `"Data/dbSNP155"`.

- logging_path:

  Directory for tidyGWAS log files. Default `"Data/logs"`.

- ...:

  Additional arguments forwarded to
  [`clean_gwas()`](https://mglev1n.github.io/gwastargets/reference/clean_gwas.md).

## Value

The file path of the written parquet file (invisibly).

## Details

Prepare GWAS summary statistics for meta-analysis

## Examples

``` r
if (FALSE) { # \dontrun{
out <- prep_gwas(
  sumstats_file = "cohort1.txt.gz",
  hm3           = hm3_snps,
  ancestry      = "EUR",
  output_path   = "Data/prepped",
  trait_type    = "binary",
  n_col         = EffectiveN
)
} # }
```
