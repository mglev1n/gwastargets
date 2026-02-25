# Summarize GWAS summary statistics

Opens one or more parquet files as an Arrow dataset, computes per-file
summary statistics (SNP count, EAF range, median SE, etc.), and appends
binary-trait columns (`n_cases`, `n_controls`, `n_effective`) when
`trait_type = "binary"`.

## Usage

``` r
summarize_sumstats(
  parquet_files,
  trait_type,
  min_mac = 100,
  se_col = SE,
  beta_col = B,
  eaf_col = EAF,
  chr_col = CHR,
  pos_col = POS_38,
  rsid_col = RSID,
  n_col = N,
  effective_n_col = EffectiveN,
  case_n_col = CaseN,
  control_n_col = ControlN
)
```

## Arguments

- parquet_files:

  Character vector of paths to parquet files.

- trait_type:

  `"binary"` or `"quantitative"` (required).

- min_mac:

  Minimum minor allele count filter. Default `100`.

- se_col:

  Bare column name for SE. Default `SE`.

- beta_col:

  Bare column name for effect size. Default `B`.

- eaf_col:

  Bare column name for effect allele frequency. Default `EAF`.

- chr_col:

  Bare column name for chromosome. Default `CHR`.

- pos_col:

  Bare column name for position. Default `POS_38`.

- rsid_col:

  Bare column name for rsid. Default `RSID`.

- n_col:

  Bare column name for total N. Default `N`.

- effective_n_col:

  Bare column name for effective N. Default `EffectiveN`.

- case_n_col:

  Bare column name for case N. Default `CaseN`.

- control_n_col:

  Bare column name for control N. Default `ControlN`.

## Value

A
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
with one row per file and summary columns.

## Details

Summarize GWAS summary statistics files

## Examples

``` r
if (FALSE) { # \dontrun{
summary <- summarize_sumstats(
  parquet_files = c("study1.parquet", "study2.parquet"),
  trait_type    = "binary"
)
} # }
```
