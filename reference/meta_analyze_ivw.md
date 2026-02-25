# Inverse-variance weighted GWAS meta-analysis

Performs a fixed-effects inverse-variance weighted (IVW) meta-analysis
across multiple cohort parquet files. Processes each chromosome via
`meta_analyze_chromosome()`, applies a minor allele count (MAC) filter,
and returns heterogeneity statistics (Q, I2) alongside the weighted
effect estimates.

## Usage

``` r
meta_analyze_ivw(
  parquet_files,
  trait_type,
  chromosomes = 1:22,
  min_mac = 100,
  se_col = SE,
  beta_col = B,
  eaf_col = EAF,
  chr_col = CHR,
  pos_col = POS_38,
  rsid_col = RSID,
  ea_col = EffectAllele,
  oa_col = OtherAllele,
  n_col = N,
  effective_n_col = EffectiveN,
  case_n_col = CaseN,
  control_n_col = ControlN
)
```

## Arguments

- parquet_files:

  Character vector of paths to prepared parquet files (output of
  [`prep_gwas()`](https://mglev1n.github.io/gwastargets/reference/prep_gwas.md)).

- trait_type:

  `"binary"` or `"quantitative"` (required).

- chromosomes:

  Integer vector of chromosomes to process. Default `1:22`.

- min_mac:

  Minimum minor allele count. Variants below this threshold are
  excluded. Default `100`.

- se_col:

  Bare column name for SE. Default `SE`.

- beta_col:

  Bare column name for effect size. Default `B`.

- eaf_col:

  Bare column name for EAF. Default `EAF`.

- chr_col:

  Bare column name for chromosome. Default `CHR`.

- pos_col:

  Bare column name for position. Default `POS_38`.

- rsid_col:

  Bare column name for rsid. Default `RSID`.

- ea_col:

  Bare column name for effect allele. Default `EffectAllele`.

- oa_col:

  Bare column name for other allele. Default `OtherAllele`.

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
with columns: `CHR`, `POS_38`, `RSID`, `EffectAllele`, `OtherAllele`,
`B`, `SE`, `p_value`, `z_score`, `EAF`, `n_contributions`, `N`,
`EffectiveN`, `Q`, `Q_df`, `Q_pval`, `I2`. Binary traits additionally
include `CaseN` and `ControlN`.

## Details

IVW meta-analysis across cohorts

## Examples

``` r
if (FALSE) { # \dontrun{
meta <- meta_analyze_ivw(
  parquet_files = c("study1.parquet", "study2.parquet"),
  trait_type    = "binary",
  chromosomes   = 1:22
)
} # }
```
