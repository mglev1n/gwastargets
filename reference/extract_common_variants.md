# Extract variants common to all GWAS studies

Opens multiple parquet files, filters to variants with minor allele
frequency above `maf_threshold`, and returns only variants present in
every file. Optionally selects one variant per megabase window per
chromosome for LD-thinned analyses.

## Usage

``` r
extract_common_variants(
  parquet_files,
  trait_type,
  maf_threshold = 0.05,
  select_per_mb_window = FALSE,
  window_size_mb = 1,
  seed = NULL,
  eaf_col = EAF,
  chr_col = CHR,
  pos_col = POS_38,
  rsid_col = RSID,
  ea_col = EffectAllele,
  oa_col = OtherAllele,
  n_col = N
)
```

## Arguments

- parquet_files:

  Character vector of paths to parquet files.

- trait_type:

  `"binary"` or `"quantitative"` (required).

- maf_threshold:

  Minimum minor allele frequency. Default `0.05`.

- select_per_mb_window:

  Logical; if `TRUE`, select one variant per `window_size_mb` Mb window
  per chromosome. Default `FALSE`.

- window_size_mb:

  Window size in megabases for thinning. Default `1`.

- seed:

  Random seed for reproducible window selection. Default `NULL`.

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

  Bare column name for N. Default `N`.

## Value

A
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
with columns `file`, chromosome, position, rsid, effect allele, other
allele, and `minor_af`. One row per variant × study combination.

## Details

Extract common variants present in all studies

## Examples

``` r
if (FALSE) { # \dontrun{
common <- extract_common_variants(
  parquet_files = c("study1.parquet", "study2.parquet"),
  trait_type    = "binary"
)
} # }
```
