# Build a column rename map from per-cohort manifest values

Constructs a named character vector mapping non-standard column names
(as they appear in the raw summary statistics file) to harmonised names.
Each `col_*` argument corresponds to a single harmonised column. Pass
the cohort-specific raw column name (a string) to override the default
mapping, or leave as `NA` to skip.

The returned vector is suitable for the `column_map` argument of
[`harmonize_sumstats_headers()`](https://mglev1n.github.io/gwastargets/reference/harmonize_sumstats_headers.md).

## Usage

``` r
build_column_map(
  col_chr = NA_character_,
  col_pos = NA_character_,
  col_rsid = NA_character_,
  col_effect_allele = NA_character_,
  col_other_allele = NA_character_,
  col_beta = NA_character_,
  col_se = NA_character_,
  col_p = NA_character_,
  col_eaf = NA_character_,
  col_n = NA_character_,
  col_n_cases = NA_character_,
  col_n_controls = NA_character_
)
```

## Arguments

- col_chr:

  Column name mapping to `CHR`.

- col_pos:

  Column name mapping to `POS`.

- col_rsid:

  Column name mapping to `RSID`.

- col_effect_allele:

  Column name mapping to `EffectAllele`.

- col_other_allele:

  Column name mapping to `OtherAllele`.

- col_beta:

  Column name mapping to `B`.

- col_se:

  Column name mapping to `SE`.

- col_p:

  Column name mapping to `P`.

- col_eaf:

  Column name mapping to `EFF_ALL_FREQ`.

- col_n:

  Column name mapping to `N`.

- col_n_cases:

  Column name mapping to `N_CASES`.

- col_n_controls:

  Column name mapping to `N_CONTROLS`.

## Value

A named character vector where names are source column names and values
are target harmonised names, or `NULL` if all arguments are `NA`.

## Examples

``` r
build_column_map(col_eaf = "MY_FREQ", col_beta = "BETA_VAL")
#>       BETA_VAL        MY_FREQ 
#>            "B" "EFF_ALL_FREQ" 
# c("MY_FREQ" = "EFF_ALL_FREQ", "BETA_VAL" = "B")

build_column_map()
#> NULL
# NULL
```
