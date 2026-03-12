# Harmonize GWAS summary statistics headers

Reads a GWAS summary statistics text file (plain or `.gz`) and renames
columns to a standardised vocabulary. Unmapped columns are left
unchanged.

## Usage

``` r
harmonize_sumstats_headers(input_file, column_map = NULL)
```

## Arguments

- input_file:

  Path to a plain-text or gzip-compressed summary statistics file
  readable by
  [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html).

- column_map:

  Optional named character vector of per-cohort column renames. Names
  are source column names in the file; values are target harmonised
  names. Applied **before** the built-in dictionary so that user
  overrides take priority. Use
  [`build_column_map()`](https://mglev1n.github.io/gwastargets/reference/build_column_map.md)
  to construct this vector from manifest `col_*` values. Default `NULL`
  (no custom mapping).

## Value

A
[`data.table::data.table()`](https://rdrr.io/pkg/data.table/man/data.table.html)
with harmonised column names.

## Details

Harmonize summary statistics column headers

## Examples

``` r
if (FALSE) { # \dontrun{
dt <- harmonize_sumstats_headers("my_gwas.txt.gz")
} # }
```
