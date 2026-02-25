# Harmonize GWAS summary statistics headers

Reads a GWAS summary statistics text file (plain or `.gz`) and renames
columns to a standardised vocabulary. Unmapped columns are left
unchanged.

## Usage

``` r
harmonize_sumstats_headers(input_file)
```

## Arguments

- input_file:

  Path to a plain-text or gzip-compressed summary statistics file
  readable by
  [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html).

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
