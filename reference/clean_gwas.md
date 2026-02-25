# Clean GWAS summary statistics with tidyGWAS

Harmonises column headers via
[`harmonize_sumstats_headers()`](https://mglev1n.github.io/gwastargets/reference/harmonize_sumstats_headers.md),
applies a minor-allele frequency pre-filter (EFF_ALL_FREQ between 0.001
and 0.999), and passes the result through
[`tidyGWAS::tidyGWAS()`](https://ararder.github.io/tidyGWAS/reference/tidyGWAS.html)
for full QC. The tidyGWAS log file is copied to `logging_path`.

## Usage

``` r
clean_gwas(sumstats_file, logging_path, ...)
```

## Arguments

- sumstats_file:

  Path to a raw GWAS summary statistics file (plain text or `.gz`).

- logging_path:

  Directory where the tidyGWAS log file will be saved.

- ...:

  Additional arguments forwarded to
  [`tidyGWAS::tidyGWAS()`](https://ararder.github.io/tidyGWAS/reference/tidyGWAS.html).

## Value

The cleaned summary statistics object returned by
[`tidyGWAS::tidyGWAS()`](https://ararder.github.io/tidyGWAS/reference/tidyGWAS.html).

## Details

Clean GWAS summary statistics

## Examples

``` r
if (FALSE) { # \dontrun{
cleaned <- clean_gwas(
  sumstats_file = "my_gwas.txt.gz",
  logging_path  = "logs/"
)
} # }
```
