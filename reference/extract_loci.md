# Extract genome-wide significant independent loci

Filters a meta-analysis data frame to genome-wide significant variants,
identifies independent loci via
[`gwasRtools::get_loci()`](https://lcpilling.github.io/gwasRtools/reference/get_loci.html),
and annotates each lead variant with the nearest gene via
[`gwasRtools::get_nearest_gene()`](https://lcpilling.github.io/gwasRtools/reference/get_nearest_gene.html).
Returns an empty
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
when no variants pass the p-value threshold.

## Usage

``` r
extract_loci(
  df,
  snp_col = RSID,
  chr_col = CHR,
  pos_col = POS_38,
  maf_col = EAF,
  beta_col = B,
  se_col = SE,
  p_col = p_value,
  p_threshold = 5e-08,
  build = 38,
  ...
)
```

## Arguments

- df:

  A data frame of meta-analysis results (e.g. output of
  [`meta_analyze_ivw()`](https://mglev1n.github.io/gwastargets/reference/meta_analyze_ivw.md)).

- snp_col:

  Bare column name for rsid. Default `RSID`.

- chr_col:

  Bare column name for chromosome. Default `CHR`.

- pos_col:

  Bare column name for position. Default `POS_38`.

- maf_col:

  Bare column name for minor allele frequency. Default `EAF`.

- beta_col:

  Bare column name for effect size. Default `B`.

- se_col:

  Bare column name for SE. Default `SE`.

- p_col:

  Bare column name for p-value. Default `p_value`.

- p_threshold:

  Genome-wide significance threshold. Default `5e-8`.

- build:

  Genome build for gene annotation (`37` or `38`). Default `38`.

- ...:

  Additional arguments forwarded to
  [`gwasRtools::get_loci()`](https://lcpilling.github.io/gwasRtools/reference/get_loci.html).

## Value

A
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
of independent lead variants with gene annotations, or an empty tibble
if no variants pass the threshold.

## Details

Extract genome-wide significant loci

## Examples

``` r
if (FALSE) { # \dontrun{
loci <- extract_loci(meta_results, p_threshold = 5e-8)
} # }
```
