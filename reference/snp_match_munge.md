# Match GWAS SNPs against a reference SNP list

Matches variants in `sumstats` against `info_snp` by rsid, effect allele
(`a1`), and other allele (`a0`). Optionally removes ambiguous SNPs (A/T
and C/G), performs strand flipping, and reverses allele coding.

## Usage

``` r
snp_match_munge(
  sumstats,
  info_snp,
  strand_flip = TRUE,
  join_by_pos = TRUE,
  remove_dups = TRUE,
  match.min.prop = 0.2,
  return_flip_and_rev = FALSE
)
```

## Arguments

- sumstats:

  A data frame of GWAS summary statistics. Must contain columns `rsid`,
  `a0`, `a1`, and `beta`.

- info_snp:

  A data frame of reference SNPs. Must contain columns `rsid`, `a0`, and
  `a1`.

- strand_flip:

  Logical; if `TRUE` (default), remove ambiguous SNPs (A/T, C/G) and
  attempt strand flipping.

- join_by_pos:

  Logical; if `TRUE`, also join by position. Default `TRUE`.

- remove_dups:

  Logical; if `TRUE` (default), remove duplicate rsid entries after
  matching.

- match.min.prop:

  Minimum proportion of variants that must match. Default `0.2`. An
  error is raised if fewer variants match.

- return_flip_and_rev:

  Logical; if `TRUE`, retain the `_FLIP_` and `_REV_` indicator columns
  in the output. Default `FALSE`.

## Value

A data frame of matched variants ordered by chromosome.

## Details

Match and munge SNPs against a reference panel

## Examples

``` r
if (FALSE) { # \dontrun{
matched <- snp_match_munge(
  sumstats = my_sumstats,
  info_snp = hm3_snps
)
} # }
```
