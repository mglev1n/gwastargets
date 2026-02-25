# Calculate MR-MEGA MDS coordinates from allele frequency data

Computes multidimensional scaling (MDS) coordinates from study-level
effect allele frequency (EAF) data. The approach follows the MR-MEGA
method (Mägi et al. 2017): a population-genetic distance matrix is
computed from EAF values, double-centred, and decomposed via SVD to
produce principal coordinates.

## Usage

``` r
calculate_mr_mega_mds(
  common_variants_data,
  pcCount = 4,
  chr_col = CHR,
  pos_col = POS_38,
  rsid_col = RSID,
  ea_col = EffectAllele,
  oa_col = OtherAllele,
  eaf_col = minor_af,
  file_col = file
)
```

## Arguments

- common_variants_data:

  A data frame (long format) with one row per variant × study
  combination. Must contain columns identified by `chr_col`, `pos_col`,
  `ea_col`, `oa_col`, `eaf_col`, and `file_col`. Typically the output of
  [`extract_common_variants()`](https://mglev1n.github.io/gwastargets/reference/extract_common_variants.md).

- pcCount:

  Number of principal coordinates to return. Default `4`.

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

- eaf_col:

  Bare column name for EAF. Default `minor_af`.

- file_col:

  Bare column name for study/file identifier. Default `file`.

## Value

A named list with elements:

- `pc_coordinates`:

  Data frame of MDS coordinates (`pcCount` columns).

- `distance_matrix`:

  Population-genetic distance matrix.

- `centered_distance_matrix`:

  Double-centred distance matrix.

- `eigenvalues`:

  Eigenvalues from SVD.

- `eaf_matrix`:

  Wide EAF matrix (variants × studies).

- `used_variants`:

  Variants with complete EAF data.

- `studies`:

  Character vector of study names.

- `summary`:

  List with `n_studies`, `n_variants_used`, `n_variants_total`,
  `pc_count`.

## Details

Calculate MR-MEGA multidimensional scaling coordinates

## References

Mägi R, et al. (2017). Trans-ethnic meta-regression of genome-wide
association studies accounting for ancestry increases power for
discovery and improves fine-mapping resolution. *Hum Mol Genet*, 26(18),
3639–3650.

## Examples

``` r
if (FALSE) { # \dontrun{
common <- extract_common_variants(files, trait_type = "binary")
mds    <- calculate_mr_mega_mds(common)
} # }
```
