# Package index

## Pipeline generation

Generate a complete `targets` pipeline from a cohort manifest.

- [`generate_gwas_meta_pipeline()`](https://mglev1n.github.io/gwastargets/reference/generate_gwas_meta_pipeline.md)
  : Generate targets pipeline code for GWAS meta-analysis

## Per-cohort preparation

Clean, harmonize, and QC individual cohort summary statistics.

- [`prep_gwas()`](https://mglev1n.github.io/gwastargets/reference/prep_gwas.md)
  : Prepare and QC GWAS summary statistics
- [`clean_gwas()`](https://mglev1n.github.io/gwastargets/reference/clean_gwas.md)
  : Clean GWAS summary statistics with tidyGWAS
- [`harmonize_sumstats_headers()`](https://mglev1n.github.io/gwastargets/reference/harmonize_sumstats_headers.md)
  : Harmonize GWAS summary statistics headers
- [`snp_match_munge()`](https://mglev1n.github.io/gwastargets/reference/snp_match_munge.md)
  : Match GWAS SNPs against a reference SNP list

## Quality control

Summarize per-cohort statistics before meta-analysis.

- [`summarize_sumstats()`](https://mglev1n.github.io/gwastargets/reference/summarize_sumstats.md)
  : Summarize GWAS summary statistics

## Meta-analysis

Fixed-effects IVW meta-analysis across cohorts.

- [`meta_analyze_ivw()`](https://mglev1n.github.io/gwastargets/reference/meta_analyze_ivw.md)
  : Inverse-variance weighted GWAS meta-analysis

## Downstream analysis

Loci extraction, common variant identification, and MR-MEGA MDS.

- [`extract_loci()`](https://mglev1n.github.io/gwastargets/reference/extract_loci.md)
  : Extract genome-wide significant independent loci
- [`extract_common_variants()`](https://mglev1n.github.io/gwastargets/reference/extract_common_variants.md)
  : Extract variants common to all GWAS studies
- [`calculate_mr_mega_mds()`](https://mglev1n.github.io/gwastargets/reference/calculate_mr_mega_mds.md)
  : Calculate MR-MEGA MDS coordinates from allele frequency data
