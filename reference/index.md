# Package index

## Pipeline generation

Generate a complete `targets` pipeline and Quarto report chunks from a
cohort manifest.

- [`generate_gwas_meta_pipeline()`](http://www.levin-lab.org/gwastargets/reference/generate_gwas_meta_pipeline.md)
  : Generate targets pipeline code for GWAS meta-analysis
- [`generate_report_chunks()`](http://www.levin-lab.org/gwastargets/reference/generate_report_chunks.md)
  : Generate Quarto report chunks for a GWAS meta-analysis trait

## Per-cohort preparation

Clean, harmonize, and QC individual cohort summary statistics.

- [`prep_gwas()`](http://www.levin-lab.org/gwastargets/reference/prep_gwas.md)
  : Prepare and QC GWAS summary statistics
- [`clean_gwas()`](http://www.levin-lab.org/gwastargets/reference/clean_gwas.md)
  : Clean GWAS summary statistics with tidyGWAS
- [`harmonize_sumstats_headers()`](http://www.levin-lab.org/gwastargets/reference/harmonize_sumstats_headers.md)
  : Harmonize GWAS summary statistics headers
- [`snp_match_munge()`](http://www.levin-lab.org/gwastargets/reference/snp_match_munge.md)
  : Match GWAS SNPs against a reference SNP list
- [`build_column_map()`](http://www.levin-lab.org/gwastargets/reference/build_column_map.md)
  : Build a column rename map from per-cohort manifest values

## Quality control

Summarize per-cohort statistics before meta-analysis.

- [`summarize_sumstats()`](http://www.levin-lab.org/gwastargets/reference/summarize_sumstats.md)
  : Summarize GWAS summary statistics

## Meta-analysis

Fixed-effects IVW meta-analysis across cohorts.

- [`meta_analyze_ivw()`](http://www.levin-lab.org/gwastargets/reference/meta_analyze_ivw.md)
  : Inverse-variance weighted GWAS meta-analysis

## Downstream analysis

Loci extraction, common variant identification, and MR-MEGA MDS.

- [`extract_loci()`](http://www.levin-lab.org/gwastargets/reference/extract_loci.md)
  : Extract genome-wide significant independent loci
- [`extract_common_variants()`](http://www.levin-lab.org/gwastargets/reference/extract_common_variants.md)
  : Extract variants common to all GWAS studies
- [`calculate_mr_mega_mds()`](http://www.levin-lab.org/gwastargets/reference/calculate_mr_mega_mds.md)
  : Calculate MR-MEGA MDS coordinates from allele frequency data
