# gwastargets

<!-- badges: start -->
[![R-CMD-check](https://github.com/mglev1n/gwastargets/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mglev1n/gwastargets/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`gwastargets` provides a suite of functions for building reproducible multi-cohort GWAS meta-analysis pipelines using the [`targets`](https://docs.ropensci.org/targets/) framework. It supports both binary (case/control) and quantitative traits and handles the full workflow from raw summary statistics through meta-analysis and downstream quality-control.

## Installation

```r
# install.packages("pak")
pak::pkg_install("mglev1n/gwastargets")
```

## Overview

The package is organized around two complementary use cases:

**1. Generating a `targets` pipeline** — given a cohort manifest and trait metadata, produce the complete `_targets.R` code for a reproducible analysis:

```r
library(gwastargets)

manifest <- data.frame(
  path     = c("/data/UKBB_EUR.txt.gz", "/data/MVP_AFR.txt.gz", "/data/BBJ_EAS.txt.gz"),
  file     = c("UKBB_EUR.txt.gz",       "MVP_AFR.txt.gz",       "BBJ_EAS.txt.gz"),
  cohort   = c("UKBB",                  "MVP",                  "BBJ"),
  ancestry = c("EUR",                   "AFR",                  "EAS"),
  study    = c("UKBB_EUR",              "MVP_AFR",              "BBJ_EAS"),
  stringsAsFactors = FALSE
)

code <- generate_gwas_meta_pipeline(
  trait       = "CAD",
  trait_type  = "binary",
  n_col       = "EffectiveN",
  manifest_df = manifest
)

cat(code)
```

**2. Using the pipeline functions directly** — each stage is available as a standalone exported function for interactive use or integration into existing workflows:

| Function | Purpose |
|---|---|
| `prep_gwas()` | Clean, harmonize, and LDSC-correct one cohort's summary statistics |
| `summarize_sumstats()` | Per-cohort QC summary (N, MAF range, median SE, precision) |
| `meta_analyze_ivw()` | Fixed-effects IVW meta-analysis across cohorts |
| `extract_loci()` | Clump significant variants into independent loci with gene annotation |
| `extract_common_variants()` | Identify variants present across all studies (for MR-MEGA) |
| `calculate_mr_mega_mds()` | Compute MR-MEGA MDS coordinates from cross-cohort effect estimates |
| `snp_match_munge()` | Match summary statistics to a reference panel and munge for LDSC |
| `harmonize_sumstats_headers()` | Standardize column names across GWAS file formats |

## Pipeline Architecture

```
Raw summary stats (per cohort)
        │
        ▼
   prep_gwas()          ← clean, match HM3, LDSC intercept correction
        │
        ├──── ldsc_h2()            per-cohort heritability (via ldscr)
        │
        ├──── summarize_sumstats() per-cohort QC table
        │
        ├──── meta_analyze_ivw()   per-ancestry IVW meta-analysis
        │         └── extract_loci()
        │
        └──── meta_analyze_ivw()   all-populations IVW meta-analysis
                  ├── extract_loci()
                  └── extract_common_variants()
                            └── calculate_mr_mega_mds()
```

## Manifest requirements

`generate_gwas_meta_pipeline()` validates the manifest before generating any code. The required columns are:

| Column | Description |
|---|---|
| `path` | Full path to the raw summary statistics file |
| `file` | Basename of the file (used to join back to per-cohort summaries) |
| `cohort` | Human-readable cohort label (e.g. `"UKBB"`) |
| `ancestry` | Ancestry label matching ldscr reference panels (e.g. `"EUR"`, `"AFR"`, `"EAS"`) |
| `study` | Unique identifier used as the `tar_map` target name — no spaces or special characters |

Additional columns are preserved in the serialized output and can be used downstream.

## Related packages

- [`ldscr`](https://github.com/mglev1n/ldscr) — LDSC heritability and genetic correlation in R
- [`tidyGWAS`](https://github.com/Ararder/tidyGWAS) — GWAS summary statistics cleaning
- [`gwasRtools`](https://github.com/lcpilling/gwasRtools) — GWAS loci extraction and annotation
- [`targets`](https://docs.ropensci.org/targets/) — reproducible pipeline toolkit
