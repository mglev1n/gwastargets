# Getting Started with gwastargets

## Overview

`gwastargets` provides functions for building reproducible multi-cohort
GWAS meta-analysis pipelines using the
[`targets`](https://docs.ropensci.org/targets/) framework. It supports
both binary (case/control) and quantitative traits.

The typical workflow is:

1.  Assemble a **manifest** data frame describing your cohort files
2.  Call
    [`generate_gwas_meta_pipeline()`](https://mglev1n.github.io/gwastargets/reference/generate_gwas_meta_pipeline.md)
    to produce `_targets.R` code
3.  Paste the output into your `_targets.R` and run
    [`targets::tar_make()`](https://docs.ropensci.org/targets/reference/tar_make.html)

The pipeline functions
([`prep_gwas()`](https://mglev1n.github.io/gwastargets/reference/prep_gwas.md),
[`meta_analyze_ivw()`](https://mglev1n.github.io/gwastargets/reference/meta_analyze_ivw.md),
etc.) are also available for interactive use outside of `targets`.

------------------------------------------------------------------------

## Required ancillary files

Two reference files must be provided before generating a pipeline:

### HapMap 3 SNP list (`hm3_path`)

Used by
[`prep_gwas()`](https://mglev1n.github.io/gwastargets/reference/prep_gwas.md)
to match variants and by LDSC for heritability and genetic correlation
estimation. The file is loaded via
[`vroom::vroom()`](https://vroom.tidyverse.org/reference/vroom.html) and
should have columns `SNP`, `A1`, `A2`. The standard LDSC HapMap 3 SNP
list (`w_hm3.snplist`) is available from the [Alkes Group LDSCORE
repository](https://alkesgroup.broadinstitute.org/LDSCORE/).

### dbSNP155 reference (`dbsnp_path`)

Used by
[`clean_gwas()`](https://mglev1n.github.io/gwastargets/reference/clean_gwas.md)
(via `tidyGWAS`) for variant QC and rsID harmonization. This should be
the path to a local dbSNP155 directory. See the [`tidyGWAS`
documentation](https://github.com/Ararder/tidyGWAS) for instructions on
downloading and preparing this reference.

------------------------------------------------------------------------

## Step 1: Build a cohort manifest

The manifest is a data frame with one row per cohort × ancestry
combination. It must contain five columns:

| Column     | Description                                                     |
|------------|-----------------------------------------------------------------|
| `path`     | Full path to the raw summary statistics file                    |
| `file`     | Basename of the file                                            |
| `cohort`   | Human-readable cohort label                                     |
| `ancestry` | Ancestry label used by `ldscr` (e.g. `"EUR"`, `"AFR"`, `"EAS"`) |
| `study`    | Unique identifier used as the `targets` target name             |

``` r
library(gwastargets)

manifest <- data.frame(
  path     = c(
    "/project/data/UKBB_EUR.txt.gz",
    "/project/data/MVP_EUR.txt.gz",
    "/project/data/BBJ_EAS.txt.gz"
  ),
  file     = c("UKBB_EUR.txt.gz", "MVP_EUR.txt.gz", "BBJ_EAS.txt.gz"),
  cohort   = c("UKBB",            "MVP",            "BBJ"),
  ancestry = c("EUR",             "EUR",            "EAS"),
  study    = c("UKBB_EUR",        "MVP_EUR",        "BBJ_EAS"),
  stringsAsFactors = FALSE
)
```

[`generate_gwas_meta_pipeline()`](https://mglev1n.github.io/gwastargets/reference/generate_gwas_meta_pipeline.md)
validates the manifest before generating any code — it will error if
required columns are missing or `study` values are duplicated, and warn
if any `path` files do not exist on disk.

------------------------------------------------------------------------

## Step 2: Generate the pipeline

``` r
code <- generate_gwas_meta_pipeline(
  trait       = "CAD",
  trait_type  = "binary",      # or "quantitative"
  n_col       = "EffectiveN",  # use "N" for quantitative traits
  manifest_df = manifest,
  hm3_path    = "/path/to/w_hm3.snplist",
  dbsnp_path  = "/path/to/dbSNP155"
)

# Review the generated code
cat(code)
```

The output is a self-contained block of R code ready to paste into
`_targets.R`. It includes:

- The manifest serialized as a
  [`tibble::tribble()`](https://tibble.tidyverse.org/reference/tribble.html)
  call
- A `tar_map()` block running per-cohort preparation and LDSC
  heritability
- LDSC genetic correlation QC targets (one per ancestry with ≥ 2
  cohorts)
- Per-ancestry IVW meta-analysis targets with Manhattan plots and loci
- All-populations IVW meta-analysis with MR-MEGA MDS
- Cohort overview table and precision plot

------------------------------------------------------------------------

## Step 3: Per-cohort preparation (`prep_gwas`)

Each cohort’s summary statistics pass through
[`prep_gwas()`](https://mglev1n.github.io/gwastargets/reference/prep_gwas.md),
which:

1.  **Cleans** the raw file via
    [`clean_gwas()`](https://mglev1n.github.io/gwastargets/reference/clean_gwas.md)
    (header harmonization + `tidyGWAS` QC)
2.  **Validates** that the expected columns are present for the declared
    `trait_type`
3.  **Matches** variants to the HapMap 3 SNP list via
    [`snp_match_munge()`](https://mglev1n.github.io/gwastargets/reference/snp_match_munge.md)
4.  **Estimates** the LDSC intercept and applies SE correction when
    intercept \> 1
5.  **Writes** a `.parquet` file with standardized columns including
    `SE_ldsc_adjusted`

``` r
# Typically called inside a targets pipeline, but works interactively too
out_path <- prep_gwas(
  sumstats_file = "/project/data/UKBB_EUR.txt.gz",
  hm3           = hm3_snps,              # data frame with SNP, A1, A2
  ancestry      = "EUR",
  output_path   = "Data/CAD/prepped_sumstats",
  trait_type    = "binary",
  n_col         = EffectiveN,            # bare name, not quoted
  dbsnp_path    = "/path/to/dbSNP155"
)
```

The output parquet file contains the original cleaned columns plus:

| Added column       | Description                                                |
|--------------------|------------------------------------------------------------|
| `intercept`        | LDSC intercept estimate                                    |
| `SE_ldsc_adjusted` | `SE * sqrt(intercept)` when intercept \> 1, otherwise `SE` |
| `Z_ldsc_adjusted`  | `B / SE_ldsc_adjusted`                                     |
| `P_ldsc_adjusted`  | Two-sided p-value from `Z_ldsc_adjusted`                   |
| `ldsc_adjustment`  | `TRUE` if intercept correction was applied                 |

------------------------------------------------------------------------

## Step 4: Cohort QC with `summarize_sumstats`

Before meta-analysis, inspect per-cohort statistics:

``` r
prepped_files <- list.files("Data/CAD/prepped_sumstats", full.names = TRUE,
                            pattern = "\\.parquet$")

summary_tbl <- summarize_sumstats(prepped_files, trait_type = "binary")
```

For binary traits, the summary includes `n_cases`, `n_controls`,
`n_effective`, and `n_total`. For quantitative traits, only `n_total` is
returned. Both include SNP counts, EAF range, and median SE — used to
build the cohort overview table and precision plot in the generated
pipeline.

------------------------------------------------------------------------

## Step 5: IVW meta-analysis

[`meta_analyze_ivw()`](https://mglev1n.github.io/gwastargets/reference/meta_analyze_ivw.md)
runs a fixed-effects inverse-variance weighted meta-analysis across all
prepped parquet files, one chromosome at a time:

``` r
meta_results <- meta_analyze_ivw(
  parquet_files = prepped_files,
  trait_type    = "binary",
  se_col        = SE_ldsc_adjusted,   # use LDSC-adjusted SE
  chromosomes   = 1:22,
  min_mac       = 100                  # minimum minor allele count filter
)
```

Output columns include: `B`, `SE`, `p_value`, `z_score`, `EAF`,
`n_contributions`, `N` (and for binary traits: `CaseN`, `ControlN`,
`EffectiveN`), plus heterogeneity statistics `Q`, `Q_df`, `Q_pval`, and
`I2`.

------------------------------------------------------------------------

## Step 6: Loci extraction

``` r
loci <- extract_loci(meta_results)
```

Calls
[`gwasRtools::get_loci()`](https://lcpilling.github.io/gwasRtools/reference/get_loci.html)
for LD-based clumping and
[`gwasRtools::get_nearest_gene()`](https://lcpilling.github.io/gwasRtools/reference/get_nearest_gene.html)
for gene annotation. Returns a tibble with one row per independent locus
and columns for lead variant, gene, and association statistics.

------------------------------------------------------------------------

## Step 7: MR-MEGA MDS

For multi-ancestry analyses, MR-MEGA requires a matrix of per-study
effect estimates at common variants.
[`extract_common_variants()`](https://mglev1n.github.io/gwastargets/reference/extract_common_variants.md)
identifies variants present in all studies above a MAF threshold, and
[`calculate_mr_mega_mds()`](https://mglev1n.github.io/gwastargets/reference/calculate_mr_mega_mds.md)
computes MDS coordinates from the resulting effect heterogeneity matrix:

``` r
common_vars <- extract_common_variants(
  parquet_files        = prepped_files,
  trait_type           = "binary",
  select_per_mb_window = TRUE   # LD-thin to ~1 variant/Mb
)

mds <- calculate_mr_mega_mds(common_vars)
# mds$pc_coordinates: study × PC matrix
# mds$distance_matrix: pairwise effect-size distance matrix
```

The MDS coordinates are passed as covariates to MR-MEGA to model
ancestry-driven effect heterogeneity.

------------------------------------------------------------------------

## Raw summary statistics format

[`prep_gwas()`](https://mglev1n.github.io/gwastargets/reference/prep_gwas.md)
accepts a wide range of column naming conventions.
[`harmonize_sumstats_headers()`](https://mglev1n.github.io/gwastargets/reference/harmonize_sumstats_headers.md)
(called internally by
[`clean_gwas()`](https://mglev1n.github.io/gwastargets/reference/clean_gwas.md))
maps common aliases to the standard names used throughout this package:

| Standard name | Common aliases accepted                |
|---------------|----------------------------------------|
| `B`           | `BETA`, `Effect`, `b`                  |
| `SE`          | `StdErr`, `standard_error`             |
| `P`           | `PVAL`, `P_VALUE`, `pvalue`, `p.value` |
| `EAF`         | `EFF_ALL_FREQ`, `A1FREQ`, `AF_Allele2` |
| `N`           | `Neff`, `NMISS`, `OBS_CT`              |
| `CaseN`       | `N_CASE`, `Ncases`                     |
| `ControlN`    | `N_CONTROL`, `Ncontrols`               |

Additional mappings can be passed via the `column_names` argument to
[`prep_gwas()`](https://mglev1n.github.io/gwastargets/reference/prep_gwas.md).
