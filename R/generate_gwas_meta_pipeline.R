# Valid ancestry labels supported by ldscr (UKB LD scores)
.VALID_ANCESTRIES <- c("AFR", "AMR", "CSA", "EAS", "EUR", "MID")

#' Generate a targets GWAS meta-analysis pipeline
#'
#' @title Generate targets pipeline code for GWAS meta-analysis
#'
#' @description
#' Generates a complete `targets` / `tarchetypes` pipeline as a character
#' string, ready to be pasted into a `_targets.R` file. The generated code
#' covers per-cohort preparation, LDSC heritability estimation, per-ancestry
#' and all-populations IVW meta-analysis, loci extraction, Manhattan plots,
#' and MR-MEGA MDS. Validates `trait_type`, `n_col`, and `manifest_df` before
#' generating any code. The manifest is serialized as a `tibble::tribble()`
#' call in the generated output so the pipeline is self-contained.
#'
#' @param trait Trait name used in target names (e.g. `"CAD"`).
#' @param trait_type `"binary"` or `"quantitative"` (required, no default).
#' @param n_col Sample size column: `"EffectiveN"` for binary traits, `"N"`
#'   for quantitative traits (required, no default).
#' @param manifest_df A data frame describing the cohort files to include.
#'   Required columns: `path` (full path to raw summary statistics file),
#'   `file` (basename of the file), `cohort` (cohort label), `ancestry`
#'   (ancestry label, e.g. `"EUR"`), `study` (unique human-readable identifier;
#'   must contain no spaces or special characters; does not need to include the
#'   ancestry). A `tar_name` column (`{study}_{ancestry}`) is automatically
#'   added and used as the `tar_map` target name so that ancestry-based target
#'   selectors work reliably. Additional columns are preserved in the serialized
#'   output.
#' @param hm3_path Path to the LDSC HapMap3 SNP list file (e.g.
#'   `w_hm3.snplist`). Required.
#' @param dbsnp_path Path to the dbSNP155 reference directory passed to
#'   [prep_gwas()]. Required.
#' @param crew_controller Name of the `crew` controller to use for
#'   resource-intensive targets (`meta_ALL` and `meta_common_variants_ALL`).
#'   If `NULL` (default), no `resources` block is added to those targets.
#' @param output_base_dir Base output directory for intermediate files.
#'   Default `"Data"`.
#'
#' @return A length-1 character string containing valid R code defining a
#'   `targets` pipeline.
#'
#' @examples
#' \dontrun{
#' manifest <- data.frame(
#'   path     = c("/data/UKBB_EUR.txt.gz", "/data/MVP_AFR.txt.gz"),
#'   file     = c("UKBB_EUR.txt.gz",       "MVP_AFR.txt.gz"),
#'   cohort   = c("UKBB",                  "MVP"),
#'   ancestry = c("EUR",                   "AFR"),
#'   study    = c("UKBB_EUR",              "MVP_AFR"),
#'   stringsAsFactors = FALSE
#' )
#' cat(generate_gwas_meta_pipeline("CAD", trait_type = "binary",
#'                                 n_col = "EffectiveN", manifest_df = manifest,
#'                                 hm3_path = "/path/to/w_hm3.snplist"))
#' }
#'
#' @importFrom cli cli_abort cli_warn
#' @importFrom glue glue
#' @export
generate_gwas_meta_pipeline <- function(trait,
                                        trait_type,
                                        n_col,
                                        manifest_df,
                                        hm3_path,
                                        dbsnp_path,
                                        crew_controller = NULL,
                                        output_base_dir = "Data") {

  # Validate trait_type
  if (missing(trait_type) || !trait_type %in% c("binary", "quantitative")) {
    cli::cli_abort(c(
      "{.arg trait_type} must be {.val binary} or {.val quantitative}",
      "i" = "No default is provided -- this must be declared explicitly"
    ))
  }

  # Validate n_col
  if (missing(n_col) || !n_col %in% c("EffectiveN", "N")) {
    cli::cli_abort(c(
      "{.arg n_col} must be {.val EffectiveN} or {.val N}",
      "i" = "Use {.val EffectiveN} for binary traits, {.val N} for quantitative traits"
    ))
  }

  # Cross-validate
  if (trait_type == "binary" && n_col == "N") {
    cli::cli_abort(c(
      "{.arg trait_type} = {.val binary} requires {.arg n_col} = {.val EffectiveN}",
      "i" = "For quantitative traits, set {.arg trait_type} = {.val quantitative}"
    ))
  }
  if (trait_type == "quantitative" && n_col == "EffectiveN") {
    cli::cli_warn(c(
      "{.arg trait_type} = {.val quantitative} with {.arg n_col} = {.val EffectiveN}",
      "i" = "This is unusual -- did you mean {.arg trait_type} = {.val binary}?"
    ))
  }

  # Validate manifest_df
  if (missing(manifest_df)) {
    cli::cli_abort(c(
      "{.arg manifest_df} must be provided",
      "i" = "Supply a data frame with columns: {.val path}, {.val file}, {.val cohort}, {.val ancestry}, {.val study}"
    ))
  }

  required_cols <- c("path", "file", "cohort", "ancestry", "study")
  missing_cols  <- setdiff(required_cols, names(manifest_df))
  if (length(missing_cols) > 0) {
    cli::cli_abort(c(
      "{.arg manifest_df} is missing required column{?s}: {.val {missing_cols}}",
      "i" = "Required columns: {.val {required_cols}}"
    ))
  }

  dup_studies <- manifest_df$study[duplicated(manifest_df$study)]
  if (length(dup_studies) > 0) {
    cli::cli_abort(c(
      "{.arg manifest_df} has duplicate {.field study} value{?s}: {.val {dup_studies}}",
      "i" = "Each study must be unique as it is used as a {.pkg targets} target name"
    ))
  }

  bad_ancestries <- setdiff(manifest_df$ancestry, .VALID_ANCESTRIES)
  if (length(bad_ancestries) > 0) {
    cli::cli_abort(c(
      "{.arg manifest_df} contains unsupported {.field ancestry} value{?s}: {.val {bad_ancestries}}",
      "i" = "Valid ancestries (ldscr UKB LD scores): {.val {(.VALID_ANCESTRIES)}}"
    ))
  }

  # Construct tar_name: always {study}_{ancestry} so ends_with(anc) selectors
  # work reliably regardless of how the user names the study column.
  manifest_df$tar_name <- paste0(manifest_df$study, "_", manifest_df$ancestry)

  missing_files <- manifest_df$path[!file.exists(manifest_df$path)]
  if (length(missing_files) > 0) {
    cli::cli_warn(c(
      "{.arg manifest_df} contains {length(missing_files)} path{?s} that do not exist",
      "i" = "Missing: {.path {missing_files}}"
    ))
  }

  # Validate hm3_path
  if (missing(hm3_path)) {
    cli::cli_abort(c(
      "{.arg hm3_path} must be provided",
      "i" = "Supply the path to the LDSC HapMap3 SNP list (e.g. {.path w_hm3.snplist})"
    ))
  }
  if (!file.exists(hm3_path)) {
    cli::cli_warn(c(
      "{.arg hm3_path} does not exist: {.path {hm3_path}}",
      "i" = "The pipeline will fail at runtime if this file is missing"
    ))
  }

  if (missing(dbsnp_path)) {
    cli::cli_abort(c(
      "{.arg dbsnp_path} must be provided",
      "i" = "Supply the path to the dbSNP155 reference directory"
    ))
  }
  if (!file.exists(dbsnp_path)) {
    cli::cli_warn(c(
      "{.arg dbsnp_path} does not exist: {.path {dbsnp_path}}",
      "i" = "The pipeline will fail at runtime if this path is missing"
    ))
  }

  # Pre-compute crew resources blocks (differ by whether they are mid-argument or final)
  if (!is.null(crew_controller)) {
    resources_block_mid <- paste0(
      '    resources = tar_resources(\n',
      '      crew = tar_resources_crew(controller = "', crew_controller, '")\n',
      '    ),\n'
    )
    resources_block_end <- paste0(
      '    resources = tar_resources(\n',
      '      crew = tar_resources_crew(controller = "', crew_controller, '")\n',
      '    )\n'
    )
  } else {
    resources_block_mid <- ""
    resources_block_end <- ""
  }

  trait_lower <- tolower(trait)
  trait_upper <- toupper(trait)

  manifest_code <- .format_manifest_tribble(manifest_df, trait_lower)

  # Cohort overview table columns differ by trait_type
  if (trait_type == "binary") {
    cohort_overview_select <- "select(cohort, ancestry, n_cases, n_controls, n_total, n_effective)"
    cohort_overview_fmt    <- "fmt_number(columns = c(n_cases, n_controls, n_total, n_effective), decimals = 0)"
  } else {
    cohort_overview_select <- "select(cohort, ancestry, n_total)"
    cohort_overview_fmt    <- "fmt_number(columns = c(n_total), decimals = 0)"
  }

  pipeline_code <- glue::glue(
    '\n## <<trait_upper>> META-ANALYSIS PIPELINE\n## trait_type: <<trait_type>>\n## n_col:      <<n_col>>\n\n<<manifest_code>>\n\n<<trait_lower>>_mapped <- tarchetypes::tar_map(\n  values = <<trait_lower>>_manifest_df,\n  names  = tar_name,\n  unlist = FALSE,\n\n  tar_file(\n    <<trait_lower>>_sumstats_raw,\n    path\n  ),\n\n  tar_file(\n    <<trait_lower>>_sumstats_prepped,\n    prep_gwas(\n      sumstats_file = <<trait_lower>>_sumstats_raw,\n      hm3           = hm3,\n      ancestry      = ancestry,\n      trait_type    = "<<trait_type>>",\n      n_col         = <<n_col>>,\n      output_path   = "<<output_base_dir>>/<<trait_upper>>/prepped_sumstats",\n      dbsnp_path    = "<<dbsnp_path>>",\n      column_names  = list(EAF = "EFF_ALL_FREQ")\n    )\n  ),\n\n  tar_target(\n    <<trait_lower>>_sumstats_munged,\n    snp_match_munge(\n      arrow::read_parquet(<<trait_lower>>_sumstats_prepped) |>\n        select(chr = CHR, rsid = RSID, a0 = OtherAllele, a1 = EffectAllele,\n               beta = Z, N = <<n_col>>, EAF),\n      info_snp       = hm3 |> select(rsid = SNP, a1 = A1, a0 = A2),\n      join_by_pos    = FALSE,\n      match.min.prop = 0\n    ) |>\n      select(SNP = rsid, N, Z = beta, A1 = a1, A2 = a0, EAF),\n    description = "<<trait_upper>> munged summary statistics",\n    format      = "parquet"\n  ),\n\n  tar_target(\n    <<trait_lower>>_ldsc_h2,\n    ldscr::ldsc_h2(\n      munged_sumstats = <<trait_lower>>_sumstats_munged,\n      ancestry        = ancestry\n    ),\n    description = "<<trait_upper>> LDSC heritability estimates"\n  )\n)\n\n# LDSC genetic correlations within each ancestry (QC)\n<<trait_lower>>_ldsc_rg_qc_targets <- purrr::map(\n  <<trait_lower>>_manifest_df |>\n    dplyr::add_count(ancestry) |> dplyr::filter(n > 1) |> dplyr::distinct(ancestry) |> dplyr::pull(ancestry),\n  function(anc) {\n    tar_map(\n      values = data.frame(ancestry = anc),\n      names  = "ancestry",\n      unlist = TRUE,\n      tarchetypes::tar_combine(\n        <<trait_lower>>_ldsc_rg_qc,\n        tarchetypes::tar_select_targets(\n          <<trait_lower>>_mapped,\n          starts_with("<<trait_lower>>_sumstats_munged") & ends_with(!!anc)\n        ),\n        command     = ldscr::ldsc_rg(munged_sumstats = list(!!!.x), ancestry = ancestry),\n        description = glue::glue("LDSC genetic correlations for {anc} ancestry")\n      ),\n      tar_target(\n        <<trait_lower>>_ldsc_rg_qc_heatmap,\n        <<trait_lower>>_ldsc_rg_qc |> plot_ldsc_rg_qc(),\n        format      = "qs",\n        description = glue::glue("LDSC rg heatmap for {anc} ancestry")\n      )\n    )\n  }\n) |> purrr::list_flatten()\n\n# Per-ancestry meta-analysis\n<<trait_lower>>_prep_meta_ancestry_targets <- purrr::map(\n  <<trait_lower>>_manifest_df |>\n    dplyr::add_count(ancestry) |> dplyr::filter(n > 1) |> dplyr::distinct(ancestry) |> dplyr::pull(ancestry),\n  function(anc) {\n    tar_map(\n      unlist = TRUE,\n      values = data.frame(ancestry = anc),\n      names  = "ancestry",\n      tarchetypes::tar_combine(\n        <<trait_lower>>_sumstats_prepped_combined,\n        tarchetypes::tar_select_targets(\n          <<trait_lower>>_mapped,\n          starts_with("<<trait_lower>>_sumstats_prepped") & ends_with(!!anc)\n        ),\n        description = glue::glue("Prepped <<trait_upper>> sumstats for {anc}"),\n        format      = "file"\n      ),\n      tar_target(\n        <<trait_lower>>_sumstats_summary,\n        summarize_sumstats(<<trait_lower>>_sumstats_prepped_combined,\n                           trait_type = "<<trait_type>>"),\n        description = glue::glue("<<trait_upper>> sumstats summary for {anc}")\n      ),\n      tar_target(\n        <<trait_lower>>_meta,\n        meta_analyze_ivw(<<trait_lower>>_sumstats_prepped_combined,\n                         trait_type  = "<<trait_type>>",\n                         se_col      = SE_ldsc_adjusted,\n                         chromosomes = 1:22,\n                         min_mac     = 100),\n        description = glue::glue("<<trait_upper>> meta-analysis for {anc}")\n      ),\n      tar_target(\n        <<trait_lower>>_meta_loci,\n        extract_loci(<<trait_lower>>_meta) |>\n          mutate(population = ancestry),\n        description = glue::glue("<<trait_upper>> meta-analysis loci for {anc}")\n      ),\n      tar_target(\n        <<trait_lower>>_meta_manhattan,\n        <<trait_lower>>_meta |>\n          dplyr::filter(p_value < 0.001) |>\n          dplyr::filter(dplyr::between(EAF, 0.05, 1 - 0.05) | p_value < 5e-8) |>\n          levinmisc::gg_manhattan_df(\n            annotation_df = if (!is.null(<<trait_lower>>_meta_loci) & nrow(<<trait_lower>>_meta_loci) > 0) {\n              <<trait_lower>>_meta_loci |> dplyr::slice_min(p_value, n = 2, with_ties = FALSE, by = CHR)\n            } else { NULL },\n            chr_col  = CHR, pos_col = POS_38, beta_col = B,\n            se_col   = SE,  pval_col = p_value, build = "hg38"\n          ),\n        description = glue::glue("<<trait_upper>> Manhattan plot for {anc}"),\n        format      = "qs"\n      )\n    )\n  }\n) |> purrr::list_flatten()\n\n# All-populations meta-analysis\n<<trait_lower>>_prep_meta_all_targets <- list(\n  tarchetypes::tar_combine(\n    <<trait_lower>>_sumstats_prepped_combined_ALL,\n    tarchetypes::tar_select_targets(\n      <<trait_lower>>_mapped,\n      starts_with("<<trait_lower>>_sumstats_prepped")\n    ),\n    description = "Prepped <<trait_upper>> sumstats for ALL populations",\n    format      = "file"\n  ),\n  tar_target(\n    <<trait_lower>>_sumstats_summary_ALL,\n    summarize_sumstats(<<trait_lower>>_sumstats_prepped_combined_ALL,\n                       trait_type = "<<trait_type>>"),\n    description = "<<trait_upper>> sumstats summary for ALL populations"\n  ),\n  tar_target(\n    <<trait_lower>>_meta_ALL,\n    meta_analyze_ivw(<<trait_lower>>_sumstats_prepped_combined_ALL,\n                     trait_type  = "<<trait_type>>",\n                     se_col      = SE_ldsc_adjusted,\n                     chromosomes = 1:22,\n                     min_mac     = 100),\n    <<resources_block_mid>>description = "<<trait_upper>> meta-analysis for ALL populations"\n  ),\n  tar_target(\n    <<trait_lower>>_meta_loci_ALL,\n    extract_loci(<<trait_lower>>_meta_ALL),\n    description = "<<trait_upper>> meta-analysis loci for ALL populations"\n  ),\n  tar_target(\n    <<trait_lower>>_meta_manhattan_ALL,\n    <<trait_lower>>_meta_ALL |>\n      dplyr::filter(p_value < 0.001) |>\n      dplyr::filter(dplyr::between(EAF, 0.05, 1 - 0.05) | p_value < 5e-8) |>\n      levinmisc::gg_manhattan_df(\n        annotation_df = if (!is.null(<<trait_lower>>_meta_loci_ALL) & nrow(<<trait_lower>>_meta_loci_ALL) > 0) {\n          <<trait_lower>>_meta_loci_ALL |> dplyr::slice_min(p_value, n = 2, with_ties = FALSE, by = CHR)\n        } else { NULL },\n        chr_col  = CHR, pos_col = POS_38, beta_col = B,\n        se_col   = SE,  pval_col = p_value, build = "hg38"\n      ),\n    description = "<<trait_upper>> Manhattan plot for ALL populations",\n    format      = "qs"\n  ),\n  tar_target(\n    <<trait_lower>>_meta_common_variants_ALL,\n    extract_common_variants(\n      parquet_files        = <<trait_lower>>_sumstats_prepped_combined_ALL,\n      trait_type           = "<<trait_type>>",\n      select_per_mb_window = TRUE\n    ),\n    description = "<<trait_upper>> common variants across all studies",\n    <<resources_block_end>>  ),\n  tar_target(\n    <<trait_lower>>_mds_ALL,\n    calculate_mr_mega_mds(<<trait_lower>>_meta_common_variants_ALL),\n    description = "<<trait_upper>> MDS for all studies"\n  )\n)\n\nlist(\n  tar_target(\n    <<trait_lower>>_manifest,\n    <<trait_lower>>_manifest_df\n  ),\n\n  <<trait_lower>>_mapped,\n\n  tar_file_read(\n    hm3,\n    "<<hm3_path>>",\n    vroom::vroom(!!.x)\n  ),\n\n  <<trait_lower>>_ldsc_rg_qc_targets,\n  <<trait_lower>>_prep_meta_ancestry_targets,\n  <<trait_lower>>_prep_meta_all_targets,\n\n  tar_target(\n    <<trait_lower>>_sumstats_summary_cleaned,\n    <<trait_lower>>_sumstats_summary_ALL |>\n      dplyr::mutate(file = basename(file)) |>\n      dplyr::mutate(file = stringr::str_to_lower(stringr::str_replace(file, ".parquet", ""))) |>\n      dplyr::left_join(<<trait_lower>>_manifest |>\n        dplyr::mutate(file = stringr::str_to_lower(stringr::str_replace(file, ".txt.gz", "")))) |>\n      dplyr::mutate(cohort = stringr::str_to_upper(cohort))\n  ),\n\n  tar_target(\n    <<trait_lower>>_cohort_overview,\n    <<trait_lower>>_sumstats_summary_cleaned |>\n      <<cohort_overview_select>> |>\n      dplyr::group_by(ancestry, cohort) |>\n      dplyr::summarize(dplyr::across(where(is.numeric), sum)) |>\n      gt::gt() |>\n      gt::<<cohort_overview_fmt>> |>\n      gt::summary_rows(\n        columns  = -cohort,\n        fns      = list(label = "Total", fn = "sum"),\n        fmt      = ~gt::fmt_number(., decimals = 0),\n        side     = "bottom"\n      ) |>\n      gt::grand_summary_rows(\n        columns = -cohort,\n        fns     = list(fn = "sum", label = "Overall"),\n        fmt     = ~gt::fmt_number(., decimals = 0)\n      ) |>\n      gt::tab_style(\n        style     = list(gt::cell_text(weight = "bold")),\n        locations = gt::cells_column_labels()\n      )\n  ),\n\n  tar_target(\n    <<trait_lower>>_precision_plot,\n    <<trait_lower>>_sumstats_summary_cleaned |>\n      ggplot2::ggplot(ggplot2::aes(x = sqrt(n_total), y = 1 / median_se, color = ancestry)) +\n      ggplot2::geom_point() +\n      ggrepel::geom_text_repel(ggplot2::aes(label = cohort)) +\n      ggplot2::labs(\n        x     = latex2exp::TeX("$\\\\\\\\sqrt{N}$"),\n        y     = latex2exp::TeX("$1/Median(SE)$"),\n        color = "Ancestry"\n      ) +\n      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4))) +\n      ggsci::scale_color_npg()\n  )\n)\n',
    .open = "<<", .close = ">>"
  )

  return(pipeline_code)
}

# Internal helper: serialize a manifest data frame as a tibble::tribble() call.
# All values are coerced to character and double-quoted.
.format_manifest_tribble <- function(manifest_df, trait_lower) {
  col_header <- paste0("  ", paste0("~", names(manifest_df), collapse = ", "))
  rows <- apply(manifest_df, 1, function(row) {
    vals <- vapply(as.list(row), function(x) paste0('"', as.character(x), '"'), character(1))
    paste0("  ", paste(vals, collapse = ", "))
  })
  paste0(
    trait_lower, "_manifest_df <- tibble::tribble(\n",
    col_header, ",\n",
    paste(rows, collapse = ",\n"), "\n",
    ")"
  )
}
