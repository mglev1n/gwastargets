test_that("returns a single character string", {
  result <- generate_report_chunks("CAD", manifest_df = make_manifest_df())
  expect_type(result, "character")
  expect_length(result, 1)
})

test_that("contains cohort overview and precision plot tar_read", {
  result <- generate_report_chunks("CAD", manifest_df = make_manifest_df())
  expect_true(grepl("tar_read(cad_cohort_overview)", result, fixed = TRUE))
  expect_true(grepl("tar_read(cad_precision_plot)", result, fixed = TRUE))
})

test_that("contains trait name in header", {
  result <- generate_report_chunks("CAD", manifest_df = make_manifest_df())
  expect_true(grepl("## CAD Results", result, fixed = TRUE))
})

test_that("manhattan chunks use knitr::include_graphics with PNG targets", {
  result <- generate_report_chunks("CAD", manifest_df = make_manifest_df())
  expect_true(grepl("knitr::include_graphics", result))
  expect_true(grepl("meta_manhattan_png_ALL", result))
  expect_true(grepl("out-width", result))
  # PNG targets used in report, not PDF
  expect_false(grepl("meta_manhattan_pdf", result))
})

test_that("ALL population manhattan is always present in tabset", {
  result <- generate_report_chunks("CAD", manifest_df = make_manifest_df())
  expect_true(grepl("#### All Populations", result, fixed = TRUE))
  expect_true(grepl("cad_meta_manhattan_png_ALL", result))
})

test_that("manhattan section uses panel-tabset", {
  result <- generate_report_chunks("CAD", manifest_df = make_manifest_df())
  expect_true(grepl("panel-tabset", result))
  expect_true(grepl("Manhattan Plots & Genome-wide Significant Loci", result, fixed = TRUE))
})

test_that("LDSC heatmap chunks generated for ancestries with 2+ cohorts", {
  # make_manifest_df has EUR and AFR each with 1 cohort -> no LDSC section
  mdf_single <- make_manifest_df()
  result_single <- generate_report_chunks("CAD", manifest_df = mdf_single)
  expect_false(grepl("ldsc_rg_qc_heatmap", result_single))
  expect_false(grepl("LDSC Genetic Correlation QC", result_single))

  # Manifest with 2 EUR cohorts -> EUR LDSC tab
  mdf_multi <- data.frame(
    path     = c("/a.txt.gz", "/b.txt.gz", "/c.txt.gz"),
    file     = c("a.txt.gz", "b.txt.gz", "c.txt.gz"),
    cohort   = c("UKBB", "BioVU", "MVP"),
    ancestry = c("EUR", "EUR", "AFR"),
    study    = c("UKBB", "BioVU", "MVP"),
    stringsAsFactors = FALSE
  )
  result_multi <- generate_report_chunks("CAD", manifest_df = mdf_multi)
  expect_true(grepl("LDSC Genetic Correlation QC", result_multi, fixed = TRUE))
  expect_true(grepl("cad_ldsc_rg_qc_heatmap_EUR", result_multi))
  # AFR has only 1 cohort -> no AFR LDSC tab
  expect_false(grepl("ldsc_rg_qc_heatmap_AFR", result_multi))
})

test_that("LDSC section uses panel-tabset when present", {
  mdf_multi <- data.frame(
    path     = c("/a.txt.gz", "/b.txt.gz"),
    file     = c("a.txt.gz", "b.txt.gz"),
    cohort   = c("UKBB", "BioVU"),
    ancestry = c("EUR", "EUR"),
    study    = c("UKBB", "BioVU"),
    stringsAsFactors = FALSE
  )
  result <- generate_report_chunks("CAD", manifest_df = mdf_multi)
  # Should have two panel-tabset markers (LDSC + Manhattan)
  expect_equal(length(gregexpr("panel-tabset", result)[[1]]), 2)
})

test_that("per-ancestry manhattan tabs present only for 2+ cohort ancestries", {
  mdf_multi <- data.frame(
    path     = c("/a.txt.gz", "/b.txt.gz", "/c.txt.gz"),
    file     = c("a.txt.gz", "b.txt.gz", "c.txt.gz"),
    cohort   = c("UKBB", "BioVU", "MVP"),
    ancestry = c("EUR", "EUR", "AFR"),
    study    = c("UKBB", "BioVU", "MVP"),
    stringsAsFactors = FALSE
  )
  result <- generate_report_chunks("CAD", manifest_df = mdf_multi)
  expect_true(grepl("#### EUR", result, fixed = TRUE))
  expect_true(grepl("cad_meta_manhattan_png_EUR", result))
  expect_false(grepl("#### AFR", result, fixed = TRUE))
})

test_that("loci chunks included by default with DT::datatable and csv export", {
  result <- generate_report_chunks("CAD", manifest_df = make_manifest_df())
  expect_true(grepl("cad_meta_loci_ALL", result))
  expect_true(grepl("DT::datatable", result, fixed = TRUE))
  expect_true(grepl("Buttons", result))
  expect_true(grepl("csv", result))
  expect_false(grepl("excel", result))
})

test_that("loci chunks excluded when include_loci = FALSE", {
  result <- generate_report_chunks("CAD", manifest_df = make_manifest_df(),
                                   include_loci = FALSE)
  expect_false(grepl("meta_loci", result))
})

test_that("per-ancestry loci chunks controlled by include_loci", {
  mdf_multi <- data.frame(
    path     = c("/a.txt.gz", "/b.txt.gz"),
    file     = c("a.txt.gz", "b.txt.gz"),
    cohort   = c("UKBB", "BioVU"),
    ancestry = c("EUR", "EUR"),
    study    = c("UKBB", "BioVU"),
    stringsAsFactors = FALSE
  )
  result_with <- generate_report_chunks("CAD", manifest_df = mdf_multi, include_loci = TRUE)
  expect_true(grepl("cad_meta_loci_EUR", result_with))

  result_without <- generate_report_chunks("CAD", manifest_df = mdf_multi, include_loci = FALSE)
  expect_false(grepl("meta_loci", result_without))
})

# --- output_file parameter ---

test_that("output_file = NULL returns visibly", {
  result <- withVisible(
    generate_report_chunks("CAD", manifest_df = make_manifest_df())
  )
  expect_true(result$visible)
})

test_that("output_file writes file, creates parent dirs, returns invisibly", {
  withr::with_tempdir({
    outpath <- file.path(getwd(), "sub", "dir", "report.qmd")
    result <- withVisible(
      generate_report_chunks("CAD", manifest_df = make_manifest_df(),
                             output_file = outpath)
    )
    expect_false(result$visible)
    expect_true(file.exists(outpath))
    # Return value matches file content
    file_content <- paste(readLines(outpath), collapse = "\n")
    expect_equal(as.character(result$value), file_content)
  })
})

test_that("errors on missing trait", {
  expect_error(
    generate_report_chunks(manifest_df = make_manifest_df()),
    "trait"
  )
})

test_that("errors on empty trait", {
  expect_error(
    generate_report_chunks("", manifest_df = make_manifest_df()),
    "trait"
  )
})

test_that("errors on missing manifest_df", {
  expect_error(
    generate_report_chunks("CAD"),
    "manifest_df"
  )
})

test_that("errors on manifest_df missing ancestry column", {
  bad_mdf <- data.frame(cohort = "UKBB", stringsAsFactors = FALSE)
  expect_error(
    generate_report_chunks("CAD", manifest_df = bad_mdf),
    "ancestry"
  )
})

test_that("errors on manifest_df missing cohort column", {
  bad_mdf <- data.frame(ancestry = "EUR", stringsAsFactors = FALSE)
  expect_error(
    generate_report_chunks("CAD", manifest_df = bad_mdf),
    "cohort"
  )
})
