test_that("returns a character string", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155")
  )
  expect_type(result, "character")
  expect_length(result, 1)
})

test_that("contains trait name in uppercase and lowercase", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("BMI", trait_type = "quantitative",
                                n_col = "N", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155")
  )
  expect_true(grepl("BMI", result))
  expect_true(grepl("bmi", result))
})

test_that("binary: output contains EffectiveN, n_cases, n_controls", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("T2D", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155")
  )
  expect_true(grepl("EffectiveN", result))
  expect_true(grepl("n_cases",   result))
  expect_true(grepl("n_controls", result))
})

test_that("quantitative: output contains n_total but NOT n_cases", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("HDL", trait_type = "quantitative",
                                n_col = "N", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155")
  )
  expect_true(grepl("n_total",  result))
  expect_false(grepl("n_cases", result))
})

test_that("generated code contains tribble with study and tar_name columns", {
  mdf <- make_manifest_df()
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = mdf,
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155")
  )
  expect_true(grepl("tibble::tribble", result))
  expect_true(grepl("~path",     result))
  expect_true(grepl("~study",    result))
  expect_true(grepl("~tar_name", result))
  # study values appear as-is
  expect_true(grepl('"UKBB"', result, fixed = TRUE))
  expect_true(grepl('"MVP"',  result, fixed = TRUE))
  # tar_name values are study_ancestry
  expect_true(grepl('"UKBB_EUR"', result, fixed = TRUE))
  expect_true(grepl('"MVP_AFR"',  result, fixed = TRUE))
})

test_that("generated code contains manifest variable assigned with trait prefix", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155")
  )
  expect_true(grepl("cad_manifest_df <- tibble::tribble", result))
})

test_that("generated code embeds hm3_path in tar_file_read call", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/my/custom/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155")
  )
  expect_true(grepl("/my/custom/w_hm3.snplist", result, fixed = TRUE))
})

test_that("generated code embeds dbsnp_path in prep_gwas call", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/my/custom/dbSNP155")
  )
  expect_true(grepl("/my/custom/dbSNP155", result, fixed = TRUE))
})

test_that("crew_controller = NULL produces no resources block", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155",
                                crew_controller = NULL)
  )
  expect_false(grepl("tar_resources_crew", result))
})

test_that("crew_controller embeds controller name in resources blocks", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155",
                                crew_controller = "my_controller")
  )
  expect_true(grepl("tar_resources_crew", result))
  expect_true(grepl("my_controller", result))
})

# --- trait_type / n_col errors (fire before hm3_path validation) ---

test_that("errors when trait_type missing", {
  expect_error(
    generate_gwas_meta_pipeline("CAD", n_col = "EffectiveN",
                                manifest_df = make_manifest_df()),
    "trait_type"
  )
})

test_that("errors when trait_type invalid", {
  expect_error(
    generate_gwas_meta_pipeline("CAD", trait_type = "continuous", n_col = "N",
                                manifest_df = make_manifest_df()),
    "trait_type"
  )
})

test_that("errors when n_col missing", {
  expect_error(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                manifest_df = make_manifest_df()),
    "n_col"
  )
})

test_that("errors when n_col invalid", {
  expect_error(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary", n_col = "BadCol",
                                manifest_df = make_manifest_df()),
    "n_col"
  )
})

test_that("errors when binary + N (should use EffectiveN)", {
  expect_error(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary", n_col = "N",
                                manifest_df = make_manifest_df()),
    "EffectiveN"
  )
})

test_that("warns when quantitative + EffectiveN", {
  expect_warning(
    generate_gwas_meta_pipeline("BMI", trait_type = "quantitative",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155"),
    "unusual"
  )
})

# --- manifest_df validation ---

test_that("errors when manifest_df missing", {
  expect_error(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary", n_col = "EffectiveN"),
    "manifest_df"
  )
})

test_that("errors when manifest_df is missing required columns", {
  bad_manifest <- make_manifest_df()[, c("path", "file", "cohort")]  # missing ancestry, study
  expect_error(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary", n_col = "EffectiveN",
                                manifest_df = bad_manifest),
    "missing"
  )
})

test_that("errors when manifest_df has duplicate study values", {
  dup_manifest <- make_manifest_df()
  dup_manifest$study <- c("UKBB", "UKBB")  # duplicate
  expect_error(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary", n_col = "EffectiveN",
                                manifest_df = dup_manifest),
    "duplicate"
  )
})

test_that("errors when ancestry values are not in the valid set", {
  bad_manifest <- make_manifest_df()
  bad_manifest$ancestry[1] <- "SAS"
  expect_error(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary", n_col = "EffectiveN",
                                manifest_df = bad_manifest),
    "unsupported"
  )
})

test_that("generated code uses names = tar_name in tar_map", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155")
  )
  expect_true(grepl("names  = tar_name", result, fixed = TRUE))
  expect_false(grepl("names  = study",   result, fixed = TRUE))
})

test_that("warns when manifest_df paths do not exist", {
  expect_warning(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary", n_col = "EffectiveN",
                                manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155"),
    "do not exist"
  )
})

# --- hm3_path validation ---

test_that("errors when hm3_path missing", {
  expect_error(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary", n_col = "EffectiveN",
                                manifest_df = make_manifest_df()),
    "hm3_path"
  )
})

test_that("warns when hm3_path does not exist", {
  expect_warning(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary", n_col = "EffectiveN",
                                manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155"),
    "hm3_path"
  )
})

# --- dbsnp_path validation ---

test_that("errors when dbsnp_path missing", {
  expect_error(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary", n_col = "EffectiveN",
                                manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist"),
    "dbsnp_path"
  )
})

test_that("warns when dbsnp_path does not exist", {
  expect_warning(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary", n_col = "EffectiveN",
                                manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155"),
    "dbsnp_path"
  )
})

# --- Manhattan PDF + PNG targets ---

test_that("generated code contains manhattan PDF and PNG targets with ggsave", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155")
  )
  expect_true(grepl("meta_manhattan_pdf", result))
  expect_true(grepl("meta_manhattan_png", result))
  expect_true(grepl("ggsave", result))
  expect_true(grepl("format      = \"file\"", result, fixed = TRUE))
})

test_that("manhattan PDF targets use output_base_dir in ggsave paths", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155",
                                output_base_dir = "MyOutput")
  )
  expect_true(grepl("MyOutput", result))
  expect_true(grepl('file.path("MyOutput"', result, fixed = TRUE))
})

test_that("custom manhattan dimensions appear in generated code", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155",
                                manhattan_width = 20, manhattan_height = 8)
  )
  expect_true(grepl("width = 20", result, fixed = TRUE))
  expect_true(grepl("height = 8", result, fixed = TRUE))
})

test_that("default manhattan dimensions appear when not specified", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155")
  )
  expect_true(grepl("width = 16", result, fixed = TRUE))
  expect_true(grepl("height = 6", result, fixed = TRUE))
})

test_that("ALL population manhattan PDF and PNG targets are present", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155")
  )
  expect_true(grepl("meta_manhattan_pdf_ALL", result))
  expect_true(grepl("manhattan_ALL.pdf", result))
  expect_true(grepl("meta_manhattan_png_ALL", result))
  expect_true(grepl("manhattan_ALL.png", result))
})

test_that("PNG targets include dpi parameter", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155")
  )
  expect_true(grepl("dpi = 300", result, fixed = TRUE))
})

test_that("custom manhattan_dpi appears in generated code", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155",
                                manhattan_dpi = 600)
  )
  expect_true(grepl("dpi = 600", result, fixed = TRUE))
})

# --- col_* column mapping in manifest ---

test_that("manifest without col_* columns: no column_map in generated code", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155")
  )
  expect_false(grepl("column_map", result))
  expect_false(grepl("build_column_map", result))
})

test_that("manifest with col_* columns: generated code includes build_column_map()", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN",
                                manifest_df = make_manifest_df_with_colmap(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155")
  )
  expect_true(grepl("column_map", result))
  expect_true(grepl("build_column_map", result))
  expect_true(grepl("col_eaf", result))
  expect_true(grepl("col_beta", result))
})

test_that("NA values in col_* serialize as unquoted NA in tribble", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN",
                                manifest_df = make_manifest_df_with_colmap(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155")
  )
  # Should have unquoted NA, not "NA"
  expect_true(grepl('"MY_FREQ"', result, fixed = TRUE))
  expect_true(grepl('"BETA_VAL"', result, fixed = TRUE))
  # The tribble should contain unquoted NA values
  # Extract the tribble section and check for NA not wrapped in quotes
  expect_false(grepl('"NA"', result, fixed = TRUE))
})

test_that("unrecognized col_* columns produce a warning", {
  mdf <- make_manifest_df()
  mdf$col_garbage <- c("FOO", "BAR")
  expect_warning(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = mdf,
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                dbsnp_path = "/nonexistent/dbSNP155"),
    "unrecognized"
  )
})

test_that("no file-existence warning when all paths exist", {
  withr::with_tempdir({
    mdf <- make_manifest_df()
    mdf$path <- c(file.path(getwd(), "c1.txt.gz"), file.path(getwd(), "c2.txt.gz"))
    file.create(mdf$path)
    hm3 <- file.path(getwd(), "w_hm3.snplist")
    file.create(hm3)
    dbsnp <- file.path(getwd(), "dbSNP155")
    dir.create(dbsnp)
    expect_no_warning(
      generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                  n_col = "EffectiveN", manifest_df = mdf,
                                  hm3_path = hm3, dbsnp_path = dbsnp)
    )
  })
})
