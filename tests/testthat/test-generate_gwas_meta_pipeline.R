test_that("returns a character string", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist")
  )
  expect_type(result, "character")
  expect_length(result, 1)
})

test_that("contains trait name in uppercase and lowercase", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("BMI", trait_type = "quantitative",
                                n_col = "N", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist")
  )
  expect_true(grepl("BMI", result))
  expect_true(grepl("bmi", result))
})

test_that("binary: output contains EffectiveN, n_cases, n_controls", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("T2D", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist")
  )
  expect_true(grepl("EffectiveN", result))
  expect_true(grepl("n_cases",   result))
  expect_true(grepl("n_controls", result))
})

test_that("quantitative: output contains n_total but NOT n_cases", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("HDL", trait_type = "quantitative",
                                n_col = "N", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist")
  )
  expect_true(grepl("n_total",  result))
  expect_false(grepl("n_cases", result))
})

test_that("generated code contains tribble with study values from manifest", {
  mdf <- make_manifest_df()
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = mdf,
                                hm3_path = "/nonexistent/w_hm3.snplist")
  )
  expect_true(grepl("tibble::tribble", result))
  expect_true(grepl("UKBB_EUR", result))
  expect_true(grepl("MVP_AFR",  result))
  expect_true(grepl("~path",    result))
  expect_true(grepl("~study",   result))
})

test_that("generated code contains manifest variable assigned with trait prefix", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist")
  )
  expect_true(grepl("cad_manifest_df <- tibble::tribble", result))
})

test_that("generated code embeds hm3_path in tar_file_read call", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/my/custom/w_hm3.snplist")
  )
  expect_true(grepl("/my/custom/w_hm3.snplist", result, fixed = TRUE))
})

test_that("crew_controller = NULL produces no resources block", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
                                crew_controller = NULL)
  )
  expect_false(grepl("tar_resources_crew", result))
})

test_that("crew_controller embeds controller name in resources blocks", {
  result <- suppressWarnings(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                n_col = "EffectiveN", manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist",
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
                                hm3_path = "/nonexistent/w_hm3.snplist"),
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
  dup_manifest$study <- c("UKBB_EUR", "UKBB_EUR")  # duplicate
  expect_error(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary", n_col = "EffectiveN",
                                manifest_df = dup_manifest),
    "duplicate"
  )
})

test_that("warns when manifest_df paths do not exist", {
  expect_warning(
    generate_gwas_meta_pipeline("CAD", trait_type = "binary", n_col = "EffectiveN",
                                manifest_df = make_manifest_df(),
                                hm3_path = "/nonexistent/w_hm3.snplist"),
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
                                hm3_path = "/nonexistent/w_hm3.snplist"),
    "hm3_path"
  )
})

test_that("no file-existence warning when all paths exist", {
  withr::with_tempdir({
    mdf <- make_manifest_df()
    mdf$path <- c(file.path(getwd(), "c1.txt.gz"), file.path(getwd(), "c2.txt.gz"))
    file.create(mdf$path)
    hm3 <- file.path(getwd(), "w_hm3.snplist")
    file.create(hm3)
    expect_no_warning(
      generate_gwas_meta_pipeline("CAD", trait_type = "binary",
                                  n_col = "EffectiveN", manifest_df = mdf,
                                  hm3_path = hm3)
    )
  })
})
