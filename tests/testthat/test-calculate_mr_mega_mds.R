make_mds_input <- function(n_variants = 50, n_studies = 4) {
  set.seed(10)
  files <- paste0("study", seq_len(n_studies), ".parquet")
  purrr::map_dfr(files, function(f) {
    data.frame(
      CHR          = rep(1L, n_variants),
      POS_38       = seq(1e6L, by = 1e4L, length.out = n_variants),
      RSID         = paste0("rs", seq_len(n_variants)),
      EffectAllele = rep("A", n_variants),
      OtherAllele  = rep("G", n_variants),
      minor_af     = runif(n_variants, 0.05, 0.5),
      file         = f,
      stringsAsFactors = FALSE
    )
  })
}

test_that("returns list with expected elements", {
  input  <- make_mds_input()
  result <- calculate_mr_mega_mds(input, pcCount = 2)

  expect_type(result, "list")
  expect_true(all(c("pc_coordinates", "distance_matrix",
                    "centered_distance_matrix", "eigenvalues",
                    "eaf_matrix", "used_variants", "studies", "summary") %in%
                    names(result)))
})

test_that("pc_coordinates has pcCount columns", {
  input  <- make_mds_input(n_variants = 50, n_studies = 4)
  result <- calculate_mr_mega_mds(input, pcCount = 3)
  expect_equal(ncol(result$pc_coordinates), 3)
})

test_that("handles n_studies < pcCount gracefully (no error)", {
  input  <- make_mds_input(n_variants = 50, n_studies = 2)
  expect_no_error(calculate_mr_mega_mds(input, pcCount = 4))
})

test_that("n_variants_used reflects complete cases only", {
  input <- make_mds_input(n_variants = 50, n_studies = 3)
  # Introduce NAs in one study for a subset of variants
  idx   <- input$file == "study1.parquet" & input$POS_38 < 1.1e6
  input$minor_af[idx] <- NA

  result <- calculate_mr_mega_mds(input, pcCount = 2)
  expect_lt(result$summary$n_variants_used, 50)
})

test_that("distance_matrix is square with correct dimnames", {
  input  <- make_mds_input(n_variants = 50, n_studies = 3)
  result <- calculate_mr_mega_mds(input, pcCount = 2)
  dm     <- result$distance_matrix
  expect_equal(nrow(dm), ncol(dm))
  expect_equal(rownames(dm), colnames(dm))
})
