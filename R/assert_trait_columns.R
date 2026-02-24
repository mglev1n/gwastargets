#' Assert trait-type-specific columns are present
#'
#' @description
#' Internal helper that validates the column schema of a summary statistics
#' dataset against the declared `trait_type`. For binary traits it errors when
#' `CaseN`, `ControlN`, or `EffectiveN` are absent. For quantitative traits it
#' warns when binary-specific columns are detected.
#'
#' @param schema_names Character vector of column names (e.g. `names(df)` or
#'   `names(arrow::schema(ds))`).
#' @param trait_type `"binary"` or `"quantitative"`.
#'
#' @return `NULL` invisibly; called for its side-effects (errors/warnings).
#'
#' @keywords internal
#' @noRd
assert_trait_columns <- function(schema_names, trait_type) {
  if (trait_type == "binary") {
    required <- c("CaseN", "ControlN", "EffectiveN")
    missing_cols <- setdiff(required, schema_names)
    if (length(missing_cols) > 0) {
      cli::cli_abort(c(
        "trait_type = {.val binary} but required columns are missing: {.field {missing_cols}}",
        "i" = "If this is a quantitative trait, set {.arg trait_type} = {.val quantitative}",
        "i" = "If this is a binary trait, check that {.field CaseN}, {.field ControlN}, and {.field EffectiveN} are present in the input data"
      ))
    }
  } else {
    # Warn rather than error: some quantitative sumstats carry these as metadata
    binary_cols_present <- intersect(c("CaseN", "ControlN"), schema_names)
    if (length(binary_cols_present) > 0) {
      cli::cli_warn(c(
        "trait_type = {.val quantitative} but binary columns found: {.field {binary_cols_present}}",
        "i" = "If this is a binary trait, set {.arg trait_type} = {.val binary}"
      ))
    }
  }
}
