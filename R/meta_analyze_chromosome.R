#' IVW meta-analysis for a single chromosome
#'
#' @description
#' Internal worker called by [meta_analyze_ivw()] for each chromosome. Filters
#' to the requested chromosome, applies MAC filtering, computes inverse-variance
#' weighted meta-analysis, and returns heterogeneity statistics (Q, I2).
#'
#' @param ds An Arrow dataset with standardised column names (as produced by
#'   [meta_analyze_ivw()]).
#' @param chr Integer chromosome number.
#' @param min_mac Minimum minor allele count.
#' @param trait_type `"binary"` or `"quantitative"`.
#' @param chr_col Bare column name for chromosome. Default `CHR`.
#' @param pos_col Bare column name for position. Default `POS`.
#' @param rsid_col Bare column name for rsid. Default `RSID`.
#' @param ea_col Bare column name for effect allele. Default `EFFECT_ALLELE`.
#' @param oa_col Bare column name for other allele. Default `OTHER_ALLELE`.
#' @param beta_col Bare column name for beta. Default `BETA`.
#' @param se_col Bare column name for SE. Default `SE`.
#' @param eaf_col Bare column name for EAF. Default `EAF`.
#' @param n_col Bare column name for N. Default `N`.
#' @param effective_n_col Bare column name for effective N. Default
#'   `EFFECTIVE_N`.
#' @param case_n_col Bare column name for case N. Default `CASE_N`.
#' @param control_n_col Bare column name for control N. Default `CONTROL_N`.
#'
#' @return A data frame of meta-analysed variants for the chromosome.
#'
#' @keywords internal
#' @noRd
meta_analyze_chromosome <- function(ds,
                                    chr,
                                    min_mac,
                                    trait_type,
                                    chr_col         = CHR,
                                    pos_col         = POS,
                                    rsid_col        = RSID,
                                    ea_col          = EFFECT_ALLELE,
                                    oa_col          = OTHER_ALLELE,
                                    beta_col        = BETA,
                                    se_col          = SE,
                                    eaf_col         = EAF,
                                    n_col           = N,
                                    effective_n_col = EFFECTIVE_N,
                                    case_n_col      = CASE_N,
                                    control_n_col   = CONTROL_N) {

  chr_col_name  <- rlang::as_name(rlang::enquo(chr_col))
  pos_col_name  <- rlang::as_name(rlang::enquo(pos_col))
  rsid_col_name <- rlang::as_name(rlang::enquo(rsid_col))
  ea_col_name   <- rlang::as_name(rlang::enquo(ea_col))
  oa_col_name   <- rlang::as_name(rlang::enquo(oa_col))
  eaf_col_name  <- rlang::as_name(rlang::enquo(eaf_col))
  n_col_name    <- rlang::as_name(rlang::enquo(n_col))

  by <- c(chr_col_name, pos_col_name, rsid_col_name, ea_col_name, oa_col_name)

  # Columns to sum in the group_by step — differs by trait_type
  if (trait_type == "binary") {
    n_cols_to_sum <- c(
      eaf_col_name,
      rlang::as_name(rlang::enquo(case_n_col)),
      rlang::as_name(rlang::enquo(control_n_col)),
      n_col_name,
      rlang::as_name(rlang::enquo(effective_n_col))
    )
  } else {
    n_cols_to_sum <- c(
      eaf_col_name,
      n_col_name,
      rlang::as_name(rlang::enquo(effective_n_col))
    )
  }

  meta_results <- ds |>
    dplyr::filter({{ chr_col }} == chr) |>
    dplyr::mutate(
      minor_af = pmin({{ eaf_col }}, 1 - {{ eaf_col }}),
      mac      = minor_af * {{ n_col }} * 2
    ) |>
    dplyr::filter(mac > min_mac) |>
    dplyr::filter(
      !is.na({{ beta_col }}),
      !is.na({{ se_col }}),
      {{ se_col }} > 0
    ) |>
    dplyr::mutate(
      W   = 1 / ({{ se_col }})^2,
      B   = {{ beta_col }} * W,
      WB2 = (B^2) / W,
      dplyr::across(dplyr::any_of(eaf_col_name), ~.x * {{ n_col }})
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(
      n_contributions = dplyr::n(),
      dplyr::across(
        dplyr::any_of(eaf_col_name),
        ~ sum(dplyr::if_else(is.na(.x), 0, {{ n_col }}), na.rm = TRUE),
        .names = "N_{.col}"
      ),
      dplyr::across(dplyr::any_of(pos_col_name), ~ min(.x)),
      WB2 = sum(WB2, na.rm = TRUE),
      dplyr::across(dplyr::any_of(c("W", "B")), sum),
      dplyr::across(dplyr::any_of(n_cols_to_sum), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::collect() |>
    dplyr::mutate(
      B_w = B,
      B   = B / W,
      SE  = 1 / sqrt(W),
      dplyr::across(
        dplyr::any_of(eaf_col_name),
        ~.x / !!rlang::sym(paste0("N_", eaf_col_name))
      ),
      CHR = chr
    ) |>
    dplyr::select(-dplyr::any_of(paste0("N_", eaf_col_name))) |>
    dplyr::mutate(
      z_score = B / SE,
      p_value = 2 * stats::pnorm(-abs(z_score)),
      Q       = WB2 - (B_w^2) / W,
      Q_df    = pmax(n_contributions - 1L, 0L),
      Q_pval  = stats::pchisq(Q, df = Q_df, lower.tail = FALSE),
      Q_pval  = dplyr::if_else(Q_df == 0, 1, Q_pval),
      I2      = dplyr::if_else(Q > 0, pmax((Q - Q_df) / Q, 0), 0)
    )

  # Final column selection differs by trait_type
  if (trait_type == "binary") {
    meta_results <- meta_results |>
      dplyr::select(
        CHR,
        POS_38       = {{ pos_col }},
        RSID         = {{ rsid_col }},
        EffectAllele = {{ ea_col }},
        OtherAllele  = {{ oa_col }},
        B, SE, p_value, z_score,
        EAF          = {{ eaf_col }},
        n_contributions,
        CaseN        = {{ case_n_col }},
        ControlN     = {{ control_n_col }},
        N            = {{ n_col }},
        EffectiveN   = {{ effective_n_col }},
        Q, Q_df, Q_pval, I2
      )
  } else {
    meta_results <- meta_results |>
      dplyr::select(
        CHR,
        POS_38       = {{ pos_col }},
        RSID         = {{ rsid_col }},
        EffectAllele = {{ ea_col }},
        OtherAllele  = {{ oa_col }},
        B, SE, p_value, z_score,
        EAF          = {{ eaf_col }},
        n_contributions,
        N            = {{ n_col }},
        EffectiveN   = {{ effective_n_col }},
        Q, Q_df, Q_pval, I2
      )
  }

  meta_results <- meta_results |> dplyr::arrange(p_value)

  if (nrow(meta_results) == 0) {
    cli::cli_alert_warning("No variants remaining for chromosome {chr}")
  }

  return(meta_results)
}
