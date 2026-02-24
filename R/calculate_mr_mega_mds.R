#' Calculate MR-MEGA multidimensional scaling coordinates
#'
#' @title Calculate MR-MEGA MDS coordinates from allele frequency data
#'
#' @description
#' Computes multidimensional scaling (MDS) coordinates from study-level effect
#' allele frequency (EAF) data. The approach follows the MR-MEGA method
#' (Mägi et al. 2017): a population-genetic distance matrix is computed from
#' EAF values, double-centred, and decomposed via SVD to produce principal
#' coordinates.
#'
#' @param common_variants_data A data frame (long format) with one row per
#'   variant × study combination. Must contain columns identified by
#'   `chr_col`, `pos_col`, `ea_col`, `oa_col`, `eaf_col`, and `file_col`.
#'   Typically the output of [extract_common_variants()].
#' @param pcCount Number of principal coordinates to return. Default `4`.
#' @param chr_col Bare column name for chromosome. Default `CHR`.
#' @param pos_col Bare column name for position. Default `POS_38`.
#' @param rsid_col Bare column name for rsid. Default `RSID`.
#' @param ea_col Bare column name for effect allele. Default `EffectAllele`.
#' @param oa_col Bare column name for other allele. Default `OtherAllele`.
#' @param eaf_col Bare column name for EAF. Default `minor_af`.
#' @param file_col Bare column name for study/file identifier. Default `file`.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{`pc_coordinates`}{Data frame of MDS coordinates (`pcCount` columns).}
#'     \item{`distance_matrix`}{Population-genetic distance matrix.}
#'     \item{`centered_distance_matrix`}{Double-centred distance matrix.}
#'     \item{`eigenvalues`}{Eigenvalues from SVD.}
#'     \item{`eaf_matrix`}{Wide EAF matrix (variants × studies).}
#'     \item{`used_variants`}{Variants with complete EAF data.}
#'     \item{`studies`}{Character vector of study names.}
#'     \item{`summary`}{List with `n_studies`, `n_variants_used`,
#'       `n_variants_total`, `pc_count`.}
#'   }
#'
#' @references
#' Mägi R, et al. (2017). Trans-ethnic meta-regression of genome-wide
#' association studies accounting for ancestry increases power for discovery
#' and improves fine-mapping resolution. *Hum Mol Genet*, 26(18), 3639–3650.
#'
#' @examples
#' \dontrun{
#' common <- extract_common_variants(files, trait_type = "binary")
#' mds    <- calculate_mr_mega_mds(common)
#' }
#'
#' @importFrom cli cli_alert_info cli_alert_success
#' @importFrom dplyr select distinct arrange
#' @importFrom rlang enquo as_name sym
#' @importFrom tidyr pivot_wider
#' @export
calculate_mr_mega_mds <- function(common_variants_data,
                                  pcCount    = 4,
                                  chr_col    = CHR,
                                  pos_col    = POS_38,
                                  rsid_col   = RSID,
                                  ea_col     = EffectAllele,
                                  oa_col     = OtherAllele,
                                  eaf_col    = minor_af,
                                  file_col   = file) {

  double_center <- function(D) {
    row_meanD <- rowMeans(D)
    mean_dist <- mean(row_meanD)
    (-0.5) * (D - row_meanD[col(D)] - row_meanD[row(D)] + mean_dist)
  }

  dist_pop <- function(X) {
    ind_X  <- (!is.na(X))
    denomi <- ind_X %*% t(ind_X)
    nomi   <- X %*% t(X)
    vec    <- apply(X, 1, function(x) sum(x * x))
    d      <- -2 * nomi + vec[col(nomi)] + vec[row(nomi)]
    diag(d) <- 0
    d / denomi
  }

  studies   <- unique(common_variants_data[[rlang::as_name(rlang::enquo(file_col))]])
  n_studies <- length(studies)
  cli::cli_alert_info("Calculating MDS for {n_studies} studies")

  eaf_wide <- common_variants_data |>
    dplyr::select({{ chr_col }}, {{ pos_col }}, {{ ea_col }}, {{ oa_col }},
                  {{ file_col }}, {{ eaf_col }}) |>
    dplyr::distinct({{ chr_col }}, {{ pos_col }}, {{ ea_col }}, {{ oa_col }},
                    {{ file_col }}, .keep_all = TRUE) |>
    tidyr::pivot_wider(
      names_from  = {{ file_col }},
      values_from = {{ eaf_col }},
      id_cols     = c({{ chr_col }}, {{ pos_col }}, {{ ea_col }}, {{ oa_col }})
    ) |>
    dplyr::arrange({{ chr_col }}, {{ pos_col }})

  eaf_matrix        <- as.matrix(dplyr::select(eaf_wide, dplyr::all_of(studies)))
  complete_variants <- complete.cases(eaf_matrix)
  eaf_matrix        <- eaf_matrix[complete_variants, ]
  used_variants     <- eaf_wide[complete_variants, ]

  n_variants_used <- nrow(eaf_matrix)
  cli::cli_alert_info("Using {n_variants_used} variants with complete EAF data across all studies")

  dist_matrix <- dist_pop(t(eaf_matrix))
  rownames(dist_matrix) <- colnames(dist_matrix) <- studies
  centered_dist <- double_center(dist_matrix)

  svd_result   <- svd(centered_dist)
  eigenvalues  <- ifelse(svd_result$d > 0, svd_result$d, 0)
  eigenvectors <- svd_result$u

  # Guard against pcCount > n_studies
  actual_pc <- min(pcCount, length(eigenvalues))
  pc_coords <- eigenvectors %*% diag(sqrt(eigenvalues))[, seq_len(actual_pc), drop = FALSE]
  rownames(pc_coords) <- studies
  colnames(pc_coords) <- paste0("PC", seq_len(actual_pc))

  variant_ids <- paste0(
    used_variants[[rlang::as_name(rlang::enquo(chr_col))]],
    ":",
    used_variants[[rlang::as_name(rlang::enquo(pos_col))]],
    "_",
    used_variants[[rlang::as_name(rlang::enquo(ea_col))]],
    "/",
    used_variants[[rlang::as_name(rlang::enquo(oa_col))]]
  )

  rownames(eaf_matrix)  <- variant_ids
  rownames(centered_dist) <- colnames(centered_dist) <- studies

  cli::cli_alert_success("MDS calculation complete")
  cli::cli_alert_info(
    "Variance explained by PC1: {round(100 * eigenvalues[1] / sum(eigenvalues), 2)}%"
  )
  if (actual_pc > 1) {
    cli::cli_alert_info(
      "Variance explained by PC2: {round(100 * eigenvalues[2] / sum(eigenvalues), 2)}%"
    )
  }

  list(
    pc_coordinates           = as.data.frame(pc_coords),
    distance_matrix          = as.data.frame(dist_matrix),
    centered_distance_matrix = as.data.frame(centered_dist),
    eigenvalues              = eigenvalues,
    eaf_matrix               = as.data.frame(eaf_matrix),
    used_variants            = used_variants,
    studies                  = studies,
    summary = list(
      n_studies        = n_studies,
      n_variants_used  = n_variants_used,
      n_variants_total = nrow(eaf_wide),
      pc_count         = actual_pc
    )
  )
}
