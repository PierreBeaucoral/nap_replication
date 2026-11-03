# ==============================================================================
# code/functions/within_share.R
# The within-country adaptation share 100 x E/S (stages 04, 07 and 13): E the
# recipient-year's adaptation-marked commitments, S its total commitments.
# Undefined (NA) when S = 0. Needs crs_positive() (functions/crs_positive.R).
# ==============================================================================

#' Within-country adaptation share, in percent
#'
#' @param commitments adaptation-marked commitments E (USD millions)
#' @param commitments_all total commitments S (USD millions)
#' @return numeric vector, 100 * E / S, NA where S is zero
within_share <- function(commitments, commitments_all) {
  ifelse(crs_positive(commitments_all), 100 * commitments / commitments_all, NA_real_)
}
