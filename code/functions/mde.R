# ==============================================================================
# code/functions/mde.R
# The single minimum-detectable-effect formula used by stages 03, 05, 13 and 14.
# ==============================================================================

#' Minimum detectable effect of a two-sided test
#'
#' (z_{1-alpha/2} + z_{power}) x SE: 2.8016 x SE at the defaults.
#'
#' @param se standard error(s) of the estimate
#' @param alpha two-sided significance level, in (0, 1)
#' @param power target power, in (0, 1)
#' @return the MDE, same length as se
mde <- function(se, alpha = 0.05, power = 0.8) {
  stopifnot(alpha > 0, alpha < 1, power > 0, power < 1)
  (qnorm(1 - alpha / 2) + qnorm(power)) * se
}
