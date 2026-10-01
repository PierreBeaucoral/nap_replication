# ==============================================================================
# code/functions/sup_t_crit.R
# Critical value of the simultaneous (uniform, sup-t) 95% band drawn on every
# Callaway-Sant'Anna event-study / cohort figure (stages 03, 04, 05, 07, 11).
# Sourced, never copied. Needs did (>= 2.5.0, via BMisc) and withr.
# ==============================================================================

#' Sup-t critical value from a multiplier bootstrap on stored influence functions
#'
#' Mirrors did 2.5.0's `mboot()` (the `crit.val` that `aggte(cband = TRUE)` /
#' `att_gt(cband = TRUE)` would use): Mammen multiplier draws on the n x K
#' influence-function matrix (`BMisc::multiplier_bootstrap`, scaled by
#' sqrt(n)), degenerate columns dropped, each column studentised by its
#' bootstrap IQR-based scale, then the 1 - alp quantile (type 1) of the
#' maximum |t| across the K columns. The fits are clustered on the unit id
#' only (`clustervars = NULL`), so, as in `mboot()`, no extra cluster sum is
#' taken. As in did, a value below the pointwise normal quantile falls back to
#' the pointwise one.
#'
#' The band is then estimate +/- crit x the SE ALREADY STORED (the one the
#' tables report); nothing is re-aggregated, because `aggte()` re-runs the
#' bootstrap and would move the SEs.
#'
#' Simultaneity is per plotted curve: the caller passes the columns of ONE
#' outcome's (or one cohort panel's) plotted coefficients, so the band covers
#' that curve jointly, not every curve in the figure jointly.
#'
#' RNG isolation: the draws run inside `withr::with_seed()`, which restores the
#' caller's `.Random.seed` on exit, so calling this function moves no later
#' random draw in the calling script.
#'
#' @param inf_func n x K influence-function matrix of the plotted coefficients
#'   (e.g. `agg$inf.function$dynamic.inf.func.e`; a sparse matrix is accepted)
#' @param se length-K vector of the stored SEs; columns with NA or ~0 SE (the
#'   normalised base period e = -1 under a universal base) are excluded
#' @param biters bootstrap replications (999, the pipeline's BITERS)
#' @param alp significance level (0.05 -> 95% band)
#' @param seed seed of the isolated draw
#' @return scalar critical value (>= qnorm(1 - alp / 2))
sup_t_crit <- function(inf_func, se, biters = 999L, alp = 0.05, seed = 1242L) {
  inf_func <- as.matrix(inf_func)
  stopifnot(ncol(inf_func) == length(se))
  keep     <- which(!is.na(se) & se > 1e-10)
  z        <- stats::qnorm(1 - alp / 2)
  if (length(keep) == 0L) return(z)
  inf_func <- inf_func[, keep, drop = FALSE]
  n        <- nrow(inf_func)

  bres <- withr::with_seed(seed,
    sqrt(n) * BMisc::multiplier_bootstrap(inf_func, biters))
  bres <- as.matrix(bres)
  ndg  <- !is.na(colSums(bres)) & colSums(bres^2) > sqrt(.Machine$double.eps) * 10
  bres <- bres[, ndg, drop = FALSE]
  if (ncol(bres) == 0L) return(z)

  b_sigma <- apply(bres, 2L, function(b) {
    b <- sort.int(b)
    nb <- length(b)
    (b[ceiling(0.75 * nb)] - b[ceiling(0.25 * nb)]) / (stats::qnorm(0.75) - stats::qnorm(0.25))
  })
  b_t  <- apply(abs(sweep(bres, 2L, b_sigma, "/")), 1L, max)
  crit <- unname(stats::quantile(b_t[is.finite(b_t)], 1 - alp, type = 1, na.rm = TRUE))
  if (!is.finite(crit) || crit < z) crit <- z  # did's fallback to pointwise
  crit
}
