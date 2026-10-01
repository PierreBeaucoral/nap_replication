# ==============================================================================
# code/functions/pretrend_test.R
# The single pre-trend test used by every stage (03, 04, 05, 07, 09, 11, 13).
# Sourced, never copied: a local copy is how the stages drifted apart before.
# Base R only.
# ==============================================================================

#' Joint Wald pre-trend test on the aggregated pre-treatment event-time ATTs
#'
#' Covariance: did (>= 2.5.0) stores the dynamic-aggregation influence function
#' at `agg_d$inf.function$dynamic.inf.func.e` (n x K, aligned with
#' `agg_d$egt`). The full covariance is `crossprod(IF) / n^2`, and the test uses
#' the whole pre-treatment block, not a diagonal approximation. Call it on an
#' ANALYTICAL fit (bstrap = FALSE): the guard below checks the covariance
#' diagonal against `agg_d$se.egt`, which a bootstrap aggregation does not
#' reproduce.
#'
#' Inverse and degrees of freedom come from ONE singular value decomposition
#' of the pre-treatment block with ONE tolerance (max(dim) * eps * d_1). The
#' numerical rank is the number of singular values above it. At full rank the
#' block is inverted with `solve()`. Below full rank, the Moore-Penrose inverse
#' is built from the same SVD, dropping exactly the directions not counted in
#' the rank, with a message; df is then the rank, because such an inverse tests
#' only the identified directions (the full lead count would overstate the df
#' and the p-value: a conservative, not a valid, test).
#'
#' Anticipation: under `att_gt(anticipation = k)` the cells e = -1, ..., -k are
#' treated by construction and e = -k-1 is the base period, so the test is
#' restricted to e < -k. (Truncating the aggregation instead, with
#' `aggte(max_e = -k-1)`, is not an option: did errors when a dynamic
#' aggregation has no post-treatment period.)
#'
#' did's own group-time pre-test (`gt_obj$W`, `gt_obj$Wpval`, computed on the
#' disaggregated ATT(g,t) cells) is returned alongside, with its restriction
#' count and, when did returns no statistic, the reason derived from the fit.
#'
#' @param agg_d dynamic AGGTEobj from `did::aggte(type = "dynamic")`; NULL (a
#'   failed aggregation) stops, so no exhibit prints an empty pre-trend cell
#' @param gt_obj the analytical `MP` object from `did::att_gt()`, or NULL
#' @param anticipation integer k of the fit (default 0)
#' @return list(stat, pval, df, n_leads, leads, ginv_used, W_did, Wpval_did,
#'   df_did, wpval_reason); statistics and p-values unrounded (callers round
#'   only for display)
compute_pretrend_test <- function(agg_d, gt_obj = NULL, anticipation = 0) {
  empty <- list(stat = NA_real_, pval = NA_real_, df = 0L, n_leads = 0L,
                leads = integer(0), ginv_used = FALSE,
                W_did = NA_real_, Wpval_did = NA_real_, df_did = NA_integer_,
                wpval_reason = NA_character_)
  if (is.null(agg_d)) stop("compute_pretrend_test(): no dynamic aggregation (the analytical ",
                           "fit or its aggregation failed)")

  keep    <- which(!is.na(agg_d$se.egt) & agg_d$se.egt > 1e-10)
  pre_pos <- which(agg_d$egt[keep] < -anticipation)
  if (length(pre_pos) == 0L) return(empty)

  pre_beta <- agg_d$att.egt[keep][pre_pos]

  IF <- agg_d$inf.function$dynamic.inf.func.e
  stopifnot(is.matrix(IF), ncol(IF) == length(agg_d$egt))
  n          <- nrow(IF)
  sigma_full <- crossprod(IF) / n^2
  # Guard: catches IF/egt column misalignment (it would silently corrupt every
  # pre-trend test) by checking the covariance diagonal against did's se.egt,
  # relative to each SE (the floor covers the zero-SE base period).
  rel_gap <- abs(sqrt(diag(sigma_full)) - agg_d$se.egt) / pmax(abs(agg_d$se.egt), 1e-10)
  stopifnot(max(rel_gap, na.rm = TRUE) < 1e-6)
  sigma_pre  <- sigma_full[keep, keep][pre_pos, pre_pos, drop = FALSE]

  sv        <- svd(sigma_pre)
  pos       <- sv$d > max(dim(sigma_pre)) * .Machine$double.eps * max(sv$d)
  df_use    <- sum(pos)
  ginv_used <- df_use < length(pre_pos)
  inv <- if (!ginv_used) solve(sigma_pre) else {
    message(sprintf(paste0("compute_pretrend_test: pre-treatment covariance is singular ",
                           "(numerical rank %d of %d leads) -- Moore-Penrose inverse on ",
                           "the identified directions, df = rank"),
                    df_use, length(pre_pos)))
    sv$v[, pos, drop = FALSE] %*% (t(sv$u[, pos, drop = FALSE]) / sv$d[pos])
  }
  W <- if (df_use == 0L) NA_real_ else as.numeric(t(pre_beta) %*% inv %*% pre_beta)

  W_did     <- NA_real_
  Wpval_did <- NA_real_
  df_did    <- NA_integer_
  if (!is.null(gt_obj)) {
    if (!is.null(gt_obj$W))     W_did     <- as.numeric(gt_obj$W)
    if (!is.null(gt_obj$Wpval)) Wpval_did <- as.numeric(gt_obj$Wpval)
    # did 2.5.0 (att_gt) drops the pre-treatment cells whose SE is NA (the
    # base-period cells) before inverting V[pre, pre], so q counts only cells
    # with a usable SE.
    df_did <- sum(gt_obj$t < gt_obj$group & !is.na(gt_obj$se))
  }

  # Why did returned no statistic, derived from the fit. did declines the
  # pre-test when (i) it never formed the analytical variance matrix, (ii)
  # there are no estimable pre-treatment cells, or (iii) rcond(preV)
  # underflows (thin subgroups: more pre-treatment cells than the subgroup's
  # recipient-level influence functions can span).
  wpval_reason <- NA_character_
  if (!is.null(gt_obj) && is.na(Wpval_did)) {
    pre_idx <- which(gt_obj$group > gt_obj$t & !is.na(gt_obj$se))  # the cells did inverts
    if (is.null(gt_obj$V_analytical)) {
      wpval_reason <- paste0("did does not form the analytical variance matrix ",
                             "for this fit, so its pre-test is unavailable by ",
                             "construction")
    } else if (length(pre_idx) == 0L) {
      wpval_reason <- "there are no estimable pre-treatment group-time cells"
    } else {
      preV <- as.matrix(gt_obj$V_analytical[pre_idx, pre_idx])
      rk   <- tryCatch({ sv <- svd(preV)$d
                         sum(sv > max(dim(preV)) * .Machine$double.eps * max(sv)) },
                       error = function(e) NA_integer_)
      wpval_reason <- sprintf(paste0("not computed: pre-treatment covariance ",
                                     "singular (%d cells, numerical rank %s)"),
                              nrow(preV),
                              if (is.na(rk)) "unavailable" else as.character(rk))
    }
  }

  list(
    stat         = W,
    pval         = pchisq(W, df = df_use, lower.tail = FALSE),
    df           = df_use,
    n_leads      = length(pre_pos),
    leads        = as.integer(agg_d$egt[keep][pre_pos]),
    ginv_used    = ginv_used,
    W_did        = W_did,
    Wpval_did    = Wpval_did,
    df_did       = as.integer(df_did),
    wpval_reason = wpval_reason
  )
}
