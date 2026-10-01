# ==============================================================================
# code/functions/read_headline_fit.R
# Reader for the headline fits that 03_main_results.R stores in output/fits/.
# Stages 04-14 report the headline specification from these fits instead of
# re-estimating it, so each specification has one ATT and one SE in the paper.
# ==============================================================================

#' Read a stored headline fit, stopping if it is missing or incomplete
#'
#' @param stem one of "adaptation", "share", "total", "nonadaptation",
#'   "disbursements"
#' @return the stored list (fields written by store_headline_fit() in
#'   03_main_results.R)
read_headline_fit <- function(stem) {
  path <- here::here("output", "fits", paste0("headline_", stem, "_dr_bs.rds"))
  if (!file.exists(path)) {
    stop("Stored headline fit not found: ", path, "\n",
         "  It is written by 03_main_results.R (run_all.R runs 03 first).")
  }
  fit <- readRDS(path)
  # Every field a later stage reads.
  needed  <- c("outcome", "spec", "gt_boot", "gt_analytical", "agg_simple",
               "agg_simple_analytic", "agg_dynamic", "agg_group", "agg_dyn_analytic",
               "ids", "pretrend", "l_vec_simple", "n_obs", "n_country")
  missing <- needed[vapply(needed, function(k) is.null(fit[[k]]), logical(1L))]
  scalar_ok <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (length(missing) > 0L || !scalar_ok(fit$att) || !scalar_ok(fit$se)) {
    stop("Stored headline fit is incomplete: ", path,
         if (length(missing) > 0L) paste0(" (missing: ", paste(missing, collapse = ", "), ")"),
         " -- rerun 03_main_results.R.")
  }
  fit
}
