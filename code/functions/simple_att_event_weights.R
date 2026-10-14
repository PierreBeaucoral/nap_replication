# ==============================================================================
# code/functions/simple_att_event_weights.R
# Event-time weights that reproduce did's simple ATT (the HonestDiD l_vec).
# Used by 03 (headline fits) and 07 (within-country share fits).
# Requires dplyr.
# ==============================================================================

#' Event-time weights implied by did's `simple` aggregation.
#'
#' did's simple ATT is the group-size-weighted average of the post-treatment
#' ATT(g,t) cells, so it equals sum_e w_e * att.egt(e) with
#'   w_e proportional to the number of treated units in the cohorts that
#'   contribute a post-treatment cell at event time e.
#' The returned vector is the l_vec that makes a HonestDiD sensitivity analysis
#' target exactly the reported simple ATT. The reconciliation against
#' the reported simple ATT is checked by the caller, not assumed.
#'
#' @param gt_obj  an `MP` object returned by did::att_gt()
#' @param panel   the estimation data frame (needs cohort_year, country_id)
#' @param egt_post numeric vector of post-treatment event times, in the order
#'   used by the dynamic aggregation
#' @return numeric vector of weights, same length as egt_post, summing to 1
simple_att_event_weights <- function(gt_obj, panel, egt_post) {
  n_by_cohort <- panel %>%
    filter(cohort_year > 0) %>%
    distinct(country_id, cohort_year) %>%
    count(cohort_year, name = "n_units")

  cells <- data.frame(
    g   = as.numeric(gt_obj$group),
    t   = as.numeric(gt_obj$t),
    att = as.numeric(gt_obj$att)
  )
  cells <- cells[cells$t >= cells$g & !is.na(cells$att), , drop = FALSE]
  # Event times are integer year differences; coerce to integer so that the
  # match below is an exact integer comparison, never a float equality test.
  cells$e <- as.integer(round(cells$t - cells$g))
  cells <- merge(cells, n_by_cohort, by.x = "g", by.y = "cohort_year",
                 all.x = TRUE, sort = FALSE)

  egt_int <- as.integer(round(egt_post))
  w <- vapply(egt_int, function(e) sum(cells$n_units[cells$e == e], na.rm = TRUE),
              numeric(1L))
  if (sum(w) <= 0) return(rep(1 / length(egt_post), length(egt_post)))
  w / sum(w)
}
