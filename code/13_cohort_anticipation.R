# ==============================================================================
# 13_cohort_anticipation.R
# Cohort-robustness battery, anticipation / event-date / placebo ladder
#, and conditioning-set sensitivity.
# Paper: Beaucoral, Goujon and Marchand (2026) — appendix to §5 (robustness)
#
# Three blocks of checks:
#   A   Cohort-robustness battery: leave-one-cohort-out over every estimation
#       cohort; drop 2024 as a cohort and as a calendar year; balanced event
#       windows (balance_e = 1, 2); the 2x2 {main, retained} x {dr, reg} grid
#       with a common SE method; and which cohorts contribute to each event
#       time.
#   B   Anticipation (1 and 2 years), an alternative event date built from the
#       NAP submission MONTH, and a placebo ladder at -3/-4/-5 years with the
#       minimum detectable effect of each placebo.
#   C   Conditioning set: baseline (2009-2012) outcome level in xformla, and a
#       common-support specification trimming recipients whose estimated
#       propensity to adopt lies outside [0.05, 0.95]; listwise-deletion losses
#       for each specification.
#
# Inputs :
#   data/processed/simple_panel_wgi.csv
#   output/fits/headline_adaptation_dr_bs.rds   (written by 03_main_results.R)
#   output/fits/headline_share_dr_bs.rds        (written by 03_main_results.R)
#   output/fits/placebo_m2_reg_bs.rds           (written by 04_robustness.R)
#   data/raw/shared_nap_data/nap_information.csv (submission dates; also carried
#                                                 in the processed panel)
# Outputs:
#   output/tables/cohort_battery/att_loco.tex
#   output/tables/cohort_battery/att_drop2024.tex
#   output/tables/cohort_battery/att_balance.tex
#   output/tables/cohort_battery/att_2x2.tex
#   output/tables/cohort_battery/cohort_contributions.tex
#   output/tables/cohort_battery/att_cohort_battery.tex   (combined summary)
#   output/tables/cohort_battery/att_conditioning.tex
#   output/tables/cohort_battery/listwise_losses_13.tex
#   output/tables/anticipation/att_anticipation.tex
#   output/tables/anticipation/att_eventdate.tex
#   output/tables/anticipation/att_placebo_ladder.tex
#
# This script NEVER touches the raw CRS tree; it reads the processed panel only.
# ==============================================================================

# ============================================================
# Paper-to-Code Naming Map
# ============================================================
# Paper Notation        | Code Name              | Description
# $Y_{it}$              | log_commits            | log(1 + adaptation commitments)
# $s_{it}$              | share_adapt            | share of global adaptation finance
# $G_i$                 | cohort_year            | NAP adoption cohort (0 = never)
# $G_i^{month}$         | cohort_year_month      | B event date: year + 1{month > 6}
# $ATT(g,t)$            | gt_obj                 | group-time ATT from att_gt()
# $\hat\theta^{simp}$   | agg_s$overall.att      | simple (group-size-weighted) ATT
# $\hat\theta^{dyn}(e)$ | agg_d$att.egt          | dynamic ATT at event time e
# $X_{it}$              | ge_est, log_population | controls
# $\bar Y_i^{0912}$     | base_outcome           | C baseline outcome level
# $\hat p_i$            | pscore                 | C adoption propensity
# ============================================================

library(data.table)
library(dplyr)
library(xtable)
library(here)
library(did)
# Version guard: SEs on unbalanced panels changed in did 2.5.0 (renv.lock pins it).
# An older install reproduces the ATTs but not the SEs / pre-trend tests.
if (utils::packageVersion("did") < "2.5.0") {
  stop(sprintf(paste0(
    "did %s is installed but the paper's results require did >= 2.5.0.\n",
    "  Older versions reproduce the ATTs but NOT the standard errors and\n",
    "  pre-trend tests on unbalanced panels (South Sudan has 14/16 years).\n",
    "  Fix: run renv::restore() from the project root (renv.lock pins did 2.5.0),\n",
    "  or install.packages(\"did\") to get >= 2.5.0, then rerun this script."),
    utils::packageVersion("did")))
}
library(parallel)   # forked evaluation guard for the segfault-prone dr cell
# NOTE: MASS is NOT attached via library() -- MASS::select() would mask
# dplyr::select(). compute_pretrend_test() calls MASS::ginv() by full namespace.

set.seed(20240601)  # global seed — local set.seed(1242) calls follow each estimator

# One bootstrap-replication constant per script.
BITERS <- 999L

# Cohort inclusion rule of the main specification (03_main_results.R).
THIN_THRESHOLD <- 5L

FITS_DIR <- here("output", "fits")

dir_tabs_cb <- here("output", "tables", "cohort_battery")
dir_tabs_an <- here("output", "tables", "anticipation")
dir.create(dir_tabs_cb, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tabs_an, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

esc_header <- function(x) {
  x <- gsub("%", "\\\\%", x)
  x <- gsub("_", "\\\\_", x)
  x <- gsub("#", "\\\\#", x)
  x
}

# ==============================================================================
# Helper: write_tex_float()
# (Identical copy to 03/04/05 — kept local per this project's no-sourced-utils
#  convention. Emits \begin{table}[H] + \adjustbox + notes minipage.)
# ==============================================================================
write_tex_float <- function(out_path, caption_title, label,
                             tabular_lines, notes_text, source_text,
                             size = "\\small") {
  inner <- tabular_lines
  is_table_open  <- grepl("^\\\\begin\\{table\\}", inner)
  is_table_close <- grepl("^\\\\end\\{table\\}", inner)
  if (any(is_table_open)) inner <- inner[!is_table_open]
  if (any(is_table_close)) inner <- inner[!is_table_close]
  inner <- inner[!grepl("^\\\\centering", inner)]
  inner <- inner[!grepl("^\\\\caption", inner)]
  inner <- inner[!grepl("^\\\\label", inner)]
  while (length(inner) > 0 && trimws(inner[1]) == "") inner <- inner[-1]
  while (length(inner) > 0 && trimws(inner[length(inner)]) == "")
    inner <- inner[-length(inner)]

  # Wrap EVERY tabular block (not only the first) in \adjustbox, so that
  # multi-panel floats (e.g. att_anticipation: pooled + dynamic panels) keep
  # each panel within \textwidth. Walk from the last block backwards so the
  # insertions do not shift the indices of blocks still to be wrapped.
  tab_starts <- which(grepl("^\\\\begin\\{tabular", inner))
  tab_ends   <- which(grepl("^\\\\end\\{tabular",   inner))
  if (length(tab_starts) > 0 && length(tab_starts) == length(tab_ends)) {
    for (b in rev(seq_along(tab_starts))) {
      tab_start <- tab_starts[b]
      tab_end   <- tab_ends[b]
      inner <- c(
        if (tab_start > 1) inner[seq_len(tab_start - 1)] else character(0),
        "\\adjustbox{max width=\\textwidth}{%",
        inner[tab_start:tab_end],
        "}",
        if (tab_end < length(inner)) inner[seq(tab_end + 1, length(inner))] else character(0)
      )
    }
  }

  lines_out <- c(
    "\\begin{table}[H]",
    "\\centering",
    paste0("\\caption{", caption_title, "}"),
    paste0("\\label{", label, "}"),
    inner,
    "\\par\\vspace{4pt}",
    "\\begin{minipage}{\\linewidth}\\footnotesize",
    paste0("Notes: ", notes_text, ".\\par"),
    paste0("Source: ", source_text, "."),
    "\\end{minipage}",
    "\\end{table}"
  )

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines_out, out_path)
  message("Saved: ", out_path)
  invisible(out_path)
}

# ==============================================================================
# Helper: compute_pretrend_test()
# (Byte-identical copy of the helper in 03/04/05/10/11 — joint Wald test on the
#  pre-treatment block of the dynamic aggregation's influence-function
#  covariance. Always called on an ANALYTICAL fit.)
# ==============================================================================
compute_pretrend_test <- function(agg_d, gt_obj = NULL, anticipation = 0) {
  # `anticipation` (default 0, i.e. the behaviour every other call site relies
  # on) excludes the event times the estimator itself treats as TREATED. Under
  # att_gt(anticipation = k) the cells e = -1, ..., -k are post-treatment by
  # construction and e = -k-1 is the normalised base period, so a pre-trend
  # test must be restricted to e < -k. Truncating the aggregation instead
  # (aggte(max_e = -k-1)) is not an option: did errors out when a dynamic
  # aggregation contains no post-treatment period.
  # Covariance: did (>= 2.5.0) stores the dynamic-aggregation influence
  # function at agg_d$inf.function$dynamic.inf.func.e (n x K, aligned with
  # agg_d$egt). The full covariance is crossprod(IF)/n^2, and the joint Wald
  # test uses the whole pre-treatment block, not a diagonal approximation
  # (which would ignore covariance across pre-treatment event-time ATTs). A
  # generalized inverse (MASS::ginv) is used, with a message, if the block is
  # singular. did's own group-time pre-test (gt_obj$W,
  # gt_obj$Wpval — computed on the disaggregated ATT(g,t) cells, independent of
  # bstrap) is captured alongside for comparison and reported as a separate
  # table row ("Pre-trend p (did Wpval)").
  keep    <- which(!is.na(agg_d$se.egt) & agg_d$se.egt > 1e-10)
  pre_pos <- which(agg_d$egt[keep] < -anticipation)
  if (length(pre_pos) == 0)
    return(list(stat = NA_real_, pval = NA_real_, df = 0L, n_leads = 0L,
                leads = integer(0), ginv_used = FALSE,
                W_did = NA_real_, Wpval_did = NA_real_, df_did = NA_integer_,
                wpval_reason = NA_character_))

  pre_beta <- agg_d$att.egt[keep][pre_pos]

  IF <- agg_d$inf.function$dynamic.inf.func.e
  stopifnot(is.matrix(IF), ncol(IF) == length(agg_d$egt))
  n          <- nrow(IF)
  sigma_full <- crossprod(IF) / n^2
  # Guard: catches IF/egt column misalignment (would silently corrupt every
  # downstream pre-trend test) by cross-checking the covariance diagonal
  # against did's own agg_d$se.egt.
  stopifnot(max(abs(sqrt(diag(sigma_full)) - agg_d$se.egt), na.rm = TRUE) < 1e-6)
  sigma_pre  <- sigma_full[keep, keep][pre_pos, pre_pos, drop = FALSE]

  # A generalized inverse is used only if the direct solve fails. When it is,
  # the statistic no longer has length(pre_pos) degrees of freedom: the correct
  # reference is the NUMERICAL RANK of the pre-treatment covariance, because a
  # Moore-Penrose inverse tests only the directions the data can identify.
  # Using the full lead count there would overstate df and understate the
  # p-value. Both the rank and the fact that a generalized inverse was used are
  # returned, so the table note can disclose them.
  # tryCatch returns BOTH the inverse and the flag, so nothing is assigned into
  # an enclosing environment (`<<-` is prohibited by this project's standards).
  inv_res <- tryCatch(
    list(inv = solve(sigma_pre), ginv = FALSE),
    error = function(e) {
      message("compute_pretrend_test: pre-treatment covariance is singular — ",
              "using MASS::ginv() generalized inverse instead of a direct solve()")
      list(inv = MASS::ginv(sigma_pre), ginv = TRUE)
    }
  )
  inv_sigma_pre <- inv_res$inv
  ginv_used     <- isTRUE(inv_res$ginv)
  W <- as.numeric(t(pre_beta) %*% inv_sigma_pre %*% pre_beta)

  df_use <- length(pre_pos)
  if (ginv_used) {
    sv <- svd(sigma_pre)$d
    df_use <- sum(sv > max(dim(sigma_pre)) * .Machine$double.eps * max(sv))
    message(sprintf(paste0("compute_pretrend_test: df set to the numerical rank ",
                           "of the pre-treatment covariance (%d of %d leads)"),
                    df_use, length(pre_pos)))
  }

  W_did     <- NA_real_
  Wpval_did <- NA_real_
  df_did    <- NA_integer_
  if (!is.null(gt_obj)) {
    if (!is.null(gt_obj$W))     W_did     <- as.numeric(gt_obj$W)
    if (!is.null(gt_obj$Wpval)) Wpval_did <- as.numeric(gt_obj$Wpval)
    # Restrictions behind did's own pre-test. This is did's OWN q: it takes
    # pre <- which(group > t) and inverts the full V[pre, pre] block, without
    # dropping cells whose standard error it has already set to NA. Using a
    # filtered count here would put two different numbers in one sentence of
    # the table note (54 cells inverted vs 50 with a usable SE), so the count
    # reported is the one did actually tests.
    df_did <- sum(gt_obj$t < gt_obj$group)
  }

  # Why did returned no statistic, derived from the fit rather than guessed.
  # did declines the pre-test when (i) it never formed the analytical variance
  # matrix (bstrap without cband), (ii) there are no estimable pre-treatment
  # cells, or (iii) rcond(preV) underflows. Case (iii) is the one that fires on
  # thin subgroups: the pre-treatment block has more cells than the subgroup's
  # recipient-level influence functions can span.
  wpval_reason <- NA_character_
  if (!is.null(gt_obj) && is.na(Wpval_did)) {
    pre_idx <- which(gt_obj$group > gt_obj$t)
    if (is.null(gt_obj$V)) {
      wpval_reason <- paste0("did does not form the analytical variance matrix ",
                             "for this fit, so its pre-test is unavailable by ",
                             "construction")
    } else if (length(pre_idx) == 0L) {
      wpval_reason <- "there are no estimable pre-treatment group-time cells"
    } else {
      preV <- as.matrix(gt_obj$V[pre_idx, pre_idx])
      rc   <- tryCatch(rcond(preV), error = function(e) NA_real_)
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
    stat         = round(W, 3),
    pval         = round(pchisq(W, df = df_use, lower.tail = FALSE), 3),
    df           = df_use,
    n_leads      = length(pre_pos),
    leads        = as.integer(agg_d$egt[keep][pre_pos]),
    ginv_used    = ginv_used,
    W_did        = if (is.na(W_did))     NA_real_ else round(W_did, 3),
    Wpval_did    = if (is.na(Wpval_did)) NA_real_ else round(Wpval_did, 3),
    df_did       = as.integer(df_did),
    wpval_reason = wpval_reason
  )
}

# canonical pre-trend wording — keep byte-identical across scripts
PRETREND_NOTE_AGG <- function(min_e, max_e, k) sprintf(
  "Pre-trend $\\chi^2$ / $p$: joint Wald test on the aggregated pre-treatment event-time coefficients ($%d \\leq e \\leq %d$; %d restrictions), using the influence-function covariance of the dynamic aggregation from the analytical (non-bootstrap) fit; a generalized inverse is used if the block is singular",
  min_e, max_e, k)
PRETREND_NOTE_DID <- "\\texttt{did} pre-test $p$: \\texttt{did}'s built-in Wald test over all pre-period $ATT(g,t)$ cells against each cohort's $g-1$ base year"
PRETREND_NOTE <- function(min_e, max_e, k) paste0(PRETREND_NOTE_AGG(min_e, max_e, k), ". ", PRETREND_NOTE_DID)

#' Minimum detectable effect at 80% power, 5% two-sided.
#' MDE = (z_{0.975} + z_{0.80}) * se = 2.8016 * se.
mde_from_se <- function(se, power = 0.80, alpha = 0.05) {
  if (is.na(se)) return(NA_real_)
  (qnorm(1 - alpha / 2) + qnorm(power)) * se
}

fmt_stars <- function(t_v) {
  if (is.na(t_v)) return("")
  if (abs(t_v) > 2.576) "***" else if (abs(t_v) > 1.960) "**" else
    if (abs(t_v) > 1.645) "*" else ""
}

# ==============================================================================
# SECTION 1. Panel construction — byte-identical to 03_main_results.R §1
# and a HARD CHECK that it reproduces the estimation sample behind Table 2.
# ==============================================================================

message("\n=== 13_cohort_anticipation.R: loading panel ===\n")

aggregated <- as.data.frame(fread(here("data", "processed", "simple_panel_wgi.csv")))

did_panel <- aggregated %>%
  mutate(country_id = as.integer(factor(recipient_name)))

first_year <- min(did_panel$year)
last_year  <- max(did_panel$year, na.rm = TRUE)

country_gname <- did_panel %>%
  group_by(recipient_name) %>%
  summarise(nap_year_c = suppressWarnings(min(nap_year, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(
    nap_year_c     = if_else(is.infinite(nap_year_c), NA_real_, nap_year_c),
    always_treated = !is.na(nap_year_c) & nap_year_c < first_year,
    cohort_year = case_when(
      is.na(nap_year_c)          ~ 0,
      nap_year_c < first_year    ~ 0,
      nap_year_c > last_year     ~ 0,
      TRUE                       ~ as.numeric(nap_year_c)
    )
  )

did_panel <- did_panel %>%
  select(-any_of(c("cohort_year", "always_treated"))) %>%
  left_join(country_gname %>% select(recipient_name, cohort_year, always_treated),
            by = "recipient_name") %>%
  mutate(log_population = log(population))

if (!"share_adapt" %in% names(did_panel)) {
  share_lookup <- aggregated %>% select(recipient_name, year, share_adapt) %>% distinct()
  did_panel <- left_join(did_panel, share_lookup, by = c("recipient_name", "year"))
}

did_panel_full <- did_panel

#' Cohort sizes (treated recipients per adoption year).
cohort_sizes_of <- function(panel) {
  panel %>%
    filter(cohort_year > 0) %>%
    distinct(recipient_name, cohort_year) %>%
    count(cohort_year, name = "n_treated") %>%
    arrange(cohort_year)
}

#' The main specification's cohort filter, exposed as a helper.
#'
#' 03_main_results.R drops adoption cohorts with fewer than `threshold` treated
#' recipients. This reproduces that rule exactly (same threshold, same
#' definition of a cohort's size) without touching 03.
#'
#' @param panel recipient-year panel with cohort_year
#' @param threshold minimum number of treated recipients per cohort
#' @return list(panel, thin_cohorts, kept_cohorts, sizes)
apply_cohort_filter <- function(panel, threshold = THIN_THRESHOLD) {
  sizes <- cohort_sizes_of(panel)
  thin  <- sizes$cohort_year[sizes$n_treated <  threshold]
  kept  <- sizes$cohort_year[sizes$n_treated >= threshold]
  list(panel        = panel %>% filter(!(cohort_year %in% thin)),
       thin_cohorts = thin,
       kept_cohorts = kept,
       sizes        = sizes)
}

cf <- apply_cohort_filter(did_panel_full)
did_panel_main <- cf$panel
thin_cohorts   <- cf$thin_cohorts
kept_cohorts   <- cf$kept_cohorts

message("Cohort sizes:")
print(cf$sizes, row.names = FALSE)
message(sprintf("Estimation cohorts (>= %d treated units): %s",
                THIN_THRESHOLD, paste(kept_cohorts, collapse = ", ")))
message(sprintf("Dropped thin cohorts: %s", paste(thin_cohorts, collapse = ", ")))
message(sprintf("Main panel: %d recipient-years, %d recipients, %d treated recipients",
                nrow(did_panel_main), n_distinct(did_panel_main$country_id),
                n_distinct(did_panel_main$country_id[did_panel_main$cohort_year > 0])))

# HARD CHECK: this panel must be the one behind Table 2. Reading the saved
# headline fit makes the comparison exact rather than a matter of comment.
path_head <- file.path(FITS_DIR, "headline_adaptation_dr_bs.rds")
if (!file.exists(path_head)) {
  stop("Saved headline fit not found: ", path_head, "\n",
       "  13_cohort_anticipation.R compares its panel against the fit behind ",
       "Table 2.\n  Fix: run `Rscript code/03_main_results.R` first.")
}
fit_head <- readRDS(path_head)
stopifnot(
  identical(as.integer(fit_head$n_obs),
            as.integer(sum(!is.na(did_panel_main$log_commits)))),
  identical(as.integer(fit_head$n_country),
            as.integer(n_distinct(did_panel_main$country_id[
              !is.na(did_panel_main$log_commits)]))),
  identical(sort(unique(did_panel_main$country_id)), fit_head$ids)
)
message(sprintf(paste0("Panel check OK: matches the stored headline fit ",
                       "(N = %d, recipients = %d); reference ATT = %.4f (SE %.4f)"),
                fit_head$n_obs, fit_head$n_country, fit_head$att, fit_head$se))

path_share <- file.path(FITS_DIR, "headline_share_dr_bs.rds")
fit_share  <- if (file.exists(path_share)) readRDS(path_share) else NULL

OUTCOMES13 <- list(
  list(var = "log_commits", label = "log(Adaptation commitments)",
       ref = fit_head),
  list(var = "share_adapt", label = "Adaptation share (\\% of global)",
       ref = fit_share)
)

# ==============================================================================
# SECTION 2. Estimation helpers
# Every fit in this script goes through run_spec(), so the specification is
# identical to the headline except for the one dimension being varied, and the
# seeding discipline (set.seed(1242) immediately before each estimator call) is
# enforced in one place.
# ==============================================================================

#' Fit the main specification on an arbitrary panel and aggregate it.
#'
#' @param panel estimation panel
#' @param yname outcome column
#' @param est_method "dr" or "reg"
#' @param anticipation anticipation periods passed to att_gt()
#' @param xformla one-sided formula of controls
#' @param balance_e balance_e passed to aggte(type = "dynamic"); NULL for none
#' @param min_e,max_e event-time window of the dynamic aggregation
#' @param gname cohort column name
#' @return list(att, se, t, stars, pretrend, agg_s, agg_d, n_obs, n_country,
#'   n_treated, gt) or NULL if estimation failed
run_spec <- function(panel, yname = "log_commits", est_method = "dr",
                     anticipation = 0, xformla = ~ ge_est + log_population,
                     balance_e = NULL, min_e = -5, max_e = Inf,
                     gname = "cohort_year", control_group = "nevertreated") {

  fit_one <- function(bstrap_flag) {
    att_gt(
      yname         = yname,
      tname         = "year",
      idname        = "country_id",
      gname         = gname,
      xformla       = xformla,
      data          = panel,
      est_method    = est_method,
      bstrap        = bstrap_flag,
      biters        = BITERS,
      cband         = FALSE,
      control_group = control_group,
      anticipation  = anticipation,
      base_period   = "universal",
      panel         = TRUE,
      allow_unbalanced_panel = TRUE
    )
  }

  # Analytical twin first (consumes no random numbers) — pre-trend test only.
  gt_an <- tryCatch(fit_one(FALSE),
                    error = function(e) { message("    att_gt (analytical) failed: ",
                                                   conditionMessage(e)); NULL })
  set.seed(1242)   # Seed rule (reproduces the published SEs): seed immediately before the estimator
  gt <- tryCatch(fit_one(TRUE),
                 error = function(e) { message("    att_gt (bootstrap) failed: ",
                                                conditionMessage(e)); NULL })
  if (is.null(gt)) return(NULL)

  agg_s <- tryCatch(aggte(gt, type = "simple", na.rm = TRUE), error = function(e) NULL)
  agg_d <- tryCatch(
    if (is.null(balance_e))
      aggte(gt, type = "dynamic", na.rm = TRUE, min_e = min_e, max_e = max_e)
    else
      aggte(gt, type = "dynamic", na.rm = TRUE, min_e = min_e, max_e = max_e,
            balance_e = balance_e),
    error = function(e) { message("    aggte(dynamic) failed: ",
                                   conditionMessage(e)); NULL })
  # Pre-trend block. With anticipation = k, event times -1, ..., -k are TREATED
  # cells, not leads: including them in a "pre-trend" test would test the
  # anticipation effect itself. The restriction is applied INSIDE
  # compute_pretrend_test() via its `anticipation` argument, not by truncating
  # the aggregation: aggte(max_e = -k-1) errors out inside did because the
  # resulting aggregation has no post-treatment period ("influence function
  # column indices exceed available columns").
  agg_d_an <- if (!is.null(gt_an)) tryCatch(
    aggte(gt_an, type = "dynamic", na.rm = TRUE, min_e = min_e, max_e = max_e),
    error = function(e) NULL) else NULL

  att <- if (!is.null(agg_s)) agg_s$overall.att else NA_real_
  se  <- if (!is.null(agg_s)) agg_s$overall.se  else NA_real_
  t_v <- if (!is.na(att) && !is.na(se) && se > 0) att / se else NA_real_
  pt  <- if (!is.null(agg_d_an))
    tryCatch(compute_pretrend_test(agg_d_an, gt_an, anticipation = anticipation),
             error = function(e) list(stat = NA_real_, pval = NA_real_, df = 0L,
                                       W_did = NA_real_, Wpval_did = NA_real_,
                                       df_did = NA_integer_))
  else list(stat = NA_real_, pval = NA_real_, df = 0L,
            W_did = NA_real_, Wpval_did = NA_real_, df_did = NA_integer_)

  # Event times that actually entered the pre-trend test, so the exhibit can
  # state the window rather than assert a formula.
  pre_leads <- if (is.null(agg_d_an)) integer(0) else {
    k_keep <- which(!is.na(agg_d_an$se.egt) & agg_d_an$se.egt > 1e-10)
    as.integer(agg_d_an$egt[k_keep][agg_d_an$egt[k_keep] < -anticipation])
  }

  # Post-submission-only average. did's SIMPLE aggregate averages every treated
  # (g, t) cell, so under anticipation = k it MECHANICALLY includes the k cells
  # before submission. This aggregation restricts attention to e >= 0 and is
  # the quantity comparable to the headline. It is computed last so that it
  # cannot disturb the bootstrap draws behind the reported ATT above; the seed
  # is set immediately before it, as for every other estimator call here.
  set.seed(1242)
  agg_e0 <- tryCatch(
    aggte(gt, type = "dynamic", na.rm = TRUE, min_e = 0, max_e = max_e),
    error = function(e) NULL)
  att_e0 <- if (!is.null(agg_e0)) agg_e0$overall.att else NA_real_
  se_e0  <- if (!is.null(agg_e0)) agg_e0$overall.se  else NA_real_

  rows_y <- panel[!is.na(panel[[yname]]), , drop = FALSE]
  list(att = att, se = se, t = t_v, stars = fmt_stars(t_v), pretrend = pt,
       agg_s = agg_s, agg_d = agg_d, agg_e0 = agg_e0, gt = gt,
       att_e0 = att_e0, se_e0 = se_e0, pre_leads = pre_leads,
       n_obs = nrow(rows_y), n_country = n_distinct(rows_y$country_id),
       n_treated = n_distinct(panel$country_id[panel$cohort_year > 0]),
       est_method = est_method, anticipation = anticipation)
}

#' Run run_spec() in a forked child so that a segfault cannot kill the script.
#'
#' att_gt(est_method = "dr") is known to segfault inside fastglm's colMax_dense
#' on thin / reduced unbalanced panels. A segfault cannot be
#' caught by tryCatch, so the risky cell of the 2x2 grid is evaluated in a
#' forked process: the child dies, mccollect() returns NULL, and the script
#' records the failure instead of aborting. The child re-seeds explicitly, so
#' the result is identical to an in-process run.
#' Falls back to an in-process call where forking is unavailable (Windows).
run_spec_forked <- function(...) {
  if (.Platform$OS.type != "unix") {
    message("    (forking unavailable on this platform — running in-process)")
    return(run_spec(...))
  }
  # Reproducibility guard: the child inherits the parent's RNG configuration
  # and then re-seeds inside run_spec(). If the generator were ever switched
  # (e.g. to L'Ecuyer-CMRG for parallel work) the forked result would not match
  # an in-process run, so the configuration is asserted rather than assumed.
  stopifnot(identical(RNGkind()[1L], "Mersenne-Twister"))
  args <- list(...)
  job <- parallel::mcparallel(do.call(run_spec, args))
  res <- parallel::mccollect(job, wait = TRUE)
  out <- res[[1L]]
  if (is.null(out) || inherits(out, "try-error")) {
    message("    forked estimation returned no result (crash or error).")
    return(NULL)
  }
  out
}

#' Estimable group-time cells and contributing treated units, per cohort.
#'
#' A cohort listed in an estimation panel does not necessarily contribute to the
#' reported ATT: `aggte(na.rm = TRUE)` silently drops cohorts whose group-time
#' cells are all missing or have a zero standard error, which is exactly what
#' happens to singleton cohorts. Reporting the INPUT treated count would then
#' overstate what the estimate rests on, so this returns both.
#'
#' @param gt an `MP` object from did::att_gt()
#' @param panel the estimation panel that produced it
#' @return data.frame(cohort, treated_in_panel, cells_pre, cells_post,
#'   contributes)
cohort_cell_report <- function(gt, panel) {
  if (is.null(gt)) return(NULL)
  cells <- data.frame(g  = as.numeric(gt$group),
                      tt = as.numeric(gt$t),
                      att = as.numeric(gt$att),
                      se  = as.numeric(gt$se))
  cells$ok <- !is.na(cells$att) & !is.na(cells$se) & cells$se > 1e-10
  n_by <- panel %>%
    filter(cohort_year > 0) %>%
    distinct(country_id, cohort_year) %>%
    count(cohort_year, name = "treated_in_panel")
  out <- lapply(seq_len(nrow(n_by)), function(i) {
    g <- n_by$cohort_year[i]
    sub <- cells[cells$g == g, , drop = FALSE]
    data.frame(cohort           = as.integer(g),
               treated_in_panel = n_by$treated_in_panel[i],
               cells_pre        = sum(sub$ok & sub$tt <  sub$g),
               cells_post       = sum(sub$ok & sub$tt >= sub$g),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  out$contributes <- out$cells_post > 0L
  out
}

#' One formatted row of a results table.
spec_row <- function(label, res) {
  if (is.null(res)) {
    return(data.frame(Specification = label, ATT = "---", SE = "---", `$t$` = "---",
                      `Pre-trend $p$` = "---", `$N$` = "---", Recipients = "---",
                      `Treated` = "---", check.names = FALSE,
                      stringsAsFactors = FALSE))
  }
  data.frame(
    Specification   = label,
    ATT             = paste0(sprintf("%.4f", res$att), res$stars),
    SE              = sprintf("(%.4f)", res$se),
    `$t$`           = sprintf("%.3f", res$t),
    `Pre-trend $p$` = if (is.na(res$pretrend$pval)) "---" else
      sprintf("%.3f", res$pretrend$pval),
    `$N$`           = format(res$n_obs, big.mark = ","),
    Recipients      = as.character(res$n_country),
    Treated         = as.character(res$n_treated),
    check.names = FALSE, stringsAsFactors = FALSE)
}

#' Print a table via xtable with the project's standard settings.
tabular_of <- function(df, label, align_str = NULL) {
  xt <- xtable(df, label = label)
  if (!is.null(align_str)) align(xt) <- align_str
  capture.output(print(xt, include.rownames = FALSE, booktabs = TRUE,
                       sanitize.text.function = identity, size = "\\small",
                       floating = FALSE))
}

# Common preamble shared by every cohort-battery table note; kept as its own
# constant so SPEC_NOTE and SPEC_NOTE_NOPT cannot drift apart (previously a
# fragile regex stripped the pre-trend clause off SPEC_NOTE; a script-editor
# hand-adding words to the pre-trend clause could then silently break the
# regex and leave a stray fragment in SPEC_NOTE_NOPT).
SPEC_NOTE_PREAMBLE <- paste0(
  "Unless noted: headline specification (CS\\,(2021) DR, never-treated controls, WGI GE + log ",
  "population, cohorts $\\geq$ ", THIN_THRESHOLD, "), multiplier-bootstrap SE (", BITERS,
  " reps, seed 1242)."
)
SRC_NOTE <- "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"
SPEC_NOTE <- paste0(SPEC_NOTE_PREAMBLE, " ", PRETREND_NOTE_AGG(-5L, -2L, 4L),
                    " (analytical twin fit)")
# Tables without a "Pre-trend p" column must not define one.
SPEC_NOTE_NOPT <- SPEC_NOTE_PREAMBLE

# Pre-sized accumulator: 4 LOCO + 2 drop-2024 + 2 balanced + 4 grid
# cells + 2 anticipation + 3 event-date + 3 placebo rungs + 2 conditioning = 22.
# add_to_battery() fills the first empty slot and grows only if that estimate is
# ever exceeded, so the common path never reallocates.
battery_summary <- vector("list", 22L)

#' Append a record to the combined battery summary.
#'
#' Returns the extended list rather than assigning into the enclosing
#' environment: `<<-` is prohibited by this project's R standards.
#'
#' @param acc current list of records
#' @param block block heading for the combined table
#' @param label specification label
#' @param res result list from run_spec() (may be NULL)
#' @return the extended list
add_to_battery <- function(acc, block, label, res) {
  slot <- which(vapply(acc, is.null, logical(1L)))
  idx  <- if (length(slot) > 0L) slot[[1L]] else length(acc) + 1L
  acc[[idx]] <- list(block = block, label = label, res = res)
  acc
}

# ==============================================================================
# SECTION 3. A.1 — leave-one-cohort-out (LOCO)
# Each estimation cohort in turn is removed from the treated group (its
# recipients leave the sample entirely; controls are untouched). With four
# cohorts of 8-11 recipients each, this is the sharpest available test of
# whether one adoption wave drives the headline result.
# ==============================================================================

message("\n=== A.1: leave-one-cohort-out ===\n")

loco_rows <- vector("list", length(kept_cohorts))   # pre-allocated
names(loco_rows) <- as.character(kept_cohorts)
for (g in kept_cohorts) {
  panel_g <- did_panel_main %>% filter(cohort_year != g)
  n_drop  <- n_distinct(did_panel_main$country_id[did_panel_main$cohort_year == g])
  message(sprintf("-- dropping cohort %d (%d treated recipients)", g, n_drop))
  res <- run_spec(panel_g, yname = "log_commits")
  lbl <- sprintf("Drop cohort %d (%d recipients)", g, n_drop)
  loco_rows[[as.character(g)]] <- spec_row(lbl, res)
  battery_summary <- add_to_battery(battery_summary, "Leave-one-cohort-out", lbl, res)
  if (!is.null(res))
    message(sprintf("   ATT = %.4f (SE %.4f, t %.3f)", res$att, res$se, res$t))
}

loco_tab <- bind_rows(
  spec_row(sprintf("Baseline: all %d cohorts (headline specification)", length(kept_cohorts)),
           list(att = fit_head$att, se = fit_head$se,
                t = fit_head$att / fit_head$se,
                stars = fmt_stars(fit_head$att / fit_head$se),
                pretrend = fit_head$pretrend, n_obs = fit_head$n_obs,
                n_country = fit_head$n_country,
                n_treated = n_distinct(did_panel_main$country_id[
                  did_panel_main$cohort_year > 0]))),
  bind_rows(loco_rows))

write_tex_float(
  file.path(dir_tabs_cb, "att_loco.tex"),
  "Leave-one-cohort-out: headline ATT dropping each adoption cohort in turn",
  "tab:att_loco",
  tabular_of(loco_tab, "tab:att_loco"),
  paste0("Outcome: log(adaptation commitments). Each row removes all recipients of one ",
         "adoption cohort; controls untouched. Baseline row read from the stored headline ",
         "fit (Table~\\ref{tab:combined_wide_main}). ", SPEC_NOTE),
  SRC_NOTE)

# ==============================================================================
# SECTION 4. A.2 — dropping 2024 as a cohort and as a calendar year
# The 2024 cohort contributes only the impact year e = 0, and 2024 is the last
# (and most provisional) CRS vintage year, so the two exclusions test different
# things: cohort composition versus data vintage.
# ==============================================================================

message("\n=== A.2: dropping 2024 (cohort / calendar year) ===\n")

panel_no2024_cohort <- did_panel_main %>% filter(cohort_year != 2024)
res_no2024_cohort   <- run_spec(panel_no2024_cohort, yname = "log_commits")

# Dropping 2024 as a CALENDAR year leaves the 2024 cohort with no post-treatment
# period; att_gt() then has no estimable ATT(g,t) for g = 2024, so that cohort
# is removed as well. We drop it explicitly to keep the message honest.
panel_no2024_year <- did_panel_main %>%
  filter(year < 2024, cohort_year != 2024)
res_no2024_year   <- run_spec(panel_no2024_year, yname = "log_commits")

drop24_tab <- bind_rows(
  spec_row("Baseline (headline specification)",
           list(att = fit_head$att, se = fit_head$se, t = fit_head$att / fit_head$se,
                stars = fmt_stars(fit_head$att / fit_head$se),
                pretrend = fit_head$pretrend, n_obs = fit_head$n_obs,
                n_country = fit_head$n_country,
                n_treated = n_distinct(did_panel_main$country_id[
                  did_panel_main$cohort_year > 0]))),
  spec_row("Drop the 2024 adoption cohort", res_no2024_cohort),
  spec_row("Drop calendar year 2024 (and the 2024 cohort)", res_no2024_year))

battery_summary <- add_to_battery(battery_summary, "Drop 2024",
                                  "Drop the 2024 adoption cohort", res_no2024_cohort)
battery_summary <- add_to_battery(battery_summary, "Drop 2024",
                                  "Drop calendar year 2024", res_no2024_year)

write_tex_float(
  file.path(dir_tabs_cb, "att_drop2024.tex"),
  "Sensitivity to the 2024 adoption cohort and to the 2024 CRS vintage year",
  "tab:att_drop2024",
  tabular_of(drop24_tab, "tab:att_drop2024"),
  paste0("Outcome: log(adaptation commitments). Row 2 drops the 2024 adoption cohort; row 3 ",
         "also drops calendar year 2024 for every recipient, which removes that cohort's ",
         "post-treatment period entirely. Baseline row read from the stored headline fit. ",
         SPEC_NOTE),
  SRC_NOTE)

# ==============================================================================
# SECTION 5. A.3 — balanced event windows (balance_e = 1, 2)
# aggte(balance_e = k) restricts the dynamic aggregation to cohorts observed for
# at least k post-treatment periods, so the event-study path is not a
# composition effect. balance_e = 1 drops the 2024 cohort from the aggregation,
# balance_e = 2 also drops 2023.
# ==============================================================================

message("\n=== A.3: balanced event windows ===\n")

balance_grid <- c(NA_integer_, 1L, 2L)
# Name-keyed result containers (<= 4 entries each): assignment is by label, so
# they are declared empty rather than pre-sized -- pre-sizing an unnamed list and
# assigning by name would append instead of fill. These are not hot loops.
balance_rows <- list()
balance_res  <- list()
for (be in balance_grid) {
  lbl <- if (is.na(be)) "Unbalanced (baseline aggregation)" else
    sprintf("balance\\_e = %d", be)
  message("-- ", lbl)
  # The unbalanced row is the STORED dynamic aggregation of the headline fit,
  # not a re-estimate: re-running it would draw a fresh multiplier bootstrap
  # and print a second standard error for a specification the paper already
  # reports. Only the balanced windows, which the stored object does not
  # contain, are estimated here.
  res <- if (is.na(be)) NULL else
    run_spec(did_panel_main, yname = "log_commits", balance_e = be)
  balance_res[[lbl]] <- res
  agg_bal <- if (is.na(be)) fit_head$agg_dynamic else
    if (!is.null(res)) res$agg_d else NULL
  # For balanced windows the object of interest is the DYNAMIC average, which
  # is what balance_e modifies; the simple ATT is unaffected by balance_e.
  dyn_att <- if (!is.null(agg_bal)) agg_bal$overall.att else NA_real_
  dyn_se  <- if (!is.null(agg_bal)) agg_bal$overall.se  else NA_real_
  dyn_t   <- if (!is.na(dyn_att) && !is.na(dyn_se) && dyn_se > 0) dyn_att / dyn_se else NA_real_
  # Event times actually estimated: e = -1 is the normalised base period (zero
  # standard error) and is excluded from the displayed window.
  cohorts_in <- if (!is.null(agg_bal)) {
    ok_e <- agg_bal$egt[!is.na(agg_bal$se.egt) & agg_bal$se.egt > 1e-10]
    paste(sprintf("%d", as.integer(sort(unique(ok_e)))), collapse = ", ")
  } else "---"
  balance_rows[[lbl]] <- data.frame(
    Specification = lbl,
    `Dynamic ATT` = if (is.na(dyn_att)) "---" else
      paste0(sprintf("%.4f", dyn_att), fmt_stars(dyn_t)),
    SE            = if (is.na(dyn_se)) "---" else sprintf("(%.4f)", dyn_se),
    `$t$`         = if (is.na(dyn_t)) "---" else sprintf("%.3f", dyn_t),
    `Event times` = cohorts_in,
    `$N$`         = if (is.na(be)) format(fit_head$n_obs, big.mark = ",") else
      if (is.null(res)) "---" else format(res$n_obs, big.mark = ","),
    check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.na(dyn_att))
    message(sprintf("   dynamic ATT = %.4f (SE %.4f)%s", dyn_att, dyn_se,
                    if (is.na(be)) "  [read from the stored headline fit]" else ""))
  if (!is.na(be)) {
    # The summary must report the same quantity as att_balance.tex: the
    # dynamic average under balance_e, not the (unchanged) simple ATT.
    res_bal <- res
    res_bal$att <- dyn_att; res_bal$se <- dyn_se; res_bal$t <- dyn_t
    p_bal <- if (is.na(dyn_t)) NA_real_ else 2 * pnorm(-abs(dyn_t))
    res_bal$stars <- if (is.na(p_bal)) "" else if (p_bal < 0.01) "***" else if (p_bal < 0.05) "**" else if (p_bal < 0.10) "*" else ""
    battery_summary <- add_to_battery(battery_summary, "Balanced event window",
                                      paste0(lbl, " (dynamic average)"), res_bal)
  }
}

write_tex_float(
  file.path(dir_tabs_cb, "att_balance.tex"),
  "Balanced event windows: dynamic ATT under \\texttt{balance\\_e} restrictions",
  "tab:att_balance",
  tabular_of(bind_rows(balance_rows), "tab:att_balance"),
  paste0("Outcome: log(adaptation commitments). Estimate: average dynamic ",
         "effect (\\texttt{aggte}, type $=$ dynamic), not the simple ATT; ",
         "\\texttt{balance\\_e} $= k$ restricts to cohorts observed $\\geq k$ post-treatment ",
         "periods ($k=1$ excludes 2024; $k=2$ also excludes 2023). Unbalanced row read from ",
         "the stored headline fit. ", SPEC_NOTE_NOPT),
  SRC_NOTE)

# ==============================================================================
# SECTION 6. A.4 — the 2x2 grid {main, retained} x {dr, reg}
# The published pair confounds two changes at once: the retained-cohort
# specification switches BOTH the cohort rule AND the estimator. This grid separates them, with the bootstrap everywhere.
# ==============================================================================

message("\n=== A.4: 2x2 cohort rule x estimator ===\n")

grid_panels <- list(
  list(key = "main",     label = "Cohorts $\\geq 5$", panel = did_panel_main),
  list(key = "retained", label = "All cohorts retained", panel = did_panel_full)
)
grid_rows <- list()   # name-keyed, 4 entries (see note in §5)
grid_notes_extra <- character(0)

for (gp in grid_panels) {
  for (em in c("dr", "reg")) {
    lbl <- sprintf("%s, est\\_method = %s", gp$label, em)
    message(sprintf("-- %s x %s", gp$key, em))
    risky <- identical(gp$key, "retained") && identical(em, "dr")
    res <- if (risky) {
      message("   (running in a forked process: att_gt(dr) can segfault on ",
              "panels containing 2-unit cohorts)")
      run_spec_forked(gp$panel, yname = "log_commits", est_method = em)
    } else {
      run_spec(gp$panel, yname = "log_commits", est_method = em)
    }
    if (risky && is.null(res))
      grid_notes_extra <- c(grid_notes_extra,
        paste0("The doubly-robust estimator could not be evaluated on the retained-cohort ",
               "panel (\\texttt{fastglm} segfaults on 2-unit cohorts); that cell ran in a ",
               "forked process and is reported as unavailable, with the outcome-regression ",
               "cell as the comparable estimate."))
    grid_rows[[lbl]] <- spec_row(lbl, res)
    battery_summary <- add_to_battery(battery_summary, "Cohort rule $\\times$ estimator", lbl, res)
    if (!is.null(res))
      message(sprintf("   ATT = %.4f (SE %.4f, t %.3f) | N = %d | treated = %d",
                      res$att, res$se, res$t, res$n_obs, res$n_treated))
  }
}

write_tex_float(
  file.path(dir_tabs_cb, "att_2x2.tex"),
  "Cohort rule and estimator, varied one at a time",
  "tab:att_2x2",
  tabular_of(bind_rows(grid_rows), "tab:att_2x2"),
  paste0("Outcome: log(adaptation commitments). This grid varies the cohort rule and ",
         "estimator separately, with common multiplier-bootstrap inference in all four ",
         "cells. ``All cohorts retained'' adds the 2015--2020 cohorts of 2--4 recipients ",
         "each. ", paste(grid_notes_extra, collapse = " "), " ", SPEC_NOTE),
  SRC_NOTE)

# ==============================================================================
# SECTION 7. A.5 — which cohorts contribute to each event time
# A dynamic ATT at event time e is an average over the cohorts that are observed
# at e. Because the panel ends in 2024, later cohorts drop out of later event
# times; this table makes the composition explicit.
# ==============================================================================

message("\n=== A.5: contributing cohorts per event time ===\n")

gt_ref <- fit_head$gt_boot
cells <- data.frame(g = as.numeric(gt_ref$group), t = as.numeric(gt_ref$t),
                    att = as.numeric(gt_ref$att), se = as.numeric(gt_ref$se))
# Integer year differences: coerced so the selection below is an exact integer
# comparison rather than a float equality test.
cells$e <- as.integer(round(cells$t - cells$g))
n_by_cohort <- did_panel_main %>%
  filter(cohort_year > 0) %>%
  distinct(country_id, cohort_year) %>%
  count(cohort_year, name = "n_units")

e_grid <- -5:3
contrib_rows <- lapply(e_grid, function(e) {
  ce <- cells[cells$e == e & !is.na(cells$att) & !is.na(cells$se) & cells$se > 1e-10, ]
  gs <- sort(unique(ce$g))
  nu <- sum(n_by_cohort$n_units[n_by_cohort$cohort_year %in% gs])
  data.frame(
    `Event time $e$`      = sprintf("%d", e),
    `Contributing cohorts` = if (length(gs) == 0L) "---" else
      paste(sprintf("%d", gs), collapse = ", "),
    `Cohorts`             = as.character(length(gs)),
    `Treated recipients`  = if (length(gs) == 0L) "0" else as.character(nu),
    check.names = FALSE, stringsAsFactors = FALSE)
})
contrib_tab <- bind_rows(contrib_rows)
print(contrib_tab, row.names = FALSE)

write_tex_float(
  file.path(dir_tabs_cb, "cohort_contributions.tex"),
  "Cohorts contributing to each event time in the headline event study",
  "tab:cohort_contributions",
  tabular_of(contrib_tab, "tab:cohort_contributions"),
  paste0("Composition of the dynamic aggregation used by Figure~\\ref{fig:did_combined_es} ",
         "and the HonestDiD analysis. A cohort contributes to event time $e$ when its ",
         "$ATT(g, g+e)$ cell is estimable in the stored headline fit (non-missing, ",
         "$se > 0$); $e = -1$ is the base period and contributes none. Counts are treated ",
         "recipients of the contributing cohorts"),
  SRC_NOTE)

# ==============================================================================
# SECTION 8. B.1 — anticipation
# att_gt(anticipation = k) treats the k years before adoption as already
# treated, i.e. it moves the base period back to g - k - 1. If donors respond to
# the drafting process rather than to submission, the headline understates the
# effect and the anticipation specifications should raise it.
# ==============================================================================

message("\n=== B.1: anticipation ===\n")

antic_grid <- c(0L, 1L, 2L)
antic_rows <- list()   # name-keyed, 3 entries (see note in §5)
antic_dyn  <- list()
antic_e0   <- list()
antic_leads <- list()
for (k in antic_grid) {
  lbl <- sprintf("Anticipation = %d year%s", k, if (k == 1L) "" else "s")
  message("-- ", lbl)
  res <- if (k == 0L) NULL else
    run_spec(did_panel_main, yname = "log_commits", anticipation = k,
             min_e = -5 - k)
  if (k == 0L) {
    # Baseline row read from the stored headline fit. Its e >= 0 average
    # is the same object as the pooled ATT up to the group-size weighting did
    # applies, and is reported for comparability with k = 1, 2.
    base_row <- spec_row(
      "Anticipation = 0 (baseline, headline specification)",
      list(att = fit_head$att, se = fit_head$se, t = fit_head$att / fit_head$se,
           stars = fmt_stars(fit_head$att / fit_head$se),
           pretrend = fit_head$pretrend, n_obs = fit_head$n_obs,
           n_country = fit_head$n_country,
           n_treated = n_distinct(did_panel_main$country_id[
             did_panel_main$cohort_year > 0])))
    set.seed(1242)
    agg_e0_base <- tryCatch(
      aggte(fit_head$gt_boot, type = "dynamic", na.rm = TRUE, min_e = 0, max_e = Inf),
      error = function(e) NULL)
    antic_e0[[lbl]] <- agg_e0_base
    # The baseline row's lead window is that of the stored headline fit.
    antic_leads[[lbl]] <- if (is.null(fit_head$agg_dyn_analytic)) integer(0) else {
      ad <- fit_head$agg_dyn_analytic
      kk <- which(!is.na(ad$se.egt) & ad$se.egt > 1e-10)
      as.integer(ad$egt[kk][ad$egt[kk] < 0])
    }
    antic_rows[[lbl]] <- base_row
    antic_dyn[[lbl]] <- fit_head$agg_dynamic
    if (!is.null(agg_e0_base))
      message(sprintf("   e >= 0 average = %.4f (SE %.4f)",
                      agg_e0_base$overall.att, agg_e0_base$overall.se))
  } else {
    antic_rows[[lbl]] <- spec_row(lbl, res)
    antic_dyn[[lbl]]  <- if (!is.null(res)) res$agg_d else NULL
    antic_e0[[lbl]]   <- if (!is.null(res)) res$agg_e0 else NULL
    antic_leads[[lbl]] <- if (!is.null(res)) res$pre_leads else integer(0)
    battery_summary <- add_to_battery(battery_summary, "Anticipation", lbl, res)
    if (!is.null(res))
      message(sprintf(paste0("   pooled ATT = %.4f (SE %.4f, t %.3f) | ",
                             "e >= 0 average = %.4f (SE %.4f) | pre-trend df = %d"),
                      res$att, res$se, res$t, res$att_e0, res$se_e0,
                      res$pretrend$df))
      message(sprintf("     pre-trend leads used: %s",
                      if (length(res$pre_leads) == 0L) "none" else
                        paste(res$pre_leads, collapse = ", ")))
  }
}

# Dynamic profiles for the three anticipation settings, on a common event-time
# grid so the table is readable.
dyn_grid <- -3:3
dyn_block <- lapply(names(antic_dyn), function(lbl) {
  ad <- antic_dyn[[lbl]]
  k_lbl <- as.integer(sub("^Anticipation = ([0-9]+).*$", "\\1", lbl))
  cells_e <- vapply(dyn_grid, function(e) {
    if (is.null(ad)) return("---")
    i <- match(e, ad$egt)
    if (is.na(i) || is.na(ad$se.egt[i]) || ad$se.egt[i] <= 1e-10) return("---")
    # Cells the estimator treats as TREATED under anticipation = k (e = -1,
    # ..., -k) are daggered: they are post-treatment by construction and must
    # not be read as leads.
    mark <- if (k_lbl > 0L && e < 0L && e >= -k_lbl) "$^{\\dagger}$" else ""
    sprintf("%.3f%s (%.3f)", ad$att.egt[i], mark, ad$se.egt[i])
  }, character(1L))
  setNames(as.data.frame(t(c(lbl, cells_e)), stringsAsFactors = FALSE),
           c("Specification", sprintf("$e = %d$", dyn_grid)))
})
dyn_tab <- bind_rows(dyn_block)

antic_tab <- bind_rows(antic_rows)

# The pooled column is did's SIMPLE aggregate, which under anticipation = k
# averages the k pre-submission cells together with the post-submission ones.
# The e >= 0 column removes that mechanical component.
antic_tab$`ATT, $e \\geq 0$` <- vapply(names(antic_rows), function(lbl) {
  a <- antic_e0[[lbl]]
  if (is.null(a) || is.na(a$overall.att)) return("---")
  sprintf("%.4f%s (%.4f)", a$overall.att,
          fmt_stars(a$overall.att / a$overall.se), a$overall.se)
}, character(1L))
# Lead window actually used by each row's pre-trend test, reported rather than
# derived from a formula (the base period shifts with k, so the window is
# k + 1 periods back, not k).
antic_tab$`Pre-trend leads` <- vapply(names(antic_rows), function(lbl) {
  lv <- antic_leads[[lbl]]
  if (is.null(lv) || length(lv) == 0L) "---" else
    sprintf("$%d \\leq e \\leq %d$ (%d)", min(lv), max(lv), length(lv))
}, character(1L))

write_tex_float(
  file.path(dir_tabs_an, "att_anticipation.tex"),
  "Anticipation: pooled and dynamic ATTs allowing donors to respond before submission",
  "tab:att_anticipation",
  c(tabular_of(antic_tab, "tab:att_anticipation"),
    "\\par\\vspace{6pt}",
    tabular_of(dyn_tab, "tab:att_anticipation_dyn")),
  paste0("Upper panel: pooled ATT and the $e \\geq 0$ average (see main text). Lower panel: ",
         "dynamic ATTs by event time, SE in parentheses. \\texttt{anticipation} $= k$ shifts the ",
         "base period to $g-k-1$. $^{\\dagger}$ marks event times treated as TREATED under that ",
         "row's $k$, not leads. Anticipation-0 row is read from the stored headline fit. ",
         SPEC_NOTE),
  SRC_NOTE)

# ==============================================================================
# SECTION 9. B.2 — event date defined by the submission MONTH
#
# nap_information.csv records the day a NAP was posted on NAP Central, and the
# processed panel carries it both as the string `date_posted` and as the
# pre-parsed `date_posted_2`. A plan posted in November cannot plausibly move
# that calendar year's commitments, so we re-date treatment as
#   G_i = year(posted) + 1{month(posted) > 6}.
#
# WHY THERE ARE THREE ROWS (2026-09-15). Re-dating is not a clean
# one-dimensional perturbation: it moves recipients ACROSS cohorts, so the
# >= 5-treated-units rule then selects a different set of cohorts, and a naive
# implementation silently (i) re-admits 2015-2020 adopters that the main
# specification excludes as thin -- they pile up into a re-dated 2020 cohort of
# five -- and (ii) recodes late-2024 adopters, whose re-dated cohort 2025 lies
# beyond the panel, into the NEVER-TREATED CONTROL GROUP. Recoding a treated
# recipient as a control is not a re-dating robustness check; it is a different
# experiment. The block therefore reports:
#   (a1) re-dating WITHIN the main treated set: the same 40 recipients, the
#        >= 5 rule re-applied to THEIR re-dated cohorts only; recipients whose
#        re-dated cohort leaves the panel are DROPPED, never recoded as
#        controls. This is the clean comparison.
#   (a2) the same 40 recipients with every re-dated cohort retained (no >= 5
#        rule), so that no treated unit is lost to the cohort filter -- the
#        natural fallback when re-dating splits a cohort below the
#        threshold.
#   (b)  the full-sample variant (all adopters re-dated), again DROPPING rather
#        than re-coding the recipients pushed beyond the panel. Its cohort set
#        differs from the main specification by construction, which the note
#        states.
# The original x re-dated cross-tabulation is exported so the reader can see
# exactly which recipients move where.
# ==============================================================================

message("\n=== B.2: event date from the submission month ===\n")

# 01_prepare_data.R already parses the posting date into `date_posted_2`; we use
# it when present. The string fallback maps English month names explicitly
# rather than calling strptime("%B"), which is LOCALE-DEPENDENT and silently
# returns NA under a non-English LC_TIME (this machine runs fr_FR).
month_names_en <- c("January", "February", "March", "April", "May", "June",
                    "July", "August", "September", "October", "November",
                    "December")

parse_posted <- function(date_str) {
  # Both orders occur in UNFCCC listings: "September 29, 2021" and
  # "29 September 2021". Take the first alphabetic run as the month name
  # wherever it sits, and the last four-digit run as the year, so neither
  # variant silently returns NA.
  m_word <- sub("^[^A-Za-z]*([A-Za-z]+).*$", "\\1", date_str)
  m_num  <- match(m_word, month_names_en)
  y_num  <- suppressWarnings(as.integer(sub("^.*?([0-9]{4})[^0-9]*$", "\\1",
                                            date_str, perl = TRUE)))
  list(month = m_num, year = y_num)
}

nap_dates <- did_panel_full %>%
  filter(!is.na(date_posted), date_posted != "") %>%
  distinct(recipient_name, date_posted, date_posted_2)

if ("date_posted_2" %in% names(nap_dates) && sum(!is.na(nap_dates$date_posted_2)) > 0L) {
  nap_dates <- nap_dates %>%
    mutate(post_month = as.integer(format(as.Date(date_posted_2), "%m")),
           post_year  = as.integer(format(as.Date(date_posted_2), "%Y")))
  message("  Submission month taken from the pre-parsed `date_posted_2` column.")
} else {
  pp <- parse_posted(nap_dates$date_posted)
  nap_dates <- nap_dates %>% mutate(post_month = pp$month, post_year = pp$year)
  message("  Submission month parsed from the `date_posted` string ",
          "(locale-independent month map).")
}

n_parsed <- sum(!is.na(nap_dates$post_month))
message(sprintf("  Submission dates parsed: %d of %d recipient records",
                n_parsed, nrow(nap_dates)))

cohort_month_ok <- n_parsed > 0L
res_month_a1 <- NULL
res_month_a2 <- NULL
res_month_b  <- NULL

if (!cohort_month_ok) {
  message("  No parseable submission month — reporting only the e >= 1 aggregation.")
} else {
  month_lookup <- nap_dates %>%
    filter(!is.na(post_month), !is.na(post_year)) %>%
    arrange(recipient_name, post_year, post_month) %>%
    group_by(recipient_name) %>%
    slice_head(n = 1L) %>%     # first (earliest) submission per recipient
    ungroup() %>%
    select(recipient_name, post_month, post_year)

  month_counts <- month_lookup %>% count(post_month, name = "n")
  message("  NAP submissions by month (all recipients with a date):")
  print(as.data.frame(month_counts), row.names = FALSE)

  # Unit-level re-dating map: original cohort, re-dated cohort, and whether the
  # re-dated cohort leaves the panel.
  redate_map <- did_panel_full %>%
    distinct(recipient_name, country_id, cohort_year) %>%
    left_join(month_lookup, by = "recipient_name") %>%
    mutate(
      g_month = case_when(
        cohort_year == 0                 ~ 0,
        is.na(post_month)                ~ cohort_year,
        post_month > 6L                  ~ cohort_year + 1,
        TRUE                             ~ cohort_year),
      beyond_panel = cohort_year > 0 & g_month > last_year,
      in_main_40   = cohort_year %in% kept_cohorts)

  n_adopters_all <- sum(redate_map$cohort_year > 0)
  n_moved        <- sum(redate_map$cohort_year > 0 & redate_map$g_month != redate_map$cohort_year)
  n_beyond       <- sum(redate_map$beyond_panel)
  message(sprintf(paste0("  Re-dated recipients: %d of %d adopters; %d pushed ",
                         "beyond the panel end (%d) and DROPPED (not recoded as controls)"),
                  n_moved, n_adopters_all, n_beyond, last_year))
  message(sprintf("  Missing submission month among adopters: %d",
                  sum(redate_map$cohort_year > 0 & is.na(redate_map$post_month))))

  #' Attach the re-dated cohort to a panel and drop beyond-panel recipients.
  #'
  #' @param panel recipient-year panel
  #' @param ids_keep treated country_ids to keep as treated (controls always kept)
  #' @return panel with `cohort_year` replaced by the re-dated cohort
  build_month_panel <- function(panel, ids_keep) {
    drop_ids <- redate_map$country_id[redate_map$beyond_panel]
    panel %>%
      filter(!(country_id %in% drop_ids)) %>%
      filter(cohort_year == 0 | country_id %in% ids_keep) %>%
      select(-cohort_year) %>%
      left_join(redate_map %>% select(country_id, cohort_year = g_month),
                by = "country_id")
  }

  # ---- (a1) re-dating within the main treated set, >= 5 rule re-applied ------
  ids_main_treated <- sort(unique(did_panel_main$country_id[did_panel_main$cohort_year > 0]))
  panel_a <- build_month_panel(did_panel_full, ids_main_treated)
  sizes_a <- cohort_sizes_of(panel_a)
  message("  (a) re-dated cohort sizes within the main treated set:")
  print(as.data.frame(sizes_a), row.names = FALSE)

  cf_a1 <- apply_cohort_filter(panel_a)
  lost_a1 <- setdiff(sizes_a$cohort_year, cf_a1$kept_cohorts)
  message(sprintf(paste0("  (a1) cohorts entering: %s | falling below the >= %d ",
                         "rule and dropped: %s"),
                  paste(cf_a1$kept_cohorts, collapse = ", "), THIN_THRESHOLD,
                  if (length(lost_a1) == 0L) "none" else paste(lost_a1, collapse = ", ")))
  res_month_a1 <- run_spec(cf_a1$panel, yname = "log_commits")
  rep_a1 <- cohort_cell_report(res_month_a1$gt, cf_a1$panel)
  if (!is.null(rep_a1)) { message("  (a1) estimable cells per cohort:")
    print(rep_a1, row.names = FALSE) }
  if (!is.null(res_month_a1))
    message(sprintf(paste0("   (a1) ATT = %.4f (SE %.4f, t %.3f) | treated in ",
                           "panel = %d | CONTRIBUTING treated = %d"),
                    res_month_a1$att, res_month_a1$se, res_month_a1$t,
                    res_month_a1$n_treated,
                    sum(rep_a1$treated_in_panel[rep_a1$contributes])))

  # ---- (a2) same recipients, every re-dated cohort retained ------------------
  # Re-dating splits original cohorts: within one adoption year some recipients
  # posted before July and some after, so a cohort that met the >= 5 rule can
  # break into two pieces that do not (the sizes are printed above and exported
  # in Table tab:cohort_redating). (a2) keeps every re-dated cohort so that no
  # treated recipient is lost to the filter. Cohorts of one or two units are
  # exactly the configuration in which att_gt(est_method = "dr") can segfault
  # inside fastglm, so this cell is evaluated in a forked process and falls back
  # to outcome regression if the doubly-robust fit cannot be produced.
  res_month_a2 <- run_spec_forked(panel_a, yname = "log_commits", est_method = "dr")
  em_a2 <- "dr"
  if (is.null(res_month_a2)) {
    message("   (a2) dr unavailable on the retained re-dated panel — using reg.")
    res_month_a2 <- run_spec(panel_a, yname = "log_commits", est_method = "reg")
    em_a2 <- "reg"
  }
  rep_a2 <- cohort_cell_report(res_month_a2$gt, panel_a)
  if (!is.null(rep_a2)) { message("  (a2) estimable cells per cohort:")
    print(rep_a2, row.names = FALSE) }
  if (!is.null(res_month_a2))
    message(sprintf(paste0("   (a2) ATT = %.4f (SE %.4f, t %.3f) | treated in ",
                           "panel = %d | CONTRIBUTING treated = %d | est_method = %s"),
                    res_month_a2$att, res_month_a2$se, res_month_a2$t,
                    res_month_a2$n_treated,
                    sum(rep_a2$treated_in_panel[rep_a2$contributes]), em_a2))
  # Reconciliation of (a1) and (a2): if the cohorts that (a1) drops contribute
  # no estimable post-treatment cell in (a2) either, aggte(na.rm = TRUE) removes
  # them and the two rows must agree to the last digit. That is asserted rather
  # than left as a coincidence for the reader to wonder about.
  a2_noncontrib <- if (is.null(rep_a2)) integer(0) else
    rep_a2$cohort[!rep_a2$contributes]
  if (!is.null(res_month_a1) && !is.null(res_month_a2)) {
    same_att <- abs(res_month_a1$att - res_month_a2$att) < 1e-8
    message(sprintf(paste0("   (a1) vs (a2): ATTs %s; cohorts in (a2) with no ",
                           "estimable post-treatment cell: %s"),
                    if (same_att) "identical" else "differ",
                    if (length(a2_noncontrib) == 0L) "none" else
                      paste(a2_noncontrib, collapse = ", ")))
    if (same_att) stopifnot(length(a2_noncontrib) > 0L)
  }

  # ---- (b) full-sample variant, beyond-panel recipients dropped -------------
  ids_all_treated <- sort(unique(redate_map$country_id[redate_map$cohort_year > 0]))
  panel_b <- build_month_panel(did_panel_full, ids_all_treated)
  cf_b <- apply_cohort_filter(panel_b)
  readmitted <- redate_map %>%
    filter(cohort_year > 0, !in_main_40, g_month %in% cf_b$kept_cohorts,
           !beyond_panel)
  message(sprintf(paste0("  (b) cohorts entering: %s | dropped as thin: %s | ",
                         "thin-cohort adopters re-admitted: %d (%s)"),
                  paste(cf_b$kept_cohorts, collapse = ", "),
                  paste(cf_b$thin_cohorts, collapse = ", "),
                  nrow(readmitted),
                  if (nrow(readmitted) == 0L) "none" else
                    paste(sort(unique(readmitted$cohort_year)), collapse = "/")))
  res_month_b <- run_spec(cf_b$panel, yname = "log_commits")
  rep_b <- cohort_cell_report(res_month_b$gt, cf_b$panel)
  if (!is.null(rep_b)) { message("  (b) estimable cells per cohort:")
    print(rep_b, row.names = FALSE) }
  if (!is.null(res_month_b))
    message(sprintf(paste0("   (b) ATT = %.4f (SE %.4f, t %.3f) | treated in ",
                           "panel = %d | CONTRIBUTING treated = %d"),
                    res_month_b$att, res_month_b$se, res_month_b$t,
                    res_month_b$n_treated,
                    sum(rep_b$treated_in_panel[rep_b$contributes])))

  battery_summary <- add_to_battery(battery_summary, "Event date",
                                    "Re-dated within the main treated set",
                                    res_month_a1)
  battery_summary <- add_to_battery(battery_summary, "Event date",
                                    "Re-dated, all cohorts retained",
                                    res_month_a2)
  battery_summary <- add_to_battery(battery_summary, "Event date",
                                    "Re-dated, full adopter sample",
                                    res_month_b)

  # ---- (c) original x re-dated cross-tabulation exhibit ----------------------
  ct <- redate_map %>%
    filter(cohort_year > 0) %>%
    mutate(g_lab = if_else(beyond_panel, "Beyond panel", as.character(g_month)))
  g_levels <- c(sort(unique(ct$g_month[!ct$beyond_panel])), NA)
  g_labs   <- c(as.character(sort(unique(ct$g_month[!ct$beyond_panel]))), "Beyond panel")
  orig_levels <- sort(unique(ct$cohort_year))

  ct_mat <- matrix("0", nrow = length(orig_levels), ncol = length(g_labs) + 2L)
  colnames(ct_mat) <- c(g_labs, "Total", "In main spec")
  for (i in seq_along(orig_levels)) {
    sub <- ct[ct$cohort_year == orig_levels[i], , drop = FALSE]
    for (j in seq_along(g_labs)) {
      n_ij <- sum(sub$g_lab == g_labs[j])
      ct_mat[i, j] <- if (n_ij == 0L) "--" else as.character(n_ij)
    }
    ct_mat[i, length(g_labs) + 1L] <- as.character(nrow(sub))
    ct_mat[i, length(g_labs) + 2L] <-
      if (orig_levels[i] %in% kept_cohorts) "Yes" else "No"
  }
  ct_tab <- cbind(`Original cohort` = as.character(orig_levels),
                  as.data.frame(ct_mat, stringsAsFactors = FALSE),
                  stringsAsFactors = FALSE)
  rownames(ct_tab) <- NULL
  print(ct_tab, row.names = FALSE)

  write_tex_float(
    file.path(dir_tabs_an, "cohort_redating_crosstab.tex"),
    "Re-dating by submission month: original against re-dated adoption cohort",
    "tab:cohort_redating",
    tabular_of(ct_tab, "tab:cohort_redating"),
    paste0("Treated recipients per (original cohort, re-dated cohort) cell. Re-dated cohort ",
           "$=$ submission year $+1$ if posted after June. ``Beyond panel'': re-dated cohort ",
           "falls in ", last_year + 1L, ", outside the window; dropped, not recoded as ",
           "controls (see main text). ``In main spec'': whether the original cohort meets ",
           "the $\\geq ", THIN_THRESHOLD, "$ rule. ", n_adopters_all,
           " adopters with a recorded date shown; ", n_moved, " change cohort"),
    SRC_NOTE)
}

# e >= 1 aggregation of the stored headline fit (excludes the impact year).
set.seed(1242)
agg_epos <- tryCatch(
  aggte(fit_head$gt_boot, type = "dynamic", na.rm = TRUE, min_e = 1, max_e = Inf),
  error = function(e) { message("  aggte(min_e = 1) failed: ",
                                 conditionMessage(e)); NULL })
if (!is.null(agg_epos))
  message(sprintf("   e >= 1 average = %.4f (SE %.4f)",
                  agg_epos$overall.att, agg_epos$overall.se))

# Pre-sized result container; filled by index, trimmed before binding.
eventdate_rows <- vector("list", 5L)
eventdate_rows[[1L]] <- (
  spec_row("Baseline: treatment dated by submission YEAR (headline specification)",
           list(att = fit_head$att, se = fit_head$se, t = fit_head$att / fit_head$se,
                stars = fmt_stars(fit_head$att / fit_head$se),
                pretrend = fit_head$pretrend, n_obs = fit_head$n_obs,
                n_country = fit_head$n_country,
                n_treated = n_distinct(did_panel_main$country_id[
                  did_panel_main$cohort_year > 0])))
)
if (cohort_month_ok) {
  eventdate_rows[[2L]] <- spec_row(
    "(a1) Re-dated, main treated set, $\\geq 5$ rule re-applied", res_month_a1)
  eventdate_rows[[3L]] <- spec_row(
    "(a2) Re-dated, main treated set, all re-dated cohorts retained", res_month_a2)
  eventdate_rows[[4L]] <- spec_row(
    "(b) Re-dated, full adopter sample, $\\geq 5$ rule re-applied", res_month_b)
}
if (!is.null(agg_epos))
  eventdate_rows[[5L]] <- data.frame(
    Specification   = "Dynamic average over $e \\geq 1$ (impact year excluded)",
    ATT             = paste0(sprintf("%.4f", agg_epos$overall.att),
                             fmt_stars(agg_epos$overall.att / agg_epos$overall.se)),
    SE              = sprintf("(%.4f)", agg_epos$overall.se),
    `$t$`           = sprintf("%.3f", agg_epos$overall.att / agg_epos$overall.se),
    `Pre-trend $p$` = if (is.na(fit_head$pretrend$pval)) "---" else
      sprintf("%.3f", fit_head$pretrend$pval),
    `$N$`           = format(fit_head$n_obs, big.mark = ","),
    Recipients      = as.character(fit_head$n_country),
    Treated         = as.character(n_distinct(did_panel_main$country_id[
      did_panel_main$cohort_year > 0])),
    check.names = FALSE, stringsAsFactors = FALSE)

eventdate_note <- paste0(
  "Outcome: log(adaptation commitments). Rows re-date treatment by posting month (post-June ",
  "$\\to$ next year); see Table~\\ref{tab:cohort_redating}. (a1)/(a2): main treated set, ",
  "$\\geq ", THIN_THRESHOLD, "$ rule re-applied vs.\\ retained. (b): full sample re-dated. ",
  "Contributing/Treated: estimable-cell counts. Last row: headline fit over $e \\geq 1$. ",
  SPEC_NOTE)

eventdate_tab <- bind_rows(eventdate_rows)

# Cohort set and CONTRIBUTING treated count per row:
# the "Treated" column of spec_row() counts recipients in the panel, which
# overstates what a row rests on when a cohort contributes no estimable cell.
ed_cohorts <- character(0)
ed_contrib <- character(0)
fmt_cohorts <- function(rep_df) {
  if (is.null(rep_df)) return("---")
  cc <- rep_df$cohort[rep_df$contributes]
  if (length(cc) == 0L) "none" else paste(sprintf("%d", cc), collapse = ", ")
}
fmt_contrib <- function(rep_df) {
  if (is.null(rep_df)) return("---")
  as.character(sum(rep_df$treated_in_panel[rep_df$contributes]))
}
rep_base <- cohort_cell_report(fit_head$gt_boot, did_panel_main)
ed_cohorts <- c(ed_cohorts, fmt_cohorts(rep_base))
ed_contrib <- c(ed_contrib, fmt_contrib(rep_base))
if (cohort_month_ok) {
  ed_cohorts <- c(ed_cohorts, fmt_cohorts(rep_a1), fmt_cohorts(rep_a2),
                  fmt_cohorts(rep_b))
  ed_contrib <- c(ed_contrib, fmt_contrib(rep_a1), fmt_contrib(rep_a2),
                  fmt_contrib(rep_b))
}
if (!is.null(agg_epos)) {
  ed_cohorts <- c(ed_cohorts, fmt_cohorts(rep_base))
  ed_contrib <- c(ed_contrib, fmt_contrib(rep_base))
}
stopifnot(length(ed_cohorts) == nrow(eventdate_tab))
eventdate_tab$`Contributing cohorts` <- ed_cohorts
eventdate_tab$Treated <- paste0(ed_contrib, " / ", eventdate_tab$Treated)

write_tex_float(
  file.path(dir_tabs_an, "att_eventdate.tex"),
  "Alternative treatment dates: submission month and the $e \\geq 1$ aggregation",
  "tab:att_eventdate",
  tabular_of(eventdate_tab, "tab:att_eventdate"),
  eventdate_note,
  SRC_NOTE)

# ==============================================================================
# SECTION 10. B.3 — placebo ladder at -3, -4 and -5 years
#
# The construction MIRRORS 04_robustness.R §25 exactly; the filters are stated
# here because the paper's text has previously described them incorrectly:
#   1. treated recipients enter with their PRE-ADOPTION observations only
#      (year < actual cohort_year), so no genuine post-treatment variation can
#      leak into the placebo;
#   2. their cohort is shifted back by k years;
#   3. never-adopters enter with all their years;
#   4. the sample is restricted to SHIFTED cohorts >= 2015 (i.e. actual cohorts
#      >= 2015 + k) -- this is the filter 04 applies, kept identical here;
#   5. cohorts with fewer than 5 treated recipients IN THE PLACEBO PANEL are
#      dropped.
# Estimator: outcome regression (as in §25), multiplier bootstrap, never-treated
# controls. For each rung we report the minimum detectable effect, because a
# placebo that cannot reject anything is not evidence of parallel trends.
# ==============================================================================

message("\n=== B.3: placebo ladder ===\n")

#' Build the placebo panel for a k-year backward shift (mirror of 04 §25).
build_placebo_panel <- function(panel_full, k) {
  p <- bind_rows(
    panel_full %>%
      filter(cohort_year > 0, year < cohort_year) %>%
      mutate(cohort_year = cohort_year - as.integer(k)),
    panel_full %>% filter(cohort_year == 0)
  ) %>%
    filter(cohort_year == 0 | cohort_year >= 2015)
  thin_p <- p %>%
    filter(cohort_year > 0) %>%
    group_by(cohort_year) %>%
    summarise(n = n_distinct(country_id), .groups = "drop") %>%
    filter(n < THIN_THRESHOLD) %>%
    pull(cohort_year)
  list(panel = p %>% filter(!(cohort_year %in% thin_p)), thin = thin_p)
}

# Pre-sized: rung -2 (read from the stored fit) plus -3, -4, -5.
placebo_rows <- vector("list", 4L)

# Rung -2 is the published placebo (Table tab:placebo). It is READ from the
# stored fit rather than re-estimated, so the ladder cannot print a second
# standard error for a specification the paper already reports.
path_plac2 <- file.path(FITS_DIR, "placebo_m2_reg_bs.rds")
if (file.exists(path_plac2)) {
  fit_p2 <- readRDS(path_plac2)
  t_p2 <- fit_p2$att / fit_p2$se
  placebo_rows[[1L]] <- data.frame(
    Placebo        = "$-2$ years (published, Table~\\ref{tab:placebo})",
    ATT            = paste0(sprintf("%.4f", fit_p2$att), fmt_stars(t_p2)),
    SE             = sprintf("(%.4f)", fit_p2$se),
    `$t$`          = sprintf("%.3f", t_p2),
    `MDE (80\\% power)` = sprintf("%.4f", mde_from_se(fit_p2$se)),
    `Cohorts`      = paste(fit_p2$cohorts, collapse = ", "),
    `Treated`      = as.character(fit_p2$n_treated),
    `$N$`          = format(fit_p2$n_obs, big.mark = ","),
    check.names = FALSE, stringsAsFactors = FALSE)
  message(sprintf("-- -2 years (from stored fit): ATT = %.4f (SE %.4f)",
                  fit_p2$att, fit_p2$se))
} else {
  message("  Stored -2-year placebo fit not found (", basename(path_plac2),
          ") — the ladder starts at -3. Run 04_robustness.R first to include it.")
}

for (k in 3:5) {
  bp <- build_placebo_panel(did_panel_full, k)
  pan <- bp$panel
  n_tr <- n_distinct(pan$country_id[pan$cohort_year > 0])
  coh  <- sort(unique(pan$cohort_year[pan$cohort_year > 0]))
  message(sprintf(paste0("-- placebo -%d: %d recipient-years, %d treated ",
                         "recipients, cohorts %s (dropped thin: %s)"),
                  k, nrow(pan), n_tr, paste(coh, collapse = ", "),
                  if (length(bp$thin) == 0L) "none" else paste(bp$thin, collapse = ", ")))
  res <- if (length(coh) < 2L) {
    message("   fewer than two placebo cohorts survive — skipped."); NULL
  } else {
    run_spec(pan, yname = "log_commits", est_method = "reg", min_e = -4)
  }
  placebo_rows[[k - 1L]] <- if (is.null(res)) data.frame(
    Placebo = sprintf("$-%d$ years", k), ATT = "---", SE = "---", `$t$` = "---",
    `MDE (80\\% power)` = "---", Cohorts = paste(coh, collapse = ", "),
    Treated = as.character(n_tr), `$N$` = format(nrow(pan), big.mark = ","),
    check.names = FALSE, stringsAsFactors = FALSE) else data.frame(
    Placebo             = sprintf("$-%d$ years", k),
    ATT                 = paste0(sprintf("%.4f", res$att), res$stars),
    SE                  = sprintf("(%.4f)", res$se),
    `$t$`               = sprintf("%.3f", res$t),
    `MDE (80\\% power)` = sprintf("%.4f", mde_from_se(res$se)),
    Cohorts             = paste(coh, collapse = ", "),
    Treated             = as.character(res$n_treated),
    `$N$`               = format(res$n_obs, big.mark = ","),
    check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(res))
    message(sprintf("   ATT = %.4f (SE %.4f, t %.3f) | MDE = %.4f",
                    res$att, res$se, res$t, mde_from_se(res$se)))
  battery_summary <- add_to_battery(battery_summary, "Placebo ladder",
                                    sprintf("Placebo -%d years", k), res)
}

write_tex_float(
  file.path(dir_tabs_an, "att_placebo_ladder.tex"),
  "Placebo ladder: fictitious adoption dates 2 to 5 years before the actual one",
  "tab:att_placebo_ladder",
  tabular_of(bind_rows(placebo_rows), "tab:att_placebo_ladder"),
  paste0("Outcome: log(adaptation commitments). Treated recipients contribute pre-adoption ",
         "years only; cohort shifted back by the stated years; sample restricted to shifted ",
         "cohorts from 2015; cohorts under ", THIN_THRESHOLD, " treated dropped. Estimator: ",
         "CS\\,(2021) OR (DR unstable on these panels), never-treated controls, multiplier-",
         "bootstrap SE (", BITERS, " reps, seed 1242). MDE $= 2.8016 \\times$ SE (80\\% power, ",
         "two-sided 5\\%). The $-2$ row is read from the stored fit behind Table~\\ref{tab:placebo}"),
  SRC_NOTE)

# ==============================================================================
# SECTION 11. Conditioning set
#   (a) add the recipient's mean 2009-2012 outcome level to xformla. The
#       balance table shows adopters and non-adopters differ most in the level
#       of adaptation finance they were already receiving (normalised
#       difference 0.537), which is exactly what a DiD cannot difference out if
#       the gap predicts DIFFERENTIAL GROWTH.
#   (b) common support: did does not expose the propensity score it uses
#       internally, so we fit the same specification's propensity model
#       ourselves -- a logit of ever-adopting on the cross-section of
#       pre-treatment covariates (mean 2009-2012 WGI government effectiveness,
#       mean 2009-2012 log population, and the baseline outcome level) -- and
#       drop recipients with a fitted propensity outside [0.05, 0.95]. Trimming
#       is at the RECIPIENT level: a recipient outside the support contributes
#       no recipient-year.
# ==============================================================================

message("\n=== C: conditioning set and common support ===\n")

BASE_YEARS <- 2009:2012

base_outcome_tab <- did_panel_full %>%
  filter(year %in% BASE_YEARS) %>%
  group_by(country_id) %>%
  summarise(base_outcome = mean(log_commits, na.rm = TRUE),
            base_ge      = mean(ge_est, na.rm = TRUE),
            base_lpop    = mean(log_population, na.rm = TRUE),
            .groups = "drop")

n_missing_base <- sum(!is.finite(base_outcome_tab$base_outcome))
message(sprintf(paste0("  Baseline (%d-%d) outcome level computed for %d ",
                       "recipients (%d with no usable baseline)"),
                min(BASE_YEARS), max(BASE_YEARS), nrow(base_outcome_tab),
                n_missing_base))

did_panel_base <- did_panel_main %>%
  left_join(base_outcome_tab, by = "country_id") %>%
  mutate(base_outcome = if_else(is.finite(base_outcome), base_outcome, NA_real_))

res_base <- run_spec(did_panel_base, yname = "log_commits",
                     xformla = ~ ge_est + log_population + base_outcome)
if (!is.null(res_base))
  message(sprintf("  Baseline-level control: ATT = %.4f (SE %.4f, t %.3f) | N = %d",
                  res_base$att, res_base$se, res_base$t, res_base$n_obs))
battery_summary <- add_to_battery(battery_summary, "Conditioning set",
                                  "Add baseline (2009--2012) outcome level", res_base)

# --- Common support ----------------------------------------------------------
cs_cross <- did_panel_main %>%
  distinct(country_id, cohort_year) %>%
  mutate(ever_treated = as.integer(cohort_year > 0)) %>%
  left_join(base_outcome_tab, by = "country_id") %>%
  filter(is.finite(base_outcome), is.finite(base_ge), is.finite(base_lpop))

message(sprintf(paste0("  Propensity model: %d recipients with complete ",
                       "pre-treatment covariates (%d treated)"),
                nrow(cs_cross), sum(cs_cross$ever_treated)))

ps_model <- glm(ever_treated ~ base_ge + base_lpop + base_outcome,
                family = binomial(link = "logit"), data = cs_cross)
# Guard the inverse link: fitted values are clamped to an open interval before
# any comparison, and the comparisons themselves use strict inequalities on
# doubles (never ==).
eps_ps <- 1e-12
cs_cross$pscore <- pmin(1 - eps_ps, pmax(eps_ps, as.numeric(fitted(ps_model))))
message(sprintf("  Propensity range: [%.4f, %.4f]; treated mean %.4f, control mean %.4f",
                min(cs_cross$pscore), max(cs_cross$pscore),
                mean(cs_cross$pscore[cs_cross$ever_treated == 1L]),
                mean(cs_cross$pscore[cs_cross$ever_treated == 0L])))

PS_LO <- 0.05
PS_HI <- 0.95
keep_ids <- cs_cross$country_id[cs_cross$pscore > PS_LO & cs_cross$pscore < PS_HI]
# Two distinct reasons a recipient can leave this specification, kept apart:
#   - no usable pre-treatment covariate window (listwise deletion BEFORE the
#     propensity model is fitted at all);
#   - a fitted propensity outside the support (the trimming itself).
dropped_missing <- setdiff(unique(did_panel_main$country_id), cs_cross$country_id)
dropped_trimmed <- setdiff(cs_cross$country_id, keep_ids)
dropped_support <- union(dropped_missing, dropped_trimmed)
message(sprintf(paste0("  Losses decomposed: %d recipients have no usable ",
                       "pre-treatment covariate window (listwise), %d are ",
                       "trimmed by the support rule"),
                length(dropped_missing), length(dropped_trimmed)))
if (length(dropped_trimmed) == 0L)
  message("  NOTE: the support restriction is NOT binding — no recipient lies ",
          "outside [", PS_LO, ", ", PS_HI, "], so this row reproduces the ",
          "headline exactly and is a disclosure, not a robustness result.")
message(sprintf(paste0("  Common support [%.2f, %.2f]: %d recipients kept, ",
                       "%d dropped (%d of them treated)"),
                PS_LO, PS_HI, length(keep_ids), length(dropped_support),
                n_distinct(did_panel_main$country_id[
                  did_panel_main$country_id %in% dropped_support &
                    did_panel_main$cohort_year > 0])))

panel_cs <- did_panel_main %>% filter(country_id %in% keep_ids)
cf_cs    <- apply_cohort_filter(panel_cs)
if (length(cf_cs$thin_cohorts) > 0L)
  message(sprintf("  Cohort rule re-applied after trimming: dropped %s",
                  paste(cf_cs$thin_cohorts, collapse = ", ")))
res_cs <- run_spec(cf_cs$panel, yname = "log_commits")
if (!is.null(res_cs))
  message(sprintf("  Common support: ATT = %.4f (SE %.4f, t %.3f) | N = %d",
                  res_cs$att, res_cs$se, res_cs$t, res_cs$n_obs))
battery_summary <- add_to_battery(battery_summary, "Conditioning set",
                                  "Common support (propensity in [0.05, 0.95])", res_cs)

cond_tab <- bind_rows(
  spec_row("Baseline: WGI GE + log population (headline specification)",
           list(att = fit_head$att, se = fit_head$se, t = fit_head$att / fit_head$se,
                stars = fmt_stars(fit_head$att / fit_head$se),
                pretrend = fit_head$pretrend, n_obs = fit_head$n_obs,
                n_country = fit_head$n_country,
                n_treated = n_distinct(did_panel_main$country_id[
                  did_panel_main$cohort_year > 0]))),
  spec_row("$+$ baseline (2009--2012) mean log adaptation commitments", res_base),
  spec_row("Common support: propensity $\\in [0.05, 0.95]$", res_cs))

write_tex_float(
  file.path(dir_tabs_cb, "att_conditioning.tex"),
  "Conditioning set and common support",
  "tab:att_conditioning",
  tabular_of(cond_tab, "tab:att_conditioning"),
  paste0("Outcome: log(adaptation commitments). Row 2 adds mean ", min(BASE_YEARS), "--",
         max(BASE_YEARS), " log commitments to the control set. Row 3 imposes common ",
         "support via a separate ever-adoption logit; propensity outside $[", PS_LO, ", ",
         PS_HI, "]$ dropped, cohort rule re-applied. ",
         "Baseline row from the stored headline fit. ", SPEC_NOTE),
  SRC_NOTE)

# --- Listwise-deletion losses for the specifications estimated here ----------
lw13 <- function(panel, yname, label, controls) {
  present <- intersect(controls, names(panel))
  miss_y <- is.na(panel[[yname]])
  miss_x <- Reduce(`|`, lapply(present, function(cc) is.na(panel[[cc]])))
  if (is.null(miss_x)) miss_x <- rep(FALSE, nrow(panel))
  keep <- !miss_y & !miss_x
  data.frame(
    Specification              = label,
    `Recipient-years (raw)`    = format(nrow(panel), big.mark = ","),
    `Dropped: outcome`         = format(sum(miss_y), big.mark = ","),
    `Dropped: controls`        = format(sum(!miss_y & miss_x), big.mark = ","),
    `Recipient-years (used)`   = format(sum(keep), big.mark = ","),
    `Recipients (raw)`         = as.character(n_distinct(panel$country_id)),
    `Recipients (used)`        = as.character(n_distinct(panel$country_id[keep])),
    check.names = FALSE, stringsAsFactors = FALSE)
}

lw_tab13 <- bind_rows(
  lw13(did_panel_main, "log_commits", "Main specification",
       c("ge_est", "log_population")),
  lw13(did_panel_base, "log_commits", "$+$ baseline outcome level",
       c("ge_est", "log_population", "base_outcome")),
  lw13(cf_cs$panel, "log_commits", "Common support (trimmed)",
       c("ge_est", "log_population")),
  lw13(did_panel_full, "log_commits", "All cohorts retained",
       c("ge_est", "log_population"))
)
print(lw_tab13, row.names = FALSE)

write_tex_float(
  file.path(dir_tabs_cb, "listwise_losses_13.tex"),
  "Listwise-deletion losses in the cohort, anticipation and conditioning specifications",
  "tab:listwise_losses_13",
  tabular_of(lw_tab13, "tab:listwise_losses_13"),
  paste0("Recipient-years and recipients lost to a missing outcome or a missing ",
         "control in each specification. \\texttt{did::att\\_gt} drops those ",
         "cells silently, so the ``used'' columns are the samples on which the ",
         "reported ATTs are computed. The baseline-outcome specification loses ",
         "additional observations for recipients with no usable ",
         min(BASE_YEARS), "--", max(BASE_YEARS), " window"),
  "OECD CRS (Rio adaptation markers); UNFCCC NAP Central; WGI; World Bank WDI")

# ==============================================================================
# SECTION 12. Combined cohort-battery summary exhibit
# One table the reader can scan: every perturbation of the cohort structure,
# the estimator and the conditioning set, against the headline.
# ==============================================================================

message("\n=== Combined battery summary ===\n")

battery_summary <- Filter(Negate(is.null), battery_summary)
bat_rows <- lapply(battery_summary, function(b) {
  r <- b$res
  data.frame(
    Block         = b$block,
    Specification = b$label,
    ATT           = if (is.null(r)) "---" else paste0(sprintf("%.4f", r$att), r$stars),
    SE            = if (is.null(r)) "---" else sprintf("(%.4f)", r$se),
    `$t$`         = if (is.null(r)) "---" else sprintf("%.3f", r$t),
    `$\\Delta$ vs headline` = if (is.null(r)) "---" else
      sprintf("%+.4f", r$att - fit_head$att),
    `Treated`     = if (is.null(r)) "---" else as.character(r$n_treated),
    `$N$`         = if (is.null(r)) "---" else format(r$n_obs, big.mark = ","),
    check.names = FALSE, stringsAsFactors = FALSE)
})

bat_tab <- bind_rows(
  data.frame(
    Block = "Headline", Specification = "Headline specification (stored fit)",
    ATT = paste0(sprintf("%.4f", fit_head$att),
                 fmt_stars(fit_head$att / fit_head$se)),
    SE = sprintf("(%.4f)", fit_head$se),
    `$t$` = sprintf("%.3f", fit_head$att / fit_head$se),
    `$\\Delta$ vs headline` = "---",
    Treated = as.character(n_distinct(did_panel_main$country_id[
      did_panel_main$cohort_year > 0])),
    `$N$` = format(fit_head$n_obs, big.mark = ","),
    check.names = FALSE, stringsAsFactors = FALSE),
  bind_rows(bat_rows))

for (r in seq_len(nrow(bat_tab)))
  message(sprintf("  %-28s %-46s ATT = %9s  SE = %9s",
                  bat_tab$Block[r], bat_tab$Specification[r],
                  bat_tab$ATT[r], bat_tab$SE[r]))

write_tex_float(
  file.path(dir_tabs_cb, "att_cohort_battery.tex"),
  "Cohort, timing and conditioning battery: every perturbation against the headline",
  "tab:att_cohort_battery",
  tabular_of(bat_tab, "tab:att_cohort_battery"),
  paste0("Outcome: log(adaptation commitments). Each row changes one thing from the headline ",
         "spec (Specification column); blocks correspond to Tables~\\ref{tab:att_loco}, ",
         "\\ref{tab:att_drop2024}, \\ref{tab:att_balance}, \\ref{tab:att_2x2}, ",
         "\\ref{tab:att_anticipation}, \\ref{tab:att_eventdate}, \\ref{tab:att_placebo_ladder}, ",
         "\\ref{tab:att_conditioning}. ``$\\Delta$ vs headline'' $=$ difference from the stored ",
         "headline ATT of ", sprintf("%.4f", fit_head$att),
         "; descriptive, not a test. ", SPEC_NOTE_NOPT),
  SRC_NOTE)

message("\n=== 13_cohort_anticipation.R: complete ===\n")
