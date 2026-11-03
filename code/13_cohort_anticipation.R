# ==============================================================================
# 13_cohort_anticipation.R
# Cohort-robustness checks, anticipation / event-date / placebo ladder,
# conditioning-set sensitivity, and the consolidated robustness table.
# Paper: Beaucoral, Goujon and Marchand (2026) — appendix to §5 (robustness)
#
# Four blocks of checks:
#   A   Cohort robustness: drop 2024 as a cohort and as a calendar year;
#       balanced event windows (aggte balance_e = 1, 2); the 2x2
#       {main, retained} x {dr, reg} grid
#       with a common SE method; and which cohorts contribute to each event
#       time.
#   B   Anticipation (1 and 2 years), an alternative event date built from the
#       NAP submission MONTH, and a placebo ladder at -3/-4/-5 years with the
#       minimum detectable effect of each placebo.
#   C   Conditioning set: baseline (2009-2012) outcome level in xformla, and a
#       common-support specification trimming recipients whose estimated
#       propensity to adopt lies outside [0.05, 0.95]; listwise-deletion losses
#       for each specification (stage log).
#   D   Timing (Section 13): an adopters-only design (never-treated dropped,
#       not-yet-treated controls) on the level, the global share and the
#       within-country share. Submission months by cohort are the lower panel
#       of the re-dating cross-tabulation (B).
#   E   tab:robust_summary: the robustness rows estimated in 04 (read from
#       output/fits/robust_rows_04.rds) and here (drop 2024, balanced event
#       windows, baseline outcome level), on log adaptation commitments, the
#       global share and the within-country share.
#
# Inputs :
#   data/processed/simple_panel_wgi.csv
#   output/fits/headline_{adaptation,total,nonadaptation}_dr_bs.rds
#                                               (written by 03_main_results.R)
#   output/fits/placebo_m2_reg_bs.rds           (written by 04_robustness.R)
#   output/fits/robust_rows_04.rds              (written by 04_robustness.R)
#   data/raw/shared_nap_data/nap_information.csv (submission dates; also carried
#                                                 in the processed panel)
# Outputs:
#   output/tables/cohort_battery/att_2x2.tex
#   output/tables/cohort_battery/cohort_contributions.tex
#   output/tables/anticipation/att_anticipation.tex
#   output/tables/anticipation/att_eventdate.tex
#   output/tables/anticipation/cohort_redating_crosstab.tex (with submission months)
#   output/tables/anticipation/att_placebo_ladder.tex
#   output/tables/timing/att_adopters_only.tex           (Section 13a)
#   output/tables/robustness/robust_summary.tex          (tab:robust_summary, Section 14)
#
# Also reads output/fits/headline_share_dr_bs.rds (03) and
# output/fits/within_share_dr_bs.rds (07) for the adopters-only comparison.
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
    "did %s is installed but the paper's results require did 2.5.0.\n",
    "  Older versions reproduce the ATTs but NOT the standard errors and\n",
    "  pre-trend tests on unbalanced panels (South Sudan has 14/16 years).\n",
    "  Fix: run renv::restore() from the project root (renv.lock pins did 2.5.0),\n",
    "  or remotes::install_version(\"did\", \"2.5.0\"), then rerun this script."),
    utils::packageVersion("did")))
} else if (utils::packageVersion("did") != "2.5.0") {
  warning(sprintf(paste0(
    "did %s is installed; the paper's results were produced with did 2.5.0 ",
    "(renv.lock). Later versions may change SEs: renv::restore() or ",
    "remotes::install_version(\"did\", \"2.5.0\")."), utils::packageVersion("did")),
    call. = FALSE)
}
library(parallel)   # forked evaluation guard for the segfault-prone dr cell
source(here("code", "functions", "pretrend_test.R"))      # compute_pretrend_test()
source(here("code", "functions", "read_headline_fit.R"))  # read_headline_fit()
source(here("code", "functions", "mde.R"))                # mde()
source(here("code", "functions", "make_country_id.R"))    # make_country_id()
source(here("code", "functions", "crs_positive.R"))       # crs_positive()
source(here("code", "functions", "table_notes.R"))       # shared table-note wording, wpval_reconciliation()
source(here("code", "functions", "within_share.R"))      # within_share()

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
dir_tabs_tm <- here("output", "tables", "timing")
dir.create(dir_tabs_tm, recursive = TRUE, showWarnings = FALSE)
dir_tabs_rs <- here("output", "tables", "robustness")
dir.create(dir_tabs_rs, recursive = TRUE, showWarnings = FALSE)
# A failed run must not leave the previous run's exhibits in place: they are
# deleted before anything is estimated, and an estimation failure stops the
# stage (the one documented exception is the forked 2x2 cell below).
unlink(list.files(c(dir_tabs_cb, dir_tabs_an, dir_tabs_tm, dir_tabs_rs), pattern = "\\.(tex|png|pdf)$",
                  full.names = TRUE))

#' Small counts as words, for table notes
num_word <- function(n) { stopifnot(n >= 1L, n <= 6L); c("one", "two", "three", "four", "five", "six")[n] }

`%||%` <- function(a, b) if (!is.null(a)) a else b

esc_header <- function(x) {
  x <- gsub("%", "\\\\%", x)
  x <- gsub("_", "\\\\_", x)
  x <- gsub("#", "\\\\#", x)
  x
}

# ==============================================================================
# Helper: write_tex_float()
# Emits \begin{table}[H] + \adjustbox around every tabular + notes minipage.
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
  # each panel within \textwidth. The wrapper lines are joined to their
  # tabular line with "\n", which writeLines() writes as separate lines.
  tab_starts <- grepl("^\\\\begin\\{tabular", inner)
  tab_ends   <- grepl("^\\\\end\\{tabular",   inner)
  if (any(tab_starts) && sum(tab_starts) == sum(tab_ends)) {
    inner[tab_starts] <- paste0("\\adjustbox{max width=\\textwidth}{%\n", inner[tab_starts])
    inner[tab_ends]   <- paste0(inner[tab_ends], "\n}")
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
  mutate(country_id = make_country_id(recipient_name))

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

did_panel_full <- did_panel %>%
  mutate(share_within = within_share(commitments, commitments_all))

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
message(paste(capture.output(print(cf$sizes, row.names = FALSE)), collapse = "\n"))
message(sprintf("Estimation cohorts (>= %d treated units): %s",
                THIN_THRESHOLD, paste(kept_cohorts, collapse = ", ")))
message(sprintf("Dropped thin cohorts: %s", paste(thin_cohorts, collapse = ", ")))
message(sprintf("Main panel: %d recipient-years, %d recipients, %d treated recipients",
                nrow(did_panel_main), n_distinct(did_panel_main$country_id),
                n_distinct(did_panel_main$country_id[did_panel_main$cohort_year > 0])))

# HARD CHECK: this panel must be the one behind Table 2. Reading the saved
# headline fit makes the comparison exact rather than a matter of comment.
fit_head <- read_headline_fit("adaptation")
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

#' A stored 03 headline fit in run_spec()'s result format (never re-estimated)
#'
#' @param fit stored fit from read_headline_fit()
#' @param panel the estimation panel it was fitted on
#' @return list with run_spec()'s fields used by the tables
headline_res <- function(fit, panel = did_panel_main) {
  stopifnot(identical(fit$ids, sort(unique(panel$country_id))))
  t_v <- fit$att / fit$se
  list(att = fit$att, se = fit$se, t = t_v, stars = fmt_stars(t_v),
       pretrend = fit$pretrend, agg_s = fit$agg_simple, agg_d = fit$agg_dynamic,
       gt = fit$gt_boot, n_obs = fit$n_obs, n_country = fit$n_country,
       n_treated = n_distinct(panel$country_id[panel$cohort_year > 0]),
       est_method = "dr", anticipation = 0)
}
res_head <- headline_res(fit_head)

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
#'   n_treated, gt); stops if an estimation step fails (run_spec_forked()
#'   returns NULL instead)
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
  # A failed fit stops here; run_spec_forked() turns that into a NULL result.
  gt_an <- tryCatch(fit_one(FALSE),
                    error = function(e) { message("    att_gt (analytical) failed: ",
                                                   conditionMessage(e)); NULL })
  set.seed(1242)   # Seed rule (reproduces the published SEs): seed immediately before the estimator
  gt <- tryCatch(fit_one(TRUE),
                 error = function(e) { message("    att_gt (bootstrap) failed: ",
                                                conditionMessage(e)); NULL })
  if (is.null(gt) || is.null(gt_an)) stop("run_spec(): att_gt() failed for ", yname)

  agg_s <- aggte(gt, type = "simple", na.rm = TRUE)
  agg_d <- if (is.null(balance_e))
    aggte(gt, type = "dynamic", na.rm = TRUE, min_e = min_e, max_e = max_e) else
    aggte(gt, type = "dynamic", na.rm = TRUE, min_e = min_e, max_e = max_e,
          balance_e = balance_e)
  # Pre-trend block. With anticipation = k, event times -1, ..., -k are TREATED
  # cells, not leads: including them in a "pre-trend" test would test the
  # anticipation effect itself. The restriction is applied INSIDE
  # compute_pretrend_test() via its `anticipation` argument, not by truncating
  # the aggregation: aggte(max_e = -k-1) errors out inside did because the
  # resulting aggregation has no post-treatment period ("influence function
  # column indices exceed available columns").
  agg_d_an <- aggte(gt_an, type = "dynamic", na.rm = TRUE, min_e = min_e, max_e = max_e)

  att <- agg_s$overall.att
  se  <- agg_s$overall.se
  t_v <- if (!is.na(att) && !is.na(se) && se > 0) att / se else NA_real_
  pt  <- compute_pretrend_test(agg_d_an, gt_an, anticipation = anticipation)

  # Event times that actually entered the pre-trend test, so the exhibit can
  # state the window rather than assert a formula.
  k_keep    <- which(!is.na(agg_d_an$se.egt) & agg_d_an$se.egt > 1e-10)
  pre_leads <- as.integer(agg_d_an$egt[k_keep][agg_d_an$egt[k_keep] < -anticipation])

  # Equal-weight average of the event-time effects e >= 0. did's SIMPLE
  # aggregate also uses only post-submission cells (t >= g) but weights them by
  # group size; the two differ only in weighting. It is computed last so that it
  # cannot disturb the bootstrap draws behind the reported ATT above; the seed
  # is set immediately before it, as for every other estimator call here.
  set.seed(1242)
  agg_e0 <- aggte(gt, type = "dynamic", na.rm = TRUE, min_e = 0, max_e = max_e)
  att_e0 <- agg_e0$overall.att
  se_e0  <- agg_e0$overall.se

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
#'
#' @param allow_null TRUE only for the forked 2x2 cell, whose failure the table
#'   note documents; any other failed estimate stops the stage.
spec_row <- function(label, res, allow_null = FALSE) {
  if (is.null(res)) {
    if (!allow_null) stop("Estimation failed for table row: ", label)
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

# Common preamble shared by every cohort-battery table note: the headline
# specification, referred to Table 2 (functions/table_notes.R).
spec_note_unless <- function(parts) paste0(
  "Unless noted: cohorts $\\geq$ ", THIN_THRESHOLD, "; ",
  as_in_main(parts, capitalize = FALSE))
# Tables without a "Pre-trend p" column use the preamble, which defines none.
SPEC_NOTE_PREAMBLE <- paste0(spec_note_unless(c("estimator", "covariates", "control group",
                                                "inference")), ".")
SRC_NOTE <- SOURCE_MAIN
# The headline test is Table 2's (same window and restriction count).
stopifnot(identical(as.integer(c(min(fit_head$pretrend$leads), max(fit_head$pretrend$leads),
                                 fit_head$pretrend$df)), unname(MAIN_PRETREND_WINDOW)))
SPEC_NOTE <- spec_note_unless(c("estimator", "covariates", "control group", "inference",
                                "pre-trend test", "significance stars"))

#' Stop unless every row's pre-trend test has Table 2's window
#'
#' SPEC_NOTE says the pre-trend test is "as in Table 2"; a table that uses it
#' must have that window (leads, restriction count) in every estimated row.
#'
#' @param rows list of list(label, res); rows with res = NULL are skipped
#' @param table table label, for the error message
assert_main_window <- function(rows, table) {
  for (r in rows) {
    pt <- r$res$pretrend
    if (is.null(pt)) next
    w <- as.integer(c(if (length(pt$leads) > 0L) range(pt$leads) else c(NA, NA), pt$df))
    if (!identical(w, unname(MAIN_PRETREND_WINDOW)))
      stop(table, ", row '", r$label, "': pre-trend window (min_e, max_e, df) = (",
           paste(w, collapse = ", "), ") differs from Table 2's (",
           paste(MAIN_PRETREND_WINDOW, collapse = ", "), "), but the note (SPEC_NOTE) ",
           "says the pre-trend test is as in Table 2; state the window with ",
           "pretrend_window_note().")
  }
  invisible(TRUE)
}


# ==============================================================================
# SECTION 4. A.2 — dropping 2024 as a cohort and as a calendar year
# The 2024 cohort contributes only the impact year e = 0, and 2024 is the last
# (and most provisional) CRS vintage year, so the two exclusions test different
# things: cohort composition versus data vintage.
# ==============================================================================

message("\n=== A.2: dropping 2024 (cohort / calendar year) ===\n")

# Dropping 2024 as a CALENDAR year leaves the 2024 cohort with no post-treatment
# period; att_gt() then has no estimable ATT(g,t) for g = 2024, so that cohort
# is removed as well. We drop it explicitly to keep the message honest.
panel_no2024_year <- did_panel_main %>%
  filter(year < 2024, cohort_year != 2024)
res_no2024_year   <- run_spec(panel_no2024_year, yname = "log_commits")


# ==============================================================================
# SECTION 5. A.3 — balanced event windows (balance_e = 1, 2)
# aggte(balance_e = k) restricts the dynamic aggregation to cohorts observed for
# at least k post-treatment periods, so the event-study path is not a
# composition effect. balance_e = 1 drops the 2024 cohort from the aggregation,
# balance_e = 2 also drops 2023.
# ==============================================================================

message("\n=== A.3: balanced event windows ===\n")

# balance_e = 1, 2 on log adaptation commitments (rows of tab:robust_summary;
# the share outcomes are added at the end of the script). The object of
# interest is the DYNAMIC average, which is what balance_e modifies.
balance_grid <- c(1L, 2L)
res_balance <- setNames(vector("list", length(balance_grid)), balance_grid)
for (be in balance_grid) {
  res_balance[[as.character(be)]] <- run_spec(did_panel_main, yname = "log_commits", balance_e = be)
  message(sprintf("-- %d post-adoption year(s): dynamic ATT = %.4f (SE %.4f)", be,
                  res_balance[[as.character(be)]]$agg_d$overall.att,
                  res_balance[[as.character(be)]]$agg_d$overall.se))
}

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
grid_rows <- vector("list", 12L)   # 3 outcomes x 2 cohort rules x 2 estimators, filled in order
grid_res  <- vector("list", 12L)   # list(label, res) per row, for assert_main_window()
i_grid    <- 0L
# The envelope outcomes are included because the retained-cohort specification
# (regression adjustment) is the only one in which they move; the grid shows
# whether that comes from the cohorts or the estimator.
grid_outcomes <- list(
  list(var = "log_commits",           tag = "Adaptation"),
  list(var = "lcommitments_all",      tag = "Total"),
  list(var = "lcommitments_nonadapt", tag = "Non-adaptation"))
# The main-rule DR cells are the headline specification: read from the fits 03
# stores (one fit, one SE), never re-estimated here.
grid_stored <- list(log_commits           = res_head,
                    lcommitments_all      = headline_res(read_headline_fit("total")),
                    lcommitments_nonadapt = headline_res(read_headline_fit("nonadaptation")))
grid_tags   <- vapply(grid_outcomes, `[[`, character(1L), "tag")
risky_fail  <- setNames(logical(length(grid_outcomes)), grid_tags)  # DR cell crashed
# Cohorts with no estimable post-treatment cell under DR, one slot per outcome.
dr_dropped_by <- setNames(vector("list", length(grid_outcomes)), grid_tags)

for (oc in grid_outcomes) {
  for (gp in grid_panels) {
    for (em in c("dr", "reg")) {
      lbl <- sprintf("%s: %s, %s", oc$tag, gp$label,
                     c(dr = "DR", reg = "OR")[[em]])
      message(sprintf("-- %s x %s x %s", oc$var, gp$key, em))
      risky <- identical(gp$key, "retained") && identical(em, "dr")
      res <- if (identical(gp$key, "main") && identical(em, "dr")) {
        grid_stored[[oc$var]]
      } else if (risky) {
        message("   (running in a forked process: att_gt(dr) can segfault on ",
                "panels containing 2-unit cohorts)")
        run_spec_forked(gp$panel, yname = oc$var, est_method = em)
      } else {
        run_spec(gp$panel, yname = oc$var, est_method = em)
      }
      if (risky && is.null(res)) risky_fail[[oc$tag]] <- TRUE
      if (risky && !is.null(res)) {
        rep_c <- cohort_cell_report(res$gt, gp$panel)
        dr_dropped_by[[oc$tag]] <- rep_c$cohort[!rep_c$contributes]
      }
      i_grid <- i_grid + 1L
      grid_rows[[i_grid]] <- spec_row(lbl, res, allow_null = risky)
      grid_res[[i_grid]]  <- list(label = lbl, res = res)
      if (!is.null(res))
        message(sprintf("   ATT = %.4f (SE %.4f, t %.3f) | N = %d | treated = %d",
                        res$att, res$se, res$t, res$n_obs, res$n_treated))
    }
  }
}
dr_dropped <- unique(unlist(dr_dropped_by))
grid_notes_extra <- vapply(names(risky_fail)[risky_fail], function(tag)
  paste0("The DR estimator could not be evaluated on the retained-cohort ",
         "panel for ", tag, " (it crashes on 2-unit cohorts); that ",
         "cell is reported as unavailable."), character(1L))
dr_dropped_txt <- if (length(dr_dropped) == 0L) "" else {
  sizes <- did_panel_full %>% filter(cohort_year %in% dr_dropped) %>%
    distinct(country_id, cohort_year) %>% count(cohort_year)
  paste0("In the DR cells with all cohorts retained, the ",
         paste(sprintf("%d cohort (%d recipients)", sizes$cohort_year, sizes$n), collapse = " and "),
         " have no estimable post-treatment cell and drop out of the aggregate; OR ",
         "estimates them. ``Treated'' counts the recipients in the estimation panel.")
}

assert_main_window(grid_res, "tab:att_2x2")
write_tex_float(
  file.path(dir_tabs_cb, "att_2x2.tex"),
  "Cohort rule and estimator, varied one at a time",
  "tab:att_2x2",
  tabular_of(bind_rows(grid_rows), "tab:att_2x2"),
  paste0("Outcomes: log adaptation, total and non-adaptation commitments. This grid varies ",
         "the cohort rule and estimator separately, with common multiplier-bootstrap inference ",
         "in all cells. ``All cohorts retained'' adds the 2015--2020 cohorts of 2--4 recipients ",
         "each. ", dr_dropped_txt, " ", paste(grid_notes_extra, collapse = " "), " ", SPEC_NOTE),
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
message(paste(capture.output(print(contrib_tab, row.names = FALSE)), collapse = "\n"))

write_tex_float(
  file.path(dir_tabs_cb, "cohort_contributions.tex"),
  "Cohorts contributing to each event time in the main event study",
  "tab:cohort_contributions",
  tabular_of(contrib_tab, "tab:cohort_contributions"),
  paste0("Composition of the dynamic aggregation used by Figure~\\ref{fig:did_combined_es} ",
         "and the HonestDiD analysis. A cohort contributes to event time $e$ when its ",
         "$ATT(g, g+e)$ cell is estimable in the main-specification estimates (non-missing, ",
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
antic_lbl  <- function(k) sprintf("Anticipation = %d year%s", k, if (k == 1L) "" else "s")
antic_lbls <- vapply(antic_grid, antic_lbl, character(1L))
# Pre-sized with the row labels as names (filled by label below).
antic_rows  <- setNames(vector("list", length(antic_grid)), antic_lbls)
antic_dyn   <- setNames(vector("list", length(antic_grid)), antic_lbls)
antic_e0    <- setNames(vector("list", length(antic_grid)), antic_lbls)
antic_leads <- setNames(vector("list", length(antic_grid)), antic_lbls)
for (k in antic_grid) {
  lbl <- antic_lbl(k)
  message("-- ", lbl)
  res <- if (k == 0L) NULL else
    run_spec(did_panel_main, yname = "log_commits", anticipation = k,
             min_e = -5 - k)
  if (k == 0L) {
    # Baseline row read from the stored headline fit. Its e >= 0 average
    # is the same object as the pooled ATT up to the group-size weighting did
    # applies, and is reported for comparability with k = 1, 2.
    base_row <- spec_row(
      "Anticipation = 0 (baseline)",
      res_head)
    # The stored dynamic aggregation (min_e = -5) has the same e >= 0 average;
    # re-running aggte() would draw a fresh bootstrap and print a second SE
    # for a statistic the balance table already reports.
    agg_e0_base <- fit_head$agg_dynamic
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

# The pooled column is did's SIMPLE aggregate: post-submission cells (t >= g)
# weighted by group size. The e >= 0 column averages the event-time effects
# with equal weights; the two differ only in weighting.
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

# Sample-size columns identical in every row are stated once in the note, and
# long headers break over two lines, so the panel prints without shrinking.
antic_const <- c(`$N$` = "observations", Recipients = "recipients", Treated = "treated")
is_const    <- vapply(names(antic_const), function(cc) length(unique(antic_tab[[cc]])) == 1L,
                      logical(1L))
antic_const_note <- if (!any(is_const)) "" else paste0(
  "Every row of the upper panel: ",
  paste(unlist(antic_tab[1L, names(antic_const)[is_const]]), antic_const[is_const],
        collapse = ", "), ". ")
antic_tab <- antic_tab[, !names(antic_tab) %in% names(antic_const)[is_const], drop = FALSE]
antic_heads <- c(`Pre-trend $p$`    = "\\makecell{Pre-trend \\\\ $p$}",
                 `ATT, $e \\geq 0$` = "\\makecell{ATT, \\\\ $e \\geq 0$}",
                 `Pre-trend leads`  = "\\makecell{Pre-trend \\\\ leads}")
stopifnot(all(names(antic_heads) %in% names(antic_tab)))
names(antic_tab)[match(names(antic_heads), names(antic_tab))] <- antic_heads

write_tex_float(
  file.path(dir_tabs_an, "att_anticipation.tex"),
  "Anticipation: pooled and dynamic ATTs allowing donors to respond before submission",
  "tab:att_anticipation",
  c(tabular_of(antic_tab, "tab:att_anticipation"),
    "\\par\\vspace{6pt}",
    tabular_of(dyn_tab, "tab:att_anticipation_dyn")),
  paste0("Upper panel: pooled ATT and the $e \\geq 0$ average, the equal-weight mean of the dynamic ATTs at $e \\geq 0$. Lower panel: ",
         "dynamic ATTs by event time, SE in parentheses. Allowing ",
         paste(vapply(antic_grid[antic_grid > 0L], num_word, character(1L)), collapse = " (or "),
         strrep(")", sum(antic_grid > 0L) - 1L), " years of anticipation shifts the base period to ",
         paste(sprintf("$g-%d$", antic_grid[antic_grid > 0L] + 1L), collapse = " ("),
         strrep(")", sum(antic_grid > 0L) - 1L), ". $^{\\dagger}$ marks event times treated as ",
         "post-treatment under that row's anticipation, not as leads. The anticipation-0 row is the ",
         "main specification (Table~\\ref{tab:combined_wide_main}). ", antic_const_note,
         SPEC_NOTE_PREAMBLE, " Pre-trend $\\chi^2$ / $p$: joint Wald test on the aggregated ",
         "pre-treatment event-time coefficients in each row's \\emph{Pre-trend leads} window, ",
         "from the analytical (non-bootstrap) fit of the same specification. ", STAR_NOTE),
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
# WHY THERE ARE THREE ROWS. Re-dating is not a clean
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
# returns NA under a non-English LC_TIME (e.g. fr_FR).
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
# Disclosure when row (a2) falls back to outcome regression (see (a2) below):
# a label suffix for its rows and one sentence for the notes of the two tables
# that report it. Both stay empty when the doubly robust fit is produced.
a2_label_suffix <- ""
a2_note <- character(0)

if (!cohort_month_ok) {
  stop("No parseable NAP submission month: the event-date table cannot be built.")
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
  message(paste(capture.output(print(as.data.frame(month_counts), row.names = FALSE)), collapse = "\n"))

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
  message(paste(capture.output(print(as.data.frame(sizes_a), row.names = FALSE)), collapse = "\n"))

  cf_a1 <- apply_cohort_filter(panel_a)
  lost_a1 <- setdiff(sizes_a$cohort_year, cf_a1$kept_cohorts)
  message(sprintf(paste0("  (a1) cohorts entering: %s | falling below the >= %d ",
                         "rule and dropped: %s"),
                  paste(cf_a1$kept_cohorts, collapse = ", "), THIN_THRESHOLD,
                  if (length(lost_a1) == 0L) "none" else paste(lost_a1, collapse = ", ")))
  res_month_a1 <- run_spec(cf_a1$panel, yname = "log_commits")
  rep_a1 <- cohort_cell_report(res_month_a1$gt, cf_a1$panel)
  if (!is.null(rep_a1)) { message("  (a1) estimable cells per cohort:")
    message(paste(capture.output(print(rep_a1, row.names = FALSE)), collapse = "\n")) }
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
    a2_label_suffix <- " (OR)"
    a2_note <- paste0("Row (a2) uses OR (multiplier-bootstrap SE): the DR ",
                      "fit could not be produced on this panel.")
  }
  rep_a2 <- cohort_cell_report(res_month_a2$gt, panel_a)
  if (!is.null(rep_a2)) { message("  (a2) estimable cells per cohort:")
    message(paste(capture.output(print(rep_a2, row.names = FALSE)), collapse = "\n")) }
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
    message(paste(capture.output(print(rep_b, row.names = FALSE)), collapse = "\n")) }
  if (!is.null(res_month_b))
    message(sprintf(paste0("   (b) ATT = %.4f (SE %.4f, t %.3f) | treated in ",
                           "panel = %d | CONTRIBUTING treated = %d"),
                    res_month_b$att, res_month_b$se, res_month_b$t,
                    res_month_b$n_treated,
                    sum(rep_b$treated_in_panel[rep_b$contributes])))


  # ---- (c) original x re-dated cross-tabulation exhibit ----------------------
  ct <- redate_map %>%
    filter(cohort_year > 0) %>%
    mutate(g_lab = if_else(beyond_panel, "Beyond panel", as.character(g_month)))
  g_levels <- c(sort(unique(ct$g_month[!ct$beyond_panel])), NA)
  g_labs   <- c(as.character(sort(unique(ct$g_month[!ct$beyond_panel]))), "Beyond panel")
  orig_levels <- sort(unique(ct$cohort_year))

  ct_mat <- matrix("0", nrow = length(orig_levels), ncol = length(g_labs) + 2L)
  colnames(ct_mat) <- c(g_labs, "Total", "In main specification")
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
  message(paste(capture.output(print(ct_tab, row.names = FALSE)), collapse = "\n"))

  # ---- (d) submission months by cohort (lower panel of the same exhibit) -----
  sub_dates <- did_panel_main %>%
    filter(cohort_year > 0) %>%
    distinct(recipient_name, cohort_year, date_posted_2) %>%
    mutate(date = as.Date(date_posted_2), month = as.integer(format(date, "%m")))
  stopifnot(!anyNA(sub_dates$date), !anyDuplicated(sub_dates$recipient_name),
            all(as.integer(format(sub_dates$date, "%Y")) == sub_dates$cohort_year))
  cohorts_tm <- sort(unique(sub_dates$cohort_year))
  month_rows <- vapply(cohorts_tm, function(g) {
    m <- as.integer(table(factor(sub_dates$month[sub_dates$cohort_year == g], levels = 1:12)))
    paste0(g, " & ", sum(m), " & ", paste(m, collapse = " & "), " & ", sum(m[10:12]), " \\\\")
  }, character(1L))
  # COP26 window (31 Oct - 13 Nov 2021) as given in the revision plan; no other
  # conference dates are hard-coded.
  cop26 <- as.Date(c("2021-10-31", "2021-11-13"))
  c21 <- sub_dates %>% filter(cohort_year == 2021) %>% arrange(date, recipient_name)
  in_cop26 <- c21$recipient_name[c21$date >= cop26[1L] & c21$date <= cop26[2L]]
  fmt_day <- function(d) paste(as.integer(format(d, "%d")), month.abb[as.integer(format(d, "%m"))])
  message(sprintf("  2021 cohort (%d): %s", nrow(c21),
                  paste(sprintf("%s %s", c21$recipient_name, format(c21$date)), collapse = "; ")))
  message(sprintf("  2021 cohort posted Oct-Dec 2021: %d; during COP26 (31 Oct-13 Nov): %d (%s)",
                  sum(c21$month >= 10L), length(in_cop26), paste(in_cop26, collapse = ", ")))
  months_lines <- c(
    "\\begin{tabular}{lcccccccccccccc}",
    "\\toprule",
    paste0("Cohort & $N$ & ", paste(month.abb, collapse = " & "), " & Oct--Dec \\\\"),
    "\\midrule",
    month_rows,
    "\\bottomrule", "\\end{tabular}"
  )

  write_tex_float(
    file.path(dir_tabs_an, "cohort_redating_crosstab.tex"),
    "Re-dating by submission month: re-dated adoption cohorts and posting months",
    "tab:cohort_redating",
    c(tabular_of(ct_tab, "tab:cohort_redating"),
      "\\par\\vspace{6pt}",
      months_lines),
    paste0("Upper panel: treated recipients per (original cohort, re-dated cohort) cell. ",
           "Re-dated cohort ",
           "$=$ submission year $+1$ if posted after June. ``Beyond panel'': re-dated cohort ",
           "falls in ", last_year + 1L, ", outside the window; dropped, not recoded as ",
           "controls. ``In main specification'': whether the original cohort meets ",
           "the $\\geq ", THIN_THRESHOLD, "$ rule. ", n_adopters_all,
           " adopters with a recorded date shown; ", n_moved, " change cohort. ",
           "Lower panel: main-sample adopters (", nrow(sub_dates), " recipients) by the ",
           "calendar month in which UNFCCC NAP Central posted their NAP, the event date of ",
           "the main specification (adoption year) and of the re-dating. ",
           "In the 2021 cohort, ", sum(c21$month >= 10L), " of ", nrow(c21),
           " NAPs were posted in October--December 2021 and ", length(in_cop26),
           " during COP26 (31 October--13 November 2021", if (length(in_cop26) > 0L)
             paste0(": ", paste(in_cop26, collapse = ", ")) else "", "). ",
           "2021 cohort posting dates (2021): ",
           paste(sprintf("%s (%s)", c21$recipient_name, fmt_day(c21$date)), collapse = ", ")),
    SRC_NOTE)
}

# Pre-sized result container; filled by index, trimmed before binding.
eventdate_rows <- vector("list", 4L)
eventdate_rows[[1L]] <- (
  spec_row("\\makecell[l]{Baseline: dated by \\\\ submission year}",
           res_head)
)
if (cohort_month_ok) {
  eventdate_rows[[2L]] <- spec_row(
    "\\makecell[l]{(a1) Main treated set, \\\\ $\\geq 5$ rule re-applied}", res_month_a1)
  eventdate_rows[[3L]] <- spec_row(
    paste0("\\makecell[l]{(a2) Main treated set, \\\\ all cohorts retained", a2_label_suffix, "}"),
    res_month_a2)
  eventdate_rows[[4L]] <- spec_row(
    "\\makecell[l]{(b) Full adopter sample, \\\\ $\\geq 5$ rule re-applied}", res_month_b)
}

assert_main_window(list(list(label = "baseline", res = res_head),
                        list(label = "(a1)", res = res_month_a1),
                        list(label = "(a2)", res = res_month_a2),
                        list(label = "(b)",  res = res_month_b)), "tab:att_eventdate")
eventdate_note <- paste0(
  "Outcome: log(adaptation commitments). Rows re-date treatment by posting month (post-June ",
  "$\\to$ next year); see Table~\\ref{tab:cohort_redating}. (a1)/(a2): main treated set, ",
  "$\\geq ", THIN_THRESHOLD, "$ rule re-applied vs.\\ retained. (b): full sample re-dated. ",
  "Contributing/Treated: estimable-cell counts. ",
  if (length(a2_note) == 0L) "" else paste0(a2_note, " "),
  SPEC_NOTE)

eventdate_tab <- bind_rows(eventdate_rows)

# Cohort set and CONTRIBUTING treated count per row:
# the "Treated" column of spec_row() counts recipients in the panel, which
# overstates what a row rests on when a cohort contributes no estimable cell.
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
# One cohort report per table row, in row order (built once, not grown).
ed_reps <- c(list(rep_base), if (cohort_month_ok) list(rep_a1, rep_a2, rep_b))
stopifnot(length(ed_reps) == nrow(eventdate_tab))
eventdate_tab$`Contributing cohorts` <- vapply(ed_reps, fmt_cohorts, character(1L))
eventdate_tab$Treated <- paste0(vapply(ed_reps, fmt_contrib, character(1L)), " / ",
                                eventdate_tab$Treated)
# Narrow layout: cohort lists over two lines (first half of the cohorts on the
# first), two-line headers.
eventdate_tab$`Contributing cohorts` <- vapply(eventdate_tab$`Contributing cohorts`, function(x) {
  cc <- strsplit(x, ", ", fixed = TRUE)[[1L]]
  if (length(cc) < 3L) return(x)
  k <- ceiling(length(cc) / 2)
  paste0("\\makecell{", paste(cc[seq_len(k)], collapse = ", "), ", \\\\ ",
         paste(cc[-seq_len(k)], collapse = ", "), "}")
}, character(1L), USE.NAMES = FALSE)
ed_heads <- c(`Pre-trend $p$`        = "\\makecell{Pre-trend \\\\ $p$}",
              `Contributing cohorts` = "\\makecell{Contributing \\\\ cohorts}")
stopifnot(all(names(ed_heads) %in% names(eventdate_tab)))
names(eventdate_tab)[match(names(ed_heads), names(eventdate_tab))] <- ed_heads

write_tex_float(
  file.path(dir_tabs_an, "att_eventdate.tex"),
  "Alternative treatment dates: re-dating adoption by submission month",
  "tab:att_eventdate",
  tabular_of(eventdate_tab, "tab:att_eventdate"),
  eventdate_note,
  SRC_NOTE)

# ==============================================================================
# SECTION 10. B.3 — placebo ladder at -3, -4 and -5 years
#
# The construction MIRRORS 04_robustness.R §25 exactly; the filters are stated
# here so that the text can be checked against them:
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

# Rung -2 is the placebo of Figure fig:placebo (04). It is READ from the
# stored fit rather than re-estimated, so the ladder cannot print a second
# standard error for a specification the paper already reports.
path_plac2 <- file.path(FITS_DIR, "placebo_m2_reg_bs.rds")
if (file.exists(path_plac2)) {
  fit_p2 <- readRDS(path_plac2)
  t_p2 <- fit_p2$att / fit_p2$se
  placebo_rows[[1L]] <- data.frame(
    Placebo        = "$-2$ years",
    ATT            = paste0(sprintf("%.4f", fit_p2$att), fmt_stars(t_p2)),
    SE             = sprintf("(%.4f)", fit_p2$se),
    `$t$`          = sprintf("%.3f", t_p2),
    `MDE (80\\% power)` = sprintf("%.4f", mde(fit_p2$se)),
    `Cohorts`      = paste(fit_p2$cohorts, collapse = ", "),
    `Treated`      = as.character(fit_p2$n_treated),
    `$N$`          = format(fit_p2$n_obs, big.mark = ","),
    check.names = FALSE, stringsAsFactors = FALSE)
  message(sprintf("-- -2 years (from stored fit): ATT = %.4f (SE %.4f)",
                  fit_p2$att, fit_p2$se))
} else {
  stop("Stored -2-year placebo fit not found: ", path_plac2,
       " (it is written by 04_robustness.R; run_all.R runs 04 before 13).")
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
  if (length(coh) < 2L) stop("Placebo -", k, " years: fewer than two placebo cohorts survive.")
  res <- run_spec(pan, yname = "log_commits", est_method = "reg", min_e = -4)
  if (is.null(res)) stop("Placebo -", k, " years: estimation failed.")
  placebo_rows[[k - 1L]] <- if (is.null(res)) data.frame(
    Placebo = sprintf("$-%d$ years", k), ATT = "---", SE = "---", `$t$` = "---",
    `MDE (80\\% power)` = "---", Cohorts = paste(coh, collapse = ", "),
    Treated = as.character(n_tr), `$N$` = format(nrow(pan), big.mark = ","),
    check.names = FALSE, stringsAsFactors = FALSE) else data.frame(
    Placebo             = sprintf("$-%d$ years", k),
    ATT                 = paste0(sprintf("%.4f", res$att), res$stars),
    SE                  = sprintf("(%.4f)", res$se),
    `$t$`               = sprintf("%.3f", res$t),
    `MDE (80\\% power)` = sprintf("%.4f", mde(res$se)),
    Cohorts             = paste(coh, collapse = ", "),
    Treated             = as.character(res$n_treated),
    `$N$`               = format(res$n_obs, big.mark = ","),
    check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(res))
    message(sprintf("   ATT = %.4f (SE %.4f, t %.3f) | MDE = %.4f",
                    res$att, res$se, res$t, mde(res$se)))
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
         "bootstrap SE (", BITERS, " reps, seed 1242). MDE $= ", sprintf("%.4f", mde(1)), " \\times$ SE (80\\% power, ",
         "two-sided 5\\%). The $-2$ row is the placebo plotted in Figure~\\ref{fig:placebo}. ", STAR_NOTE),
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
# Not added to the combined battery: the support restriction is not binding,
# so the estimate is the headline; the log above records it.

# --- Listwise-deletion losses for the specifications estimated here ----------
# Printed to the stage log (the text states them).
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
  lw13(did_panel_base, "log_commits", "+ baseline outcome level",
       c("ge_est", "log_population", "base_outcome")),
  lw13(cf_cs$panel, "log_commits", "Common support (trimmed)",
       c("ge_est", "log_population")),
  lw13(did_panel_full, "log_commits", "All cohorts retained",
       c("ge_est", "log_population"))
)
for (r in seq_len(nrow(lw_tab13)))
  message(sprintf("  Listwise losses, %-28s raw %6s -> used %6s  (outcome %s, controls %s; recipients %s -> %s)",
                  lw_tab13$Specification[r], lw_tab13$`Recipient-years (raw)`[r],
                  lw_tab13$`Recipient-years (used)`[r], lw_tab13$`Dropped: outcome`[r],
                  lw_tab13$`Dropped: controls`[r], lw_tab13$`Recipients (raw)`[r],
                  lw_tab13$`Recipients (used)`[r]))

# ==============================================================================
# SECTION 13. Timing investigation
# (a) Adopters-only design: the 40 main-sample adopters only, never-treated
#     recipients dropped, not-yet-treated controls; otherwise the headline
#     specification (WGI GE + log population, universal base period,
#     multiplier bootstrap, seed 1242). With no never-treated unit, did uses
#     the last cohort (2024) as the comparison group and drops calendar years
#     from 2024 on, so ATT(g, t) is identified for the 2021-2023 cohorts up to
#     2023 only (read from the fit below, not assumed). Estimator by the
#     project rule: doubly robust only with >= 40 treated units in the
#     estimated cohorts, else outcome regression; the other estimator is
#     reported beside it. Pre-trend test: the four-lead full-covariance Wald of
#     compute_pretrend_test() on an analytical twin.
# (b) Submission months by cohort, from the NAP Central posting dates carried
#     in the panel (data/raw/shared_nap_data/nap_information.csv, parsed by 01).
# ==============================================================================

message("\n=== Section 13: timing investigation ===\n")

panel_adopt <- did_panel_main %>%
  filter(cohort_year > 0) %>%
  mutate(share_within = if_else(crs_positive(commitments_all),
                                100 * commitments / commitments_all, NA_real_))
n_adopt <- n_distinct(panel_adopt$country_id)
stopifnot(n_adopt == n_distinct(did_panel_main$country_id[did_panel_main$cohort_year > 0]))

path_within <- file.path(FITS_DIR, "within_share_dr_bs.rds")
if (!file.exists(path_within))
  stop("Stored within-share fit not found: ", path_within, " (written by 07; run 07 before 13).")
adopt_outcomes <- list(
  list(var = "log_commits",  label = "log(Adaptation commitments)",
       main_att = fit_head$att),
  list(var = "share_adapt",  label = "Adaptation share (\\% of global)",
       main_att = read_headline_fit("share")$att),
  list(var = "share_within", label = "Within-country adaptation share (\\%)",
       main_att = readRDS(path_within)$att)
)

#' Adopters-only CS fit for one outcome and estimator: bootstrap fit for the
#' ATT/SE (aggregated inside the same seed), analytical twin for the pre-trend.
fit_adopters_only <- function(y, em) {
  args <- list(yname = y, tname = "year", idname = "country_id", gname = "cohort_year",
               xformla = ~ ge_est + log_population, data = panel_adopt, est_method = em,
               cband = FALSE, control_group = "notyettreated", anticipation = 0,
               base_period = "universal", panel = TRUE, allow_unbalanced_panel = TRUE)
  gt_an  <- suppressMessages(do.call(att_gt, c(args, list(bstrap = FALSE))))
  agg_da <- aggte(gt_an, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf)
  pt     <- compute_pretrend_test(agg_da, gt_an)
  bt <- withr::with_seed(1242L, {
    gt_b <- suppressMessages(do.call(att_gt, c(args, list(bstrap = TRUE, biters = BITERS))))
    list(gt = gt_b, s = aggte(gt_b, type = "simple", na.rm = TRUE),
         d = aggte(gt_b, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf))
  })
  groups <- sort(unique(bt$gt$group))
  att <- bt$s$overall.att; se <- bt$s$overall.se
  if (!is.finite(att) || !is.finite(se) || !is.finite(pt$stat))
    stop("Adopters-only fit failed for ", y, " (", em, ")")
  list(var = y, em = em, att = att, se = se, p = 2 * pnorm(-abs(att / se)),
       pt = pt, groups = groups, years = range(bt$gt$t),
       n_treated = n_distinct(panel_adopt$country_id[panel_adopt$cohort_year %in% groups]),
       n_obs = sum(!is.na(panel_adopt[[y]]) & panel_adopt$year <= max(bt$gt$t)),
       egt = bt$d$egt, att_egt = bt$d$att.egt, se_egt = bt$d$se.egt)
}

adopt_fits <- lapply(adopt_outcomes, function(oc)
  lapply(c(reg = "reg", dr = "dr"), function(em) fit_adopters_only(oc$var, em)))
names(adopt_fits) <- vapply(adopt_outcomes, `[[`, character(1L), "var")

f0 <- adopt_fits[[1L]]$reg
em_rule <- if (f0$n_treated >= 40L) "dr" else "reg"
control_only <- setdiff(sort(unique(panel_adopt$cohort_year)), f0$groups)
message(sprintf(paste0("  Adopters-only: %d adopters; estimated cohorts %s (%d treated units); ",
                       "comparison-only cohort(s) %s; years %d-%d; rule estimator: %s"),
                n_adopt, paste(f0$groups, collapse = ", "), f0$n_treated,
                paste(control_only, collapse = ", "), f0$years[1L], f0$years[2L], em_rule))
for (v in names(adopt_fits)) for (em in c("reg", "dr")) {
  f <- adopt_fits[[v]][[em]]
  message(sprintf(paste0("  %-12s %-3s ATT = %.4f (SE %.4f, p = %.3f) | pre-trend chi2(%d) = %.3f, ",
                         "p = %.3f | dynamic: %s"),
                  v, em, f$att, f$se, f$p, f$pt$df, f$pt$stat, f$pt$pval,
                  paste(sprintf("e=%d %.4f (%.4f)", as.integer(f$egt), f$att_egt, f$se_egt),
                        collapse = "; ")))
}

em_lab <- c(reg = "OR", dr = "DR")
adopt_rows <- unlist(lapply(seq_along(adopt_outcomes), function(k) {
  oc <- adopt_outcomes[[k]]
  ems <- c(em_rule, setdiff(c("reg", "dr"), em_rule))
  vapply(seq_along(ems), function(j) {
    f <- adopt_fits[[oc$var]][[ems[j]]]
    paste0(if (j == 1L) oc$label else "", " & ", em_lab[[ems[j]]], " & ",
           sprintf("%.4f", f$att), fmt_stars(f$att / f$se), " & (", sprintf("%.4f", f$se), ") & ",
           sprintf("%.3f", f$p), " & ", sprintf("%.3f", f$pt$stat), " & ",
           sprintf("%.3f", f$pt$pval), " & ", f$n_treated, " & ",
           format(f$n_obs, big.mark = ","), " & ",
           if (j == 1L) sprintf("%.4f", oc$main_att) else "", " \\\\")
  }, character(1L))
}))
pt_leads <- f0$pt$leads
adopt_lines <- c(
  "\\begin{tabular}{llcccccccc}",
  "\\toprule",
  paste0("Outcome & Estimator & ATT & SE & $p$ & Pre-trend $\\chi^2(", f0$pt$df, ")$ & ",
         "Pre-trend $p$ & Treated & $N$ & Main ATT \\\\"),
  "\\midrule",
  adopt_rows,
  "\\bottomrule", "\\end{tabular}"
)
post_e <- sort(unique(as.integer(f0$egt[f0$egt >= 0])))
write_tex_float(
  file.path(dir_tabs_tm, "att_adopters_only.tex"),
  "Adopters-only design: not-yet-treated comparisons among the NAP adopters",
  "tab:att_adopters_only",
  adopt_lines,
  paste0("CS\\,(2021) on the ", n_adopt, " main-sample adopters only (cohorts ",
         paste(range(panel_adopt$cohort_year), collapse = "--"), "); never-treated recipients ",
         "dropped; not-yet-treated controls; WGI GE + log population; universal base period ",
         "($g-1$). With no never-treated unit, \\texttt{did} uses the ",
         paste(control_only, collapse = ", "), " cohort only as a comparison group and drops ",
         "calendar years after ", f0$years[2L], ", so the estimated cohorts are ",
         paste(f0$groups, collapse = ", "), " (``Treated'' units) and the post-treatment ",
         "event times are $e=", min(post_e), ",\\dots,", max(post_e), "$. ",
         "Estimator: the first row of each outcome uses ", em_lab[[em_rule]],
         ", which our rule selects (DR only with at least 40 treated units); the ",
         "second row reports ", em_lab[[setdiff(c("reg", "dr"), em_rule)]],
         " for comparison. ATT/SE: multiplier bootstrap (", BITERS, " reps, clustered by recipient, ",
         "seed 1242). ", as_in_main("pre-trend test"), ". ",
         pretrend_window_note(min(pt_leads), max(pt_leads), f0$pt$df),
         "Main ATT: the main specification (never-treated controls, all ",
         n_distinct(did_panel_main$country_id), " recipients)",
         if (any(grepl("\\*", adopt_rows))) paste0(". ", STAR_NOTE) else ""),
  SRC_NOTE)


# ==============================================================================
# SECTION 14. tab:robust_summary -- one table for the sample, conditioning,
# estimator and outcome-unit checks, on log adaptation commitments, the global
# share and the within-country share (Panel A), and on log total and
# non-adaptation commitments where they are estimated (Panel B).
# Rows estimated in 04 are read from output/fits/robust_rows_04.rds (cells as
# formatted there); the main row is read from the fits 03 and 07 store; the
# drop-2024, balanced-window and baseline-level rows are estimated here. These
# fits run last in the script, so they move no draw of an earlier exhibit.
# ==============================================================================

message("\n=== Section 14: tab:robust_summary ===\n")

path_rr <- file.path(FITS_DIR, "robust_rows_04.rds")
if (!file.exists(path_rr))
  stop("Robustness rows not found: ", path_rr, " (written by 04_robustness.R; run_all.R runs 04 before 13).")
rr04 <- readRDS(path_rr)

rs_vars_a <- c("log_commits", "share_adapt", "share_within")
rs_vars_b <- c("lcommitments_all", "lcommitments_nonadapt")

#' One table cell from an estimate, its SE and a pre-trend p (NA -> "---")
rs_cell <- function(att, se, p = NA_real_) {
  list(att_fmt = paste0(sprintf("%.4f", att), fmt_stars(att / se)),
       se_fmt = sprintf("(%.4f)", se),
       p_fmt = if (is.na(p)) "---" else sprintf("%.3f", p), att = att, se = se)
}
#' Cell from a run_spec() result: simple ATT, or the dynamic average for the
#' balanced-window rows (whose pre-trend test is the main fit's: "---")
rs_cell_res <- function(res, dynamic = FALSE) {
  if (dynamic) rs_cell(res$agg_d$overall.att, res$agg_d$overall.se)
  else rs_cell(res$att, res$se, res$pretrend$pval)
}

# Main specification: the stored fits behind Table 2 and tab:within_share_wide.
fit_within_main <- readRDS(file.path(FITS_DIR, "within_share_dr_bs.rds"))
main_fits <- list(log_commits = fit_head, share_adapt = read_headline_fit("share"),
                  share_within = fit_within_main,
                  lcommitments_all = read_headline_fit("total"),
                  lcommitments_nonadapt = read_headline_fit("nonadaptation"))
rs_main <- lapply(main_fits, function(f) rs_cell(f$att, f$se, f$pretrend$pval))
# Unbalanced dynamic average (all post-adoption years, equal weights over
# e >= 0) of the same stored fits: the reference of the balanced-window rows.
rs_unbal <- lapply(main_fits[rs_vars_a], function(f)
  rs_cell(f$agg_dynamic$overall.att, f$agg_dynamic$overall.se))
stopifnot("Unbalanced dynamic average differs from tab:att_anticipation's ATT, e >= 0" =
            identical(sprintf("%.4f (%.4f)", rs_unbal$log_commits$att, rs_unbal$log_commits$se),
                      sprintf("%.4f (%.4f)", agg_e0_base$overall.att, agg_e0_base$overall.se)))

# Rows estimated here, for the two share outcomes (log adaptation: above).
rs_drop24 <- list(log_commits = rs_cell_res(res_no2024_year))
rs_bal    <- lapply(res_balance, function(r) list(log_commits = rs_cell_res(r, dynamic = TRUE)))
rs_base   <- list(log_commits = rs_cell_res(res_base))
for (y in c("share_adapt", "share_within")) {
  message("-- ", y, ": drop 2024, balanced windows, baseline level")
  rs_drop24[[y]] <- rs_cell_res(run_spec(panel_no2024_year, yname = y))
  for (be in names(res_balance))
    rs_bal[[be]][[y]] <- rs_cell_res(run_spec(did_panel_main, yname = y, balance_e = as.integer(be)),
                                     dynamic = TRUE)
  base_y <- did_panel_full %>%
    filter(year %in% BASE_YEARS) %>%
    group_by(country_id) %>%
    summarise(base_outcome = mean(.data[[y]], na.rm = TRUE), .groups = "drop") %>%
    mutate(base_outcome = if_else(is.finite(base_outcome), base_outcome, NA_real_))
  rs_base[[y]] <- rs_cell_res(run_spec(did_panel_main %>% left_join(base_y, by = "country_id"),
                                       yname = y,
                                       xformla = ~ ge_est + log_population + base_outcome))
}

n_bal <- unname(rr04$balanced_n[["recipients"]])
rs_rows <- list(
  list(label = "Main specification (Table~\\ref{tab:combined_wide_main})", cells = rs_main),
  list(label = "Balanced panel, covariates at $g-1$",  cells = rr04$balanced),
  list(label = "Panel from 2010",                        cells = rr04$panel2010),
  list(label = "Excluding India",                        cells = rr04$noindia),
  list(label = "Not-yet-treated control group",          cells = rr04$notyettreated),
  list(label = "Baseline (2009--2012) outcome level as covariate", cells = rs_base),
  list(label = "de Chaisemartin--D'Haultf{\\oe}uille estimator", cells = rr04$dcdh),
  list(label = "Outcome in raw USD, $\\log(1+x)$",         cells = rr04$rawusd),
  list(label = "PPML in levels",                           cells = rr04$ppml),
  list(label = "Drop the 2024 cohort and calendar year 2024", cells = rs_drop24),
  list(label = "Unbalanced, all post-adoption years (dynamic average)", cells = rs_unbal),
  list(label = "Balanced event window, 1 post-adoption year",  cells = rs_bal[["1"]]),
  list(label = "Balanced event window, 2 post-adoption years", cells = rs_bal[["2"]])
)
# Every row carries log adaptation commitments; the two log-outcome cells
# that are not meaningful for a share are "---".
stopifnot(all(vapply(rs_rows, function(r) !is.null(r$cells$log_commits), logical(1L))))

#' Two table lines (estimates and pre-trend p; SEs) for one row
#'
#' @param r row: list(label, cells)
#' @param vars outcomes, one ATT/p column pair each ("---" if the cell is absent)
#' @param pad empty cells appended (Panel B fills the 7-column tabular)
rs_lines <- function(r, vars, pad = 0L) {
  est <- unlist(lapply(vars, function(v) {
    ce <- r$cells[[v]]; if (is.null(ce)) c("---", "---") else c(ce$att_fmt, ce$p_fmt) }))
  ses <- unlist(lapply(vars, function(v) {
    ce <- r$cells[[v]]; c(if (is.null(ce)) "" else ce$se_fmt, "") }))
  fill <- rep("", pad)
  paste0(c(paste(c(r$label, est, fill), collapse = " & "),
           paste(c("", ses, fill), collapse = " & ")), " \\\\")
}
rows_b <- Filter(function(r) !is.null(r$cells$lcommitments_all), rs_rows)
col_head <- function(labs) paste0(" & ", paste(sprintf("\\multicolumn{2}{c}{%s}", labs),
                                               collapse = " & "))
rs_tabular <- c(
  "\\begin{tabular}{lcccccc}",
  "\\toprule",
  "\\multicolumn{7}{l}{\\textit{Panel A: adaptation outcomes}} \\\\",
  paste0(col_head(c("log(Adaptation commitments)", "Adaptation share (\\% of global)",
                    "Within-country share (pp)")), " \\\\"),
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}",
  "Specification & ATT & Pre-trend $p$ & ATT & Pre-trend $p$ & ATT & Pre-trend $p$ \\\\",
  "\\midrule",
  unlist(lapply(rs_rows, rs_lines, vars = rs_vars_a)),
  "\\midrule",
  "\\multicolumn{7}{l}{\\textit{Panel B: envelope outcomes}} \\\\",
  paste0(col_head(c("log(Total commitments)", "log(Non-adaptation commitments)")), " & & \\\\"),
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
  "Specification & ATT & Pre-trend $p$ & ATT & Pre-trend $p$ & & \\\\",
  "\\midrule",
  unlist(lapply(rows_b, rs_lines, vars = rs_vars_b, pad = 2L)),
  "\\bottomrule", "\\end{tabular}"
)

for (r in rs_rows)
  message(sprintf("  %-55s %s", r$label, paste(vapply(names(r$cells), function(v)
    sprintf("%s %s %s p=%s", v, r$cells[[v]]$att_fmt, r$cells[[v]]$se_fmt, r$cells[[v]]$p_fmt),
    character(1L)), collapse = " | ")))

write_tex_float(
  file.path(dir_tabs_rs, "robust_summary.tex"),
  "Robustness of the main estimates: sample, comparison group, conditioning, estimator and outcome units",
  "tab:robust_summary",
  rs_tabular,
  paste0(
    as_in_main(c("estimator", "covariates", "control group", "inference",
                 "pre-trend test", "significance stars")),
    "; SE in parentheses below each ATT. Each row changes one element of the main ",
    "specification: balanced panel, the ", n_bal, " recipients observed in every year, so ",
    "covariates enter at $g-1$; panel from 2010, 2009 dropped; excluding India, a ",
    "never-treated recipient; not-yet-treated control group; baseline outcome level, the ",
    "recipient's ", min(BASE_YEARS), "--", max(BASE_YEARS), " mean of the outcome added to the ",
    "covariates; de Chaisemartin--D'Haultf{\\oe}uille (2026) estimator, average total effect ",
    "(Table~\\ref{tab:dcdh}); raw USD, outcome $\\log(1+\\text{commitments in USD})$ and PPML, ",
    "Sun--Abraham Poisson in levels (both Table~\\ref{tab:units_zeros}); drop 2024, the 2024 ",
    "cohort and calendar year 2024; unbalanced dynamic average, the equal-weight mean of the ",
    "main fit's dynamic ATTs over all post-adoption years (as in Table~\\ref{tab:att_anticipation}); ",
    "balanced event window, the same average restricted to the cohorts observed for at least ",
    "one (two) post-adoption years, to be compared with the unbalanced dynamic average, not ",
    "with the main simple ATT. Within-country share: ",
    "$100 \\times E/S$, as in Table~\\ref{tab:within_share_wide}. ``---'' in an ATT cell: ",
    "not meaningful (a share has no raw-USD or levels version); in a pre-trend cell: no ",
    "comparable four-lead test (de Chaisemartin--D'Haultf{\\oe}uille and PPML), or the ",
    "test of the main fit (the three dynamic-average rows, which re-aggregate it). Panel B: rows ",
    "in which total and non-adaptation commitments are estimated"),
  SRC_NOTE)

message("\n=== 13_cohort_anticipation.R: complete ===\n")
