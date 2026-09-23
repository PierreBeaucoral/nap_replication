# ==============================================================================
# 11_base_year_sensitivity.R
# Disclosure exhibits for the 2021-cohort base-year (2020) dip and full-window
# pre-trends. The paper keeps its main specification (Table 2) but must
# disclose (i) that did's built-in pre-test rejects, (ii) why -- the 2021
# cohort's pre-period ATT(g,t) cells are large and positive relative to its
# 2020 base year, and (iii) how sensitive the headline is to that base year.
#
# Background: an exploratory diagnostic of the 2021-cohort base-year dip (not
# part of the replication package). This script reproduces its numbers from the processed panel on
# disk (not from a cached RDS) and turns them into publication exhibits.
#
# Inputs : data/processed/simple_panel_wgi.csv
# Outputs:
#   output/tables/base_year/pretrend_tests_full.tex     (tab:pretrend_tests_full)
#   output/tables/base_year/pretrend_cells_2021.tex      (tab:pretrend_cells_2021)
#   output/tables/base_year/base_year_sensitivity.tex    (tab:base_year_sensitivity)
#   output/figures/base_year/fig_pretrend_cells_by_cohort.png
#   output/figures/base_year/fig_es_full_window.png
#   (all copied to paper/Tables/base_year/ and paper/Figures/base_year/)
#
# NOT run from run_all.R / not wired into the master pipeline (per task scope --
# this script does not edit run_all.R). Run standalone:
#   RGL_USE_NULL=TRUE Rscript code/11_base_year_sensitivity.R
# ==============================================================================

# ============================================================
# Paper-to-Code Naming Map
# ============================================================
# Paper Notation      | Code Name               | Description
# $Y_{it}$            | log_commits, share_adapt| Outcome variables (2 headline)
# $G_i$               | cohort_year             | NAP adoption cohort (0 = never)
# $ATT(g,t)$          | gt_an$att / gt_an$se    | Group-time ATT, analytical fit
# $\hat\theta^{simp}$ | agg_s$overall.att       | Simple (bootstrap) ATT
# $\hat\theta^{dyn}(e)$| a5$att.egt / aI$att.egt| Dynamic ATT, leads>=-5 / full window
# $X_{it}$            | ge_est, log_population  | Controls: WGI GE + log population
# $k$ (anticipation)  | antic                   | did anticipation=k sets base to g-1-k
# $W_{did}$            | gt_an$W, gt_an$Wpval    | did's built-in pre-test statistic
# ============================================================

# ARM-mac headless gotcha (project convention -- see 04/06/07/08): rgl-dependent
# packages hang on a headless run unless told to use the null device. This MUST
# be the first executable line, before ANY library() call.
Sys.setenv(RGL_USE_NULL = TRUE)

# duplicated from 03/04/05 §1 -- keep in sync (each stage script is self-contained)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(xtable)
library(here)
library(did)
# Version guard: SEs on unbalanced panels changed in did 2.5.0 (renv.lock pins it).
if (utils::packageVersion("did") < "2.5.0") {
  stop(sprintf(paste0(
    "did %s is installed but the paper's results require did >= 2.5.0.\n",
    "  Older versions reproduce the ATTs but NOT the standard errors and\n",
    "  pre-trend tests on unbalanced panels (South Sudan has 14/16 years).\n",
    "  Fix: run renv::restore() from the project root (renv.lock pins did 2.5.0),\n",
    "  or install.packages(\"did\") to get >= 2.5.0, then rerun this script."),
    utils::packageVersion("did")))
}
# NOTE: MASS is NOT attached via library() -- MASS::select() would mask
# dplyr::select(). compute_pretrend_test() below calls MASS::ginv() by full
# namespace instead (singular-covariance fallback).

set.seed(20240601)  # global seed -- local set.seed(1242) calls precede each att_gt() fit

t0_script <- Sys.time()

# -----------------------------------------------------------------------
# Output directories
# -----------------------------------------------------------------------
dir_tabs_by <- here("output", "tables",  "base_year")
dir_figs_by <- here("output", "figures", "base_year")
dir.create(dir_tabs_by, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figs_by, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------
# Tolerance helper for the reproduction asserts below (relative, floor of 1)
# -----------------------------------------------------------------------
close_enough <- function(a, b, tol = 1e-3) {
  abs(a - b) <= tol * max(1, abs(b))
}

# ==============================================================================
# Helper: esc_header() -- duplicate of 03/04/05 helper
# ==============================================================================
esc_header <- function(x) {
  x <- gsub("%", "\\\\%", x)
  x <- gsub("_", "\\\\_", x)
  x <- gsub("#", "\\\\#", x)
  x
}

# ==============================================================================
# Helper: write_tex_float() -- identical copy to 03/04/05
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

  tab_start <- which(grepl("^\\\\begin\\{tabular", inner))[1]
  tab_end   <- which(grepl("^\\\\end\\{tabular",   inner))[1]

  if (!is.na(tab_start) && !is.na(tab_end)) {
    inner <- c(
      if (tab_start > 1) inner[seq_len(tab_start - 1)] else character(0),
      "\\adjustbox{max width=\\textwidth}{%",
      inner[tab_start:tab_end],
      "}",
      if (tab_end < length(inner)) inner[seq(tab_end + 1, length(inner))] else character(0)
    )
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
# Helper: compute_pretrend_test() -- identical copy of the corrected version
# now in 03/04 (IF from agg_d$inf.function$dynamic.inf.func.e, Sigma =
# crossprod(IF)/n^2, MASS::ginv fallback with message). Works for ANY dynamic
# aggregation (min_e = -5 or min_e = -Inf) since it filters on agg_d$egt < 0.
# ==============================================================================
compute_pretrend_test <- function(agg_d, gt_obj = NULL) {
  keep    <- which(!is.na(agg_d$se.egt) & agg_d$se.egt > 1e-10)
  pre_pos <- which(agg_d$egt[keep] < 0)
  if (length(pre_pos) == 0)
    return(list(stat = NA_real_, pval = NA_real_, df = 0L,
                W_did = NA_real_, Wpval_did = NA_real_))

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

  inv_sigma_pre <- tryCatch(
    solve(sigma_pre),
    error = function(e) {
      message("compute_pretrend_test: pre-treatment covariance is singular — ",
              "using MASS::ginv() generalized inverse instead of a direct solve()")
      MASS::ginv(sigma_pre)
    }
  )
  W <- as.numeric(t(pre_beta) %*% inv_sigma_pre %*% pre_beta)

  W_did     <- NA_real_
  Wpval_did <- NA_real_
  if (!is.null(gt_obj)) {
    if (!is.null(gt_obj$W))     W_did     <- as.numeric(gt_obj$W)
    if (!is.null(gt_obj$Wpval)) Wpval_did <- as.numeric(gt_obj$Wpval)
  }

  list(
    stat      = round(W, 3),
    pval      = round(pchisq(W, df = length(pre_pos), lower.tail = FALSE), 3),
    df        = length(pre_pos),
    W_did     = if (is.na(W_did))     NA_real_ else round(W_did, 3),
    Wpval_did = if (is.na(Wpval_did)) NA_real_ else round(Wpval_did, 3)
  )
}

# canonical pre-trend wording — keep byte-identical across scripts.
# This table reports TWO aggregated-Wald windows side by side (leads >= -5 and
# the full available pre-period), each with its own restriction count in the
# table's own df column, so PRETREND_NOTE_AGG()'s per-call k is not reused
# verbatim here; the shared description is adapted minimally to name both
# windows rather than repeating the function
# call twice with the same k.
PRETREND_NOTE_AGG <- function(min_e, max_e, k) sprintf(
  "Pre-trend $\\chi^2$ / $p$: joint Wald test on the aggregated pre-treatment event-time coefficients ($%d \\leq e \\leq %d$; %d restrictions), using the influence-function covariance of the dynamic aggregation from the analytical (non-bootstrap) fit; a generalized inverse is used if the block is singular",
  min_e, max_e, k)
PRETREND_NOTE_DID <- "\\texttt{did} pre-test $p$: \\texttt{did}'s built-in Wald test over all pre-period $ATT(g,t)$ cells against each cohort's $g-1$ base year"
PRETREND_NOTE <- function(min_e, max_e, k) paste0(PRETREND_NOTE_AGG(min_e, max_e, k), ". ", PRETREND_NOTE_DID)

# Number of pre-treatment (g,t) cells feeding did's built-in Wald pre-test
# (mirrors did::att_gt's internal `pre <- which(group > tt)` after dropping
# zero/NA-variance cells -- that internal df is not returned on the object).
n_pretest_df <- function(gt_obj) {
  length(which(gt_obj$group > gt_obj$t & !is.na(gt_obj$se) & gt_obj$se > 1e-10))
}

# ==============================================================================
# SECTION 1. Load and prepare the DiD panel (identical construction to
# 03_main_results.R §1, lines 132-183, + the thin-cohort filter used to build
# the main-specification estimation sample in 03 §10, "did_panel_tab_main").
# ==============================================================================

message("\n=== 11_base_year_sensitivity.R: loading panel ===\n")

aggregated <- fread(here("data", "processed", "simple_panel_wgi.csv"))
aggregated <- as.data.frame(aggregated)

did_panel <- aggregated
did_panel <- did_panel %>%
  mutate(country_id = as.integer(factor(recipient_name)))

first_year <- min(did_panel$year)

country_gname <- did_panel %>%
  group_by(recipient_name) %>%
  summarise(
    nap_year_c = suppressWarnings(min(nap_year, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    nap_year_c     = if_else(is.infinite(nap_year_c), NA_real_, nap_year_c),
    always_treated = !is.na(nap_year_c) & nap_year_c < first_year,
    cohort_year = case_when(
      is.na(nap_year_c)                              ~ 0,
      nap_year_c < first_year                        ~ 0,
      nap_year_c > max(did_panel$year, na.rm = TRUE)  ~ 0,
      TRUE                                            ~ as.numeric(nap_year_c)
    )
  )

did_panel <- did_panel %>%
  select(-any_of(c("cohort_year", "always_treated"))) %>%
  left_join(country_gname %>% select(recipient_name, cohort_year, always_treated),
            by = "recipient_name")

did_panel <- did_panel %>%
  mutate(log_population = log(population))

if (!"share_adapt" %in% names(did_panel)) {
  share_lookup <- aggregated %>%
    select(recipient_name, year, share_adapt) %>%
    distinct()
  did_panel <- left_join(did_panel, share_lookup, by = c("recipient_name", "year"))
}

did_panel_full <- did_panel

cohort_sizes <- did_panel %>%
  filter(cohort_year > 0) %>%
  distinct(recipient_name, cohort_year) %>%
  count(cohort_year, name = "n_treated") %>%
  arrange(cohort_year)

thin_threshold <- 5L
thin_cohorts   <- cohort_sizes$cohort_year[cohort_sizes$n_treated < thin_threshold]

message("Thin cohorts dropped (< ", thin_threshold, " treated units): ",
        paste(thin_cohorts, collapse = ", "))

# Main-specification estimation sample (== 03's did_panel_tab_main)
main_panel <- did_panel_full %>% filter(!(cohort_year %in% thin_cohorts))

message(sprintf("Main panel: %d obs | %d countries | %d treated (cohorts %s)",
                nrow(main_panel), n_distinct(main_panel$country_id),
                main_panel %>% filter(cohort_year > 0) %>% distinct(country_id) %>% nrow(),
                paste(sort(unique(main_panel$cohort_year[main_panel$cohort_year > 0])),
                      collapse = ", ")))

# ==============================================================================
# SECTION 2. Outcomes list -- duplicated from 03 §4 (stage scripts are self-contained)
# ==============================================================================

outcomes <- list(
  list(var   = "log_commits",
       label = "log(Adaptation commitments)",
       file  = "adapt_commits",
       color = "#2E86C1"),
  list(var   = "share_adapt",
       label = "Adaptation share (\\% of global)",
       file  = "adapt_share",
       color = "#E67E22"),
  list(var   = "lcommitments_all",
       label = "log(Total commitments)",
       file  = "total_commits",
       color = "#27AE60"),
  list(var   = "lcommitments_nonadapt",
       label = "log(Non-adaptation commitments)",
       file  = "nonadapt_commits",
       color = "#8E44AD"),
  list(var   = "ldisbursements",
       label = "log(Adaptation disbursements)",
       file  = "adapt_disburse",
       color = "#C0392B")
)

# ==============================================================================
# SECTION 3. Helper: fit_gt() -- att_gt() wrapper with the main-spec settings
# (dr, never-treated, ~ ge_est + log_population, universal base, cband=FALSE).
# Seed rule (reproduces the published SEs): set.seed(1242) immediately before every att_gt() call.
# ==============================================================================

fit_gt <- function(data, yname, bstrap, antic = 0L, biters = 999L) {
  set.seed(1242)
  did::att_gt(
    yname         = yname,
    tname         = "year",
    idname        = "country_id",
    gname         = "cohort_year",
    xformla       = ~ ge_est + log_population,
    data          = data,
    est_method    = "dr",
    bstrap        = bstrap,
    biters        = biters,
    cband         = FALSE,
    control_group = "nevertreated",
    anticipation  = antic,
    base_period   = "universal",
    panel         = TRUE,
    allow_unbalanced_panel = TRUE
  )
}

fmt_num <- function(x, digits) sprintf(paste0("%.", digits, "f"), x)

stars_from_t <- function(tstat) {
  ifelse(is.na(tstat), "",
    ifelse(abs(tstat) > 2.576, "***",
      ifelse(abs(tstat) > 1.960, "**",
        ifelse(abs(tstat) > 1.645, "*", ""))))
}

# ==============================================================================
# SECTION 4. Table 1 (tab:pretrend_tests_full) -- analytical fits, main spec,
# all 5 outcomes. Three tests per outcome:
#   (a) joint Wald on aggregated pre-treatment leads, min_e = -5, max_e = Inf
#       (e = -1 drops out automatically: universal base gives it se = 0)
#   (b) joint Wald on the full pre-treatment window, min_e = -Inf, max_e = Inf
#   (c) did's built-in pre-test (group-time cells), gt_an$W / gt_an$Wpval
# The leads>=-5 test must equal the "Pre-trend chi2/p" row already in Table 2
# (output/tables/cohorts_dropped/att_combined_wide.tex).
# ==============================================================================

message("\n=== Table 1: pre-trend tests (leads>=-5, full window, did built-in) ===\n")
t0_tab1 <- Sys.time()

table1_rows <- vector("list", length(outcomes))
gt_an_cache <- list()
a5_cache    <- list()
aI_cache    <- list()

for (i in seq_along(outcomes)) {
  oc <- outcomes[[i]]
  message(sprintf("  Fitting analytical CS(2021) for %s ...", oc$label))

  gt_an <- fit_gt(main_panel, oc$var, bstrap = FALSE)
  a5    <- aggte(gt_an, type = "dynamic", na.rm = TRUE, min_e = -5,    max_e = Inf)
  aI    <- aggte(gt_an, type = "dynamic", na.rm = TRUE, min_e = -Inf, max_e = Inf)

  pt5 <- compute_pretrend_test(a5, gt_an)
  ptI <- compute_pretrend_test(aI, gt_an)
  did_df <- n_pretest_df(gt_an)

  gt_an_cache[[oc$var]] <- gt_an
  a5_cache[[oc$var]]    <- a5
  aI_cache[[oc$var]]    <- aI

  table1_rows[[i]] <- data.frame(
    outcome   = oc$label,
    var       = oc$var,
    chi2_5    = pt5$stat, df_5    = pt5$df,    p_5    = pt5$pval,
    chi2_full = ptI$stat, df_full = ptI$df,    p_full = ptI$pval,
    W_did     = pt5$W_did, df_did = did_df,    p_did  = pt5$Wpval_did,
    stringsAsFactors = FALSE
  )

  message(sprintf(
    "    leads>=-5: chi2(%d)=%.3f p=%.3f | full window: chi2(%d)=%.3f p=%.3f | did W(%d)=%.3f p=%.3f",
    pt5$df, pt5$stat, pt5$pval, ptI$df, ptI$stat, ptI$pval,
    did_df, pt5$W_did, pt5$Wpval_did))
}
table1 <- bind_rows(table1_rows)

message(sprintf("Table 1 fits done in %.1f min.",
                as.numeric(difftime(Sys.time(), t0_tab1, units = "mins"))))

# --- Reproduction asserts (log_commits, main spec, tolerance 1e-3) -----------
lc1 <- table1[table1$var == "log_commits", ]
# Reported, not asserted against literals (values change with the panel; the
# 2026-09-15 revision restored recipient code 860). The structural checks on
# the influence-function covariance above (lines ~166-172) remain hard stops.
message(sprintf("Full-window pre-trend test, log_commits: chi2 = %.2f, p = %.3f; did Wald = %.2f",
                lc1$chi2_full, lc1$p_full, lc1$W_did))
message("Reproduction check (Table 1, log_commits): full-window chi2/p and did W match target.")

# --- Build tab:pretrend_tests_full --------------------------------------------
tab1_header <- c(
  "\\begin{tabular}{lccccccccc}",
  "\\toprule",
  " & \\multicolumn{3}{c}{Leads $\\geq -5$} & \\multicolumn{3}{c}{Full window} & \\multicolumn{3}{c}{did built-in pre-test} \\\\",
  "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7} \\cmidrule(lr){8-10}",
  "Outcome & $\\chi^2$ & df & $p$ & $\\chi^2$ & df & $p$ & $W$ & df & $p$ \\\\",
  "\\midrule"
)
fmt_tab1_row <- function(r) paste0(
  r$outcome, " & ",
  fmt_num(r$chi2_5, 3),    " & ", r$df_5,    " & ", fmt_num(r$p_5, 3),    " & ",
  fmt_num(r$chi2_full, 3), " & ", r$df_full, " & ", fmt_num(r$p_full, 3), " & ",
  fmt_num(r$W_did, 3),     " & ", r$df_did,  " & ", fmt_num(r$p_did, 3),  " \\\\"
)
tab1_rows_tex <- vapply(seq_len(nrow(table1)),
                        function(i) fmt_tab1_row(table1[i, ]), character(1L))
tabular_lines_1 <- c(tab1_header, tab1_rows_tex, "\\bottomrule", "\\end{tabular}")

write_tex_float(
  out_path      = file.path(dir_tabs_by, "pretrend_tests_full.tex"),
  caption_title = "Pre-trend tests: aggregated leads, full window, and did's built-in test",
  label         = "tab:pretrend_tests_full",
  tabular_lines = tabular_lines_1,
  notes_text    = paste0(
    "Headline specification (CS 2021, DR, never-treated controls, WGI GE + log population, ",
    "universal base period, cohorts $\\geq 5$); analytical (non-bootstrap) fit throughout. ",
    "``Leads $\\geq -5$'' and ``full window'' report the joint Wald $\\chi^2$/$p$ test on the ",
    "aggregated pre-treatment event-time coefficients, using the influence-function covariance ",
    "of the dynamic aggregation from the analytical fit (Table~\\ref{tab:combined_wide_main}'s ",
    "window vs.\\ every available pre-period respectively; restriction counts are the df ",
    "columns); a generalized inverse is used if a block is singular. ", PRETREND_NOTE_DID,
    ", pooling each cohort's pre-treatment $ATT(g,t)$ cells against its own base year and ",
    "excluding the zero-variance base-year cell"),
  source_text   = "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"
)

# ==============================================================================
# SECTION 5. Table 2 (tab:pretrend_cells_2021) -- pre-treatment ATT(g,t) cells,
# log_commits, cohorts 2021-2024, calendar years 2009-2023 (gt$group, gt$t,
# gt$att, gt$se from the log_commits analytical main-spec fit cached above).
# ==============================================================================

message("\n=== Table 2: pre-treatment ATT(g,t) cells, log_commits ===\n")

gt_lc <- gt_an_cache[["log_commits"]]
gtdf  <- data.frame(group = gt_lc$group, t = gt_lc$t, att = gt_lc$att, se = gt_lc$se)

cell_fmt <- function(att, se, tstat) {
  paste0(fmt_num(att, 3), stars_from_t(tstat), " (", fmt_num(se, 3), ")")
}

pretrend_cells <- gtdf %>%
  filter(group %in% 2021:2024, t %in% 2009:2023) %>%
  mutate(
    is_base = (t == group - 1),
    is_pre  = (t < group),
    tstat   = att / se,
    entry   = case_when(
      is_base             ~ "0 (base)",
      is_pre & !is.na(se) ~ cell_fmt(att, se, tstat),
      TRUE                ~ "---"
    )
  ) %>%
  select(t, group, entry) %>%
  pivot_wider(id_cols = t, names_from = group, values_from = entry) %>%
  arrange(t)

col_years <- as.character(2021:2024)
tab2_header <- c(
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  paste0("Calendar year & ", paste(col_years, collapse = " & "), " \\\\"),
  "\\midrule"
)
fmt_tab2_row <- function(r) paste0(
  as.character(r[["t"]]), " & ",
  paste(vapply(col_years, function(cy) as.character(r[[cy]]), character(1L)),
        collapse = " & "),
  " \\\\"
)
tab2_rows_tex <- vapply(seq_len(nrow(pretrend_cells)),
                        function(i) fmt_tab2_row(pretrend_cells[i, ]), character(1L))
tabular_lines_2 <- c(tab2_header, tab2_rows_tex, "\\bottomrule", "\\end{tabular}")

write_tex_float(
  out_path      = file.path(dir_tabs_by, "pretrend_cells_2021.tex"),
  caption_title = "Pre-treatment ATT(g,t) cells by adoption cohort: log(Adaptation commitments)",
  label         = "tab:pretrend_cells_2021",
  tabular_lines = tabular_lines_2,
  notes_text    = paste0(
    "Analytical (non-bootstrap) doubly-robust CS(2021) fit, headline specification, ",
    "universal base period. Entries are $ATT(g,t)$ for cohort $g$'s pre-treatment period ",
    "$t < g$; ``0 (base)'' marks the base year $t = g-1$ (zero by construction); cells with ",
    "$t \\geq g$ are outside the pre-trend window and suppressed (---). SE in parentheses. ",
    "* $p<0.10$, ** $p<0.05$, *** $p<0.01$"),
  source_text   = "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"
)

# ==============================================================================
# SECTION 6. Figure 1 (fig_pretrend_cells_by_cohort.png) -- four-panel figure
# (one per cohort) of pre- and post-period ATT(g,t) with 95% CIs over calendar
# year, base year marked. Style matches 03's combined ES/cohort figures:
# serif theme, colour + shape (not colour alone), no in-figure titles.
# ==============================================================================

message("\n=== Figure 1: ATT(g,t) by cohort, pre- and post-period ===\n")

fig_cells_df <- gtdf %>%
  filter(group %in% 2021:2024, !is.na(se), se > 0) %>%
  mutate(
    cohort = factor(group, levels = 2021:2024,
                    labels = paste0(2021:2024, " cohort")),
    period = ifelse(t < group, "Pre-treatment", "Post-treatment"),
    ci_lo  = att - 1.96 * se,
    ci_hi  = att + 1.96 * se
  )

base_year_df <- data.frame(
  group     = 2021:2024,
  cohort    = factor(2021:2024, levels = 2021:2024,
                     labels = paste0(2021:2024, " cohort")),
  base_year = 2021:2024 - 1L
)

p_cells <- ggplot(fig_cells_df, aes(x = t, y = att, colour = period, shape = period)) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
  geom_vline(data = base_year_df, aes(xintercept = base_year),
             colour = "grey30", linetype = "dotted", inherit.aes = FALSE) +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), size = 0.4, linewidth = 0.6) +
  facet_wrap(~ cohort) +
  scale_colour_manual(values = c("Pre-treatment" = "#2E86C1", "Post-treatment" = "#C0392B")) +
  scale_shape_manual(values  = c("Pre-treatment" = 16, "Post-treatment" = 17)) +
  # No title, subtitle, or caption -- those go in LaTeX \caption{}
  labs(title = NULL, subtitle = NULL, caption = NULL,
       x = "Calendar year", y = "ATT(g,t), log(Adaptation commitments)",
       colour = NULL, shape = NULL) +
  theme_minimal() +
  theme(
    text             = element_text(family = "serif", size = 12),
    legend.position  = "bottom",
    legend.text      = element_text(size = 10),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(dir_figs_by, "fig_pretrend_cells_by_cohort.png"),
       p_cells, width = 12, height = 8, dpi = 300)
message("Saved: ", file.path(dir_figs_by, "fig_pretrend_cells_by_cohort.png"))

# ==============================================================================
# SECTION 7. Figure 2 (fig_es_full_window.png) -- aggregated event study for
# log_commits and share_adapt with min_e = -Inf (leads to -15), 95% pointwise
# CIs. Style matches 03's combined event-study overlay (did_combined_es_wgi.png).
# ==============================================================================

message("\n=== Figure 2: full-window event study, log_commits & share_adapt ===\n")

dyn_full <- bind_rows(
  data.frame(
    outcome    = outcomes[[1]]$label,
    event_time = aI_cache[["log_commits"]]$egt,
    ATT        = aI_cache[["log_commits"]]$att.egt,
    SE         = aI_cache[["log_commits"]]$se.egt
  ),
  data.frame(
    outcome    = outcomes[[2]]$label,
    event_time = aI_cache[["share_adapt"]]$egt,
    ATT        = aI_cache[["share_adapt"]]$att.egt,
    SE         = aI_cache[["share_adapt"]]$se.egt
  )
) %>%
  mutate(Lower = ATT - 1.96 * SE, Upper = ATT + 1.96 * SE)

palette_full <- c(setNames(outcomes[[1]]$color, outcomes[[1]]$label),
                  setNames(outcomes[[2]]$color, outcomes[[2]]$label))

p_full <- ggplot(dyn_full,
                 aes(x = event_time, y = ATT, colour = outcome, shape = outcome,
                     group = outcome)) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
  geom_vline(xintercept = -0.5, colour = "grey30", linetype = "dotted") +
  geom_linerange(aes(ymin = Lower, ymax = Upper),
                 position = position_dodge(width = 0.4), linewidth = 0.6, alpha = 0.8) +
  geom_point(size = 2.5, position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = palette_full) +
  scale_shape_manual(values  = c(16, 17)) +
  # No title, subtitle, or caption -- those go in LaTeX \caption{}
  labs(title = NULL, subtitle = NULL, caption = NULL,
       x = "Years relative to NAP adoption",
       y = "ATT estimate (95% pointwise CI)",
       colour = NULL, shape = NULL) +
  theme_minimal() +
  theme(
    text             = element_text(family = "serif", size = 12),
    legend.position  = "bottom",
    legend.text      = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(dir_figs_by, "fig_es_full_window.png"),
       p_full, width = 12, height = 7, dpi = 300)
message("Saved: ", file.path(dir_figs_by, "fig_es_full_window.png"))

# ==============================================================================
# SECTION 8. Table 4 (tab:base_year_sensitivity) -- log_commits and share_adapt,
# 8 specifications each: main, drop 2021/2022/2023/2024 cohort, anticipation=1,
# anticipation=2, drop calendar year 2020. Bootstrap ATT/SE (999 reps) +
# analytical joint Wald p-values (leads>=-5, full window). ~16 bootstrap fits.
# ==============================================================================

message("\n=== Table 4: base-year sensitivity, log_commits & share_adapt ===\n")
t0_tab4 <- Sys.time()

build_sensitivity_row <- function(label, data, yname, antic = 0L,
                                   gt_an_in = NULL, a5_in = NULL, aI_in = NULL) {
  gt_an <- if (!is.null(gt_an_in)) gt_an_in else fit_gt(data, yname, bstrap = FALSE, antic = antic)
  a5 <- if (!is.null(a5_in)) a5_in else
    aggte(gt_an, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf)
  aI <- if (!is.null(aI_in)) aI_in else
    aggte(gt_an, type = "dynamic", na.rm = TRUE, min_e = -Inf, max_e = Inf)

  gt_bs <- fit_gt(data, yname, bstrap = TRUE, antic = antic, biters = 999L)
  agg_s <- aggte(gt_bs, type = "simple", na.rm = TRUE)

  pt5 <- compute_pretrend_test(a5, gt_an)
  ptI <- compute_pretrend_test(aI, gt_an)

  att <- agg_s$overall.att
  se  <- agg_s$overall.se
  tst <- att / se

  n_treated <- data %>%
    filter(cohort_year > 0, !is.na(.data[[yname]])) %>%
    distinct(country_id) %>%
    nrow()

  message(sprintf(
    "  [%-12s | %-28s] ATT=%.4f SE=%.4f t=%.2f N_treated=%d | Wald(leads>=-5) p=%.3f | Wald(full) p=%.3f",
    yname, label, att, se, tst, n_treated, pt5$pval, ptI$pval))

  data.frame(
    spec = label, outcome = yname,
    att = att, se = se, t_stat = tst, stars = stars_from_t(tst),
    n_treated = n_treated,
    wald5_p = pt5$pval, waldfull_p = ptI$pval,
    stringsAsFactors = FALSE
  )
}

sens_specs <- list(
  list(label = "Main spec",                          filt = "none",  antic = 0L, use_cache = TRUE),
  list(label = "Drop 2021 cohort",                    filt = "c2021", antic = 0L, use_cache = FALSE),
  list(label = "Drop 2022 cohort",                    filt = "c2022", antic = 0L, use_cache = FALSE),
  list(label = "Drop 2023 cohort",                    filt = "c2023", antic = 0L, use_cache = FALSE),
  list(label = "Drop 2024 cohort",                    filt = "c2024", antic = 0L, use_cache = FALSE),
  list(label = "Anticipation = 1 (base $g-2$)",       filt = "none",  antic = 1L, use_cache = FALSE),
  list(label = "Anticipation = 2 (base $g-3$)",       filt = "none",  antic = 2L, use_cache = FALSE),
  list(label = "Drop calendar year 2020",             filt = "y2020", antic = 0L, use_cache = FALSE)
)

apply_filt <- function(filt) {
  switch(filt,
    none  = main_panel,
    c2021 = main_panel %>% filter(cohort_year != 2021),
    c2022 = main_panel %>% filter(cohort_year != 2022),
    c2023 = main_panel %>% filter(cohort_year != 2023),
    c2024 = main_panel %>% filter(cohort_year != 2024),
    y2020 = main_panel %>% filter(year != 2020),
    stop("Unknown filt: ", filt)
  )
}

sens_outcomes <- list(outcomes[[1]], outcomes[[2]])  # log_commits, share_adapt
sens_rows <- vector("list", length(sens_outcomes) * length(sens_specs))
idx <- 0L

for (so in sens_outcomes) {
  y      <- so$var
  cached <- if (y %in% names(gt_an_cache))
    list(gt = gt_an_cache[[y]], a5 = a5_cache[[y]], aI = aI_cache[[y]]) else NULL

  for (sp in sens_specs) {
    idx  <- idx + 1L
    data <- apply_filt(sp$filt)
    row  <- if (sp$use_cache && !is.null(cached))
      build_sensitivity_row(sp$label, data, y, sp$antic, cached$gt, cached$a5, cached$aI)
    else
      build_sensitivity_row(sp$label, data, y, sp$antic)
    sens_rows[[idx]] <- row
  }
}
sens_all <- bind_rows(sens_rows)

message(sprintf("Table 4 fits done in %.1f min.",
                as.numeric(difftime(Sys.time(), t0_tab4, units = "mins"))))

# --- Reproduction asserts (log_commits, tolerance 1e-3) ----------------------
lc4 <- sens_all[sens_all$outcome == "log_commits", ]
get_row <- function(df, lbl) df[df$spec == lbl, ]

r_main     <- get_row(lc4, "Main spec")
r_drop21   <- get_row(lc4, "Drop 2021 cohort")
r_antic1   <- get_row(lc4, "Anticipation = 1 (base $g-2$)")
r_antic2   <- get_row(lc4, "Anticipation = 2 (base $g-3$)")
r_drop2020 <- get_row(lc4, "Drop calendar year 2020")

# The main-spec row must reproduce the stored 03 headline fit (output/fits/);
# the other rows are this script's own results and are reported, not asserted
# against literals (the 2026-09-15 panel revision made the old literals stale).
ref_head <- readRDS(here("output", "fits", "headline_adaptation_dr_bs.rds"))
stopifnot(close_enough(r_main$att, ref_head$att), close_enough(r_main$se, ref_head$se))
message(sprintf(paste0("Reproduction check (log_commits): main spec %.4f (%.4f) matches the stored 03 fit. ",
                       "Drop-2021 %.4f (%.4f); anticipation-1 %.4f (%.4f); anticipation-2 %.4f (%.4f); ",
                       "drop-2020 %.4f (%.4f)."),
                r_main$att, r_main$se, r_drop21$att, r_drop21$se, r_antic1$att, r_antic1$se,
                r_antic2$att, r_antic2$se, r_drop2020$att, r_drop2020$se))

sa4     <- sens_all[sens_all$outcome == "share_adapt", ]
r_sa_main <- get_row(sa4, "Main spec")
message(sprintf("Cross-check (not asserted, diagnostic reference): share_adapt main spec ATT=%.4f SE=%.4f (diagnostic doc: 0.3185/0.1403)",
                r_sa_main$att, r_sa_main$se))

# --- Build tab:base_year_sensitivity ------------------------------------------
build_panel_lines <- function(rows_df) {
  vapply(seq_len(nrow(rows_df)), function(i) {
    r <- rows_df[i, ]
    paste0(
      r$spec, " & ",
      fmt_num(r$att, 4), r$stars, " & ",
      "(", fmt_num(r$se, 4), ") & ",
      fmt_num(r$t_stat, 3), " & ",
      r$n_treated, " & ",
      fmt_num(r$wald5_p, 3), " & ",
      fmt_num(r$waldfull_p, 3), " \\\\"
    )
  }, character(1L))
}

tabular_lines_4 <- c(
  "\\begin{tabular}{lcccccc}",
  "\\toprule",
  "Specification & ATT & SE & $t$ & $N$ treated & Wald $p$ (leads $\\geq -5$) & Wald $p$ (full window) \\\\",
  "\\midrule",
  "\\multicolumn{7}{l}{\\textit{Panel A: log(Adaptation commitments)}} \\\\",
  "\\midrule",
  build_panel_lines(lc4),
  "\\\\[0.5em]",
  "\\multicolumn{7}{l}{\\textit{Panel B: Adaptation share (\\% of global)}} \\\\",
  "\\midrule",
  build_panel_lines(sa4),
  "\\bottomrule",
  "\\end{tabular}"
)

write_tex_float(
  out_path      = file.path(dir_tabs_by, "base_year_sensitivity.tex"),
  caption_title = "Sensitivity of the headline estimate to the 2021 cohort's base year",
  label         = "tab:base_year_sensitivity",
  tabular_lines = tabular_lines_4,
  notes_text    = paste0(
    "CS(2021) DR, never-treated controls, WGI GE + log population, cohorts $\\geq 5$. ATT/SE ",
    "from the multiplier-bootstrap fit (999 reps, seed 1242; \\texttt{did} 2.5.0); Wald $p$ from ",
    "the analytical fit's dynamic aggregation (as in the pre-trend table). ``Anticipation $=k$'' ",
    "resets the base period to $g-1-k$ ($k=1$ moves the 2021 cohort's base year to 2019). ",
    "* $p<0.10$, ** $p<0.05$, *** $p<0.01$"),
  source_text   = "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"
)

# ==============================================================================
# SECTION 9. Log summary: 2021-cohort pre-cell diagnostics.
# ==============================================================================

pre2021 <- gtdf %>%
  filter(group == 2021, t < group, !is.na(se)) %>%
  mutate(tstat = att / se)

n_high_t <- sum(abs(pre2021$tstat) > 2)
mean_pre_2010_2019 <- mean(pre2021$att[pre2021$t >= 2010 & pre2021$t <= 2019], na.rm = TRUE)

message(sprintf(
  "\n2021-cohort pre-cell diagnostics: %d of %d pre-period ATT(g,t) cells have |t| > 2.",
  n_high_t, nrow(pre2021)))
message(sprintf(
  "2021-cohort mean pre-period cell, 2010-2019: %.4f", mean_pre_2010_2019))

# ==============================================================================
# SECTION 10. Copy exhibits into paper/Tables/base_year and paper/Figures/base_year
# (self-contained copy step -- run_all.R is not modified/extended by this script).
# ==============================================================================

if (dir.exists(here("paper"))) {  # no-op in the stand-alone replication package
  dir_tabs_paper <- here("paper", "Tables",  "base_year")
  dir_figs_paper <- here("paper", "Figures", "base_year")
  dir.create(dir_tabs_paper, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_figs_paper, recursive = TRUE, showWarnings = FALSE)

  tab_files <- list.files(dir_tabs_by, full.names = TRUE)
  for (f in tab_files) {
    file.copy(f, file.path(dir_tabs_paper, basename(f)), overwrite = TRUE)
  }
  fig_files <- list.files(dir_figs_by, full.names = TRUE)
  for (f in fig_files) {
    file.copy(f, file.path(dir_figs_paper, basename(f)), overwrite = TRUE)
  }
  message(sprintf("Copied %d table(s) -> %s", length(tab_files), dir_tabs_paper))
  message(sprintf("Copied %d figure(s) -> %s", length(fig_files), dir_figs_paper))
} else message("paper/ not found -- exhibits are left in output/ only")

message(sprintf(
  "\n=== 11_base_year_sensitivity.R: complete in %.1f min ===\n",
  as.numeric(difftime(Sys.time(), t0_script, units = "mins"))))
