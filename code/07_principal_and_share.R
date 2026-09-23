# ==============================================================================
# 07_principal_and_share.R
# Exhibits for three checks: principal-marker-only outcome (Part A), within-
# country adaptation share (Part B), and global-share reconciliation
# arithmetic (Part C).
#
# Inputs : data/processed/simple_panel_wgi.csv
# Outputs:
#   output/tables/principal_share/principal_wide.tex        (tab:principal_wide)
#   output/tables/principal_share/principal_retained.tex     (tab:principal_retained)
#   output/tables/principal_share/within_share_wide.tex      (tab:within_share_wide)
#   output/tables/principal_share/share_reconciliation.tex   (tab:share_reconciliation)
#   output/figures/principal_share/fig_principal_es.png
#   output/figures/principal_share/fig_within_share_es.png
#   output/figures/principal_share/fig_share_reconciliation.png
#   (mirrored to paper/Tables/principal_share/, paper/Figures/principal_share/)
#
# This script is self-contained: it reloads and rebuilds the DiD panel from
# data/processed/simple_panel_wgi.csv exactly as 03_main_results.R /
# 04_robustness.R do (no sourcing: each stage script is self-contained, and
# duplication is the convention in this project; see 04_robustness.R §1).
# ==============================================================================

# ============================================================
# Paper-to-Code Naming Map
# ============================================================
# Paper Notation        | Code Name                | Description
# $Y_{it}$               | log_commits               | Headline: log1p(principal + significant Rio-marker commitments)
# $Y_{it}^{principal}$   | lcommitments_principal    | log1p(principal-only, Rio marker = 2)
# $Y_{it}^{share}$       | share_adapt               | Country i's adaptation commitments as % of the 144-country pool, year t
# $S_{it}$               | share_within              | 100 * commitments / commitments_all (% of country's own total ODA that is adaptation-marked)
# $\text{logit}(S_{it})$ | logit_within              | log((commitments+0.5)/(commitments_all-commitments+0.5))
# $G_i$                  | cohort_year               | Year of first NAP adoption (0 = never)
# $ATT(g,t)$             | gt_obj / gt_main           | Group-time ATT from att_gt()
# $\hat\theta^{simp}$    | agg_s$overall.att          | Calendar-time simple ATT
# $\hat\theta^{dyn}(e)$  | agg_d$att.egt              | Dynamic ATT by event time e
# $X_{it}$                | ge_est, log_population    | Controls: WGI GE + log population
# $W$ (joint pre-trend)   | compute_pretrend_wald()$stat | Corrected joint Wald using dynamic.inf.func.e
# $p_{did}$                | gt_main$Wpval             | did package's own internal pre-test p-value
# $\Delta = \theta_{head}-\theta_{princ}$ | diff  | Headline-minus-principal ATT difference (Part A)
# ============================================================

# Defensive: some transitive dependencies of did/ggplot2 on this machine pull
# in rgl, which hangs headless runs unless told to use the null device (see
# 04_robustness.R). This script
# does not load DIDmultiplegtDYN itself, but setting this first is cheap
# insurance and keeps the convention consistent across code/*.R.
Sys.setenv(RGL_USE_NULL = TRUE)

library(data.table)
library(dplyr)
library(ggplot2)
library(xtable)
library(here)
library(did)
# Version guard: SEs on unbalanced panels changed in did 2.5.0 (renv.lock pins it).
# An older install reproduces the ATTs but not the SEs / pre-trend tests.
# (Identical guard to 03_main_results.R l.35-46.)
if (utils::packageVersion("did") < "2.5.0") {
  stop(sprintf(paste0(
    "did %s is installed but the paper's results require did >= 2.5.0.\n",
    "  Older versions reproduce the ATTs but NOT the standard errors and\n",
    "  pre-trend tests on unbalanced panels (South Sudan has 14/16 years).\n",
    "  Fix: run renv::restore() from the project root (renv.lock pins did 2.5.0),\n",
    "  or install.packages(\"did\") to get >= 2.5.0, then rerun this script."),
    utils::packageVersion("did")))
}

set.seed(20240601)  # global seed — local set.seed(1242) calls follow each estimator

# ------------------------------------------------------------------------
# Output directories
# ------------------------------------------------------------------------
dir.create(here("output", "tables",  "principal_share"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "figures", "principal_share"), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------
# Escape LaTeX-special characters in TABLE COLUMN HEADERS only.
# (Identical helper to 03/04/05.)
# ------------------------------------------------------------------------
esc_header <- function(x) {
  x <- gsub("%", "\\\\%", x)
  x <- gsub("_", "\\\\_", x)
  x <- gsub("#", "\\\\#", x)
  x
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ==============================================================================
# Helper: write_tex_float()
# Identical to 03/04/05_...: wraps a bare tabular block in the standard
# complete float (\begin{table}[H] ... \adjustbox ... minipage notes/source).
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

# ------------------------------------------------------------------------
# Helper: copy_to_paper()
# Mirrors run_all.R's final assembly step (output/ -> paper/), scoped to a
# single file, exactly as 05_heterogeneity.R §25 does for het_difference_tests.tex.
# This script owns only the principal_share/ subtree, so it copies its own
# outputs directly rather than depending on run_all.R (not edited by this task).
# ------------------------------------------------------------------------
copy_to_paper <- function(out_path, kind = c("Tables", "Figures")) {
  kind <- match.arg(kind)
  if (!dir.exists(here("paper"))) {  # stand-alone replication package: no paper/
    message("paper/ not found -- exhibits are left in output/ only")
    return(invisible(NULL))
  }
  rel  <- sub(paste0("^.*output/(tables|figures)/principal_share/"), "", out_path)
  dest <- here("paper", kind, "principal_share", rel)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  file.copy(out_path, dest, overwrite = TRUE)
  message("Mirrored to: ", dest)
  invisible(dest)
}

# ==============================================================================
# SECTION 1. Load and prepare the DiD panel
# (identical to 03_main_results.R §1 — duplicated because each stage script
#  is self-contained; do not diverge without a matching note there.)
# ==============================================================================

message("\n=== 07_principal_and_share.R: loading panel ===\n")

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
      is.na(nap_year_c)                                        ~ 0,
      nap_year_c < first_year                                  ~ 0,
      nap_year_c > max(did_panel$year, na.rm = TRUE)           ~ 0,
      TRUE                                                     ~ as.numeric(nap_year_c)
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

cohort_sizes <- did_panel_full %>%
  filter(cohort_year > 0) %>%
  distinct(recipient_name, cohort_year) %>%
  count(cohort_year, name = "n_treated") %>%
  arrange(cohort_year)

thin_threshold <- 5L
thin_cohorts   <- cohort_sizes$cohort_year[cohort_sizes$n_treated < thin_threshold]

message(sprintf("Full panel: %d rows, %d countries, %d-%d.",
                nrow(did_panel_full), n_distinct(did_panel_full$recipient_name),
                min(did_panel_full$year), max(did_panel_full$year)))
message(sprintf("Ever-adopters (cohort_year > 0): %d.",
                n_distinct(did_panel_full$recipient_name[did_panel_full$cohort_year > 0])))
message(sprintf("thin_cohorts (< %d treated units, dropped from main spec): %s",
                thin_threshold, paste(thin_cohorts, collapse = ", ")))

# Main-spec estimation panel: thin cohorts dropped (identical construction to
# 03_main_results.R §10 did_panel_tab_main).
did_panel_tab_main <- did_panel_full %>% filter(!(cohort_year %in% thin_cohorts))
message(sprintf("Main-spec panel (thin cohorts dropped): %d rows, %d countries, %d treated countries.",
                nrow(did_panel_tab_main), n_distinct(did_panel_tab_main$country_id),
                n_distinct(did_panel_tab_main$country_id[did_panel_tab_main$cohort_year > 0])))

# ==============================================================================
# SECTION 2. Construct outcome variables for Parts A and B
# ==============================================================================

# --- Part Lcommitments_principal already built in 01_prepare_data.R §18 as
#     log1p(commitments_principal) (Rio marker = 2 only) — CONFIRMED log1p, not
#     log(x): see code/01_prepare_data.R l.933. log_commits is likewise log1p
#     (l.673): log_commits = log1p(commitments), commitments = principal (2) +
#     significant (1). Both series therefore share the same log1p functional
#     form, so the (exp(ATT)-1) back-transform used below (identical to
#     03_main_results.R's raw_var_map convention) applies to both unchanged.
stopifnot(
  "commitments_principal" %in% names(did_panel_full),
  "lcommitments_principal" %in% names(did_panel_full),
  "commitments_all" %in% names(did_panel_full)
)

# --- Part Within-country adaptation share ---
# share_within: % of a country's OWN total reported commitments that carry an
#   adaptation Rio marker (1 or 2). NA where commitments_all == 0 (no reported
#   commitments of any kind that year — a genuinely undefined share, not a
#   zero share).
# logit_within: log((commitments + 0.5)/(commitments_all - commitments + 0.5)).
#   The 0.5 offset (a standard continuity correction) keeps the logit defined
#   at the two closed endpoints (commitments == 0 and commitments == commitments_all,
#   the latter guaranteed by 01_prepare_data.R's stopifnot(commitments_all >=
#   commitments)). commitments_all >= commitments always, so the log() argument
#   is algebraically positive even when commitments_all == 0 (it reduces to
#   log(0.5/0.5) = 0) -- but that cell is not a genuine zero logit, it is the
#   same "no reported commitments of any kind that year" case that makes
#   share_within undefined. logit_within is therefore set to NA under the
#   identical commitments_all == 0 condition, not left at its spurious 0.
n_within_undefined <- sum(did_panel_full$commitments_all == 0)
message(sprintf(
  "share_within: %d / %d country-year cells (%.2f%%) have commitments_all == 0 and are set to NA.",
  n_within_undefined, nrow(did_panel_full), 100 * n_within_undefined / nrow(did_panel_full)))

did_panel_full <- did_panel_full %>%
  mutate(
    share_within = if_else(commitments_all > 0, 100 * commitments / commitments_all, NA_real_),
    logit_within = if_else(
      commitments_all > 0,
      log((commitments + 0.5) / (commitments_all - commitments + 0.5)),
      NA_real_
    )
  )
did_panel_tab_main <- did_panel_tab_main %>%
  mutate(
    share_within = if_else(commitments_all > 0, 100 * commitments / commitments_all, NA_real_),
    logit_within = if_else(
      commitments_all > 0,
      log((commitments + 0.5) / (commitments_all - commitments + 0.5)),
      NA_real_
    )
  )

pretreat_share_within <- did_panel_tab_main %>%
  filter(cohort_year > 0, year < cohort_year, !is.na(share_within)) %>%
  pull(share_within)
message("\n--- Distribution of pre-treatment within-country share (share_within, %) ---")
message(sprintf(
  "  Mean = %.2f | Median = %.2f | p90 = %.2f | Max = %.2f | N = %d",
  mean(pretreat_share_within, na.rm = TRUE), stats::median(pretreat_share_within, na.rm = TRUE),
  stats::quantile(pretreat_share_within, probs = 0.90, na.rm = TRUE, names = FALSE),
  max(pretreat_share_within, na.rm = TRUE), length(pretreat_share_within)))
message(sprintf("  Share of pre-treatment cells with share_within < 50%%: %.1f%%",
                100 * mean(pretreat_share_within < 50, na.rm = TRUE)))

# ==============================================================================
# SECTION 3. Corrected joint pre-trend Wald test
#
# 03/04/05_*.R's compute_pretrend_test() never locates the dynamic-aggregation
# influence function (it looks for inf.function$inf.func.egt / $inf.func /
# $inffunc / $inf.func.egt, none of which did >= 2.5.0 populates) and silently
# falls back to a diagonal covariance approximation. The correct covariance
# uses inf.function$dynamic.inf.func.e (with a dimension check); did's own
# att_gt()$Wpval is an alternative. This script implements BOTH and reports
# them side by side. (03/04/05 were corrected the same way on 2026-09-14.)
# ==============================================================================

# canonical pre-trend wording — keep byte-identical across scripts.
# This script's compute_pretrend_wald() has NO generalized-inverse fallback (it
# returns NA on a singular pre-treatment covariance instead of using
# MASS::ginv(), unlike 03/04/05/09/13's compute_pretrend_test()), so the
# generalized-inverse clause in the other scripts' PRETREND_NOTE_AGG is
# dropped here rather than claimed falsely.
PRETREND_NOTE_AGG <- function(min_e, max_e, k) sprintf(
  "Pre-trend $\\chi^2$ / $p$: joint Wald test on the aggregated pre-treatment event-time coefficients ($%d \\leq e \\leq %d$; %d restrictions), using the influence-function covariance of the dynamic aggregation from the analytical (non-bootstrap) fit",
  min_e, max_e, k)
PRETREND_NOTE_DID <- "\\texttt{did} pre-test $p$: \\texttt{did}'s built-in Wald test over all pre-period $ATT(g,t)$ cells against each cohort's $g-1$ base year"
PRETREND_NOTE <- function(min_e, max_e, k) paste0(PRETREND_NOTE_AGG(min_e, max_e, k), ". ", PRETREND_NOTE_DID)

compute_pretrend_wald <- function(agg_d) {
  if (is.null(agg_d)) return(list(stat = NA_real_, pval = NA_real_, df = 0L))

  keep    <- which(!is.na(agg_d$se.egt) & agg_d$se.egt > 1e-10)
  pre_pos <- which(agg_d$egt[keep] < 0)
  if (length(pre_pos) == 0L) return(list(stat = NA_real_, pval = NA_real_, df = 0L))

  pre_beta <- agg_d$att.egt[keep][pre_pos]
  IF       <- agg_d$inf.function$dynamic.inf.func.e

  if (is.null(IF) || ncol(IF) != length(agg_d$egt)) {
    warning("compute_pretrend_wald(): dynamic.inf.func.e missing or the wrong ",
            "shape (expected ncol == length(egt)) -- falling back to the ",
            "diagonal-covariance approximation for this test only.")
    pre_se    <- agg_d$se.egt[keep][pre_pos]
    sigma_pre <- diag(pre_se^2, nrow = length(pre_pos))
  } else {
    n          <- nrow(IF)
    # Var-cov of the full event-time ATT vector from the multiplier/analytical
    # influence function: Sigma = E[IF IF'] / n (matches the did package's own
    # internal construction of se.egt = sqrt(diag(Sigma))).
    sigma_full <- crossprod(IF) / n^2
    sigma_pre  <- sigma_full[keep, keep][pre_pos, pre_pos, drop = FALSE]
  }

  W <- tryCatch(
    as.numeric(t(pre_beta) %*% solve(sigma_pre) %*% pre_beta),
    error = function(e) {
      message("  compute_pretrend_wald(): sigma_pre is singular -- ", conditionMessage(e))
      NA_real_
    }
  )
  list(
    stat = if (is.na(W)) NA_real_ else round(W, 3),
    pval = if (is.na(W)) NA_real_ else round(pchisq(W, df = length(pre_pos), lower.tail = FALSE), 3),
    df   = length(pre_pos)
  )
}

# ==============================================================================
# SECTION 4. Shared estimation helper: fit_main_spec()
#
# Reproduces the main-spec att_gt() call from 03_main_results.R §9-10 exactly:
#   xformla = ~ ge_est + log_population, est_method = "dr" (or "reg" for the
#   retained-cohorts robustness spec), control_group = "nevertreated",
#   anticipation = 0, base_period = "universal", panel = TRUE,
#   allow_unbalanced_panel = TRUE. A separate analytical fit (bstrap = FALSE)
#   feeds the pre-trend test, exactly as in 03_main_results.R's make_wide_table()
#   (so bootstrap SE noise never touches the pre-trend statistic). When the
#   bootstrap fit is not requested (bstrap_flag = FALSE, e.g. the retained-
#   cohorts spec), the analytical fit IS the main fit — no duplicate
#   computation, matching 03_main_results.R's retain_thin branch.
#   Seed rule (reproduces the published SEs): set.seed(1242) immediately before the bootstrap fit.
# ==============================================================================

fit_main_spec <- function(did_panel_in, yname, est_method = "dr",
                          bstrap_flag = TRUE, biters = 999L) {

  gt_analytical <- tryCatch(
    att_gt(
      yname = yname, tname = "year", idname = "country_id", gname = "cohort_year",
      xformla = ~ ge_est + log_population, data = did_panel_in,
      est_method = est_method, bstrap = FALSE, cband = FALSE,
      control_group = "nevertreated", anticipation = 0,
      base_period = "universal", panel = TRUE, allow_unbalanced_panel = TRUE
    ),
    error = function(e) {
      message("  att_gt (analytical) failed for ", yname, ": ", conditionMessage(e)); NULL
    }
  )

  if (bstrap_flag) {
    set.seed(1242L)  # Seed rule (reproduces the published SEs): immediately before the bootstrap fit
    gt_main <- tryCatch(
      att_gt(
        yname = yname, tname = "year", idname = "country_id", gname = "cohort_year",
        xformla = ~ ge_est + log_population, data = did_panel_in,
        est_method = est_method, bstrap = TRUE, biters = biters, cband = FALSE,
        control_group = "nevertreated", anticipation = 0,
        base_period = "universal", panel = TRUE, allow_unbalanced_panel = TRUE
      ),
      error = function(e) {
        message("  att_gt (bootstrap) failed for ", yname, ": ", conditionMessage(e)); NULL
      }
    )
  } else {
    gt_main <- gt_analytical  # reuse — no duplicate computation
  }
  if (is.null(gt_main)) return(NULL)

  agg_s <- tryCatch(aggte(gt_main, type = "simple", na.rm = TRUE), error = function(e) NULL)
  agg_d <- tryCatch(
    aggte(gt_main, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf),
    error = function(e) NULL
  )
  agg_d_analytical <- if (!is.null(gt_analytical)) tryCatch(
    aggte(gt_analytical, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf),
    error = function(e) NULL
  ) else NULL

  list(gt_main = gt_main, gt_analytical = gt_analytical,
       agg_s = agg_s, agg_d = agg_d, agg_d_analytical = agg_d_analytical)
}

# ------------------------------------------------------------------------
# Helper: summarize_fit() -- formats one fit_main_spec() output into the
# common row set used by every wide table below.
# ------------------------------------------------------------------------
summarize_fit <- function(fit, did_panel_in, raw_var, is_share = FALSE,
                          estimator_label = "CS (2021), doubly-robust",
                          bootstrap_label = "Yes (multiplier, 999 reps)") {

  na_out <- list(att_fmt = "---", se_fmt = "---", t_fmt = "---",
                mean_pre_fmt = "---", effect_fmt = "---",
                n_obs = "---", n_country = "---", n_treated_fmt = "---",
                zero_fmt = "---", pt_stat_fmt = "---", pt_pval_fmt = "---",
                wpval_fmt = "---", estimator_label = estimator_label,
                bootstrap_label = bootstrap_label,
                att_num = NA_real_, se_num = NA_real_, agg_s = NULL)
  if (is.null(fit) || is.null(fit$agg_s)) return(na_out)

  att <- fit$agg_s$overall.att
  se  <- fit$agg_s$overall.se
  t_v <- if (!is.na(se) && se > 1e-12) att / se else NA_real_
  stars <- if (is.na(t_v)) "" else
    if (abs(t_v) > 2.576) "***" else if (abs(t_v) > 1.960) "**" else
    if (abs(t_v) > 1.645) "*"   else ""

  pt    <- compute_pretrend_wald(fit$agg_d_analytical)
  wpval <- fit$gt_main$Wpval %||% NA_real_

  n_obs     <- nrow(did_panel_in)
  n_country <- n_distinct(did_panel_in$country_id)
  n_treated <- did_panel_in %>% filter(cohort_year > 0) %>% distinct(country_id) %>% nrow()

  pre_rows <- did_panel_in %>% filter(cohort_year > 0, year < cohort_year, !is.na(.data[[raw_var]]))
  mean_pre <- mean(pre_rows[[raw_var]], na.rm = TRUE)

  # Exact `== 0` is safe here: raw_var is a sum of CRS commitment rows built
  # with if_else(is.na(.), 0, .) in 01_prepare_data.R -- a "no rows matched"
  # cell is a literal 0.0, not a near-zero float from cancellation, so this is
  # not the float-equality pattern the numerical-discipline rule warns against.
  zero_share <- 100 * mean(did_panel_in[[raw_var]] == 0, na.rm = TRUE)

  mean_pre_fmt <- if (is_share) sprintf("%.2f pp", mean_pre) else sprintf("%.1f", mean_pre)
  # Implied USD effect only where the ATT is significant at 5% (|t| > 1.960),
  # matching the suppression rule of the main-spec tables (03/04/05).
  is_sig5      <- !is.na(att) && !is.na(se) && se > 0 && abs(att / se) > 1.960
  effect_fmt   <- if (is_share || is.na(mean_pre) || !is_sig5) "---" else
    sprintf("%.1f", (exp(att) - 1) * mean_pre)

  list(
    att_fmt       = paste0(sprintf("%.4f", att), stars),
    se_fmt        = sprintf("(%.4f)", se),
    t_fmt         = sprintf("%.3f", t_v),
    mean_pre_fmt  = mean_pre_fmt,
    effect_fmt    = effect_fmt,
    n_obs         = format(n_obs, big.mark = ","),
    n_country     = as.character(n_country),
    n_treated_fmt = as.character(n_treated),
    zero_fmt      = sprintf("%.1f", zero_share),
    pt_stat_fmt   = if (is.na(pt$stat)) "---" else sprintf("%.3f", pt$stat),
    pt_pval_fmt   = if (is.na(pt$pval)) "---" else sprintf("%.3f", pt$pval),
    wpval_fmt     = if (is.na(wpval)) "---" else sprintf("%.3f", wpval),
    estimator_label = estimator_label,
    bootstrap_label = bootstrap_label,
    att_num       = att,
    se_num        = se,
    agg_s         = fit$agg_s
  )
}

# ------------------------------------------------------------------------
# Helper: build_wide_table() -- assembles the common row block (metrics x
# outcome columns) into a bare tabular, via xtable, matching 03/04/05's style.
# ------------------------------------------------------------------------
build_wide_table <- function(stats_list, col_labels, extra_row_lines = character(0)) {

  row_labels <- c(
    "ATT", "SE", "$t$-statistic",
    "Mean (pre-treat.)", "Implied effect (USD M, exp(ATT)$-$1 $\\times$ pre-treat.\\ mean)",
    "Observations", "Countries", "N treated countries",
    "Share of zero/undefined cells (\\%)",
    "Pre-trend $\\chi^2$", "Pre-trend $p$",
    "\\texttt{did} pre-test $p$",
    "\\midrule Estimator", "Control group", "Bootstrap SE"
  )

  tab_wide <- data.frame(` ` = row_labels, check.names = FALSE, stringsAsFactors = FALSE)
  for (i in seq_along(stats_list)) {
    s <- stats_list[[i]]
    col <- c(
      s$att_fmt, s$se_fmt, s$t_fmt, s$mean_pre_fmt, s$effect_fmt,
      s$n_obs, s$n_country, s$n_treated_fmt, s$zero_fmt,
      s$pt_stat_fmt, s$pt_pval_fmt, s$wpval_fmt,
      s$estimator_label, "Never-treated", s$bootstrap_label
    )
    tab_wide[[esc_header(col_labels[i])]] <- col
  }

  xtab <- xtable(tab_wide)
  align(xtab) <- paste0("ll", paste(rep("c", length(col_labels)), collapse = ""))

  raw_lines <- capture.output(
    print(xtab, include.rownames = FALSE, booktabs = TRUE,
          sanitize.text.function = identity, size = "\\small", floating = FALSE)
  )

  if (length(extra_row_lines) > 0) {
    bottom_idx <- which(grepl("^\\s*\\\\bottomrule", raw_lines))[1]
    stopifnot(!is.na(bottom_idx))
    raw_lines <- append(raw_lines, c("\\midrule", extra_row_lines), after = bottom_idx - 1L)
  }
  raw_lines
}

message("\n=== Section 4 helpers ready (fit_main_spec / summarize_fit / build_wide_table) ===\n")

# ==============================================================================
# SECTION 5. Reference fits (main spec, did_panel_tab_main): headline and
# adaptation-share outcome. These feed Part A (headline column + difference
# test), Part C (reconciliation), and the exact-reproduction check.
# ==============================================================================

message("\n=== Section 5: reference fits (log_commits, share_adapt) ===\n")

fit_headline <- fit_main_spec(did_panel_tab_main, "log_commits",
                              est_method = "dr", bstrap_flag = TRUE, biters = 999L)
fit_share    <- fit_main_spec(did_panel_tab_main, "share_adapt",
                              est_method = "dr", bstrap_flag = TRUE, biters = 999L)

stopifnot(!is.null(fit_headline), !is.null(fit_headline$agg_s),
         !is.null(fit_share),    !is.null(fit_share$agg_s))

message(sprintf("  log_commits : ATT = %.4f  SE = %.4f", fit_headline$agg_s$overall.att, fit_headline$agg_s$overall.se))
message(sprintf("  share_adapt : ATT = %.4f  SE = %.4f", fit_share$agg_s$overall.att,    fit_share$agg_s$overall.se))

# HARD REQUIREMENT: the headline column must reproduce Table 2
# (output/tables/cohorts_dropped/att_combined_wide.tex) exactly, since this
# script re-fits the same main spec on the same panel with the same seed.
# If it does not, something about the panel construction or estimator call
# has silently diverged from 03_main_results.R -- stop rather than publish a
# table that looks like Table 2 but is not.
# Reference values are read from the fits 03_main_results.R stored in
# output/fits/ (one fit, one SE), never hardcoded: a panel revision (e.g. the
# 2026-09-15 restoration of recipient code 860) must not trip this check.
ref_head  <- readRDS(here("output", "fits", "headline_adaptation_dr_bs.rds"))
ref_share <- readRDS(here("output", "fits", "headline_share_dr_bs.rds"))
tol_att <- 1e-3
tol_se  <- 2e-3
att_ok <- abs(fit_headline$agg_s$overall.att - ref_head$att)  < tol_att &&
  abs(fit_headline$agg_s$overall.se  - ref_head$se)   < tol_se &&
  abs(fit_share$agg_s$overall.att    - ref_share$att) < tol_att &&
  abs(fit_share$agg_s$overall.se     - ref_share$se)  < tol_se
if (!att_ok) {
  stop(sprintf(paste0(
    "Headline reproduction check FAILED.\n",
    "  log_commits : got ATT = %.4f, SE = %.4f (stored 03 fit: ", sprintf("%.4f / %.4f", ref_head$att, ref_head$se), ")\n",
    "  share_adapt : got ATT = %.4f, SE = %.4f (stored 03 fit: ", sprintf("%.4f / %.4f", ref_share$att, ref_share$se), ")\n",
    "  Panel construction or att_gt() call has diverged from 03_main_results.R -- ",
    "investigate before trusting any table produced by this script."),
    fit_headline$agg_s$overall.att, fit_headline$agg_s$overall.se,
    fit_share$agg_s$overall.att,    fit_share$agg_s$overall.se))
}
message("  Headline reproduction check: PASSED (matches Table 2 to within tolerance).")

# ==============================================================================
# SECTION 6. Part A -- principal-marker-only outcome, full Table-2 analogue
# ==============================================================================

message("\n=== Section 6 (Part A): principal-marker-only outcome ===\n")

fit_principal <- fit_main_spec(did_panel_tab_main, "lcommitments_principal",
                               est_method = "dr", bstrap_flag = TRUE, biters = 999L)
stopifnot(!is.null(fit_principal), !is.null(fit_principal$agg_s))
message(sprintf("  lcommitments_principal : ATT = %.4f  SE = %.4f",
                fit_principal$agg_s$overall.att, fit_principal$agg_s$overall.se))

stats_headline  <- summarize_fit(fit_headline,  did_panel_tab_main, raw_var = "commitments")
stats_principal <- summarize_fit(fit_principal, did_panel_tab_main, raw_var = "commitments_principal")

# --- Difference row: headline minus principal ---
# The two outcomes are estimated on the identical panel (did_panel_tab_main)
# with no missing yname values (checked in Section 2), so the simple-ATT
# influence-function vectors are aligned unit-for-unit and
# Var(theta_head - theta_princ) = Var(IF_head - IF_princ) / n is the correct
# (non-independence) SE. If for any reason the two IF vectors do not align
# (different length), fall back to the independence formula sqrt(se_a^2+se_b^2)
# with an explicit caveat: because principal-only commitments are a strict
# subset of headline (principal+significant) commitments, the two series are
# mechanically positively correlated, so ignoring that covariance likely
# overstates -- but does not certainly overstate, hence "conservative,
# direction not separately verified" -- the true SE of the difference.
# --- Empirical verification of the unit-level-IF SE scaling convention -----
# Before trusting sqrt(var(IFa - IFb) / n_units) below, confirm empirically
# which scaling convention did's own overall.se uses: sd(IF)/sqrt(n) (Bessel-
# corrected sample variance, stats::var()) or sqrt(mean(IF^2)/n) (the
# uncentered second-moment convention used by compute_pretrend_test's
# crossprod(IF)/n^2 in 03/04/05/10/11). Uses the headline's ANALYTICAL
# (bstrap = FALSE) fit so the check is not contaminated by multiplier-
# bootstrap resampling noise -- the bootstrap SE (0.1245) is expected to
# differ from either analytical formula.
stopifnot(!is.null(fit_headline$gt_analytical))
agg_s_head_analytical <- aggte(fit_headline$gt_analytical, type = "simple", na.rm = TRUE)
IF_head <- agg_s_head_analytical$inf.function$simple.att
stopifnot(!is.null(IF_head))
n_head      <- length(IF_head)
se_var_form <- sqrt(stats::var(IF_head) / n_head)
se_msq_form <- sqrt(mean(IF_head^2) / n_head)
se_official <- agg_s_head_analytical$overall.se
message(sprintf(paste0(
  "  SE convention check (headline, ANALYTICAL fit, bstrap = FALSE): official overall.se ",
  "= %.6f | sqrt(var(IF)/n) = %.6f | sqrt(mean(IF^2)/n) = %.6f"),
  se_official, se_var_form, se_msq_form))
match_var <- abs(se_var_form - se_official) < 1e-6
match_msq <- abs(se_msq_form - se_official) < 1e-6
stopifnot(match_var || match_msq)
message(sprintf("  -> matches the %s convention within 1e-6",
                if (match_var && match_msq) "var(IF)/n AND mean(IF^2)/n"
                else if (match_var) "sqrt(var(IF)/n)"
                else "sqrt(mean(IF^2)/n)"))

compute_att_difference <- function(fit_a, fit_b) {
  att_a <- fit_a$agg_s$overall.att; att_b <- fit_b$agg_s$overall.att
  se_a  <- fit_a$agg_s$overall.se;  se_b  <- fit_b$agg_s$overall.se
  diff  <- att_a - att_b
  IFa   <- fit_a$agg_s$inf.function$simple.att
  IFb   <- fit_b$agg_s$inf.function$simple.att
  if (!is.null(IFa) && !is.null(IFb) && length(IFa) == length(IFb)) {
    n_units <- length(IFa)
    # did's convention is the uncentred mean-square form (verified above: it
    # reproduces overall.se exactly; the centred var() form is ~0.4% off).
    se_diff <- sqrt(mean((IFa - IFb)^2) / n_units)
    method  <- paste0("aligned unit-level influence functions (n = ", n_units, ")")
  } else {
    se_diff <- sqrt(se_a^2 + se_b^2)
    method  <- paste0("independence approximation (IF vectors did not align; ",
                      "outcomes overlap mechanically -- principal is a subset ",
                      "of headline -- so this SE is conservative in direction ",
                      "but its magnitude relative to the true SE is not separately verified)")
  }
  z    <- if (se_diff > 1e-12) diff / se_diff else NA_real_
  pval <- if (is.na(z)) NA_real_ else 2 * pnorm(-abs(z))
  list(diff = diff, se = se_diff, z = z, pval = pval, method = method)
}

diff_hp <- compute_att_difference(fit_headline, fit_principal)
message(sprintf("  Difference (headline - principal) = %.4f, SE = %.4f, z = %.3f, p = %.3f [%s]",
                diff_hp$diff, diff_hp$se, diff_hp$z, diff_hp$pval, diff_hp$method))

extra_lines_r2 <- c(
  sprintf("\\multicolumn{3}{l}{Difference (headline $-$ principal) = %.4f, SE = %.4f, $z$ = %.3f, $p$ = %.3f} \\\\",
          diff_hp$diff, diff_hp$se, diff_hp$z, diff_hp$pval),
  sprintf("\\multicolumn{3}{l}{\\footnotesize Difference SE method: %s.} \\\\", diff_hp$method)
)

tab_lines_r2 <- build_wide_table(
  stats_list = list(stats_headline, stats_principal),
  col_labels = c("Principal + significant (headline)", "Principal only"),
  extra_row_lines = extra_lines_r2
)

notes_r2 <- paste0(
 "CS\\,(2021) DR, never-treated controls; multiplier-bootstrap SE (999 reps, seed 1242). ",
  "``Headline'' = log1p(commitments), Rio marker $\\in\\{1,2\\}$; ``Principal only'' = ",
  "log1p(commitments\\_principal), Rio marker $=2$ only. ",
  PRETREND_NOTE(-5L, -2L, 4L), ". ",
  "Implied effects: back-transform on pre-treatment mean; suppressed if insig.\\ at 5\\%. ",
  "* $p<0.10$, ** $p<0.05$, *** $p<0.01$"
)
source_r2 <- "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"

out_path_r2 <- here("output", "tables", "principal_share", "principal_wide.tex")
write_tex_float(
  out_path_r2,
  "Effect of NAP adoption on adaptation finance: headline vs.\\ principal-marker-only outcome",
  "tab:principal_wide", tab_lines_r2, notes_r2, source_r2
)
copy_to_paper(out_path_r2, "Tables")

# --- fig_principal_es.png: combined event study, headline vs principal ------
dyn_headline  <- data.frame(outcome = "Principal + significant (headline)",
                            event_time = fit_headline$agg_d$egt,
                            ATT = fit_headline$agg_d$att.egt, SE = fit_headline$agg_d$se.egt,
                            Lower = fit_headline$agg_d$att.egt - fit_headline$agg_d$crit.val.egt * fit_headline$agg_d$se.egt,
                            Upper = fit_headline$agg_d$att.egt + fit_headline$agg_d$crit.val.egt * fit_headline$agg_d$se.egt)
dyn_principal <- data.frame(outcome = "Principal only",
                            event_time = fit_principal$agg_d$egt,
                            ATT = fit_principal$agg_d$att.egt, SE = fit_principal$agg_d$se.egt,
                            Lower = fit_principal$agg_d$att.egt - fit_principal$agg_d$crit.val.egt * fit_principal$agg_d$se.egt,
                            Upper = fit_principal$agg_d$att.egt + fit_principal$agg_d$crit.val.egt * fit_principal$agg_d$se.egt)
dyn_r2 <- bind_rows(dyn_headline, dyn_principal) %>%
  mutate(outcome = factor(outcome, levels = c("Principal + significant (headline)", "Principal only")))

p_r2 <- ggplot(dyn_r2, aes(x = event_time, y = ATT, colour = outcome, shape = outcome, group = outcome)) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
  geom_vline(xintercept = -0.5, colour = "grey30", linetype = "dotted") +
  geom_linerange(aes(ymin = Lower, ymax = Upper), position = position_dodge(width = 0.4),
                linewidth = 0.6, alpha = 0.8) +
  geom_point(size = 2.5, position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = c("Principal + significant (headline)" = "#2E86C1", "Principal only" = "#C0392B")) +
  scale_shape_manual(values = c("Principal + significant (headline)" = 16, "Principal only" = 17)) +
  labs(title = NULL, subtitle = NULL, caption = NULL,
      x = "Years relative to NAP adoption", y = "ATT estimate (log points)",
      colour = NULL, shape = NULL) +
  theme_minimal() +
  theme(text = element_text(family = "serif", size = 12), legend.position = "bottom",
       legend.text = element_text(size = 10), panel.grid.minor = element_blank())

fig_path_r2 <- here("output", "figures", "principal_share", "fig_principal_es.png")
ggsave(fig_path_r2, p_r2, width = 12, height = 7, dpi = 300)
message("Saved: ", fig_path_r2)
copy_to_paper(fig_path_r2, "Figures")

# ------------------------------------------------------------------------
# tab:principal_retained -- principal-only under the retained-cohorts spec
# (all cohorts kept, est_method = "reg", bstrap = FALSE), headline shown
# alongside for direct comparison: is the cohort-inclusion attenuation
# confined to the significant-marker component of the headline outcome?
# ------------------------------------------------------------------------

message("\n=== Section 6b: retained-cohorts spec (Part A) ===\n")

fit_headline_retained  <- fit_main_spec(did_panel_full, "log_commits",
                                        est_method = "reg", bstrap_flag = FALSE)
fit_principal_retained <- fit_main_spec(did_panel_full, "lcommitments_principal",
                                        est_method = "reg", bstrap_flag = FALSE)
stopifnot(!is.null(fit_headline_retained), !is.null(fit_principal_retained))
message(sprintf("  log_commits (retained)             : ATT = %.4f  SE = %.4f",
                fit_headline_retained$agg_s$overall.att, fit_headline_retained$agg_s$overall.se))
message(sprintf("  lcommitments_principal (retained)  : ATT = %.4f  SE = %.4f",
                fit_principal_retained$agg_s$overall.att, fit_principal_retained$agg_s$overall.se))

stats_headline_retained  <- summarize_fit(fit_headline_retained, did_panel_full,
                                          raw_var = "commitments",
                                          estimator_label = "CS (2021), outcome regression",
                                          bootstrap_label = "No (analytical SE)")
stats_principal_retained <- summarize_fit(fit_principal_retained, did_panel_full,
                                          raw_var = "commitments_principal",
                                          estimator_label = "CS (2021), outcome regression",
                                          bootstrap_label = "No (analytical SE)")

tab_lines_r2b <- build_wide_table(
  stats_list = list(stats_headline_retained, stats_principal_retained),
  col_labels = c("Headline (retained cohorts)", "Principal only (retained cohorts)")
)

notes_r2b <- paste0(
 "CS\\,(2021) OR, never-treated controls; analytical (IF) SE; all cohorts retained ",
  "(2015--2024, incl.\\ $<$ ", thin_threshold, " treated units). Pre-trend statistics: see ",
  "notes to Table~\\ref{tab:principal_wide}. \\texttt{did}'s pre-test not computed: ",
  "pre-treatment covariance singular with all cohorts retained. ",
  "Implied effects: back-transform on pre-treatment mean; suppressed if insig.\\ at 5\\%. ",
  "* $p<0.10$, ** $p<0.05$, *** $p<0.01$"
)
out_path_r2b <- here("output", "tables", "principal_share", "principal_retained.tex")
write_tex_float(
  out_path_r2b,
  "Retained-cohorts robustness: headline vs.\\ principal-marker-only outcome",
  "tab:principal_retained", tab_lines_r2b, notes_r2b, source_r2
)
copy_to_paper(out_path_r2b, "Tables")

message("\n=== Section 6 complete ===\n")

# ==============================================================================
# SECTION 7. Part B -- within-country adaptation share
# share_within and logit_within, main-spec settings, did_panel_tab_main.
# ==============================================================================

message("\n=== Section 7 (Part B): within-country adaptation share ===\n")

fit_share_within <- fit_main_spec(did_panel_tab_main, "share_within",
                                  est_method = "dr", bstrap_flag = TRUE, biters = 999L)
fit_logit_within <- fit_main_spec(did_panel_tab_main, "logit_within",
                                  est_method = "dr", bstrap_flag = TRUE, biters = 999L)
stopifnot(!is.null(fit_share_within), !is.null(fit_share_within$agg_s),
         !is.null(fit_logit_within), !is.null(fit_logit_within$agg_s))

message(sprintf("  share_within : ATT = %.4f pp  SE = %.4f", fit_share_within$agg_s$overall.att, fit_share_within$agg_s$overall.se))
message(sprintf("  logit_within : ATT = %.4f     SE = %.4f", fit_logit_within$agg_s$overall.att, fit_logit_within$agg_s$overall.se))

stats_share_within <- summarize_fit(fit_share_within, did_panel_tab_main,
                                    raw_var = "share_within", is_share = TRUE)
stats_logit_within <- summarize_fit(fit_logit_within, did_panel_tab_main,
                                    raw_var = "share_within", is_share = TRUE)

tab_lines_r4 <- build_wide_table(
  stats_list = list(stats_share_within, stats_logit_within),
  col_labels = c("Within-country share (%, level)", "Within-country share (logit)")
)

notes_r4 <- paste0(
  "CS\\,(2021) DR, never-treated controls; multiplier-bootstrap SE (999 reps, seed 1242). ",
  "share\\_within $=100 \\times$ commitments/commitments\\_all; logit\\_within ",
  "$= \\log[(\\text{commitments}+0.5)/(\\text{commitments\\_all}-\\text{commitments}+0.5)]$; ",
  "both NA if commitments\\_all $=0$. Mean/zero-share rows use share\\_within (pp) levels. ",
  "Pre-trend statistics: see Table~\\ref{tab:principal_wide} notes. ",
  "* $p<0.10$, ** $p<0.05$, *** $p<0.01$"
)
source_r4 <- "OECD CRS (Rio adaptation markers)"

out_path_r4 <- here("output", "tables", "principal_share", "within_share_wide.tex")
write_tex_float(
  out_path_r4,
  "Effect of NAP adoption on the within-country adaptation share of reported commitments",
  "tab:within_share_wide", tab_lines_r4, notes_r4, source_r4
)
copy_to_paper(out_path_r4, "Tables")

# --- fig_within_share_es.png: faceted (scales differ, level % vs logit) -----
dyn_share_within <- data.frame(outcome = "Level (%)",
                               event_time = fit_share_within$agg_d$egt,
                               ATT = fit_share_within$agg_d$att.egt, SE = fit_share_within$agg_d$se.egt,
                               Lower = fit_share_within$agg_d$att.egt - fit_share_within$agg_d$crit.val.egt * fit_share_within$agg_d$se.egt,
                               Upper = fit_share_within$agg_d$att.egt + fit_share_within$agg_d$crit.val.egt * fit_share_within$agg_d$se.egt)
dyn_logit_within <- data.frame(outcome = "Logit",
                               event_time = fit_logit_within$agg_d$egt,
                               ATT = fit_logit_within$agg_d$att.egt, SE = fit_logit_within$agg_d$se.egt,
                               Lower = fit_logit_within$agg_d$att.egt - fit_logit_within$agg_d$crit.val.egt * fit_logit_within$agg_d$se.egt,
                               Upper = fit_logit_within$agg_d$att.egt + fit_logit_within$agg_d$crit.val.egt * fit_logit_within$agg_d$se.egt)
dyn_r4 <- bind_rows(dyn_share_within, dyn_logit_within) %>%
  mutate(outcome = factor(outcome, levels = c("Level (%)", "Logit")))

p_r4 <- ggplot(dyn_r4, aes(x = event_time, y = ATT)) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
  geom_vline(xintercept = -0.5, colour = "grey30", linetype = "dotted") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, fill = "#2E86C1") +
  geom_line(colour = "#2E86C1", linewidth = 0.8) +
  geom_point(colour = "#2E86C1", shape = 16, size = 2.2) +
  facet_wrap(~ outcome, scales = "free_y") +
  labs(title = NULL, subtitle = NULL, caption = NULL,
      x = "Years relative to NAP adoption", y = "ATT estimate") +
  theme_minimal() +
  theme(text = element_text(family = "serif", size = 12),
       strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())

fig_path_r4 <- here("output", "figures", "principal_share", "fig_within_share_es.png")
ggsave(fig_path_r4, p_r4, width = 12, height = 6, dpi = 300)
message("Saved: ", fig_path_r4)
copy_to_paper(fig_path_r4, "Figures")

message("\n=== Section 7 complete ===\n")

# ==============================================================================
# SECTION 8. Part C -- global share reconciliation
# ==============================================================================

message("\n=== Section 8 (Part C): share_adapt construction + reconciliation ===\n")

# --- C.1: exact construction of share_adapt (as in 01_prepare_data.R) --------
# code/01_prepare_data.R §18 (l.966-972):
#   setDT(aggregated)
#   aggregated[, global_adapt_t := sum(commitments, na.rm = TRUE), by = year]
#   aggregated[, share_adapt := fifelse(global_adapt_t > 0,
#                                        commitments / global_adapt_t * 100, NA_real_)]
# `aggregated` at that point is the SAME object later written verbatim to
# data/processed/simple_panel_wgi.csv (l.980) -- i.e. this script's
# did_panel_full. There is no filter() on `aggregated` between the listwise
# deletion on ge_est + population (§15, l.835-838, which fixes the panel at
# 144 recipient countries) and the share_adapt computation (§18) -- confirmed
# by grep (no `filter(` on the bare `aggregated` object in that range; the
# filter() calls that DO appear there operate on the raw CRS donor-level
# lists used to build commitments_all / commitments_principal / commitments_oda
# etc., not on `aggregated` itself). So:
#   DENOMINATOR = sum(commitments) over the 144-country panel that has
#     non-missing WGI governance-effectiveness AND non-missing population data
#     (the listwise-deletion sample), by calendar year -- NOT the full universe
#     of all OECD CRS recipient countries (smaller/data-poor recipients that
#     lack WGI or population coverage are excluded before this point), and NOT
#     restricted to the 126-country / 40-adopter DiD ESTIMATION sample (that
#     further restriction -- dropping thin cohorts -- happens only downstream,
#     in 03/04/05/07's did_panel construction).
n_pool_by_year <- did_panel_full %>% group_by(year) %>%
  summarise(n_countries_pool = n_distinct(recipient_name), .groups = "drop")
message("--- share_adapt denominator: number of countries in the pool, by year ---")
for (i in seq_len(nrow(n_pool_by_year)))
  message(sprintf("  %d : %d countries", n_pool_by_year$year[i], n_pool_by_year$n_countries_pool[i]))

# Internal consistency check: recompute share_adapt from raw commitments and
# compare to the column shipped in simple_panel_wgi.csv (should match to
# floating precision if the description above is correct).
pool_denom <- did_panel_full %>% group_by(year) %>%
  summarise(global_adapt_t = sum(commitments, na.rm = TRUE), .groups = "drop")
share_check <- did_panel_full %>% left_join(pool_denom, by = "year") %>%
  mutate(share_recomputed = if_else(global_adapt_t > 0, 100 * commitments / global_adapt_t, NA_real_))
max_diff_share <- max(abs(share_check$share_recomputed - share_check$share_adapt), na.rm = TRUE)
message(sprintf("  Recomputed share_adapt vs shipped column: max abs diff = %.8f pp (expect ~0).", max_diff_share))
if (max_diff_share > 0.01) {
  warning("share_adapt reconstruction differs from the shipped column by > 0.01pp -- ",
         "the denominator description above may not exactly match 01_prepare_data.R's ",
         "current logic; verify before citing the C.1 finding.")
}

pretreat_share_adapt_treated <- did_panel_full %>%
  filter(cohort_year > 0, year < cohort_year) %>% pull(share_adapt)
message(sprintf("  Pre-treatment mean share_adapt among treated units: %.3f pp (N = %d).",
                mean(pretreat_share_adapt_treated, na.rm = TRUE), length(pretreat_share_adapt_treated)))

# --- C.2: reconciliation exhibit --------------------------------------------
# (a) observed pooled share held by the 58 ever-adopters (any cohort_year > 0,
#     including thin cohorts not in the main estimation sample)
# (b) observed pooled share held by the 40 main-sample adopters (cohort not thin)
# (c) counterfactual pooled share for the 40 main-sample adopters implied by the
#     main-spec share ATT: for each treated country-year, subtract the
#     event-time dynamic ATT (fit_share$agg_d, from Section 5 -- same panel,
#     same seed as Table 2) from the observed share, leaving pre-treatment
#     and out-of-window cells unadjusted.
adopters58_by_year <- did_panel_full %>% filter(cohort_year > 0) %>%
  group_by(year) %>% summarise(share_a_observed = sum(share_adapt, na.rm = TRUE), .groups = "drop")

att_by_e <- data.frame(event_time = fit_share$agg_d$egt, att = fit_share$agg_d$att.egt)

adopters40_cf <- did_panel_tab_main %>% filter(cohort_year > 0) %>%
  mutate(event_time = year - cohort_year) %>%
  left_join(att_by_e, by = "event_time") %>%
  mutate(att_adj = if_else(event_time >= 0 & !is.na(att), att, 0))

n_treated_cells_r7 <- sum(adopters40_cf$event_time >= 0 & !is.na(adopters40_cf$att))
message(sprintf("  Treated country-year cells (40-adopter sample, event_time >= 0, ATT available): %d",
                n_treated_cells_r7))

recon_by_year <- adopters40_cf %>%
  group_by(year) %>%
  summarise(
    share_b_observed       = sum(share_adapt, na.rm = TRUE),
    total_att_adj          = sum(att_adj, na.rm = TRUE),
    share_c_counterfactual = share_b_observed - total_att_adj,
    .groups = "drop"
  ) %>%
  left_join(adopters58_by_year, by = "year") %>%
  arrange(year)

reallocation_correct <- sum(adopters40_cf$att_adj, na.rm = TRUE)
reallocation_naive    <- n_distinct(adopters40_cf$recipient_name) * fit_share$agg_s$overall.att

message(sprintf(
  "  Naive arithmetic  : %d adopters x ATT %.4f pp = %.2f pp",
  n_distinct(adopters40_cf$recipient_name), fit_share$agg_s$overall.att, reallocation_naive))
message(sprintf(
  "  Correct arithmetic: sum of event-time ATT over %d treated cells = %.2f pp",
  n_treated_cells_r7, reallocation_correct))
message(sprintf(paste0(
  "  The naive figure treats the pooled (cell-averaged) simple ATT as if it applied once ",
  "per COUNTRY; it is applied %d times (once per treated cell), and most adopters contribute ",
  "more than one post-adoption year, so naive != correct whenever n_treated_cells != n_adopters."),
  n_treated_cells_r7))

# The treated-cell-weighted "total reallocation" above (reallocation_correct)
# sums att_adj over all n_treated_cells_r7 treated country-years pooled across
# the whole 2021-2024 panel; it is a cumulative, multi-year quantity, not a
# single-year magnitude, and is not bounded by 100 pp (a sum of per-cell
# percentage-point ATTs across years/countries has no such bound). For a
# single-year magnitude, report the per-year reallocation below instead: the
# sum of ATT-adjusted cells within each calendar year (recon_by_year$total_att_adj).
message("  Per-year reallocation (sum of ATT-adjusted treated cells within each calendar year):")
for (yr_i in seq_len(nrow(recon_by_year))) {
  message(sprintf("    %d: %.2f pp", recon_by_year$year[yr_i], recon_by_year$total_att_adj[yr_i]))
}

obs_change_b <- recon_by_year$share_b_observed[recon_by_year$year == max(recon_by_year$year)] -
  recon_by_year$share_b_observed[recon_by_year$year == min(recon_by_year$year)]
cf_change_c  <- recon_by_year$share_c_counterfactual[recon_by_year$year == max(recon_by_year$year)] -
  recon_by_year$share_c_counterfactual[recon_by_year$year == min(recon_by_year$year)]
message(sprintf(
  "  Observed change in the 40-adopter pooled share, %d -> %d: %.2f pp (implied counterfactual change: %.2f pp).",
  min(recon_by_year$year), max(recon_by_year$year), obs_change_b, cf_change_c))

# --- tab:share_reconciliation -----------------------------------------------
recon_tab <- recon_by_year %>%
  transmute(
    Year = year,
    `Ever-adopters (58), observed (\\%)`             = sprintf("%.2f", share_a_observed),
    `Main-sample adopters (40), observed (\\%)`       = sprintf("%.2f", share_b_observed),
    `Main-sample adopters (40), counterfactual (\\%)` = sprintf("%.2f", share_c_counterfactual),
    `ATT-implied reallocation, same year (pp)`        = sprintf("%.2f", total_att_adj)
  )
xtab_r7 <- xtable(recon_tab)
align(xtab_r7) <- "llcccc"
raw_lines_r7 <- capture.output(
  print(xtab_r7, include.rownames = FALSE, booktabs = TRUE,
       sanitize.text.function = identity, size = "\\small", floating = FALSE)
)
footer_r7 <- c(
  sprintf("\\multicolumn{4}{l}{Total reallocation, correct (treated-cell-weighted): %.2f pp over %d treated cells} \\\\",
         reallocation_correct, n_treated_cells_r7),
  sprintf("\\multicolumn{4}{l}{Total reallocation, naive ($n_{\\text{adopters}} \\times$ simple ATT): %d $\\times$ %.4f = %.2f pp} \\\\",
         n_distinct(adopters40_cf$recipient_name), fit_share$agg_s$overall.att, reallocation_naive),
  sprintf("\\multicolumn{4}{l}{Observed change in the 40-adopter share, %d$\\to$%d: %.2f pp; implied counterfactual change: %.2f pp} \\\\",
         min(recon_by_year$year), max(recon_by_year$year), obs_change_b, cf_change_c)
)
bottom_idx_r7 <- which(grepl("^\\s*\\\\bottomrule", raw_lines_r7))[1]
stopifnot(!is.na(bottom_idx_r7))
raw_lines_r7 <- append(raw_lines_r7, c("\\midrule", footer_r7), after = bottom_idx_r7 - 1L)

notes_r7 <- paste0(
  "share\\_adapt denominator: see main-text footnote (145-country listwise-deletion panel). ",
  "Counterfactual column subtracts the main-spec event-time ATT, cell by cell, from each ",
  "treated country's observed share. Naive $n\\times$ATT is not the correct arithmetic; ",
  "treated-cell-weighted total sums att\\_adj over all ", n_treated_cells_r7,
  " treated country-years (cumulative, not bounded by 100 pp). Last column: single-year ",
  "sum of ATT-adjusted cells"
)
source_r7 <- "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"

out_path_r7 <- here("output", "tables", "principal_share", "share_reconciliation.tex")
write_tex_float(
  out_path_r7,
  "Global adaptation-finance share: observed vs.\\ ATT-implied counterfactual reconciliation",
  "tab:share_reconciliation", raw_lines_r7, notes_r7, source_r7
)
copy_to_paper(out_path_r7, "Tables")

# --- fig_share_reconciliation.png -------------------------------------------
# (built with dplyr::bind_rows() rather than tidyr::pivot_longer() to avoid an
# extra dependency not already loaded by this script)
recon_long <- bind_rows(
  recon_by_year %>% transmute(year, value = share_a_observed, series = "Ever-adopters (58), observed"),
  recon_by_year %>% transmute(year, value = share_b_observed, series = "Main-sample adopters (40), observed"),
  recon_by_year %>% transmute(year, value = share_c_counterfactual, series = "Main-sample adopters (40), counterfactual")
) %>%
  mutate(series = factor(series, levels = c(
    "Ever-adopters (58), observed",
    "Main-sample adopters (40), observed",
    "Main-sample adopters (40), counterfactual"
  )))

p_r7 <- ggplot(recon_long, aes(x = year, y = value, colour = series, linetype = series, shape = series)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_colour_manual(values = c(
    "Ever-adopters (58), observed"               = "#2E86C1",
    "Main-sample adopters (40), observed"         = "#27AE60",
    "Main-sample adopters (40), counterfactual"   = "#C0392B"
  )) +
  scale_linetype_manual(values = c(
    "Ever-adopters (58), observed"               = "solid",
    "Main-sample adopters (40), observed"         = "solid",
    "Main-sample adopters (40), counterfactual"   = "dashed"
  )) +
  scale_shape_manual(values = c(
    "Ever-adopters (58), observed"               = 16,
    "Main-sample adopters (40), observed"         = 17,
    "Main-sample adopters (40), counterfactual"   = 15
  )) +
  scale_x_continuous(breaks = min(recon_long$year):max(recon_long$year)) +
  labs(title = NULL, subtitle = NULL, caption = NULL,
      x = "Year", y = "Share of the global adaptation-commitment pool (%)",
      colour = NULL, linetype = NULL, shape = NULL) +
  theme_minimal() +
  theme(text = element_text(family = "serif", size = 12),
       axis.text.x = element_text(angle = 45, hjust = 1),
       legend.position = "bottom", legend.text = element_text(size = 9),
       panel.grid.minor = element_blank())

fig_path_r7 <- here("output", "figures", "principal_share", "fig_share_reconciliation.png")
ggsave(fig_path_r7, p_r7, width = 12, height = 7, dpi = 300)
message("Saved: ", fig_path_r7)
copy_to_paper(fig_path_r7, "Figures")

message("\n=== Section 8 complete ===\n")

# ==============================================================================
# SECTION 9. Done
# ==============================================================================

message("\n=== 07_principal_and_share.R: complete ===\n")
