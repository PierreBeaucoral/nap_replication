# ==============================================================================
# 05_heterogeneity.R
# Heterogeneity analyses: donor type, LDC status, governance, income group.
# Paper: Beaucoral, Goujon and Marchand (2026) — §6 heterogeneity (§21-24)
#
# Inputs:
#   data/processed/simple_panel_wgi.csv
#
# Outputs (all under output/):
#   tables/heterogeneity/donor_type/att_combined_wide.tex   (§21)
#   tables/heterogeneity/ldc/het_ldc_wide.tex              (§22)
#   tables/heterogeneity/governance/het_gov_wide.tex       (§23)
#   tables/heterogeneity/income_group/het_income_wide.tex  (§24)
#   tables/heterogeneity/het_capacity.tex                  (§24b, compact H4 table)
#   figures/heterogeneity/fig_het_es_panel.png             (§24b, 2x2 of the four event studies)
#   tables/heterogeneity/het_difference_tests.tex          (§25)
# ==============================================================================

# ============================================================
# Paper-to-Code Naming Map
# ============================================================
# Paper Notation      | Code Name               | Description
# $Y_{it}$            | log_commits, etc.       | Outcome variables
# $G_i$               | cohort_year             | NAP adoption cohort (0=never)
# $ATT(g,t)$          | gt_obj                  | Group-time ATT from att_gt()
# $\hat\theta^{simp}$ | agg_s$overall.att       | Simple ATT
# $\hat\theta^{dyn}$  | agg_d$att.egt           | Dynamic ATT by event time e
# $X_{it}$            | ge_est, log_population  | Controls: WGI GE + log pop
# $\Delta = \theta_a-\theta_b$ | diff                  | Subgroup ATT difference
# $se(\Delta)$        | se_diff                 | sqrt(se_a^2 + se_b^2)
# $W$                 | joint$stat              | Wald chi2, equal income ATTs
# ============================================================

# duplicated from 03/04 §1 — keep in sync
library(data.table)
library(dplyr)
library(ggplot2)
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
library(countrycode)
# FY2013 (July 2012) World Bank income classification is read from the
# World Bank's OGHIST workbook (data/raw/oghist/OGHIST.xlsx); readxl is
# therefore a top-level dependency of this script.
library(readxl)
library(cowplot)  # 2x2 heterogeneity event-study panel (pinned in renv.lock)
source(here("code", "functions", "pretrend_test.R"))      # compute_pretrend_test()
source(here("code", "functions", "two_line_head.R"))     # two_line_head()
source(here("code", "functions", "read_headline_fit.R"))  # read_headline_fit()
source(here("code", "functions", "mde.R"))                # mde()
source(here("code", "functions", "make_country_id.R"))    # make_country_id()
source(here("code", "functions", "make_wide_table.R"))    # make_wide_table()
source(here("code", "functions", "sup_t_crit.R"))         # sup_t_crit()

set.seed(20240601)  # global seed — local set.seed(1242) calls follow each estimator

# -----------------------------------------------------------------------
# One fit, one SE: single bootstrap-replication constant, and the
# multiplier bootstrap for every heterogeneity cell.
#
# The inference method does not depend on the number of treated units: thin
# cells are where analytical influence-function SEs are most
# anti-conservative, and a 40-unit cut-off is a convention, not an estimator
# requirement. All reported subgroup ATTs and SEs come from a
# multiplier-bootstrap fit clustered by recipient country (did clusters the
# multiplier bootstrap on `idname` by construction).
# Analytical twins are retained ONLY where an influence function is required:
# the pre-trend Wald test and the correlated-sample contrasts of §25a.
# -----------------------------------------------------------------------
BITERS <- 999L

# Saved fits written by 03_main_results.R / 04_robustness.R.
FITS_DIR <- here("output", "fits")

# A failed run must not leave the previous run's exhibits in place: every
# exhibit this stage owns is deleted before anything is re-estimated.
unlink(list.files(here("output", c("tables", "figures"), "heterogeneity"),
                  pattern = "\\.(tex|png|pdf)$", recursive = TRUE, full.names = TRUE))

# canonical pre-trend wording — keep byte-identical across scripts
PRETREND_NOTE_AGG <- function(min_e, max_e, k) sprintf(
  "Pre-trend $\\chi^2$ / $p$: joint Wald test on the aggregated pre-treatment event-time coefficients ($%d \\leq e \\leq %d$; %d restrictions), using the influence-function covariance of the dynamic aggregation from the analytical (non-bootstrap) fit; a generalized inverse is used if the block is singular",
  min_e, max_e, k)
PRETREND_NOTE_DID <- "\\texttt{did} pre-test $p$: \\texttt{did}'s built-in Wald test over all pre-period $ATT(g,t)$ cells against each cohort's $g-1$ base year"
PRETREND_NOTE <- function(min_e, max_e, k) paste0(PRETREND_NOTE_AGG(min_e, max_e, k), ". ", PRETREND_NOTE_DID)

#' Render the two pre-trend columns honestly, from the fit that produced them.
#'
#' The tables report two pre-trend statistics that often disagree. They are
#' different tests of different objects, and the note must say so WITH the
#' numbers of the specification it sits under -- nothing here is hardcoded.
#'
#' The lead window is taken from the window actually used, "rejects" is
#' stated only when did returns a statistic, and no mechanism is asserted.
#'
#' @param pre_egt integer vector of pre-treatment event times that entered the
#'   aggregated Wald test (from compute_pretrend_test()$leads)
#' @param df_did restrictions in did's built-in pre-test, i.e. the number of
#'   estimable pre-treatment group-time cells
#' @param n_clusters recipients the influence function is clustered on
#' @param wpval_did did's pre-test p-value(s); NA where it returned none
#' @param pval_wald aggregated Wald p-value(s), same length as wpval_did
#' @param labels optional column labels, same length, for per-column disclosure
#' @param wpval_reason machine-derived reason did returned no statistic
#' @param ginv_used TRUE if the aggregated test used a generalized inverse
#' @param df degrees of freedom of the aggregated test(s) (compute_pretrend_test()$df),
#'   one per column; below full rank they are fewer than the leads
#' @return a character string for the table note
wpval_reconciliation <- function(pre_egt, df_did, n_clusters,
                                 wpval_did = NA_real_, pval_wald = NA_real_,
                                 labels = NULL, wpval_reason = NA_character_,
                                 ginv_used = FALSE, df = length(pre_egt)) {
  # Compact pre-trend disclosure for table notes (byte-identical
  # in 03/04/05). The averaging mechanism and the interpretation of each
  # result are stated once in the main text; the note keeps the restriction
  # counts and the per-column verdicts, which the text does not repeat table
  # by table.

  # Column labels are interpolated into LaTeX prose, so they must be escaped
  # here: "Adaptation share (% of global)" would otherwise comment out the rest
  # of the note.
  esc_lab <- function(x) {
    x <- gsub("\\\\", "", x)
    x <- gsub("([%#&_])", "\\\\\\1", x)
    x
  }
  if (!is.null(labels)) labels <- esc_lab(labels)

  # --- lead window, rendered from the vector actually used -----------------
  n_lead <- length(pre_egt)
  lead_list <- if (n_lead == 0L) "none" else
    paste0("$e = ", paste(sprintf("%d", as.integer(sort(pre_egt))),
                          collapse = ", "), "$")

  head_txt <- if (n_lead == 0L) {
    "No aggregated pre-trend test: no pre-treatment event time has a usable SE. "
  } else {
    # canonical wording: the lead window itself is disclosed via
    # min/max of pre_egt, not re-derived from the (possibly non-contiguous)
    # lead_list string built above.
    # Restrictions = the test's df (the rank actually inverted), not the lead
    # count; if the df differ across columns the lead count is quoted and the
    # singular columns are flagged just below.
    df_u <- unique(as.integer(df[!is.na(df)]))
    paste0(PRETREND_NOTE_AGG(min(as.integer(pre_egt)), max(as.integer(pre_egt)),
                             if (length(df_u) == 1L) df_u else n_lead),
           if (isTRUE(any(ginv_used)))
             paste0("; the block was singular",
                    if (!is.null(labels) && length(labels) == length(ginv_used))
                      paste0(" for ", paste(labels[ginv_used], collapse = ", ")) else "",
                    " and a generalized inverse was used")
           else "",
           ". ")
  }

  wp <- suppressWarnings(as.numeric(wpval_did))
  pw <- suppressWarnings(as.numeric(pval_wald))
  if (length(wp) == 0L) wp <- NA_real_
  if (length(pw) == 0L) pw <- NA_real_

  # Counts may differ across columns (subgroup tables); render one number when
  # they agree and the range when they do not, never a single column's value
  # presented as if it were the table's.
  fmt_count <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0L) return(NA_character_)
    if (length(unique(x)) == 1L) as.character(x[[1L]]) else
      paste0(min(x), "--", max(x))
  }
  df_did_txt <- fmt_count(df_did)
  clust_txt  <- fmt_count(n_clusters)
  # canonical wording: PRETREND_NOTE_DID states what the test IS;
  # the parenthetical keeps the per-column restriction/recipient counts, which
  # are data, not wording, and must not be lost.
  did_head <- paste0(
    PRETREND_NOTE_DID, " (",
    if (is.na(df_did_txt)) "one restriction per pre-treatment cell" else
      paste0(df_did_txt, " restrictions"),
    if (is.na(clust_txt)) "" else paste0(", ", clust_txt, " recipients"),
    ")")

  if (all(is.na(wp))) {
    rr <- wpval_reason[!is.na(wpval_reason)]
    sing_re <- "^not computed: pre-treatment covariance singular \\((\\d+) cells, numerical rank (\\d+)\\)$"
    reason_txt <- if (length(rr) == 0L) "no statistic returned"
      else if (length(unique(rr)) == 1L) sub("^not computed: ", "", unique(rr)[[1L]])
      else if (!is.null(labels) && length(wpval_reason) == length(labels) &&
               all(grepl(sing_re, wpval_reason)))
        paste0("pre-treatment covariance singular (",
               sub(sing_re, "\\1", wpval_reason[[1L]]), " cells; rank ",
               paste(sprintf("%s %s", labels, sub(sing_re, "\\2", wpval_reason)),
                     collapse = ", "), ")")
      else if (!is.null(labels) && length(wpval_reason) == length(labels))
        paste(sprintf("%s, %s", labels,
                      ifelse(is.na(wpval_reason), "no reason recorded",
                             sub("^not computed: ", "", wpval_reason))),
              collapse = "; ")
      else "reason varies by column (see text)"
    did_txt <- paste0(did_head, " not computed: ", reason_txt,
                      ". ")
  } else {
    both <- !is.na(wp) & !is.na(pw)
    rej_wald <- both & pw < 0.05
    rej_did  <- both & wp < 0.05
    per_col <- ""
    if (length(pw) == 1L && both[[1L]]) {
      per_col <- sprintf(": aggregate $p = %.3f$, \\texttt{did} $p = %.3f$", pw[[1L]], wp[[1L]])
    } else if (!is.null(labels) && length(labels) == length(pw) && length(pw) > 1L) {
      nm_rej <- labels[which(rej_wald)]
      wald_clause <- if (length(nm_rej) == 0L) "aggregate rejects (5\\%) in no column"
        else if (all(rej_wald[both])) "aggregate rejects (5\\%) in every column"
        else paste0("aggregate rejects (5\\%) for ", paste(nm_rej, collapse = ", "), " only")
      did_clause <- if (all(rej_did[both])) "\\texttt{did} rejects in every column"
        else if (!any(rej_did[both])) "\\texttt{did} rejects in none"
        else paste0("\\texttt{did} rejects for ", paste(labels[which(rej_did)], collapse = ", "))
      per_col <- paste0(": ", wald_clause, "; ", did_clause)
    }
    did_txt <- paste0(did_head, per_col, ". ")
  }
  paste0(head_txt, did_txt)
}

# -----------------------------------------------------------------------
# Null-coalescing operator
# -----------------------------------------------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------------------------------------------------
# Escape LaTeX-special characters in TABLE COLUMN HEADERS only.
# (Duplicate of helper in 03_main_results.R — each stage script is self-contained)
# -----------------------------------------------------------------------
esc_header <- function(x) {
  x <- gsub("%", "\\\\%", x)
  x <- gsub("_", "\\\\_", x)
  x <- gsub("#", "\\\\#", x)
  x
}

# ==============================================================================
# Helper: write_tex_float()
# (Identical copy to 04_robustness.R — each stage script is self-contained)
# Wraps a bare tabular block in the project's standard complete float.
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
# SECTION 1. Load and rebuild the DiD panel
# (duplicated from 03_main_results.R §1 — keep in sync)
# ==============================================================================

message("\n=== 05_heterogeneity.R: loading panel ===\n")

aggregated <- fread(here("data", "processed", "simple_panel_wgi.csv"))
aggregated  <- as.data.frame(aggregated)

did_panel <- aggregated %>%
  mutate(country_id = make_country_id(recipient_name))

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
            by = "recipient_name") %>%
  mutate(log_population = log(population))

if (!"share_adapt" %in% names(did_panel)) {
  share_lookup <- aggregated %>%
    select(recipient_name, year, share_adapt) %>%
    distinct()
  did_panel <- left_join(did_panel, share_lookup, by = c("recipient_name", "year"))
}

did_panel_full <- did_panel

# Cohort sizes and thin-cohort threshold
cohort_sizes <- did_panel_full %>%
  filter(cohort_year > 0) %>%
  distinct(recipient_name, cohort_year) %>%
  count(cohort_year, name = "n_treated") %>%
  arrange(cohort_year)

thin_threshold <- 5L
thin_cohorts   <- cohort_sizes$cohort_year[cohort_sizes$n_treated < thin_threshold]

message(sprintf("thin_cohorts: %s", paste(thin_cohorts, collapse = ", ")))

# ==============================================================================
# SECTION 2. Outcomes list (identical to 03_main_results.R §4)
# (duplicated from 03/04 §1 — keep in sync)
# ==============================================================================

outcomes <- list(
  list(var   = "log_commits",
       label = "log(Adaptation commitments)",
       file  = "adapt_commits",
       color = "#2E86C1"),
  list(var   = "share_adapt",
       label = "Adaptation share (% of global)",
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
# SECTION 3. Shared helpers
# (make_wide_table() is sourced from code/functions/make_wide_table.R and
#  compute_pretrend_test() from code/functions/pretrend_test.R)
# ==============================================================================

# --- Difference between two subgroup ATTs ---
#' Wald test that two subgroup ATTs are equal.
#'
#' The two subgroups are disjoint sets of countries and each subgroup ATT is
#' estimated on its own sample with an influence function clustered by country,
#' so the two estimators are asymptotically independent and
#'   Var(theta_a - theta_b) = Var(theta_a) + Var(theta_b).
#'
#' @param stats_a list with att_num, se_num (low-capacity subgroup)
#' @param stats_b list with att_num, se_num (high-capacity subgroup)
#' @return named list with diff, se, z, pval
att_difference_test <- function(stats_a, stats_b) {
  stopifnot(is.list(stats_a), is.list(stats_b))
  att_a <- stats_a$att_num
  att_b <- stats_b$att_num
  se_a  <- stats_a$se_num
  se_b  <- stats_b$se_num
  if (any(is.na(c(att_a, att_b, se_a, se_b))))
    return(list(diff = NA_real_, se = NA_real_, z = NA_real_, pval = NA_real_))

  diff    <- att_a - att_b
  se_diff <- sqrt(se_a^2 + se_b^2)          # independence across disjoint samples
  z       <- if (se_diff > 1e-12) diff / se_diff else NA_real_
  pval    <- if (is.na(z)) NA_real_ else 2 * pnorm(-abs(z))
  list(diff = diff, se = se_diff, z = z, pval = pval)
}

# --- Joint Wald test that k subgroup ATTs are all equal ---
#' @param theta numeric vector of subgroup ATTs
#' @param se_vec numeric vector of subgroup SEs (same order)
#' @param ref_index integer index of the reference cell for the contrasts
#' @return named list with stat, df, pval
att_joint_wald <- function(theta, se_vec, ref_index) {
  stopifnot(length(theta) == length(se_vec), length(theta) >= 2L)
  k <- length(theta)
  if (any(is.na(theta)) || any(is.na(se_vec)))
    return(list(stat = NA_real_, df = NA_integer_, pval = NA_real_))

  rows <- setdiff(seq_len(k), as.integer(ref_index))
  R    <- matrix(0, nrow = length(rows), ncol = k)   # pre-allocated contrast matrix
  for (i in seq_along(rows)) {
    R[i, rows[i]]            <- 1
    R[i, as.integer(ref_index)] <- -1
  }
  # Diagonal covariance: disjoint subsamples => zero cross-subgroup covariance.
  V     <- diag(se_vec^2, nrow = k)
  r_th  <- R %*% theta
  r_v_r <- R %*% V %*% t(R)
  W <- as.numeric(t(r_th) %*% solve(r_v_r) %*% r_th)  # a singular contrast covariance stops
  list(stat = W,
       df   = nrow(R),
       pval = if (is.na(W)) NA_real_ else pchisq(W, df = nrow(R), lower.tail = FALSE))
}

# --- Contrast between two ATTs estimated on the SAME units -------------
#' Difference between two ATTs with correlated-sample inference.
#'
#' Ported verbatim (logic and scaling convention) from
#' 07_principal_and_share.R::compute_att_difference. Do not "simplify"
#' (e.g. dividing by n instead of sqrt(n)): the scaling is sqrt(mean((IFa - IFb)^2) / n),
#' i.e. the UNCENTRED mean-square convention that did itself uses (verified
#' below against overall.se), not sqrt(var(IFa - IFb)) / n and not
#' sqrt(var(IFa - IFb)).
#'
#' Applies when the two ATTs are estimated on the same recipient-year panel
#' with different outcome variables (donor-type decompositions, adaptation vs
#' mitigation): the unit-level influence functions are then aligned
#' unit-for-unit and Var(theta_a - theta_b) = E[(IFa - IFb)^2] / n, which
#' subtracts the (positive) covariance the independence formula ignores.
#'
#' @param att_a,att_b point estimates
#' @param IFa,IFb unit-level influence-function vectors from
#'   aggte(type = "simple")$inf.function$simple.att of the ANALYTICAL fits
#' @param se_a,se_b reported standard errors (used only for the fallback)
#' @param ids_a,ids_b unit identifiers the influence-function rows refer to
#' @return list(diff, se, z, pval, method, n_units)
compute_att_difference <- function(att_a, att_b, IFa, IFb, se_a, se_b,
                                   ids_a = NULL, ids_b = NULL) {
  diff <- att_a - att_b
  have_if <- !is.null(IFa) && !is.null(IFb)
  same_units <- have_if && length(IFa) == length(IFb) &&
    (is.null(ids_a) || is.null(ids_b) || identical(ids_a, ids_b))

  # If the two fits used different unit sets, restrict to their intersection
  # (reported explicitly in the method string): the resulting SE treats the
  # intersection as the estimation sample, which is exact when the unit sets
  # coincide and an approximation otherwise.
  if (have_if && !same_units && !is.null(ids_a) && !is.null(ids_b) &&
      length(IFa) == length(ids_a) && length(IFb) == length(ids_b)) {
    common <- intersect(ids_a, ids_b)
    if (length(common) >= 2L) {
      IFa <- IFa[match(common, ids_a)]
      IFb <- IFb[match(common, ids_b)]
      n_units <- length(common)
      se_diff <- sqrt(mean((IFa - IFb)^2) / n_units)
      z    <- if (se_diff > 1e-12) diff / se_diff else NA_real_
      return(list(diff = diff, se = se_diff, z = z,
                  pval = if (is.na(z)) NA_real_ else 2 * pnorm(-abs(z)),
                  method = paste0("influence functions aligned on the ",
                                  "intersection of the two estimation samples ",
                                  "(n = ", n_units, " of ", length(ids_a), " and ",
                                  length(ids_b), " units)"),
                  n_units = n_units))
    }
  }

  if (same_units) {
    n_units <- length(IFa)
    se_diff <- sqrt(mean((IFa - IFb)^2) / n_units)
    method  <- paste0("aligned unit-level influence functions (n = ", n_units, ")")
  } else {
    n_units <- NA_integer_
    se_diff <- sqrt(se_a^2 + se_b^2)
    method  <- paste0("independence approximation (influence functions did not ",
                      "align unit-for-unit); the two outcomes are measured on ",
                      "the same recipient-years and are positively correlated, ",
                      "so this SE is conservative in direction but its magnitude ",
                      "relative to the true SE is not separately verified")
  }
  z    <- if (se_diff > 1e-12) diff / se_diff else NA_real_
  pval <- if (is.na(z)) NA_real_ else 2 * pnorm(-abs(z))
  list(diff = diff, se = se_diff, z = z, pval = pval, method = method,
       n_units = n_units)
}

#' Joint Wald test on a set of correlated ATTs using aligned influence functions.
#'
#' @param theta numeric vector of ATTs
#' @param IFmat n x k matrix of aligned unit-level influence functions
#' @param ref_index index of the reference cell
#' @return list(stat, df, pval)
att_joint_wald_if <- function(theta, IFmat, ref_index) {
  stopifnot(is.matrix(IFmat), ncol(IFmat) == length(theta), length(theta) >= 2L)
  k <- length(theta)
  n <- nrow(IFmat)
  rows <- setdiff(seq_len(k), as.integer(ref_index))
  R <- matrix(0, nrow = length(rows), ncol = k)   # pre-allocated contrast matrix
  for (i in seq_along(rows)) {
    R[i, rows[i]]               <- 1
    R[i, as.integer(ref_index)] <- -1
  }
  # Uncentred mean-square convention, matching compute_att_difference().
  V <- crossprod(IFmat) / n^2
  r_th  <- R %*% theta
  r_v_r <- R %*% V %*% t(R)
  W <- as.numeric(t(r_th) %*% solve(r_v_r) %*% r_th)  # a singular contrast covariance stops
  list(stat = W, df = nrow(R),
       pval = if (is.na(W)) NA_real_ else pchisq(W, df = nrow(R), lower.tail = FALSE))
}

# ==============================================================================
# SECTION 4. Helper: make_het_wide_table()
# Shared across §22-24.
# Estimator and inference labels are derived from the em_g/bstrap actually
# used; every table note defines the significance stars.
# Seed rule (reproduces the published SEs): set.seed(1242) before het loop.
#
# Returns (invisibly) the per-group stat list. Each element carries both the
# formatted strings used to build the wide table and the unrounded att_num /
# se_num / n_treated used by §25 for the subgroup difference tests, so the
# difference tests re-use exactly the fits that feed the published tables.
# ==============================================================================

make_het_wide_table <- function(groups, outcome_var, raw_var,
                                dir_tabs, tex_label, caption_txt,
                                thin_cohorts_vec = thin_cohorts) {

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # One slot per subgroup; a subgroup that cannot be estimated stops the run.
  grp_labels <- vapply(groups, `[[`, character(1L), "label")
  col_stats  <- setNames(vector("list", length(groups)), grp_labels)
  em_used    <- setNames(vector("list", length(groups)), grp_labels)
  bs_used    <- setNames(vector("list", length(groups)), grp_labels)

  # Seed rule (reproduces the published SEs): seed before het loop
  set.seed(1242)

  for (grp in groups) {
    lbl <- grp$label
    message(sprintf("  [het_wide] %s: %s", lbl, outcome_var))

    panel_g <- grp$panel %>% filter(!(cohort_year %in% thin_cohorts_vec))

    if (n_distinct(panel_g$cohort_year[panel_g$cohort_year > 0]) < 2)
      stop("Subgroup ", lbl, " has fewer than two treated cohorts: no column.")

    n_treated_g <- n_distinct(panel_g$country_id[panel_g$cohort_year > 0])
    # est_method still follows the sample-size rule (doubly robust needs enough
    # treated units for the propensity-score leg to be stable, and att_gt(dr)
    # is known to segfault inside fastglm on thin unbalanced subpanels). The
    # BOOTSTRAP switch, by contrast, is gone: bootstrap is always on.
    em_g <- if (n_treated_g >= 40L) "dr" else "reg"
    message(sprintf("    N treated = %d  ->  est_method = %s  |  bstrap = TRUE (%d reps)",
                    n_treated_g, em_g, BITERS))
    em_used[[lbl]] <- em_g
    bs_used[[lbl]] <- TRUE

    # Helper so the analytical twin and the bootstrap fit cannot drift apart.
    run_gt_g <- function(bstrap_flag) {
      att_gt(
        yname         = outcome_var,
        tname         = "year",
        idname        = "country_id",
        gname         = "cohort_year",
        xformla       = ~ ge_est + log_population,
        data          = panel_g,
        est_method    = em_g,
        bstrap        = bstrap_flag,
        biters        = BITERS,
        cband         = FALSE,
        control_group = "notyettreated",
        anticipation  = 0,
        base_period   = "universal",
        panel         = TRUE,
        allow_unbalanced_panel = TRUE
      )
    }

    # Analytical twin first (consumes no random numbers): supplies the
    # influence function for the pre-trend Wald test and for the §25a contrasts.
    gt_g_analytical <- tryCatch(run_gt_g(FALSE),
      error = function(e) { message("  att_gt (analytical) failed: ",
                                     conditionMessage(e)); NULL })

    # Reported ATT and SE: multiplier bootstrap, clustered by recipient.
    set.seed(1242)   # Seed rule (reproduces the published SEs): seed immediately before the estimator
    gt_g <- tryCatch(run_gt_g(TRUE),
      error = function(e) { message("  att_gt (bootstrap) failed: ",
                                     conditionMessage(e)); NULL })
    if (is.null(gt_g) || is.null(gt_g_analytical)) stop("att_gt() failed for subgroup ", lbl)

    # Pre-trend test uses the analytical fit (as in 03/04): its aggregation
    # reports the influence-function SEs the test's covariance reproduces.
    agg_d <- aggte(gt_g_analytical, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf)
    agg_s_analytical <- aggte(gt_g_analytical, type = "simple", na.rm = TRUE)
    agg_s <- aggte(gt_g, type = "simple", na.rm = TRUE)

    att  <- if (!is.null(agg_s)) agg_s$overall.att else NA_real_
    se   <- if (!is.null(agg_s)) agg_s$overall.se  else NA_real_
    t_v  <- if (!is.na(att) && !is.na(se) && se > 0) att / se else NA_real_
    stars <- if (is.na(t_v)) "" else
      if (abs(t_v) > 2.576) "***" else if (abs(t_v) > 1.960) "**" else
      if (abs(t_v) > 1.645) "*"   else ""

    # gt_g_analytical (not gt_g) is passed so that did's own pre-test and the
    # influence-function covariance are never bootstrap-contaminated -- the
    # convention used at every other call site in this project.
    pt <- compute_pretrend_test(agg_d, gt_g_analytical)

    rows_g    <- panel_g[!is.na(panel_g[[outcome_var]]), ]
    n_obs     <- nrow(rows_g)
    n_country <- length(unique(rows_g$country_id))

    # The naive back-transform exp(ATT)-1 x pre-treatment mean is not
    # additive across outcomes/subgroups, so it is reported only when the
    # underlying ATT is significant at 5% (|t| > 1.960); otherwise "---".
    is_sig5 <- !is.na(t_v) && abs(t_v) > 1.960

    if (!is.na(raw_var) && raw_var %in% names(panel_g)) {
      pre_rows <- panel_g %>%
        filter(cohort_year > 0, year < cohort_year, !is.na(.data[[raw_var]]))
      mean_pre <- mean(pre_rows[[raw_var]], na.rm = TRUE)
      mean_pre_fmt  <- sprintf("%.1f", mean_pre)
      implied_fmt <- if (!is.na(att) && is_sig5) {
        sprintf("%.1f", (exp(att) - 1) * mean_pre)
      } else {
        "---"
      }
    } else {
      mean_pre_fmt <- "---"
      implied_fmt  <- "---"
    }

    col_stats[[lbl]] <- list(
      att_fmt      = paste0(sprintf("%.4f", att), stars),
      se_fmt       = sprintf("(%.4f)", se),
      t_fmt        = sprintf("%.3f", t_v),
      mean_pre_fmt = mean_pre_fmt,
      implied_fmt  = implied_fmt,
      n_obs        = format(n_obs,  big.mark = ","),
      n_country    = as.character(n_country),
      pt_stat      = sprintf("%.3f", pt$stat),
      pt_pval      = sprintf("%.3f", pt$pval),
      pt_df        = pt$df,
      pt_df_did    = pt$df_did,
      pt_leads     = pt$leads,
      pt_ginv      = pt$ginv_used,
      pt_reason    = pt$wpval_reason,
      pt_pval_num  = pt$pval,
      pt_wpval_num = pt$Wpval_did,
      pt_wpval_did = if (is.na(pt$Wpval_did)) "---" else sprintf("%.3f", pt$Wpval_did),
      # Unrounded values for the §25 difference tests (the wide table selects
      # its cells by name).
      att_num      = att,
      se_num       = se,
      n_treated    = n_treated_g,
      # Analytical influence function of the same subgroup fit, for the
      # correlated-sample contrasts; se_analytic lets the note report how far
      # the bootstrap SE sits from the analytical one.
      inf_func     = agg_s_analytical$inf.function$simple.att,
      se_analytic  = agg_s_analytical$overall.se,
      att_analytic = agg_s_analytical$overall.att,
      ids          = sort(unique(panel_g$country_id)),
      est_method   = em_g,
      # Reported as its own table row: the split cells are small and the text
      # cites the per-cell treated counts, so they must be verifiable.
      n_treated_fmt = as.character(n_treated_g)
    )
    message(sprintf("    ATT = %s  pre-trend chi2(%d)=%.3f p=%.3f  did Wpval=%s",
                    col_stats[[lbl]]$att_fmt, pt$df, pt$stat, pt$pval,
                    col_stats[[lbl]]$pt_wpval_did))
    if (is.na(pt$Wpval_did))
      message("    did pre-test unavailable for this subgroup (singular pooled ",
              "pre-treatment covariance); the aggregated 4-lead Wald test is reported.")
    message(sprintf("    SE: bootstrap = %.4f | analytical = %.4f (ratio %.3f)",
                    se, col_stats[[lbl]]$se_analytic,
                    se / col_stats[[lbl]]$se_analytic))
  }

  if (all(vapply(col_stats, is.null, logical(1L))))
    stop("No subgroup results: ", tex_label, " not written.")

  # Estimator label from the em_g actually used per group.
  unique_em <- unique(unlist(em_used))
  estimator_label <- if (length(unique_em) == 1L) {
    if (unique_em == "dr") "CS (2021) DR" else "CS (2021) regression adjustment"
  } else {
    "CS (2021) DR or regression adjustment (by subgroup N)"
  }
  # Derived from what was actually used, not hard-coded (every cell is
  # bootstrapped, so this is TRUE for every cell; the check keeps the label
  # honest if a cell ever is not).
  bootstrap_label <- if (all(unlist(bs_used))) paste0("Yes (multiplier, ", BITERS, " reps)") else
    "Mixed (see text)"

  row_labels <- c(
    "ATT", "SE", "$t$-statistic",
    "Mean (pre-treat., USD M)",
    "Implied effect (USD M, exp(ATT)$-$1 $\\times$ pre-treat.\\ mean)",
    "Observations", "Countries", "N treated countries",
    "Pre-trend $\\chi^2$", "Pre-trend $p$", "\\texttt{did} pre-test $p$",
    "\\midrule Estimator", "Control group", "Bootstrap SE"
  )

  tab_wide <- data.frame(` ` = row_labels, check.names = FALSE, stringsAsFactors = FALSE)

  for (lbl in names(col_stats)) {
    s <- col_stats[[lbl]]
    stopifnot(!is.null(s))
    col <- c(
      s$att_fmt, s$se_fmt, s$t_fmt,
      s$mean_pre_fmt, s$implied_fmt,
      s$n_obs, s$n_country, s$n_treated_fmt,
      s$pt_stat, s$pt_pval, s$pt_wpval_did,
      estimator_label, "Not-yet-treated", bootstrap_label
    )
    tab_wide[[esc_header(lbl)]] <- col
  }

  xtab <- xtable(tab_wide, label = tex_label)
  n_grp <- length(col_stats)
  align(xtab) <- paste0("ll", paste(rep("c", n_grp), collapse = ""))

  raw_lines <- capture.output(
    print(xtab, include.rownames = FALSE, booktabs = TRUE,
          sanitize.text.function = identity, size = "\\small",
          floating = FALSE)
  )

  # Per-column inputs for the pre-trend note:
  # the lead window actually used, both p-values, the cell and cluster counts
  # and did's own machine-derived reason -- per subgroup, not borrowed from the
  # first column.
  het_cells <- Filter(Negate(is.null), col_stats)
  het_labels   <- names(het_cells)
  het_leads    <- if (length(het_cells) == 0L) integer(0) else het_cells[[1L]]$pt_leads
  het_df_did   <- vapply(het_cells, function(s)
    if (is.null(s$pt_df_did)) NA_integer_ else as.integer(s$pt_df_did), integer(1L))
  het_nclust   <- vapply(het_cells, function(s) as.integer(s$n_country), integer(1L))
  het_wpval    <- vapply(het_cells, function(s)
    if (is.null(s$pt_wpval_num)) NA_real_ else s$pt_wpval_num, numeric(1L))
  het_pwald    <- vapply(het_cells, function(s)
    if (is.null(s$pt_pval_num)) NA_real_ else s$pt_pval_num, numeric(1L))
  het_reason   <- vapply(het_cells, function(s)
    if (is.null(s$pt_reason)) NA_character_ else s$pt_reason, character(1L))
  het_ginv     <- vapply(het_cells, function(s) isTRUE(s$pt_ginv), logical(1L))
  het_df       <- vapply(het_cells, function(s) as.integer(s$pt_df), integer(1L))

  # Star definitions end the note.
  notes_txt <- paste0(
    caption_txt,
    " ATT/SE: multiplier-bootstrap (", BITERS,
    " reps, seed 1242); pre-trend: separate analytical fit. ",
    wpval_reconciliation(pre_egt = het_leads, df_did = het_df_did,
                         n_clusters = het_nclust, wpval_did = het_wpval,
                         pval_wald = het_pwald, labels = het_labels,
                         wpval_reason = het_reason, ginv_used = het_ginv,
                         df = het_df),
    "Implied effects: back-transform on pre-treatment mean; suppressed if insig.\\ at 5\\%. ",
    "* $p<0.10$, ** $p<0.05$, *** $p<0.01$"
  )
  source_txt <- "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"

  dir.create(dir_tabs, recursive = TRUE, showWarnings = FALSE)
  # Strip "tab:" prefix from tex_label for the filename — the LaTeX
  # \label{} keeps the full "tab:het_..." prefix, but filenames must not
  # contain colons (breaks on Windows and confuses latexmk includes).
  out_path <- file.path(dir_tabs, paste0(sub("^tab:", "", tex_label), ".tex"))

  # Caption title is the brief part of caption_txt (before the method detail)
  # For het tables we split at the first period; fallback to full text.
  cap_short <- strsplit(caption_txt, "\\.")[[1]][1]
  if (is.na(cap_short) || nchar(cap_short) == 0) cap_short <- caption_txt

  write_tex_float(out_path, cap_short, tex_label, raw_lines, notes_txt, source_txt)

  invisible(col_stats)
}

# ==============================================================================
# SECTION 5. §21 Donor-type heterogeneity
# Seed rule (reproduces the published SEs): set.seed(1242) before donor loop.
# Note: donor_type uses DR + multiplier bootstrap in BOTH the event-study
#   figure and make_wide_table; make_wide_table additionally fits an analytical
#   twin whose influence function feeds the pre-trend test and the §25a contrasts.
# ==============================================================================

message("\n=== Section 21: Heterogeneity — Donor type ===\n")

dir_tabs_h1 <- here("output", "tables",  "heterogeneity", "donor_type")
dir.create(here("output", "figures", "heterogeneity"), recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tabs_h1, recursive = TRUE, showWarnings = FALSE)

outcomes_donor <- list(
  list(var = "log_commits_dac",   label = "DAC bilateral",  color = "#2166ac"),
  list(var = "log_commits_multi", label = "Multilateral",   color = "#d6604d"),
  list(var = "log_commits_other", label = "Other donors",   color = "#4dac26")
)

did_panel_h1 <- did_panel_full %>% filter(!(cohort_year %in% thin_cohorts))

# --- Event-study curves (panel (a) of the 2x2 figure, §24b) ---
h1_dyn <- setNames(vector("list", length(outcomes_donor)),
                   vapply(outcomes_donor, `[[`, character(1L), "var"))
# Seed set immediately before the estimator (reproduces the published SE)
set.seed(1242)
for (oc in outcomes_donor) {
  if (!oc$var %in% names(did_panel_h1)) stop("Donor-type outcome missing from the panel: ", oc$var)
  gt_h1 <- tryCatch(
    att_gt(yname = oc$var, tname = "year", idname = "country_id", gname = "cohort_year",
           xformla = ~ ge_est + log_population, data = did_panel_h1,
           est_method = "dr", bstrap = TRUE, biters = BITERS, cband = FALSE,
           control_group = "nevertreated", anticipation = 0,
           base_period = "universal", panel = TRUE, allow_unbalanced_panel = TRUE),
    error = function(e) NULL)
  if (is.null(gt_h1)) stop("Donor-type event study: att_gt() failed for ", oc$var)
  agg_d <- tryCatch(aggte(gt_h1, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf),
                    error = function(e) NULL)
  if (is.null(agg_d)) stop("Donor-type event study: aggte() failed for ", oc$var)
  # Simultaneous (sup-t) 95% band, uniform over this curve's event times (sup_t_crit()).
  cv <- sup_t_crit(agg_d$inf.function$dynamic.inf.func.e, agg_d$se.egt, biters = BITERS)
  message(sprintf("  Donor-type ES sup-t crit (%s): %.4f", oc$label, cv))
  h1_dyn[[oc$var]] <- data.frame(
    outcome    = oc$label, color = oc$color,
    event_time = agg_d$egt, ATT = agg_d$att.egt, SE = agg_d$se.egt,
    Lower      = agg_d$att.egt - cv * agg_d$se.egt,
    Upper      = agg_d$att.egt + cv * agg_d$se.egt,
    stringsAsFactors = FALSE)
}
if (!all(vapply(h1_dyn, is.null, logical(1L)))) {
  h1_df <- bind_rows(h1_dyn) %>%
    filter(!is.na(SE), SE > 1e-10) %>%
    mutate(outcome = factor(outcome, levels = vapply(outcomes_donor, `[[`, character(1L), "label")))
  palette_h1 <- setNames(vapply(outcomes_donor, `[[`, character(1L), "color"),
                          vapply(outcomes_donor, `[[`, character(1L), "label"))
  p_h1 <- ggplot(h1_df, aes(x = event_time, y = ATT,
                              colour = outcome, shape = outcome, group = outcome)) +
    geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = -0.5, colour = "grey60", linetype = "dotted", linewidth = 0.4) +
    geom_linerange(aes(ymin = Lower, ymax = Upper),
                   position = position_dodge(0.4), linewidth = 0.6, alpha = 0.8) +
    geom_point(size = 2.5, position = position_dodge(0.4)) +
    scale_colour_manual(values = palette_h1) +
    scale_shape_manual(values = c(16, 17, 15)) +
    # No title, subtitle, or caption — those go in LaTeX \caption{}
    labs(title = NULL, subtitle = NULL, caption = NULL,
         x = "Event time (years relative to NAP adoption)",
         y = "ATT (log commitments)", colour = NULL, shape = NULL) +
    theme_minimal() +
    theme(text             = element_text(family = "serif", size = 11),
          legend.position  = "bottom", panel.grid.minor = element_blank())
}

# --- Wide table: uses make_wide_table (not het version) ---
# The return value is captured so that §25a can contrast the three donor
# ATTs using the influence functions of these very fits (no re-estimation).
donor_stats <- make_wide_table(
  did_panel_in  = did_panel_h1,
  retain_thin   = FALSE,
  outcomes      = outcomes_donor,
  dir_tabs      = dir_tabs_h1,
  tex_label     = "tab:het_donor_wide",
  caption_spec  = "heterogeneity by donor type; cohorts $\\geq 5$ units, DR estimator"
)

message("\n=== Section 21 complete ===\n")

# ==============================================================================
# SECTION 6. §22 LDC vs. non-LDC heterogeneity
# Seed rule (reproduces the published SEs): set.seed(1242) before LDC loop.
# ==============================================================================

message("\n=== Section 22: Heterogeneity — LDC vs. non-LDC ===\n")

dir_tabs_h2 <- here("output", "tables",  "heterogeneity", "ldc")
dir.create(dir_tabs_h2, recursive = TRUE, showWarnings = FALSE)

# PRE-TREATMENT LDC VINTAGE.
# The split is a proxy for baseline capacity, so it must be measured before any
# recipient in the sample adopts a NAP (first estimation cohort: 2021). The UN
# DESA/CDP list as of 1 January 2024 would encode a decade of graduations that
# are themselves outcomes of the development process the split is meant to
# condition on. The classification used is therefore the UN list
# as of 2013 (49 economies), the first full year of the panel's pre-treatment
# window; the 2024 list is retained below only to document what changes.
# Sources:
#   https://www.un.org/development/desa/dpad/least-developed-country-category.html
#   UN CDP: South Sudan added December 2012 (hence an LDC in 2013).
# Graduations AFTER 2013 and therefore still LDC in the 2013 vintage:
#   Samoa (2014), Equatorial Guinea (2017), Vanuatu (2020), Bhutan (2023),
#   Sao Tome and Principe (December 2024).
# Graduations BEFORE 2013 and therefore non-LDC in both vintages:
#   Botswana (1994), Cape Verde (2007), Maldives (2011).
ldc_iso3_2013 <- c(
  "AFG", "AGO", "BGD", "BEN", "BTN", "BFA", "BDI", "KHM", "CAF", "TCD",
  "COM", "COD", "DJI", "GNQ", "ERI", "ETH", "GMB", "GIN", "GNB", "HTI",
  "KIR", "LAO", "LSO", "LBR", "MDG", "MWI", "MLI", "MRT", "MOZ", "MMR",
  "NPL", "NER", "RWA", "WSM", "STP", "SEN", "SLE", "SLB", "SOM", "SSD",
  "SDN", "TLS", "TGO", "TUV", "UGA", "TZA", "VUT", "YEM", "ZMB"
)
stopifnot(length(ldc_iso3_2013) == 49L, !anyDuplicated(ldc_iso3_2013))

# 2024 vintage of the UN LDC list, used only for the disclosure below.
ldc_iso3_2024 <- c(
  "AFG", "AGO", "BGD", "BEN", "BFA", "BDI", "KHM", "CAF", "TCD", "COM",
  "COD", "DJI", "ERI", "ETH", "GMB", "GIN", "GNB", "HTI", "KIR", "LAO",
  "LSO", "LBR", "MDG", "MWI", "MLI", "MRT", "MOZ", "MMR", "NPL", "NER",
  "RWA", "SEN", "SLE", "SLB", "SOM", "SSD", "SDN", "STP", "TLS", "TGO",
  "TUV", "UGA", "TZA", "YEM", "ZMB"
)
stopifnot(length(ldc_iso3_2024) == 45L, !anyDuplicated(ldc_iso3_2024))

ldc_iso3 <- ldc_iso3_2013

ldc_switch <- setdiff(ldc_iso3_2013, ldc_iso3_2024)
message(sprintf(paste0("  LDC vintage: UN list as of 2013 (%d economies). ",
                       "Reclassified from non-LDC (2024 list) to LDC (2013 list): %s"),
                length(ldc_iso3_2013), paste(ldc_switch, collapse = ", ")))
in_panel_switch <- intersect(ldc_switch, unique(did_panel_full$recipient_iso))
message(sprintf("  Of those, present in the estimation panel: %s",
                if (length(in_panel_switch) == 0L) "none" else
                  paste(sort(in_panel_switch), collapse = ", ")))
stopifnot(length(setdiff(ldc_iso3_2024, ldc_iso3_2013)) == 0L)

did_panel_full <- did_panel_full %>%
  mutate(is_ldc = as.integer(recipient_iso %in% ldc_iso3))

# Cross-check against the legacy NAP-Central flag: every legacy-flagged LDC
# must also be on the 2013 UN list (the reverse cannot hold — see above).
if ("ldc_sids" %in% names(did_panel_full)) {
  legacy_ldc <- unique(did_panel_full$recipient_iso[
    grepl("LDC", did_panel_full$ldc_sids, fixed = TRUE)])
  if (length(setdiff(legacy_ldc, ldc_iso3)) > 0L) {
    warning("Legacy ldc_sids flags not on the UN list: ",
            paste(setdiff(legacy_ldc, ldc_iso3), collapse = ", "))
  }
}

message(sprintf("  LDC countries: %d | Non-LDC: %d (UN list, ISO3 match)",
  n_distinct(did_panel_full$recipient_name[did_panel_full$is_ldc == 1L]),
  n_distinct(did_panel_full$recipient_name[did_panel_full$is_ldc == 0L])))

ldc_groups <- list(
  list(label = "LDC",     panel = did_panel_full %>% filter(is_ldc == 1L)),
  list(label = "Non-LDC", panel = did_panel_full %>% filter(is_ldc == 0L))
)

# Wide table: Seed rule (reproduces the published SEs): set.seed(1242) inside make_het_wide_table
het_stats_ldc <- make_het_wide_table(
  groups      = ldc_groups,
  outcome_var = "log_commits",
  raw_var     = "commitments",
  dir_tabs    = dir_tabs_h2,
  tex_label   = "tab:het_ldc_wide",
  caption_txt = sprintf(paste0(
    "Heterogeneity by LDC status: ATT on log(adaptation commitments). CS ",
    "(2021), not-yet-treated controls. LDC status: pre-treatment vintage, UN ",
    "2013 list (%d economies); reclassifications vs.\\ 2024 in text."),
    length(ldc_iso3_2013))
)

# Event-study curves (one panel of the 2x2 figure, §24b)
h2_dyn <- setNames(vector("list", length(ldc_groups)),
                   vapply(ldc_groups, `[[`, character(1L), "label"))
# Seed set immediately before the estimator (reproduces the published SE)
set.seed(1242)
for (grp in ldc_groups) {
  panel_g <- grp$panel %>% filter(!(cohort_year %in% thin_cohorts))
  if (n_distinct(panel_g$cohort_year[panel_g$cohort_year > 0]) < 2)
    stop("Event study: subgroup ", grp$label, " has fewer than two treated cohorts")
  n_tr_h2 <- n_distinct(panel_g$country_id[panel_g$cohort_year > 0])
  em_h2   <- if (n_tr_h2 >= 40L) "dr" else "reg"
  bs_h2   <- TRUE   # Bootstrap everywhere (no n >= 40 switch)
  gt_h2 <- tryCatch(
    att_gt(yname = "log_commits", tname = "year", idname = "country_id",
           gname = "cohort_year", xformla = ~ ge_est + log_population, data = panel_g,
           est_method = em_h2, bstrap = bs_h2, biters = BITERS, cband = FALSE,
           control_group = "notyettreated", anticipation = 0,
           base_period = "universal", panel = TRUE, allow_unbalanced_panel = TRUE),
    error = function(e) NULL)
  if (is.null(gt_h2)) stop("Event study: att_gt() failed for subgroup ", grp$label)
  agg_d2 <- tryCatch(aggte(gt_h2, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf),
                     error = function(e) NULL)
  if (is.null(agg_d2)) stop("Event study: aggte() failed for subgroup ", grp$label)
  # Simultaneous (sup-t) 95% band, uniform over this curve's event times (sup_t_crit()).
  cv2    <- sup_t_crit(agg_d2$inf.function$dynamic.inf.func.e, agg_d2$se.egt, biters = BITERS)
  message(sprintf("  LDC ES sup-t crit (%s): %.4f", grp$label, cv2))
  color2 <- if (grp$label == "LDC") "#e08214" else "#542788"
  h2_dyn[[grp$label]] <- data.frame(
    group = grp$label, color = color2,
    event_time = agg_d2$egt, ATT = agg_d2$att.egt, SE = agg_d2$se.egt,
    Lower = agg_d2$att.egt - cv2 * agg_d2$se.egt,
    Upper = agg_d2$att.egt + cv2 * agg_d2$se.egt,
    stringsAsFactors = FALSE)
}
if (!all(vapply(h2_dyn, is.null, logical(1L)))) {
  h2_df <- bind_rows(h2_dyn) %>%
    filter(!is.na(SE), SE > 1e-10) %>%
    mutate(group = factor(group, levels = c("LDC", "Non-LDC")))
  palette_h2 <- c(LDC = "#e08214", `Non-LDC` = "#542788")
  p_h2 <- ggplot(h2_df, aes(x = event_time, y = ATT,
                              colour = group, shape = group, group = group)) +
    geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = -0.5, colour = "grey60", linetype = "dotted", linewidth = 0.4) +
    geom_linerange(aes(ymin = Lower, ymax = Upper),
                   position = position_dodge(0.4), linewidth = 0.6, alpha = 0.8) +
    geom_point(size = 2.5, position = position_dodge(0.4)) +
    scale_colour_manual(values = palette_h2) +
    scale_shape_manual(values = c(16, 17)) +
    # No title, subtitle, or caption — those go in LaTeX \caption{}
    labs(title = NULL, subtitle = NULL, caption = NULL,
         x = "Event time", y = "ATT — log(adaptation commitments)",
         colour = NULL, shape = NULL) +
    theme_minimal() +
    theme(text             = element_text(family = "serif", size = 11),
          legend.position  = "bottom", panel.grid.minor = element_blank())
}

message("\n=== Section 22 complete ===\n")

# ==============================================================================
# SECTION 7. §23 Governance heterogeneity (high vs. low WGI GE)
# Seed rule (reproduces the published SEs): set.seed(1242) before gov loop.
# ==============================================================================

message("\n=== Section 23: Heterogeneity — High vs. low governance ===\n")

dir_tabs_h3 <- here("output", "tables",  "heterogeneity", "governance")
dir.create(dir_tabs_h3, recursive = TRUE, showWarnings = FALSE)

ge_baseline <- did_panel_full %>%
  filter(cohort_year > 0) %>%
  group_by(country_id, recipient_name) %>%
  summarise(ge_base = mean(ge_est[year < cohort_year], na.rm = TRUE), .groups = "drop")

ge_never <- did_panel_full %>%
  filter(cohort_year == 0) %>%
  group_by(country_id, recipient_name) %>%
  summarise(ge_base = mean(ge_est, na.rm = TRUE), .groups = "drop")

ge_all  <- bind_rows(ge_baseline, ge_never)
ge_med  <- median(ge_all$ge_base, na.rm = TRUE)
ge_all  <- ge_all %>% mutate(hi_gov = as.integer(ge_base >= ge_med))
message(sprintf("  GE median split: %.3f  |  hi_gov N=%d  lo_gov N=%d",
                ge_med,
                sum(ge_all$hi_gov == 1L, na.rm = TRUE),
                sum(ge_all$hi_gov == 0L, na.rm = TRUE)))

did_panel_gov <- did_panel_full %>%
  left_join(ge_all %>% select(country_id, ge_base, hi_gov), by = "country_id")

gov_groups <- list(
  list(label = "High governance", panel = did_panel_gov %>% filter(hi_gov == 1L)),
  list(label = "Low governance",  panel = did_panel_gov %>% filter(hi_gov == 0L))
)

# Wide table: Seed rule (reproduces the published SEs): set.seed(1242) inside make_het_wide_table
het_stats_gov <- make_het_wide_table(
  groups      = gov_groups,
  outcome_var = "log_commits",
  raw_var     = "commitments",
  dir_tabs    = dir_tabs_h3,
  tex_label   = "tab:het_gov_wide",
  caption_txt = "Heterogeneity by governance level: ATT on log(adaptation commitments). CS (2021), not-yet-treated controls; median WGI GE split at baseline."
)

# Event-study curves (one panel of the 2x2 figure, §24b)
h3_dyn <- setNames(vector("list", length(gov_groups)),
                   vapply(gov_groups, `[[`, character(1L), "label"))
# Seed set immediately before the estimator (reproduces the published SE)
set.seed(1242)
for (grp in gov_groups) {
  panel_g <- grp$panel %>% filter(!(cohort_year %in% thin_cohorts))
  if (n_distinct(panel_g$cohort_year[panel_g$cohort_year > 0]) < 2)
    stop("Event study: subgroup ", grp$label, " has fewer than two treated cohorts")
  n_tr_h3 <- n_distinct(panel_g$country_id[panel_g$cohort_year > 0])
  em_h3   <- if (n_tr_h3 >= 40L) "dr" else "reg"
  bs_h3   <- TRUE   # Bootstrap everywhere (no n >= 40 switch)
  gt_h3 <- tryCatch(
    att_gt(yname = "log_commits", tname = "year", idname = "country_id",
           gname = "cohort_year", xformla = ~ ge_est + log_population, data = panel_g,
           est_method = em_h3, bstrap = bs_h3, biters = BITERS, cband = FALSE,
           control_group = "notyettreated", anticipation = 0,
           base_period = "universal", panel = TRUE, allow_unbalanced_panel = TRUE),
    error = function(e) NULL)
  if (is.null(gt_h3)) stop("Event study: att_gt() failed for subgroup ", grp$label)
  agg_d3 <- tryCatch(aggte(gt_h3, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf),
                     error = function(e) NULL)
  if (is.null(agg_d3)) stop("Event study: aggte() failed for subgroup ", grp$label)
  # Simultaneous (sup-t) 95% band, uniform over this curve's event times (sup_t_crit()).
  cv3    <- sup_t_crit(agg_d3$inf.function$dynamic.inf.func.e, agg_d3$se.egt, biters = BITERS)
  message(sprintf("  Governance ES sup-t crit (%s): %.4f", grp$label, cv3))
  # Event-time coefficients behind the governance figure, quoted in the text.
  message(sprintf("  Governance event study [%s]: %s", grp$label,
                  paste(sprintf("e=%d %.4f (SE %.4f)", as.integer(agg_d3$egt),
                                agg_d3$att.egt, agg_d3$se.egt), collapse = "; ")))
  color3 <- if (grp$label == "High governance") "#1b7837" else "#762a83"
  h3_dyn[[grp$label]] <- data.frame(
    group = grp$label, color = color3,
    event_time = agg_d3$egt, ATT = agg_d3$att.egt, SE = agg_d3$se.egt,
    Lower = agg_d3$att.egt - cv3 * agg_d3$se.egt,
    Upper = agg_d3$att.egt + cv3 * agg_d3$se.egt,
    stringsAsFactors = FALSE)
}
if (!all(vapply(h3_dyn, is.null, logical(1L)))) {
  h3_df <- bind_rows(h3_dyn) %>%
    filter(!is.na(SE), SE > 1e-10) %>%
    mutate(group = factor(group, levels = c("High governance", "Low governance")))
  palette_h3 <- c(`High governance` = "#1b7837", `Low governance` = "#762a83")
  p_h3 <- ggplot(h3_df, aes(x = event_time, y = ATT,
                              colour = group, shape = group, group = group)) +
    geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = -0.5, colour = "grey60", linetype = "dotted", linewidth = 0.4) +
    geom_linerange(aes(ymin = Lower, ymax = Upper),
                   position = position_dodge(0.4), linewidth = 0.6, alpha = 0.8) +
    geom_point(size = 2.5, position = position_dodge(0.4)) +
    scale_colour_manual(values = palette_h3) +
    scale_shape_manual(values = c(16, 17)) +
    # No title, subtitle, or caption — those go in LaTeX \caption{}
    labs(title = NULL, subtitle = NULL, caption = NULL,
         x = "Event time (years relative to NAP adoption)",
         y = "ATT — log(adaptation commitments)", colour = NULL, shape = NULL) +
    theme_minimal() +
    theme(text             = element_text(family = "serif", size = 11),
          legend.position  = "bottom", panel.grid.minor = element_blank())
}

message("\n=== Section 23 complete ===\n")

# ==============================================================================
# SECTION 8. §24 Income group heterogeneity
# Seed rule (reproduces the published SEs): set.seed(1242) before income loop.
# ==============================================================================

message("\n=== Section 24: Heterogeneity — Income group ===\n")

het_stats_income <- NULL   # filled only if the income split is estimable

dir_tabs_h4 <- here("output", "tables",  "heterogeneity", "income_group")
dir.create(dir_tabs_h4, recursive = TRUE, showWarnings = FALSE)

# PRE-TREATMENT INCOME VINTAGE (FY2013 = July 2012 classification).
# Income classification joined on ISO3 (recipient_iso <-> iso3c); a name-based
# join would drop 18 countries whose CRS names differ from WDI names (e.g.,
# "China (People's Republic of)", "Democratic Republic of the Congo", "Egypt",
# "Yemen", "Turkiye"). The VINTAGE also matters. The classification
# shipped with the installed WDI package is a CURRENT cross-section: a recipient
# that moved from lower-middle to upper-middle income during the estimation
# window would be assigned its post-treatment group, so the split would condition
# partly on an outcome. We therefore use the World Bank's own historical file
# (OGHIST.xlsx, "Country Analytical History" sheet) and read the FY2013 column
# -- the classification announced in July 2012, i.e. strictly before the first
# NAP cohort in the estimation sample (2021).
# The file is downloaded once to data/raw/oghist/ and cached; if neither the
# cache nor the download is available the script STOPS (no fallback to the WDI
# snapshot, which is post-treatment; see read_income_classification()).
OGHIST_URL  <- paste0("https://datacatalogfiles.worldbank.org/ddh-published/",
                      "0037712/DR0090754/OGHIST.xlsx")
oghist_path <- here("data", "raw", "oghist", "OGHIST.xlsx")
dir.create(dirname(oghist_path), recursive = TRUE, showWarnings = FALSE)

# SHA-256 of the vintage behind the published exhibits, recorded in the
# README (Data availability and provenance). A silently different OGHIST would change the income
# split without changing anything visible, so the checksum is verified on every
# run: a mismatch is a warning with both digests, not a stop, because the World
# Bank does reissue the file and the authors must decide whether to adopt the new
# vintage and update the README checksum.
OGHIST_SHA256 <- "17eb9e67b2eaf7d489ceb303f39f53697ce5b6238ed13c0b11c0846c33024d7b"

if (!file.exists(oghist_path)) {
  message("  OGHIST cache missing — downloading from ", OGHIST_URL)
  ok_dl <- tryCatch({
    utils::download.file(OGHIST_URL, destfile = oghist_path, mode = "wb",
                         quiet = TRUE)
    file.exists(oghist_path) && file.size(oghist_path) > 10000
  }, error = function(e) { message("  download failed: ", conditionMessage(e)); FALSE },
     warning = function(w) { message("  download warning: ", conditionMessage(w)); FALSE })
  if (!isTRUE(ok_dl) && file.exists(oghist_path)) unlink(oghist_path)
}

if (file.exists(oghist_path)) {
  oghist_sha <- tryCatch(
    as.character(tools::sha256sum(oghist_path)),
    error = function(e) NA_character_)
  if (is.na(oghist_sha)) {
    warning("Could not compute the SHA-256 of ", oghist_path)
  } else if (!identical(oghist_sha, OGHIST_SHA256)) {
    warning("OGHIST.xlsx checksum mismatch.\n",
            "  expected (README, Data availability): ", OGHIST_SHA256, "\n",
            "  found on disk:                   ", oghist_sha, "\n",
            "  The World Bank has reissued the workbook. Verify the FY13 ",
            "column, then update the README checksum and OGHIST_SHA256 here.")
  } else {
    message("  OGHIST.xlsx SHA-256 verified against the published checksum (README).")
  }
}

#' Read the FY2013 World Bank income classification.
#'
#' HARD STOP, no silent fallback. The FY2013 vintage
#' is not a nicety: the split is a proxy for pre-treatment capacity, and 53 of
#' the panel's recipients are classified differently under the current World
#' Bank vintage, which would silently condition the heterogeneity analysis on
#' post-treatment information. Falling back to the current WDI snapshot would
#' change the reported ATTs while the table note still claimed a pre-treatment
#' classification unless every downstream caption were rebuilt, so a missing or
#' unreadable OGHIST file stops the script with instructions instead.
#' The file is cached under data/raw/oghist/ and documented in the README
#' (Data availability and provenance).
#'
#' @param path path to the cached OGHIST.xlsx
#' @return list(lookup = data.frame, vintage = character)
read_income_classification <- function(path) {
  if (!file.exists(path)) {
    stop("World Bank OGHIST workbook not found: ", path, "\n",
         "  The income-group split requires the FY2013 (July 2012) vintage; ",
         "the current WDI snapshot is a post-treatment classification and is ",
         "NOT an acceptable substitute.\n",
         "  Fix: download it once (the script attempts this automatically when ",
         "the cache is absent):\n",
         "    curl -L -o data/raw/oghist/OGHIST.xlsx \\\n",
         "      ", OGHIST_URL, "\n",
         "  See README.md, section Data availability and provenance.")
  }
  raw_og <- as.data.frame(read_excel(path,
                                     sheet = "Country Analytical History",
                                     col_names = FALSE, .name_repair = "minimal"))
  # Row 5 of the sheet holds the Bank fiscal-year headers (FY89, FY90, ...);
  # country rows begin below the four threshold rows and carry an ISO3 code in
  # column 1 and the country name in column 2.
  fy_row <- as.character(unlist(raw_og[5L, ]))
  j_fy13 <- which(fy_row == "FY13")
  if (length(j_fy13) != 1L)
    stop("FY13 column not found in ", path, " (found ", length(j_fy13),
         " matches in the fiscal-year header row). The workbook layout has ",
         "changed; update read_income_classification().")
  body_og <- raw_og[-seq_len(11L), , drop = FALSE]
  iso_og  <- as.character(body_og[[1L]])
  keep_og <- !is.na(iso_og) & nchar(iso_og) == 3L
  if (sum(keep_og) < 150L)
    stop("Only ", sum(keep_og), " economies parsed from ", path,
         "; expected at least 150. Aborting rather than estimating the income ",
         "split on a truncated classification.")
  code_map <- c(L = "Low income", LM = "Lower middle income",
                UM = "Upper middle income", H = "High income")
  lk <- data.frame(
    recipient_iso = iso_og[keep_og],
    og_code       = as.character(body_og[[j_fy13]])[keep_og],
    stringsAsFactors = FALSE)
  lk$income_group <- unname(code_map[lk$og_code])
  lk$income_group[is.na(lk$income_group)] <- "Not classified"
  list(lookup  = distinct(lk[, c("recipient_iso", "income_group")]),
       vintage = paste0("the World Bank historical income classification, FY2013 ",
                        "(July 2012) column"))
}

income_res    <- read_income_classification(oghist_path)
income_lookup <- income_res$lookup
income_vintage <- income_res$vintage
message("  Income classification vintage: ", income_vintage)
if (!is.null(income_lookup)) print(table(income_lookup$income_group))

# Disclosure: which panel recipients change income group between the two
# vintages. Reported in the console and summarised in the table note.
income_switchers <- character(0)
# read_income_classification() stops unless the FY2013 column exists, so a
# non-NULL lookup is always the FY2013 vintage.
if (!is.null(income_lookup)) {
  wdi_now <- local({
    cty <- WDI::WDI_data$country
    inc_col <- names(cty)[grepl("income", names(cty), ignore.case = TRUE)][1]
    iso_col <- names(cty)[grepl("^iso3c$", names(cty), ignore.case = TRUE)][1]
    cty %>% select(all_of(c(iso_col, inc_col))) %>%
      rename(recipient_iso = 1, income_now = 2) %>% distinct()
  })
  if (!is.null(wdi_now)) {
    panel_iso <- unique(did_panel_full$recipient_iso)
    cmp <- income_lookup %>%
      filter(recipient_iso %in% panel_iso) %>%
      inner_join(wdi_now, by = "recipient_iso") %>%
      filter(!is.na(income_now), income_group != income_now)
    income_switchers <- cmp$recipient_iso
    message(sprintf("  Income vintage: %d of %d panel recipients change group ",
                    nrow(cmp), length(panel_iso)),
            "between FY2013 and the current WDI snapshot.")
    if (nrow(cmp) > 0L)
      for (r in seq_len(nrow(cmp)))
        message(sprintf("    %-5s FY2013: %-20s -> current: %s",
                        cmp$recipient_iso[r], cmp$income_group[r], cmp$income_now[r]))
  }
}

if (!is.null(income_lookup)) {
  did_panel_inc <- did_panel_full %>%
    left_join(income_lookup, by = "recipient_iso")

  # Coverage audit: every panel country must carry an income group; the split
  # then excludes High income / Not classified by design (documented in the
  # table note), never by silent join failure.
  inc_audit <- did_panel_inc %>%
    distinct(recipient_name, income_group) %>%
    count(income_group, name = "n_countries") %>%
    mutate(income_group = ifelse(is.na(income_group), "Unmatched (join failure)",
                                 income_group))
  message(sprintf("  Income-group coverage (all %d panel countries):",
                  n_distinct(did_panel_inc$recipient_name)))
  for (r in seq_len(nrow(inc_audit))) {
    message(sprintf("    %-22s %d", inc_audit$income_group[r],
                    inc_audit$n_countries[r]))
  }
  n_highinc <- sum(inc_audit$n_countries[inc_audit$income_group == "High income"])
  n_notclass <- sum(inc_audit$n_countries[
    inc_audit$income_group == "Not classified"])
  n_unmatched <- sum(is.na(did_panel_inc %>%
                             distinct(recipient_name, income_group) %>%
                             pull(income_group)))
  if (n_unmatched > 0L) {
    warning(sprintf("Income join left %d countries unmatched — check ISO3 codes",
                    n_unmatched))
  }

  income_levels <- c("Low income", "Lower middle income", "Upper middle income")
  inc_colors    <- c("Low income" = "#d73027",
                     "Lower middle income" = "#fc8d59",
                     "Upper middle income" = "#4575b4")

  inc_groups <- lapply(income_levels, function(lvl) {
    list(label = lvl,
         panel = did_panel_inc %>% filter(income_group == lvl))
  })
  inc_groups <- inc_groups[vapply(inc_groups, function(g)
    n_distinct(g$panel$cohort_year[g$panel$cohort_year > 0]) >= 2L,
    logical(1L))]

  if (length(inc_groups) >= 2L) {
    # Wide table: Seed rule (reproduces the published SEs): set.seed(1242) inside make_het_wide_table
    het_stats_income <- make_het_wide_table(
      groups      = inc_groups,
      outcome_var = "log_commits",
      raw_var     = "commitments",
      dir_tabs    = dir_tabs_h4,
      tex_label   = "tab:het_income_wide",
      caption_txt = sprintf(paste0(
        "Heterogeneity by World Bank income group: ATT on log(adaptation commitments). CS ",
        "(2021), not-yet-treated controls. Income groups: %s, held fixed."),
        income_vintage)
    )

    # Event-study curves (one panel of the 2x2 figure, §24b)
    h4_dyn <- setNames(vector("list", length(inc_groups)),
                       vapply(inc_groups, `[[`, character(1L), "label"))
    # Seed set immediately before the estimator (reproduces the published SE)
    set.seed(1242)
    for (grp in inc_groups) {
      panel_g <- grp$panel %>% filter(!(cohort_year %in% thin_cohorts))
      if (n_distinct(panel_g$cohort_year[panel_g$cohort_year > 0]) < 2L)
        stop("Event study: subgroup ", grp$label, " has fewer than two treated cohorts")
      n_tr_h4 <- n_distinct(panel_g$country_id[panel_g$cohort_year > 0])
      em_h4   <- if (n_tr_h4 >= 40L) "dr" else "reg"
      bs_h4   <- TRUE   # Bootstrap everywhere (no n >= 40 switch)
      gt_h4 <- tryCatch(
        att_gt(yname = "log_commits", tname = "year", idname = "country_id",
               gname = "cohort_year", xformla = ~ ge_est + log_population, data = panel_g,
               est_method = em_h4, bstrap = bs_h4, biters = BITERS, cband = FALSE,
               control_group = "notyettreated", anticipation = 0,
               base_period = "universal", panel = TRUE, allow_unbalanced_panel = TRUE),
        error = function(e) NULL)
      if (is.null(gt_h4)) stop("Event study: att_gt() failed for subgroup ", grp$label)
      agg_d4 <- tryCatch(aggte(gt_h4, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf),
                         error = function(e) NULL)
      if (is.null(agg_d4)) stop("Event study: aggte() failed for subgroup ", grp$label)
    # Simultaneous (sup-t) 95% band, uniform over this curve's event times (sup_t_crit()).
      cv4 <- sup_t_crit(agg_d4$inf.function$dynamic.inf.func.e, agg_d4$se.egt, biters = BITERS)
      message(sprintf("  Income ES sup-t crit (%s): %.4f", grp$label, cv4))
      h4_dyn[[grp$label]] <- data.frame(
        group = grp$label, color = inc_colors[grp$label],
        event_time = agg_d4$egt, ATT = agg_d4$att.egt, SE = agg_d4$se.egt,
        Lower = agg_d4$att.egt - cv4 * agg_d4$se.egt,
        Upper = agg_d4$att.egt + cv4 * agg_d4$se.egt,
        stringsAsFactors = FALSE)
    }
    if (!all(vapply(h4_dyn, is.null, logical(1L)))) {
      h4_df <- bind_rows(h4_dyn) %>%
        filter(!is.na(SE), SE > 1e-10) %>%
        mutate(group = factor(group, levels = income_levels))
      palette_h4 <- inc_colors[levels(h4_df$group)]
      p_h4 <- ggplot(h4_df, aes(x = event_time, y = ATT,
                                  colour = group, shape = group, group = group)) +
        geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed", linewidth = 0.4) +
        geom_vline(xintercept = -0.5, colour = "grey60", linetype = "dotted", linewidth = 0.4) +
        geom_linerange(aes(ymin = Lower, ymax = Upper),
                       position = position_dodge(0.4), linewidth = 0.6, alpha = 0.8) +
        geom_point(size = 2.5, position = position_dodge(0.4)) +
        scale_colour_manual(values = palette_h4) +
        scale_shape_manual(values = c(16, 17, 15)) +
        # No title, subtitle, or caption — those go in LaTeX \caption{}
        labs(title = NULL, subtitle = NULL, caption = NULL,
             x = "Event time", y = "ATT — log(adaptation commitments)",
             colour = NULL, shape = NULL) +
        theme_minimal() +
        theme(text             = element_text(family = "serif", size = 11),
              legend.position  = "bottom", panel.grid.minor = element_blank())
    }
  } else {
    stop("Too few income groups with sufficient cohorts: no income event study.")
  }
} else {
  stop("Income group data unavailable: Section 24 not run.")
}

message("\n=== Section 24 complete ===\n")

# ==============================================================================
# SECTION 8a. §24b Compact capacity table and 2x2 event-study panel
#
# Layout only: nothing is re-estimated. het_capacity.tex re-uses the formatted
# cells that make_het_wide_table() wrote into the three wide tables (ATT with
# stars, SE, N treated, four-lead pre-trend p), so every value is identical to
# Tables tab:het_gov_wide, tab:het_ldc_wide and tab:het_income_wide. The 2x2
# figure re-uses the four plotted objects saved above (p_h1, p_h3, p_h2, p_h4),
# with harmonised axis titles and text size only.
# ==============================================================================

message("\n=== Section 24b: compact capacity table and ES panel ===\n")

cap_panels <- list(
  list(title = "Panel A: Governance (WGI Government Effectiveness, median split at baseline)",
       stats = het_stats_gov, cells = c("High governance", "Low governance")),
  list(title = sprintf("Panel B: LDC status (UN list as of 2013, %d economies)",
                       length(ldc_iso3_2013)),
       stats = het_stats_ldc, cells = c("LDC", "Non-LDC")),
  list(title = "Panel C: World Bank income group (FY2013 classification, held fixed)",
       stats = het_stats_income,
       cells = c("Low income", "Lower middle income", "Upper middle income"))
)
cap_cells <- unlist(lapply(cap_panels, function(pn) pn$stats[pn$cells]), recursive = FALSE)
if (any(vapply(cap_cells, is.null, logical(1L))))
  stop("het_capacity: a subgroup fit is missing.")
# One lead window and restriction count for all seven cells, so the note can
# state it once (it is the same window as in the three wide tables).
cap_leads <- unique(lapply(cap_cells, `[[`, "pt_leads"))
cap_df    <- unique(vapply(cap_cells, function(s) as.integer(s$pt_df), integer(1L)))
stopifnot(length(cap_leads) == 1L, length(cap_df) == 1L)
cap_em    <- unique(vapply(cap_cells, `[[`, character(1L), "est_method"))
cap_ntr   <- vapply(cap_cells, function(s) as.integer(s$n_treated), integer(1L))

cap_lines <- c(
  "\\begingroup\\small",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  sprintf("Subgroup & ATT & SE & $N$ treated & Pre-trend $p$ (%d leads) \\\\", length(cap_leads[[1L]])),
  "\\midrule",
  unlist(lapply(seq_along(cap_panels), function(k) {
    pn <- cap_panels[[k]]
    c(if (k > 1L) "\\addlinespace" else character(0),
      sprintf("\\multicolumn{5}{l}{\\textit{%s}} \\\\", pn$title),
      vapply(pn$cells, function(cl) {
        st <- pn$stats[[cl]]
        sprintf("\\quad %s & %s & %s & %s & %s \\\\", cl, st$att_fmt, st$se_fmt,
                st$n_treated_fmt, st$pt_pval)
      }, character(1L)))
  })),
  "\\bottomrule",
  "\\end{tabular}",
  "\\endgroup"
)

cap_est_txt <- if (identical(cap_em, "reg")) {
  sprintf(paste0("CS (2021) regression adjustment (doubly robust only for a subgroup with at ",
                 "least 40 treated recipients; none here, %d--%d per cell)"),
          min(cap_ntr), max(cap_ntr))
} else {
  "CS (2021), doubly robust for subgroups with at least 40 treated recipients, regression adjustment otherwise"
}
cap_notes <- paste0(
  "ATT on log(adaptation commitments), estimated separately on each subgroup; full columns in ",
  "Tables~\\ref{tab:het_gov_wide}, \\ref{tab:het_ldc_wide} and \\ref{tab:het_income_wide}. ",
  cap_est_txt, ", not-yet-treated controls, WGI GE + log population. ",
  "ATT/SE: multiplier-bootstrap (", BITERS, " reps, clustered by recipient, seed 1242); ",
  "pre-trend: separate analytical fit. ",
  PRETREND_NOTE_AGG(min(cap_leads[[1L]]), max(cap_leads[[1L]]), cap_df), ". ",
  "Panel A: median of pre-treatment mean WGI GE. Panel B: pre-treatment vintage of the UN LDC ",
  "list. Panel C: ", income_vintage, "; high-income and unclassified recipients excluded. ",
  "* $p<0.10$, ** $p<0.05$, *** $p<0.01$"
)
write_tex_float(
  file.path(here("output", "tables", "heterogeneity"), "het_capacity.tex"),
  "NAP effect by baseline capacity: governance, LDC status and income group",
  "tab:het_capacity", cap_lines, cap_notes,
  paste0("OECD CRS (Rio adaptation markers); UNFCCC NAP Central; WGI; UN LDC list; ",
         "World Bank income classification"))

# --- 2x2 panel: donor type, governance, LDC, income -------------------------
es_plots <- list(p_h1, p_h3, p_h2, p_h4)
if (!all(vapply(es_plots, inherits, logical(1L), what = "ggplot")))
  stop("ES panel: one of the four heterogeneity event-study plots is missing.")
es_plots <- lapply(es_plots, function(pl) pl +
  labs(x = "Event time (years relative to NAP adoption)", y = "ATT (log points)") +
  guides(colour = guide_legend(nrow = 1L), shape = guide_legend(nrow = 1L)) +
  theme(text = element_text(family = "serif", size = 12),
        legend.text = element_text(size = 11),
        legend.key.width = unit(10, "pt"),
        legend.key.spacing.x = unit(6, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(20, 6, 4, 6)))
p_het_panel <- cowplot::plot_grid(
  plotlist = es_plots, ncol = 2L, align = "hv",
  labels = c("(a) Donor type", "(b) Governance", "(c) LDC status", "(d) Income group"),
  label_size = 12, label_fontfamily = "serif", label_fontface = "plain",
  hjust = 0, label_x = 0.02, label_y = 0.995)
out_panel <- here("output", "figures", "heterogeneity", "fig_het_es_panel.png")
ggsave(out_panel, p_het_panel, width = 10, height = 8, dpi = 300, bg = "white")
message("Saved: ", out_panel)

# ==============================================================================
# SECTION 8b. §25a Correlated-sample contrasts
#
# Two families of contrast that the H4 machinery above cannot handle, because
# they compare ATTs estimated on the SAME recipient-years rather than on
# disjoint subsamples, so sqrt(se_a^2 + se_b^2) is wrong (it ignores a
# positive covariance and is therefore anti-conservative in the usual case):
#
#   H3 (donor type): DAC bilateral, multilateral and other-donor commitments
#     are a decomposition of the same recipient-year flows. Contrasts use the
#     unit-level influence functions of the very fits behind
#     Table~\ref{tab:het_donor_wide} (analytical twins of the reported
#     bootstrap fits), so Var(theta_a - theta_b) = E[(IFa - IFb)^2]/n.
#
#   Mitigation falsification: the headline adaptation ATT versus the
#     mitigation ATT of §26 of 04_robustness.R. Both fits are read from
#     output/fits/ (One fit, one SE) and their influence functions are
#     aligned on the common unit set, which is asserted, not assumed.
#
# Nothing is re-estimated here and nothing is stochastic: no seed is set.
# ==============================================================================

message("\n=== Section 25a: correlated-sample contrasts ===\n")

# Pre-sized result container: three donor-type pairwise contrasts
# plus the mitigation falsification; NULL slots are dropped before use.
extra_records <- vector("list", 4L)

# --- H3: donor-type contrasts -------------------------------------------------
donor_keys  <- c(DAC = "log_commits_dac", Multilateral = "log_commits_multi",
                 Other = "log_commits_other")
have_donor <- !is.null(donor_stats) && all(donor_keys %in% names(donor_stats)) &&
  all(vapply(donor_keys, function(k) !is.null(donor_stats[[k]]$inf_func),
             logical(1L)))

if (!have_donor) {
  stop("Donor-type fits or their influence functions are unavailable: no H3 contrasts.")
} else {
  ids_donor <- donor_stats[[donor_keys[["DAC"]]]]$ids
  if_len <- vapply(donor_keys, function(k) length(donor_stats[[k]]$inf_func),
                   integer(1L))
  id_len <- vapply(donor_keys, function(k) length(donor_stats[[k]]$ids),
                   integer(1L))
  # One influence-function entry per estimation unit, and the same units in all
  # three donor fits.
  stopifnot(length(unique(if_len)) == 1L, length(unique(id_len)) == 1L,
            if_len[[1L]] == length(ids_donor),
            all(vapply(donor_keys,
                       function(k) identical(donor_stats[[k]]$ids, ids_donor),
                       logical(1L))))
  message(sprintf("  H3: influence functions aligned across %d recipients ",
                  length(ids_donor)),
          "(identical estimation panel for the three donor outcomes).")

  # Sanity check on the SE convention: the analytical overall.se of each donor
  # fit must be reproduced by sqrt(mean(IF^2)/n). This is the same hard stop
  # 07_principal_and_share.R uses; it catches an IF/SE scaling mismatch before
  # any contrast is reported.
  for (nm in names(donor_keys)) {
    st <- donor_stats[[donor_keys[[nm]]]]
    se_from_if <- sqrt(mean(st$inf_func^2) / length(st$inf_func))
    message(sprintf("    %-12s analytical SE = %.6f | sqrt(mean(IF^2)/n) = %.6f",
                    nm, st$se_analytic, se_from_if))
    stopifnot(abs(se_from_if - st$se_analytic) < 1e-6)
  }

  donor_pairs <- list(
    list(label = "DAC bilateral $-$ multilateral", a = "DAC", b = "Multilateral"),
    list(label = "DAC bilateral $-$ other donors",  a = "DAC", b = "Other"),
    list(label = "Multilateral $-$ other donors",   a = "Multilateral", b = "Other")
  )
  for (i_dp in seq_along(donor_pairs)) {
    dp <- donor_pairs[[i_dp]]
    sa <- donor_stats[[donor_keys[[dp$a]]]]
    sb <- donor_stats[[donor_keys[[dp$b]]]]
    res <- compute_att_difference(
      att_a = sa$att_num, att_b = sb$att_num,
      IFa = sa$inf_func, IFb = sb$inf_func,
      se_a = sa$se_num, se_b = sb$se_num,
      ids_a = sa$ids, ids_b = sb$ids)
    extra_records[[i_dp]] <- list(
      split = "Donor type (H3)", label = dp$label, res = res,
      family = "H3")
    message(sprintf("  %-34s diff = %8.4f  SE = %7.4f  z = %7.3f  p = %.3f  [%s]",
                    dp$label, res$diff, res$se, res$z, res$pval, res$method))
  }

  # Joint test that the three donor-type ATTs are equal (2 df), using the same
  # aligned influence functions (NOT a diagonal covariance).
  IF_donor <- cbind(donor_stats[[donor_keys[["DAC"]]]]$inf_func,
                    donor_stats[[donor_keys[["Multilateral"]]]]$inf_func,
                    donor_stats[[donor_keys[["Other"]]]]$inf_func)
  theta_donor <- c(donor_stats[[donor_keys[["DAC"]]]]$att_analytic,
                   donor_stats[[donor_keys[["Multilateral"]]]]$att_analytic,
                   donor_stats[[donor_keys[["Other"]]]]$att_analytic)
  donor_joint <- att_joint_wald_if(theta_donor, IF_donor, ref_index = 3L)
  message(sprintf("  Joint (H3): chi2(%d) = %.3f  p = %.3f",
                  donor_joint$df, donor_joint$stat, donor_joint$pval))
}

# --- Mitigation falsification contrast ---------------------------------------
# The headline fit is stored by 03 (read_headline_fit() stops if it is missing
# or incomplete), the mitigation fit by 04. A missing fit stops the stage.
path_mit   <- file.path(FITS_DIR, "mitigation_dr_bs.rds")
mit_contrast <- NULL

if (!file.exists(path_mit)) {
  stop("Saved mitigation fit not found: ", path_mit, ". Run 04 before 05.")
} else {
  fit_adapt <- read_headline_fit("adaptation")
  fit_mit   <- readRDS(path_mit)
  IF_adapt  <- fit_adapt$agg_simple_analytic$inf.function$simple.att
  IF_mit    <- fit_mit$agg_simple_analytic$inf.function$simple.att

  # Hard checks before any contrast is formed: each influence-function vector
  # must have exactly one entry per unit of its own estimation sample, and the
  # two samples must be the same units. Without this the difference is formed
  # across misaligned rows and the SE is meaningless.
  if (!is.null(IF_adapt))
    stopifnot(length(IF_adapt) == length(fit_adapt$ids))
  if (!is.null(IF_mit))
    stopifnot(length(IF_mit) == length(fit_mit$ids))
  if (is.null(IF_adapt) || is.null(IF_mit)) {
    stop("A saved fit carries no simple-aggregation influence function: ",
         "no mitigation contrast.")
  } else {
    message(sprintf("  Mitigation contrast: n_adapt = %d units, n_mit = %d units, ",
                    length(fit_adapt$ids), length(fit_mit$ids)),
            sprintf("identical unit sets: %s",
                    identical(fit_adapt$ids, fit_mit$ids)))
    mit_contrast <- compute_att_difference(
      att_a = fit_adapt$att, att_b = fit_mit$att,
      IFa = IF_adapt, IFb = IF_mit,
      se_a = fit_adapt$se, se_b = fit_mit$se,
      ids_a = fit_adapt$ids, ids_b = fit_mit$ids)
    extra_records[[4L]] <- list(
      split = "Falsification", label = "Adaptation $-$ mitigation commitments",
      res = mit_contrast, family = "Mitigation")
    message(sprintf(paste0("  %-34s diff = %8.4f  SE = %7.4f  z = %7.3f  ",
                           "p = %.3f  [%s]"),
                    "Adaptation - mitigation", mit_contrast$diff, mit_contrast$se,
                    mit_contrast$z, mit_contrast$pval, mit_contrast$method))
    message(sprintf("    inputs: adaptation ATT = %.4f (SE %.4f) | mitigation ATT = %.4f (SE %.4f)",
                    fit_adapt$att, fit_adapt$se, fit_mit$att, fit_mit$se))
  }
}

# Rows appended to tab:het_difftests below (existing rows are untouched).
fmt_num_x <- function(x, digits = 4L)
  if (is.na(x)) "---" else sprintf(paste0("%.", digits, "f"), x)
fmt_p_x <- function(x) {
  if (is.na(x)) return("---")
  if (x < 0.001) "$<$0.001" else sprintf("%.3f", x)
}

#' Format a block of contrast records as table rows.
rows_from_records <- function(recs) {
  if (length(recs) == 0L) return(NULL)
  data.frame(
    Split      = vapply(recs, `[[`, character(1L), "split"),
    Contrast   = vapply(recs, `[[`, character(1L), "label"),
    Difference = vapply(recs, function(r) fmt_num_x(r$res$diff), character(1L)),
    SE         = vapply(recs, function(r) fmt_num_x(r$res$se),   character(1L)),
    z          = vapply(recs, function(r) fmt_num_x(r$res$z, 3L), character(1L)),
    p          = vapply(recs, function(r) fmt_p_x(r$res$pval),   character(1L)),
    stringsAsFactors = FALSE)
}

extra_records <- Filter(Negate(is.null), extra_records)
is_h3 <- vapply(extra_records, function(r) identical(r$family, "H3"), logical(1L))
joint_h3_row <- if (exists("donor_joint") && !is.null(donor_joint) &&
                    !is.na(donor_joint$stat)) data.frame(
  Split      = "Donor type (H3)",
  Contrast   = paste0("Joint: all three donor-type ATTs equal ",
                      "($\\chi^2 = ", fmt_num_x(donor_joint$stat, 3L), "$, ",
                      donor_joint$df, " df)"),
  Difference = "---", SE = "---", z = "---",
  p          = fmt_p_x(donor_joint$pval),
  stringsAsFactors = FALSE) else NULL

# Order: the three H3 pairwise contrasts, their joint test, then the
# falsification contrast.
extra_rows <- bind_rows(rows_from_records(extra_records[is_h3]),
                        joint_h3_row,
                        rows_from_records(extra_records[!is_h3]))
if (!is.null(extra_rows) && nrow(extra_rows) == 0L) extra_rows <- NULL

# ==============================================================================
# SECTION 9. §25 Formal difference tests between subgroup ATTs (H4)
#
# H4: the NAP effect is larger where baseline capacity is weaker. The three
# splits above (LDC status, governance median, World Bank income group) are read
# jointly as one test of H4, so the paper needs the contrasts themselves, not
# just the subgroup ATTs.
#
# No re-estimation: the ATTs and SEs below are the ones returned by
# make_het_wide_table(), i.e. exactly the fits behind Tables tab:het_ldc_wide,
# tab:het_gov_wide and tab:het_income_wide (CS 2021, est_method "reg" unless a
# subgroup has >= 40 treated units, multiplier-bootstrap SE clustered by
# recipient, not-yet-treated controls, cohorts with < 5 treated units dropped,
# xformla ~ ge_est + log_population, outcome log_commits).
#
# Inference: subgroups are disjoint sets of countries, each clustered by
# country, so the two subgroup estimators are asymptotically independent and
# se(diff) = sqrt(se_a^2 + se_b^2), with the bootstrap SEs. The joint
# income test uses the same argument through a diagonal covariance matrix.
# No new seed is set here: nothing in this section is stochastic.
# ==============================================================================

message("\n=== Section 25: Formal subgroup difference tests (H4) ===\n")

dir_tabs_h5 <- here("output", "tables", "heterogeneity")

have_ldc <- !is.null(het_stats_ldc) &&
  all(c("LDC", "Non-LDC") %in% names(het_stats_ldc))
have_gov <- !is.null(het_stats_gov) &&
  all(c("Low governance", "High governance") %in% names(het_stats_gov))
income_cells <- c("Low income", "Lower middle income", "Upper middle income")
have_inc <- !is.null(het_stats_income) &&
  all(income_cells %in% names(het_stats_income))

if (!have_ldc || !have_gov || !have_inc) {
  stop("Missing subgroup fits (LDC=", have_ldc, " gov=", have_gov,
       " income=", have_inc, "): difference-test table not written.")
} else {

  # --- Per-group inputs (unrounded, straight from the published fits) ---
  used_stats <- c(
    het_stats_ldc[c("LDC", "Non-LDC")],
    het_stats_gov[c("Low governance", "High governance")],
    het_stats_income[income_cells]
  )
  for (lbl in names(used_stats)) {
    message(sprintf("  input | %-20s ATT = %8.5f  SE = %7.5f  N treated = %d",
                    lbl, used_stats[[lbl]]$att_num, used_stats[[lbl]]$se_num,
                    used_stats[[lbl]]$n_treated))
  }

  # --- Pairwise contrasts: low capacity minus high capacity ---
  contrast_spec <- list(
    list(split = "LDC status",
         label = "LDC $-$ Non-LDC",
         a = het_stats_ldc[["LDC"]],          b = het_stats_ldc[["Non-LDC"]]),
    list(split = "Governance (WGI GE, median split)",
         label = "Low governance $-$ High governance",
         a = het_stats_gov[["Low governance"]],
         b = het_stats_gov[["High governance"]]),
    list(split = "Income group",
         label = "Low income $-$ Upper middle income",
         a = het_stats_income[["Low income"]],
         b = het_stats_income[["Upper middle income"]]),
    list(split = "Income group",
         label = "Lower middle income $-$ Upper middle income",
         a = het_stats_income[["Lower middle income"]],
         b = het_stats_income[["Upper middle income"]])
  )

  contrast_res <- lapply(contrast_spec, function(cs)
    att_difference_test(cs$a, cs$b))

  # --- Joint test: the three income-cell ATTs are all equal (2 df) ---
  inc_theta <- vapply(income_cells,
                      function(l) het_stats_income[[l]]$att_num, numeric(1L))
  inc_se    <- vapply(income_cells,
                      function(l) het_stats_income[[l]]$se_num,  numeric(1L))
  inc_joint <- att_joint_wald(inc_theta, inc_se,
                              ref_index = which(income_cells == "Upper middle income"))

  fmt_num <- function(x, digits = 4L)
    if (is.na(x)) "---" else sprintf(paste0("%.", digits, "f"), x)
  fmt_p <- function(x) {
    if (is.na(x)) return("---")
    if (x < 0.001) "$<$0.001" else sprintf("%.3f", x)
  }

  diff_tab <- data.frame(
    Split        = vapply(contrast_spec, `[[`, character(1L), "split"),
    Contrast     = vapply(contrast_spec, `[[`, character(1L), "label"),
    Difference   = vapply(contrast_res, function(r) fmt_num(r$diff), character(1L)),
    SE           = vapply(contrast_res, function(r) fmt_num(r$se),   character(1L)),
    z            = vapply(contrast_res, function(r) fmt_num(r$z, 3L), character(1L)),
    p            = vapply(contrast_res, function(r) fmt_p(r$pval),   character(1L)),
    stringsAsFactors = FALSE
  )

  joint_row <- data.frame(
    Split      = "Income group",
    Contrast   = paste0("Joint: all three income ATTs equal (",
                        "$\\chi^2 = ", fmt_num(inc_joint$stat, 3L),
                        "$, ", inc_joint$df, " df)"),
    Difference = "---",
    SE         = "---",
    z          = "---",
    p          = fmt_p(inc_joint$pval),
    stringsAsFactors = FALSE
  )
  diff_tab <- rbind(diff_tab, joint_row)

  # Append the correlated-sample contrasts computed in §25a below the four H4
  # rows and the joint income row.
  if (!is.null(extra_rows)) {
    names(extra_rows) <- names(diff_tab)
    diff_tab <- rbind(diff_tab, extra_rows)
  }

  # MDE column: minimum detectable difference at 80% power (code/functions/mde.R);
  # joint Wald rows have no single contrast SE, hence "---".
  fmt_mde <- function(se) if (is.na(se)) "---" else sprintf("%.4f", mde(se))
  diff_tab$MDE <- c(
    vapply(contrast_res, function(r) fmt_mde(r$se), character(1L)),
    "---",
    if (is.null(extra_rows)) character(0) else c(
      vapply(extra_records[is_h3], function(r) fmt_mde(r$res$se), character(1L)),
      if (is.null(joint_h3_row)) character(0) else "---",
      vapply(extra_records[!is_h3], function(r) fmt_mde(r$res$se), character(1L)))
  )
  stopifnot(length(diff_tab$MDE) == nrow(diff_tab))
  mde_reached <- c(
    vapply(contrast_res, function(r) abs(r$diff) >= mde(r$se), logical(1L)),
    vapply(extra_records, function(r) abs(r$res$diff) >= mde(r$res$se), logical(1L)))

  for (r in seq_len(nrow(diff_tab))) {
    message(sprintf("  %-36s diff = %8s  SE = %7s  z = %7s  p = %s",
                    diff_tab$Contrast[r], diff_tab$Difference[r],
                    diff_tab$SE[r], diff_tab$z[r], diff_tab$p[r]))
  }

  names(diff_tab) <- c("Split", "Contrast", "Difference", "SE", "$z$", "$p$-value",
                       "MDE (80\\% power)")

  xtab_h5 <- xtable(diff_tab, label = "tab:het_difftests")
  align(xtab_h5) <- "lllccccc"

  raw_lines_h5 <- capture.output(
    print(xtab_h5, include.rownames = FALSE, booktabs = TRUE,
          sanitize.text.function = identity, size = "\\small",
          floating = FALSE)
  )

  # Estimator label derived from the est_method rule actually applied (dr if a
  # subgroup has >= 40 treated units, otherwise regression adjustment).
  n_tr_used <- vapply(used_stats, function(s) as.integer(s$n_treated), integer(1L))
  est_label_h5 <- if (all(n_tr_used < 40L)) {
    "CS\\,(2021), reg.\\ adjustment"
  } else {
    "CS\\,(2021), DR/reg.\\ (by subgroup $N$)"
  }

  notes_h5 <- paste0(
    "Each row tests equality of two ATTs from the paired heterogeneity table ",
    "(not re-estimated). H4 rows: ", est_label_h5, " on log(adaptation commitments); H3 ",
    "(donor-type) rows: headline specification (DR) on the donor-group commitments; ",
    "falsification row: adaptation minus mitigation. H4 differences = low- minus ",
    "high-capacity (positive: larger effect where capacity weaker). H4 rows (LDC/governance/income): SE $=\\sqrt{se_a^2+se_b^2}$ on ",
    "multiplier-bootstrap SEs (", BITERS, " reps, seed 1242), disjoint subsamples. ",
    "H3/mitigation rows: aligned analytical-IF SEs (not bootstrap), same recipient-years. ",
    "$z$ vs.\\ standard normal; joint stats vs.\\ $\\chi^2$; $p$ two-sided. ",
    "MDE: minimum detectable difference, $(z_{0.975}+z_{0.80}) \\times SE = ",
    sprintf("%.4f", mde(1)), " \\times SE$, the smallest true difference a two-sided 5\\% ",
    "test rejects with 80\\% probability; none for the joint tests. ",
    if (all(!mde_reached)) "No estimated difference reaches its MDE" else
      sprintf("%d of %d estimated differences reach their MDE", sum(mde_reached),
              length(mde_reached))
  )
  source_h5 <- paste0("OECD CRS (Rio adaptation markers); UNFCCC NAP Central; ",
                      "WGI; World Bank income classification")

  out_path_h5 <- file.path(dir_tabs_h5, "het_difference_tests.tex")
  write_tex_float(
    out_path_h5,
    paste0("Formal tests of differences in the NAP effect across subgroups, donor types ",
           "and outcomes (H3, H4)"),
    "tab:het_difftests", raw_lines_h5, notes_h5, source_h5
  )

  # Mirror run_all.R's assembly step for this single file: output/tables/... ->
  # paper/Tables/... preserving the relative path. write_tex_float already emits
  # the \begin{table}[H] + \adjustbox form that fix_result_tables() enforces.
  # No-op when paper/ is absent (e.g. in the stand-alone replication package).
  if (dir.exists(here("paper"))) {
    paper_path_h5 <- here("paper", "Tables", "heterogeneity", "het_difference_tests.tex")
    dir.create(dirname(paper_path_h5), recursive = TRUE, showWarnings = FALSE)
    file.copy(out_path_h5, paper_path_h5, overwrite = TRUE)
    message("Saved: ", paper_path_h5)
  } else message("paper/ not found -- exhibits are left in output/ only")
}

message("\n=== Section 25 complete ===\n")

# ==============================================================================
# SECTION 11. §27 Zero shares by donor type
#
# The donor-type results are estimated on log1p outcomes, so how often each
# donor group commits nothing at all to a recipient in a year determines how
# much of the donor-type contrast is an extensive-margin phenomenon. This is a
# pure description of the estimation panel (cohorts with >= 5 treated units):
# no estimation, no randomness, no seed.
# ==============================================================================

message("\n=== Section 27: zero shares by donor type ===\n")

zero_groups <- list(
  list(label = "Adopters, pre-adoption years",
       rows  = did_panel_h1$cohort_year > 0 & did_panel_h1$year < did_panel_h1$cohort_year),
  list(label = "Adopters, post-adoption years",
       rows  = did_panel_h1$cohort_year > 0 & did_panel_h1$year >= did_panel_h1$cohort_year),
  list(label = "Never-adopters, all years",
       rows  = did_panel_h1$cohort_year == 0),
  # Time-matched benchmark: the post-adoption row above covers 2021-2024 only,
  # so the comparable never-adopter figure is the same calendar window.
  list(label = "Never-adopters, 2021--2024 only",
       rows  = did_panel_h1$cohort_year == 0 & did_panel_h1$year >= 2021)
)
zero_vars <- c("DAC bilateral" = "commits_dac",
               "Multilateral"  = "commits_multi",
               "Other donors"  = "commits_other",
               "Any donor"     = "commitments")

missing_zero_vars <- setdiff(zero_vars, names(did_panel_h1))
if (length(missing_zero_vars) > 0L) {
  stop("Columns missing from the panel (", paste(missing_zero_vars, collapse = ", "),
       "): zero-share table not written.")
} else {
  zero_mat <- matrix("---", nrow = length(zero_groups), ncol = length(zero_vars) + 2L)
  colnames(zero_mat) <- c(names(zero_vars), "Recipient-years", "Recipients")
  for (i in seq_along(zero_groups)) {
    sub <- did_panel_h1[zero_groups[[i]]$rows, , drop = FALSE]
    for (j in seq_along(zero_vars)) {
      v <- sub[[zero_vars[[j]]]]
      v <- v[!is.na(v)]
      zero_mat[i, j] <- if (length(v) == 0L) "---" else
        sprintf("%.1f", 100 * mean(v <= 0))
    }
    zero_mat[i, length(zero_vars) + 1L] <- format(nrow(sub), big.mark = ",")
    zero_mat[i, length(zero_vars) + 2L] <- as.character(n_distinct(sub$country_id))
    message(sprintf("  %-32s DAC %5s%%  Multi %5s%%  Other %5s%%  Any %5s%%  (N = %s)",
                    zero_groups[[i]]$label, zero_mat[i, 1L], zero_mat[i, 2L],
                    zero_mat[i, 3L], zero_mat[i, 4L],
                    zero_mat[i, length(zero_vars) + 1L]))
  }

  zero_tab <- cbind(Sample = vapply(zero_groups, `[[`, character(1L), "label"),
                    as.data.frame(zero_mat, stringsAsFactors = FALSE),
                    stringsAsFactors = FALSE)
  rownames(zero_tab) <- NULL
  names(zero_tab)[seq_along(zero_vars) + 1L] <-
    paste0(names(zero_vars), " (\\%)")

  xtab_zero <- xtable(zero_tab, label = "tab:zero_shares_donor")
  align(xtab_zero) <- "llcccccc"
  raw_zero <- capture.output(
    print(xtab_zero, include.rownames = FALSE, booktabs = TRUE,
          sanitize.text.function = identity, size = "\\small", floating = FALSE))

  notes_zero <- paste0(
    "Each cell: \\% of recipient-years where the donor group commits nothing under an ",
    "adaptation Rio marker. Sample: main estimation panel (cohorts $\\geq 5$, 2021--2024; ",
    format(nrow(did_panel_h1), big.mark = ","), " recipient-years, ",
    n_distinct(did_panel_h1$country_id), " recipients). Outcomes are $\\log(1+x)$ ",
    "transforms, so a high zero share reflects the extensive margin. Pre-/post-adoption ",
    "rows are not time-comparable; see text"
  )

  out_path_zero <- file.path(dir_tabs_h1, "zero_shares.tex")
  write_tex_float(
    out_path_zero,
    "Share of recipient-years with zero adaptation-marked commitments, by donor type",
    "tab:zero_shares_donor", raw_zero, notes_zero,
    "OECD CRS (Rio adaptation markers); UNFCCC NAP Central")
}

# ==============================================================================
# SECTION 12. §28 Listwise-deletion losses
#
# att_gt() silently drops recipient-years with a missing outcome or a missing
# control, so every specification in this script estimates on a slightly
# different sample. Those losses are printed to the stage log (the text states
# them): recipient-years lost to a missing outcome, to a missing control (WGI
# government effectiveness or log population), and in total, per specification.
# ==============================================================================

message("\n=== Section 28: listwise-deletion losses ===\n")

#' Count recipient-years and recipients lost to missing outcome/controls.
#'
#' @param panel data frame with country_id and the columns below
#' @param outcome_var name of the outcome column
#' @param label specification label
#' @param controls character vector of control column names
#' @return one-row data frame
listwise_losses <- function(panel, outcome_var, label,
                            controls = c("ge_est", "log_population")) {
  present_controls <- intersect(controls, names(panel))
  miss_y <- is.na(panel[[outcome_var]])
  miss_x <- Reduce(`|`, lapply(present_controls, function(cc) is.na(panel[[cc]])))
  if (is.null(miss_x)) miss_x <- rep(FALSE, nrow(panel))
  keep <- !miss_y & !miss_x
  data.frame(
    Specification   = label,
    `Recipient-years (raw)`  = format(nrow(panel), big.mark = ","),
    `Dropped: outcome`       = format(sum(miss_y), big.mark = ","),
    `Dropped: controls`      = format(sum(!miss_y & miss_x), big.mark = ","),
    `Recipient-years (used)` = format(sum(keep), big.mark = ","),
    `Recipients (raw)`       = as.character(n_distinct(panel$country_id)),
    `Recipients (used)`      = as.character(n_distinct(panel$country_id[keep])),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

# Every split section must have run: a missing one stops the stage rather than
# silently dropping its rows from the count.
split_objs <- c("ldc_groups", "gov_groups", "inc_groups")
if (!all(vapply(split_objs, exists, logical(1L))))
  stop("Listwise-loss count: missing split object(s) ",
       paste(split_objs[!vapply(split_objs, exists, logical(1L))], collapse = ", "))
split_specs <- function(groups, prefix) lapply(groups, function(g) list(
  panel = g$panel %>% filter(!(cohort_year %in% thin_cohorts)),
  y = "log_commits", lbl = paste0(prefix, g$label)))
lw_specs <- c(
  list(
    list(panel = did_panel_h1, y = "log_commits",
         lbl = "Main panel: log(adaptation commitments)"),
    list(panel = did_panel_h1, y = "log_commits_dac",   lbl = "Donor type: DAC bilateral"),
    list(panel = did_panel_h1, y = "log_commits_multi", lbl = "Donor type: multilateral"),
    list(panel = did_panel_h1, y = "log_commits_other", lbl = "Donor type: other donors")
  ),
  split_specs(ldc_groups, "LDC split: "),
  split_specs(gov_groups, "Governance split: "),
  split_specs(inc_groups, "Income split: ")
)

lw_tab <- do.call(rbind, lapply(lw_specs, function(sp) {
  if (!sp$y %in% names(sp$panel)) stop("Listwise-loss count: ", sp$y, " missing for ", sp$lbl)
  listwise_losses(sp$panel, sp$y, sp$lbl)
}))

for (r in seq_len(nrow(lw_tab)))
  message(sprintf("  %-44s raw %7s -> used %7s  (outcome %s, controls %s)",
                  lw_tab$Specification[r], lw_tab$`Recipient-years (raw)`[r],
                  lw_tab$`Recipient-years (used)`[r],
                  lw_tab$`Dropped: outcome`[r], lw_tab$`Dropped: controls`[r]))

n_lost <- sum(as.integer(gsub(",", "", lw_tab$`Recipient-years (raw)`, fixed = TRUE)) -
              as.integer(gsub(",", "", lw_tab$`Recipient-years (used)`, fixed = TRUE)))
message(sprintf("  Recipient-years lost to listwise deletion, summed over the %d specifications: %d",
                nrow(lw_tab), n_lost))

message("\n=== 05_heterogeneity.R: complete ===\n")
