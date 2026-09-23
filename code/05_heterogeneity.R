# ==============================================================================
# 05_heterogeneity.R
# Heterogeneity analyses: donor type, LDC status, governance, income group.
# Paper: Beaucoral, Goujon and Marchand (2026) — §6 heterogeneity (§21-24)
#
# Inputs:
#   data/processed/simple_panel_wgi.csv
#
# Outputs (all under output/):
#   figures/heterogeneity/donor_type/did_donor_type_es.png  (§21)
#   tables/heterogeneity/donor_type/att_combined_wide.tex   (§21)
#   figures/heterogeneity/ldc/did_ldc_es.png               (§22)
#   tables/heterogeneity/ldc/het_ldc_wide.tex              (§22)
#   figures/heterogeneity/governance/did_governance_es.png (§23)
#   tables/heterogeneity/governance/het_gov_wide.tex       (§23)
#   figures/heterogeneity/income_group/did_income_es.png   (§24)
#   tables/heterogeneity/income_group/het_income_wide.tex  (§24)
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
# NOTE: MASS is NOT attached via library() -- MASS::select() would mask
# dplyr::select() used throughout this script. compute_pretrend_test() below
# calls MASS::ginv() by full namespace instead (singular-covariance fallback).

set.seed(20240601)  # global seed — local set.seed(1242) calls follow each estimator

# -----------------------------------------------------------------------
# One fit, one SE: single bootstrap-replication constant, and the
# multiplier bootstrap is now ON for every heterogeneity cell.
#
# What changed: this script used to select the inference method
# with `bstrap <- n_treated >= 40`, so every subgroup table was estimated with
# analytical influence-function standard errors -- precisely in the thin cells
# where analytical SEs are most anti-conservative. The 40-unit rule is a
# convention, not an estimator requirement. All reported subgroup ATTs and SEs
# now come from a multiplier-bootstrap fit clustered by recipient country
# (did clusters the multiplier bootstrap on `idname` by construction).
# Analytical twins are retained ONLY where an influence function is required:
# the pre-trend Wald test and the correlated-sample contrasts of §26.
# -----------------------------------------------------------------------
BITERS <- 999L

# Saved fits written by 03_main_results.R / 04_robustness.R.
FITS_DIR <- here("output", "fits")

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
#' @return a character string for the table note
wpval_reconciliation <- function(pre_egt, df_did, n_clusters,
                                 wpval_did = NA_real_, pval_wald = NA_real_,
                                 labels = NULL, wpval_reason = NA_character_,
                                 ginv_used = FALSE) {
  # Compact pre-trend disclosure for table notes (byte-identical
  # in 03/04/05). The averaging mechanism and the interpretation of each
  # result are stated once in the main text; the note keeps the restriction
  # counts and the per-column verdicts, which the text does not repeat table
  # by table.

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
    paste0(PRETREND_NOTE_AGG(min(as.integer(pre_egt)), max(as.integer(pre_egt)), n_lead),
           if (isTRUE(any(ginv_used)))
             "; for this column the block was singular and a generalized inverse was used"
           else "",
           ". ")
  }

  # Column labels are interpolated into LaTeX prose, so they must be escaped
  # here: "Adaptation share (% of global)" would otherwise comment out the rest
  # of the note.
  esc_lab <- function(x) {
    x <- gsub("\\\\", "", x)
    x <- gsub("([%#&_])", "\\\\\\1", x)
    x
  }
  if (!is.null(labels)) labels <- esc_lab(labels)

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
# (compute_pretrend_test, make_wide_table, raw_var_map — duplicated from
#  03_main_results.R; each stage script is self-contained, no sourced utils file)
# ==============================================================================

# --- Wald pre-trend test ---
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
  W <- tryCatch(as.numeric(t(r_th) %*% solve(r_v_r) %*% r_th),
                error = function(e) NA_real_)
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
  W <- tryCatch(as.numeric(t(r_th) %*% solve(r_v_r) %*% r_th),
                error = function(e) NA_real_)
  list(stat = W, df = nrow(R),
       pval = if (is.na(W)) NA_real_ else pchisq(W, df = nrow(R), lower.tail = FALSE))
}

#' Minimum detectable effect at 80% power and a 5% two-sided level.
#'
#' MDE = (z_{1-alpha/2} + z_{power}) * se = (1.96 + 0.84) * se at the defaults.
#' @param se standard error of the estimate or contrast
#' @param power target power; alpha two-sided significance level
mde_from_se <- function(se, power = 0.80, alpha = 0.05) {
  if (is.na(se)) return(NA_real_)
  (qnorm(1 - alpha / 2) + qnorm(power)) * se
}

# --- raw-variable map ---
raw_var_map <- c(
  log_commits           = "commitments",
  share_adapt           = "share_adapt",
  lcommitments_all      = "commitments_all",
  lcommitments_nonadapt = "commitments_nonadapt",
  ldisbursements        = "disbursements",
  log_commits_dac       = "commits_dac",
  log_commits_multi     = "commits_multi",
  log_commits_other     = "commits_other"
)

# --- Combined wide table (donor-type version) ---
# Used by §21 only. Seed rule (reproduces the published SEs): set.seed(1242) inside before outcomes loop.
make_wide_table <- function(did_panel_in, retain_thin, outcomes,
                             dir_tabs, tex_label, caption_spec) {

  use_dr <- if (retain_thin) "reg" else "dr"

  table_stats <- list()
  # Seed rule (reproduces the published SEs): seed immediately before estimator loop
  set.seed(1242)

  for (oc in outcomes) {
    message(sprintf("  Table stats [%s]: %s", caption_spec, oc$label))

    if (!retain_thin) {
      # Main spec (cohorts >= 5): doubly-robust + multiplier-bootstrap SEs.
      # A separate analytical fit is used exclusively for the pre-trend Wald test
      # so that the pre-trend test is computed from an analytical
      # influence-function covariance rather than from bootstrap SEs. (compute_pretrend_test() builds the full
      # covariance from inf.function$dynamic.inf.func.e.)
      gt_tab_analytical <- tryCatch(
        att_gt(
          yname         = oc$var,
          tname         = "year",
          idname        = "country_id",
          gname         = "cohort_year",
          xformla       = ~ ge_est + log_population,
          data          = did_panel_in,
          est_method    = use_dr,
          bstrap        = FALSE,   # analytical IF — feeds pre-trend test only
          cband         = FALSE,
          control_group = "nevertreated",
          anticipation  = 0,
          base_period   = "universal",
          panel         = TRUE,
          allow_unbalanced_panel = TRUE
        ),
        error = function(e) NULL
      )
      # Bootstrap fit — reported ATT and SE
      set.seed(1242)
      gt_tab <- tryCatch(
        att_gt(
          yname         = oc$var,
          tname         = "year",
          idname        = "country_id",
          gname         = "cohort_year",
          xformla       = ~ ge_est + log_population,
          data          = did_panel_in,
          est_method    = use_dr,
          bstrap        = TRUE,    # multiplier-bootstrap SE (DECISION 1)
          biters        = BITERS,
          cband         = FALSE,
          control_group = "nevertreated",
          anticipation  = 0,
          base_period   = "universal",
          panel         = TRUE,
          allow_unbalanced_panel = TRUE
        ),
        error = function(e) NULL
      )
    } else {
      # Thin-cohort robustness spec: outcome regression + analytical SEs.
      gt_tab_analytical <- tryCatch(
        att_gt(
          yname         = oc$var,
          tname         = "year",
          idname        = "country_id",
          gname         = "cohort_year",
          xformla       = ~ ge_est + log_population,
          data          = did_panel_in,
          est_method    = use_dr,
          bstrap        = FALSE,   # analytical IF for pre-trend test
          cband         = FALSE,
          control_group = "nevertreated",
          anticipation  = 0,
          base_period   = "universal",
          panel         = TRUE,
          allow_unbalanced_panel = TRUE
        ),
        error = function(e) NULL
      )
      gt_tab <- gt_tab_analytical  # same fit used for ATT and pre-trend
    }
    if (is.null(gt_tab)) next

    # Pre-trend Wald test always uses the analytical fit: did returns no
    # analytical variance matrix for a bootstrap fit, and the influence-function
    # covariance must come from the same fit as the coefficients.
    agg_dyn_analytical <- if (!is.null(gt_tab_analytical)) tryCatch(
      aggte(gt_tab_analytical, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf),
      error = function(e) NULL
    ) else NULL
    agg_dyn <- agg_dyn_analytical  # alias used for compute_pretrend_test below

    agg_s <- tryCatch(
      aggte(gt_tab, type = "simple", na.rm = TRUE),
      error = function(e) NULL
    )
    # Analytical simple aggregation: influence function for the §26 contrasts
    #. Consumes no random numbers, so it cannot disturb the bootstrap
    # draws behind the reported SE above.
    agg_s_analytical <- if (!is.null(gt_tab_analytical)) tryCatch(
      aggte(gt_tab_analytical, type = "simple", na.rm = TRUE),
      error = function(e) NULL
    ) else NULL

    att <- if (!is.null(agg_s)) agg_s$overall.att else NA_real_
    se  <- if (!is.null(agg_s)) agg_s$overall.se  else NA_real_
    t_v <- if (!is.na(att) && !is.na(se) && se > 0) att / se else NA_real_
    stars <- if (is.na(t_v)) "" else
      if (abs(t_v) > 2.576) "***" else
      if (abs(t_v) > 1.960) "**"  else
      if (abs(t_v) > 1.645) "*"   else ""

    pt <- if (!is.null(agg_dyn)) compute_pretrend_test(agg_dyn, gt_tab_analytical) else
      list(stat = NA_real_, pval = NA_real_, df = 0L,
           W_did = NA_real_, Wpval_did = NA_real_,
           df_did = NA_integer_)

    rows_oc   <- did_panel_in[!is.na(did_panel_in[[oc$var]]), ]
    n_obs     <- nrow(rows_oc)
    n_country <- length(unique(rows_oc$country_id))

    raw_var  <- raw_var_map[oc$var]
    is_share <- (oc$var == "share_adapt")

    if (!is.na(raw_var) && raw_var %in% names(did_panel_in)) {
      pre_rows <- did_panel_in %>%
        filter(cohort_year > 0, year < cohort_year,
               !is.na(.data[[raw_var]]))
      mean_pre <- mean(pre_rows[[raw_var]], na.rm = TRUE)
    } else {
      mean_pre <- NA_real_
    }

    # The naive back-transform exp(ATT)-1 x pre-treatment mean is not
    # additive across outcomes (it violates the paper's own additivity check),
    # so it is reported only when the underlying ATT is significant at 5%
    # (|t| > 1.960); otherwise the cell is suppressed ("---").
    is_sig5 <- !is.na(t_v) && abs(t_v) > 1.960
    if (!is_share && !is.na(att) && !is.na(mean_pre) && is_sig5) {
      implied_usd <- (exp(att) - 1) * mean_pre
    } else {
      implied_usd <- NA_real_
    }

    mean_pre_fmt <- if (is_share) {
      sprintf("%.2f pp", mean_pre)
    } else if (!is.na(mean_pre)) {
      sprintf("%.1f", mean_pre)
    } else {
      "---"
    }

    implied_fmt <- if (is_share) {
      "---"
    } else if (!is.na(implied_usd)) {
      sprintf("%.1f", implied_usd)
    } else {
      "---"
    }

    table_stats[[oc$var]] <- list(
      att_fmt      = paste0(sprintf("%.4f", att), stars),
      se_fmt       = sprintf("(%.4f)", se),
      t_fmt        = sprintf("%.3f", t_v),
      mean_pre_fmt = mean_pre_fmt,
      implied_fmt  = implied_fmt,
      n_obs        = format(n_obs, big.mark = ","),
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
      pt_label     = oc$label,
      pt_wpval_did = if (is.na(pt$Wpval_did)) "---" else sprintf("%.3f", pt$Wpval_did),
      # Unrounded values and the unit-level influence function of the
      # ANALYTICAL twin, used by §26 to contrast donor-type ATTs estimated on
      # the same recipient-years. Adding these fields leaves the published
      # table untouched (cells are selected by name below).
      att_num      = att,
      se_num       = se,
      inf_func     = if (!is.null(agg_s_analytical))
        agg_s_analytical$inf.function$simple.att else NULL,
      se_analytic  = if (!is.null(agg_s_analytical))
        agg_s_analytical$overall.se else NA_real_,
      att_analytic = if (!is.null(agg_s_analytical))
        agg_s_analytical$overall.att else NA_real_,
      ids          = sort(unique(did_panel_in$country_id))
    )
    message(sprintf("    ATT = %s  |  Pre-trend chi2(%d) = %.3f  p = %.3f  |  did Wpval = %s",
                    table_stats[[oc$var]]$att_fmt, pt$df, pt$stat, pt$pval,
                    table_stats[[oc$var]]$pt_wpval_did))
  }

  row_labels <- c(
    "ATT",
    "SE",
    "$t$-statistic",
    "Mean (pre-treat., USD M)",
    "Implied effect (USD M, exp(ATT)$-$1 $\\times$ pre-treat.\\ mean)",
    "Observations",
    "Countries",
    "Pre-trend $\\chi^2$",
    "Pre-trend $p$",
    "\\texttt{did} pre-test $p$",
    "\\midrule Estimator",
    "Control group",
    "Bootstrap SE"
  )

  col_short <- sapply(outcomes, `[[`, "label")
  tab_wide  <- data.frame(` ` = row_labels, check.names = FALSE, stringsAsFactors = FALSE)

  # Row label reflects the actual inference method: bootstrap for main spec
  # (retain_thin = FALSE), analytical for thin-cohort robustness spec.
  # DECISION 1: main spec uses multiplier-bootstrap; thin-cohort spec uses analytical.
  bootstrap_row_val <- if (!retain_thin)
    paste0("Yes (multiplier, ", BITERS, " reps)") else "No (analytical SE)"

  for (i in seq_along(outcomes)) {
    oc  <- outcomes[[i]]
    s   <- table_stats[[oc$var]]
    col <- if (is.null(s)) rep("---", length(row_labels)) else c(
      s$att_fmt, s$se_fmt, s$t_fmt,
      s$mean_pre_fmt, s$implied_fmt,
      s$n_obs, s$n_country,
      s$pt_stat, s$pt_pval, s$pt_wpval_did,
      "CS (2021)", "Never-treated", bootstrap_row_val
    )
    tab_wide[[esc_header(col_short[i])]] <- col
  }

  xtab <- xtable(tab_wide, label = tex_label)
  align(xtab) <- paste0("ll", paste(rep("c", length(col_short)), collapse = ""))

  raw_lines <- capture.output(
    print(xtab, include.rownames = FALSE, booktabs = TRUE,
          sanitize.text.function = identity, size = "\\small",
          floating = FALSE)
  )

  # Degrees of freedom behind the two pre-trend columns, taken from the first
  # outcome that produced a test (they are identical across outcomes here, as
  # all five are estimated on the same panel).
  first_stats <- Filter(Negate(is.null), table_stats)
  first_pt_df       <- if (length(first_stats) == 0L) NA_integer_ else first_stats[[1L]]$pt_df
  first_pt_df_did   <- if (length(first_stats) == 0L) NA_integer_ else first_stats[[1L]]$pt_df_did
  first_n_country   <- if (length(first_stats) == 0L) NA_integer_ else
    as.integer(first_stats[[1L]]$n_country)
  # Per-column inputs for the pre-trend note: the lead window actually used,
  # both p-values for every outcome, whether a generalized inverse was needed,
  # and did's own reason when it returned no statistic. Nothing is hardcoded.
  first_pt_leads    <- if (length(first_stats) == 0L) integer(0) else
    first_stats[[1L]]$pt_leads
  col_wpval <- vapply(first_stats, function(s)
    if (is.null(s$pt_wpval_num)) NA_real_ else s$pt_wpval_num, numeric(1L))
  col_pwald <- vapply(first_stats, function(s)
    if (is.null(s$pt_pval_num)) NA_real_ else s$pt_pval_num, numeric(1L))
  col_labels_pt <- vapply(first_stats, function(s) s$pt_label, character(1L))
  col_ginv  <- vapply(first_stats, function(s) isTRUE(s$pt_ginv), logical(1L))
  first_reason <- {
    rs <- unlist(lapply(first_stats, function(s) s$pt_reason))
    rs <- rs[!is.na(rs)]
    if (length(rs) == 0L) NA_character_ else rs[[1L]]
  }

  est_method_label <- if (use_dr == "dr") "DR" else "OR"
  # DECISION 1: main spec (retain_thin=FALSE) uses multiplier-bootstrap SEs;
  # thin-cohort robustness spec (retain_thin=TRUE) uses analytical SEs.
  se_label <- if (!retain_thin)
    paste0("multiplier-bootstrap SE (", BITERS, " reps), seed 1242") else
    "analytical (IF) SE"

  notes_txt <- paste0(
    "CS\\,(2021) ", est_method_label,
    ", never-treated controls; ", se_label, ". ",
    wpval_reconciliation(pre_egt = first_pt_leads, df_did = first_pt_df_did,
                         n_clusters = first_n_country, wpval_did = col_wpval,
                         pval_wald = col_pwald, labels = col_labels_pt,
                         wpval_reason = first_reason, ginv_used = col_ginv),
    "Implied effects: back-transform on pre-treatment mean; suppressed if insig.\\ at 5\\%. ",
    "* $p<0.10$, ** $p<0.05$, *** $p<0.01$"
  )
  source_txt <- "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"

  cap_title <- paste0(
    "Effect of NAP adoption on climate finance: simple ATT across outcomes (",
    caption_spec, ")"
  )

  out_path <- file.path(dir_tabs, "att_combined_wide.tex")
  write_tex_float(out_path, cap_title, tex_label, raw_lines, notes_txt, source_txt)

  invisible(table_stats)
}

# ==============================================================================
# SECTION 4. Helper: make_het_wide_table()
# Shared across §22-24.
# FIX §2a: derive estimator/inference labels from actual em_g/bstrap.
# FIX §2b: every table note defines significance stars.
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

  col_stats  <- list()
  em_used    <- list()
  bs_used    <- list()

  # Seed rule (reproduces the published SEs): seed before het loop
  set.seed(1242)

  for (grp in groups) {
    lbl <- grp$label
    message(sprintf("  [het_wide] %s: %s", lbl, outcome_var))

    panel_g <- grp$panel %>% filter(!(cohort_year %in% thin_cohorts_vec))

    if (n_distinct(panel_g$cohort_year[panel_g$cohort_year > 0]) < 2) {
      message("  Too few treated cohorts — skipping"); col_stats[[lbl]] <- NULL; next
    }

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
    # influence function for the pre-trend Wald test and for the §26 contrasts.
    gt_g_analytical <- tryCatch(run_gt_g(FALSE),
      error = function(e) { message("  att_gt (analytical) failed: ",
                                     conditionMessage(e)); NULL })

    # Reported ATT and SE: multiplier bootstrap, clustered by recipient.
    set.seed(1242)   # Seed rule (reproduces the published SEs): seed immediately before the estimator
    gt_g <- tryCatch(run_gt_g(TRUE),
      error = function(e) { message("  att_gt (bootstrap) failed: ",
                                     conditionMessage(e)); NULL })
    if (is.null(gt_g)) { col_stats[[lbl]] <- NULL; next }

    # Pre-trend test uses the analytical fit (as in 03/04): did returns no
    # analytical variance matrix for a bootstrap fit.
    agg_d <- if (!is.null(gt_g_analytical))
      tryCatch(aggte(gt_g_analytical, type = "dynamic", na.rm = TRUE,
                     min_e = -5, max_e = Inf), error = function(e) NULL) else NULL
    agg_s_analytical <- if (!is.null(gt_g_analytical))
      tryCatch(aggte(gt_g_analytical, type = "simple", na.rm = TRUE),
               error = function(e) NULL) else NULL
    agg_s <- tryCatch(aggte(gt_g, type = "simple",  na.rm = TRUE), error = function(e) NULL)

    att  <- if (!is.null(agg_s)) agg_s$overall.att else NA_real_
    se   <- if (!is.null(agg_s)) agg_s$overall.se  else NA_real_
    t_v  <- if (!is.na(att) && !is.na(se) && se > 0) att / se else NA_real_
    stars <- if (is.na(t_v)) "" else
      if (abs(t_v) > 2.576) "***" else if (abs(t_v) > 1.960) "**" else
      if (abs(t_v) > 1.645) "*"   else ""

    # gt_g_analytical (not gt_g) is passed so that did's own pre-test and the
    # influence-function covariance are never bootstrap-contaminated -- the
    # convention used at every other call site in this project.
    # Empty pre-trend result, used when the aggregation or the test is
    # unavailable; declared once so the two branches cannot drift apart.
    pt_empty <- list(stat = NA_real_, pval = NA_real_, df = 0L, n_leads = 0L,
                     leads = integer(0), ginv_used = FALSE, W_did = NA_real_,
                     Wpval_did = NA_real_, df_did = NA_integer_,
                     wpval_reason = NA_character_)
    pt <- if (!is.null(agg_d))
      tryCatch(compute_pretrend_test(agg_d, gt_g_analytical),
               error = function(e) pt_empty)
    else pt_empty

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
      # Unrounded values for §25 difference tests. The wide table below selects
      # its cells by name, so adding these fields leaves every published
      # heterogeneity table numerically and textually unchanged.
      att_num      = att,
      se_num       = se,
      n_treated    = n_treated_g,
      # Analytical influence function of the same subgroup fit, kept for
      # correlated-sample contrasts; se_analytic lets the note report how far
      # the bootstrap SE sits from the analytical one it replaced.
      inf_func     = if (!is.null(agg_s_analytical))
        agg_s_analytical$inf.function$simple.att else NULL,
      se_analytic  = if (!is.null(agg_s_analytical))
        agg_s_analytical$overall.se else NA_real_,
      att_analytic = if (!is.null(agg_s_analytical))
        agg_s_analytical$overall.att else NA_real_,
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

  if (length(col_stats) == 0) {
    message("  No results — table not saved."); return(invisible(NULL))
  }

  # FIX §2a: derive estimator label from actual em_g used per group.
  unique_em <- unique(unlist(em_used))
  estimator_label <- if (length(unique_em) == 1L) {
    if (unique_em == "dr") "CS (2021) DR" else "CS (2021) regression adjustment"
  } else {
    "CS (2021) DR or regression adjustment (by subgroup N)"
  }
  # Derived from what was actually used, not hard-coded (the
  # n >= 40 switch was removed, so this is TRUE for every cell; the check keeps the
  # label honest if that ever changes again).
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
    col <- if (is.null(s)) rep("---", length(row_labels)) else c(
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

  # FIX §2b: append star definition to notes
  notes_txt <- paste0(
    caption_txt,
    " ATT/SE: multiplier-bootstrap (", BITERS,
    " reps, seed 1242); pre-trend: separate analytical fit. ",
    wpval_reconciliation(pre_egt = het_leads, df_did = het_df_did,
                         n_clusters = het_nclust, wpval_did = het_wpval,
                         pval_wald = het_pwald, labels = het_labels,
                         wpval_reason = het_reason, ginv_used = het_ginv),
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

dir_figs_h1 <- here("output", "figures", "heterogeneity", "donor_type")
dir_tabs_h1 <- here("output", "tables",  "heterogeneity", "donor_type")
dir.create(dir_figs_h1, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tabs_h1, recursive = TRUE, showWarnings = FALSE)

outcomes_donor <- list(
  list(var = "log_commits_dac",   label = "DAC bilateral",  color = "#2166ac"),
  list(var = "log_commits_multi", label = "Multilateral",   color = "#d6604d"),
  list(var = "log_commits_other", label = "Other donors",   color = "#4dac26")
)

did_panel_h1 <- did_panel_full %>% filter(!(cohort_year %in% thin_cohorts))

# --- Event-study figure ---
h1_dyn <- list()
# Seed set immediately before the estimator (reproduces the published SE)
set.seed(1242)
for (oc in outcomes_donor) {
  if (!oc$var %in% names(did_panel_h1)) { message("  Skipping ", oc$var); next }
  gt_h1 <- tryCatch(
    att_gt(yname = oc$var, tname = "year", idname = "country_id", gname = "cohort_year",
           xformla = ~ ge_est + log_population, data = did_panel_h1,
           est_method = "dr", bstrap = TRUE, biters = BITERS, cband = FALSE,
           control_group = "nevertreated", anticipation = 0,
           base_period = "universal", panel = TRUE, allow_unbalanced_panel = TRUE),
    error = function(e) NULL)
  if (is.null(gt_h1)) next
  agg_d <- tryCatch(aggte(gt_h1, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf),
                    error = function(e) NULL)
  if (is.null(agg_d)) next
  cv <- agg_d$crit.val.egt
  h1_dyn[[oc$var]] <- data.frame(
    outcome    = oc$label, color = oc$color,
    event_time = agg_d$egt, ATT = agg_d$att.egt, SE = agg_d$se.egt,
    Lower      = agg_d$att.egt - cv * agg_d$se.egt,
    Upper      = agg_d$att.egt + cv * agg_d$se.egt,
    stringsAsFactors = FALSE)
}
if (length(h1_dyn) > 0) {
  h1_df <- bind_rows(h1_dyn) %>%
    filter(!is.na(SE), SE > 1e-10) %>%
    mutate(outcome = factor(outcome, levels = sapply(outcomes_donor, `[[`, "label")))
  palette_h1 <- setNames(sapply(outcomes_donor, `[[`, "color"),
                          sapply(outcomes_donor, `[[`, "label"))
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
  ggsave(file.path(dir_figs_h1, "did_donor_type_es.png"), p_h1, width = 10, height = 5, dpi = 300)
  message("Saved: ", file.path(dir_figs_h1, "did_donor_type_es.png"))
}

# --- Wide table: uses make_wide_table (not het version) ---
# The return value is captured so that §26 can contrast the three donor
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

dir_figs_h2 <- here("output", "figures", "heterogeneity", "ldc")
dir_tabs_h2 <- here("output", "tables",  "heterogeneity", "ldc")
dir.create(dir_figs_h2, recursive = TRUE, showWarnings = FALSE)
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

# Previous (2024) vintage, kept for the disclosure below only.
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

# Event-study figure
h2_dyn <- list()
# Seed set immediately before the estimator (reproduces the published SE)
set.seed(1242)
for (grp in ldc_groups) {
  panel_g <- grp$panel %>% filter(!(cohort_year %in% thin_cohorts))
  if (n_distinct(panel_g$cohort_year[panel_g$cohort_year > 0]) < 2) next
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
  if (is.null(gt_h2)) next
  agg_d2 <- tryCatch(aggte(gt_h2, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf),
                     error = function(e) NULL)
  if (is.null(agg_d2)) next
  cv2    <- agg_d2$crit.val.egt
  color2 <- if (grp$label == "LDC") "#e08214" else "#542788"
  h2_dyn[[grp$label]] <- data.frame(
    group = grp$label, color = color2,
    event_time = agg_d2$egt, ATT = agg_d2$att.egt, SE = agg_d2$se.egt,
    Lower = agg_d2$att.egt - cv2 * agg_d2$se.egt,
    Upper = agg_d2$att.egt + cv2 * agg_d2$se.egt,
    stringsAsFactors = FALSE)
}
if (length(h2_dyn) > 0) {
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
  ggsave(file.path(dir_figs_h2, "did_ldc_es.png"), p_h2, width = 10, height = 5, dpi = 300)
  message("Saved: ", file.path(dir_figs_h2, "did_ldc_es.png"))
}

message("\n=== Section 22 complete ===\n")

# ==============================================================================
# SECTION 7. §23 Governance heterogeneity (high vs. low WGI GE)
# Seed rule (reproduces the published SEs): set.seed(1242) before gov loop.
# ==============================================================================

message("\n=== Section 23: Heterogeneity — High vs. low governance ===\n")

dir_figs_h3 <- here("output", "figures", "heterogeneity", "governance")
dir_tabs_h3 <- here("output", "tables",  "heterogeneity", "governance")
dir.create(dir_figs_h3, recursive = TRUE, showWarnings = FALSE)
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

# Event-study figure
h3_dyn <- list()
# Seed set immediately before the estimator (reproduces the published SE)
set.seed(1242)
for (grp in gov_groups) {
  panel_g <- grp$panel %>% filter(!(cohort_year %in% thin_cohorts))
  if (n_distinct(panel_g$cohort_year[panel_g$cohort_year > 0]) < 2) next
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
  if (is.null(gt_h3)) next
  agg_d3 <- tryCatch(aggte(gt_h3, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf),
                     error = function(e) NULL)
  if (is.null(agg_d3)) next
  cv3    <- agg_d3$crit.val.egt
  color3 <- if (grp$label == "High governance") "#1b7837" else "#762a83"
  h3_dyn[[grp$label]] <- data.frame(
    group = grp$label, color = color3,
    event_time = agg_d3$egt, ATT = agg_d3$att.egt, SE = agg_d3$se.egt,
    Lower = agg_d3$att.egt - cv3 * agg_d3$se.egt,
    Upper = agg_d3$att.egt + cv3 * agg_d3$se.egt,
    stringsAsFactors = FALSE)
}
if (length(h3_dyn) > 0) {
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
  ggsave(file.path(dir_figs_h3, "did_governance_es.png"), p_h3, width = 10, height = 5, dpi = 300)
  message("Saved: ", file.path(dir_figs_h3, "did_governance_es.png"))
}

message("\n=== Section 23 complete ===\n")

# ==============================================================================
# SECTION 8. §24 Income group heterogeneity
# Seed rule (reproduces the published SEs): set.seed(1242) before income loop.
# ==============================================================================

message("\n=== Section 24: Heterogeneity — Income group ===\n")

het_stats_income <- NULL   # filled only if the income split is estimable

dir_figs_h4 <- here("output", "figures", "heterogeneity", "income_group")
dir_tabs_h4 <- here("output", "tables",  "heterogeneity", "income_group")
dir.create(dir_figs_h4, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tabs_h4, recursive = TRUE, showWarnings = FALSE)

# PRE-TREATMENT INCOME VINTAGE (FY2013 = July 2012 classification).
# Income classification joined on ISO3 (recipient_iso <-> iso3c); a name-based
# join would drop 18 countries whose CRS names differ from WDI names (e.g.,
# "China (People's Republic of)", "Democratic Republic of the Congo", "Egypt",
# "Yemen", "Turkiye"). The VINTAGE also matters. The classification
# shipped with the installed WDI package is a CURRENT cross-section: a recipient
# that moved from lower-middle to upper-middle income during the estimation
# window was being assigned its post-treatment group, so the split conditioned
# partly on an outcome. We now use the World Bank's own historical file
# (OGHIST.xlsx, "Country Analytical History" sheet) and read the FY2013 column
# -- the classification announced in July 2012, i.e. strictly before the first
# NAP cohort in the estimation sample (2021).
# The file is downloaded once to data/raw/oghist/ and cached; if neither the
# cache nor the download is available the script falls back to the WDI snapshot
# and says so in the table note (income_vintage below is the single source of
# truth for that wording).
OGHIST_URL  <- paste0("https://datacatalogfiles.worldbank.org/ddh-published/",
                      "0037712/DR0090754/OGHIST.xlsx")
oghist_path <- here("data", "raw", "oghist", "OGHIST.xlsx")
dir.create(dirname(oghist_path), recursive = TRUE, showWarnings = FALSE)

# SHA-256 of the vintage behind the published exhibits, recorded in
# DATA_AVAILABILITY.md. A silently different OGHIST would change the income
# split without changing anything visible, so the checksum is verified on every
# run: a mismatch is a warning with both digests, not a stop, because the World
# Bank does reissue the file and the authors must decide whether to adopt the new
# vintage and update DATA_AVAILABILITY.md.
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
            "  expected (DATA_AVAILABILITY.md): ", OGHIST_SHA256, "\n",
            "  found on disk:                   ", oghist_sha, "\n",
            "  The World Bank has reissued the workbook. Verify the FY13 ",
            "column, then update DATA_AVAILABILITY.md and OGHIST_SHA256 here.")
  } else {
    message("  OGHIST.xlsx SHA-256 verified against DATA_AVAILABILITY.md.")
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
#' unreadable OGHIST file now stops the script with instructions instead.
#' The file is cached under data/raw/oghist/ and documented in
#' DATA_AVAILABILITY.md.
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
         "  See DATA_AVAILABILITY.md.")
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
       vintage = paste0("World Bank FY2013 (July 2012) analytical ",
                        "classification, OGHIST.xlsx"))
}

income_res    <- read_income_classification(oghist_path)
income_lookup <- income_res$lookup
income_vintage <- income_res$vintage
message("  Income classification vintage: ", income_vintage)
if (!is.null(income_lookup)) print(table(income_lookup$income_group))

# Disclosure: which panel recipients change income group between the two
# vintages. Reported in the console and summarised in the table note.
income_switchers <- character(0)
if (!is.null(income_lookup) && grepl("^World Bank FY2013", income_vintage %||% "")) {
  wdi_now <- tryCatch({
    cty <- WDI::WDI_data$country
    inc_col <- names(cty)[grepl("income", names(cty), ignore.case = TRUE)][1]
    iso_col <- names(cty)[grepl("^iso3c$", names(cty), ignore.case = TRUE)][1]
    cty %>% select(all_of(c(iso_col, inc_col))) %>%
      rename(recipient_iso = 1, income_now = 2) %>% distinct()
  }, error = function(e) NULL)
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
  message("  Income-group coverage (all 144 panel countries):")
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

    # Event-study figure
    h4_dyn <- list()
    # Seed set immediately before the estimator (reproduces the published SE)
    set.seed(1242)
    for (grp in inc_groups) {
      panel_g <- grp$panel %>% filter(!(cohort_year %in% thin_cohorts))
      if (n_distinct(panel_g$cohort_year[panel_g$cohort_year > 0]) < 2L) next
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
      if (is.null(gt_h4)) next
      agg_d4 <- tryCatch(aggte(gt_h4, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf),
                         error = function(e) NULL)
      if (is.null(agg_d4)) next
      cv4 <- agg_d4$crit.val.egt
      h4_dyn[[grp$label]] <- data.frame(
        group = grp$label, color = inc_colors[grp$label],
        event_time = agg_d4$egt, ATT = agg_d4$att.egt, SE = agg_d4$se.egt,
        Lower = agg_d4$att.egt - cv4 * agg_d4$se.egt,
        Upper = agg_d4$att.egt + cv4 * agg_d4$se.egt,
        stringsAsFactors = FALSE)
    }
    if (length(h4_dyn) > 0) {
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
      ggsave(file.path(dir_figs_h4, "did_income_es.png"), p_h4, width = 10, height = 5, dpi = 300)
      message("Saved: ", file.path(dir_figs_h4, "did_income_es.png"))
    }
  } else {
    message("  Too few income groups with sufficient cohorts — skipping figure")
  }
} else {
  message("  Income group data unavailable — skipping Section 24")
}

message("\n=== Section 24 complete ===\n")

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
  message("  Donor-type fits or their influence functions are unavailable — ",
          "H3 contrasts skipped.")
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
path_adapt <- file.path(FITS_DIR, "headline_adaptation_dr_bs.rds")
path_mit   <- file.path(FITS_DIR, "mitigation_dr_bs.rds")
mit_contrast <- NULL

if (!file.exists(path_adapt) || !file.exists(path_mit)) {
  message("  Saved adaptation and/or mitigation fit not found (",
          basename(path_adapt), ", ", basename(path_mit), ") — ",
          "mitigation contrast skipped. Run 03 then 04 before 05.")
} else {
  fit_adapt <- readRDS(path_adapt)
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
    message("  One of the saved fits carries no simple-aggregation influence ",
            "function — mitigation contrast skipped.")
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
# subgroup has >= 40 treated units, analytical influence-function SE,
# not-yet-treated controls, cohorts with < 5 treated units dropped,
# xformla ~ ge_est + log_population, outcome log_commits).
#
# Inference: subgroups are disjoint sets of countries and each subgroup's
# influence function is clustered by country, so the two subgroup estimators are
# asymptotically independent and se(diff) = sqrt(se_a^2 + se_b^2). The joint
# income test uses the same argument through a diagonal covariance matrix.
# No new seed is set here: nothing in this section is stochastic.
# ==============================================================================

message("\n=== Section 25: Formal subgroup difference tests (H4) ===\n")

# Filled inside the else-branch below; used by §26 (MDE table). Assigned at
# script top level, so the branch writes into the global environment.
h4_records <- NULL

dir_tabs_h5 <- here("output", "tables", "heterogeneity")

have_ldc <- !is.null(het_stats_ldc) &&
  all(c("LDC", "Non-LDC") %in% names(het_stats_ldc))
have_gov <- !is.null(het_stats_gov) &&
  all(c("Low governance", "High governance") %in% names(het_stats_gov))
income_cells <- c("Low income", "Lower middle income", "Upper middle income")
have_inc <- !is.null(het_stats_income) &&
  all(income_cells %in% names(het_stats_income))

if (!have_ldc || !have_gov || !have_inc) {
  message("  Missing subgroup fits (LDC=", have_ldc, " gov=", have_gov,
          " income=", have_inc, ") — difference-test table not written.")
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

  # Append the correlated-sample contrasts computed in §25a. The four H4
  # rows and the joint income row above are built exactly as before; only new
  # rows are added.
  if (!is.null(extra_rows)) {
    names(extra_rows) <- names(diff_tab)
    diff_tab <- rbind(diff_tab, extra_rows)
  }

  # Records for the MDE table in §26.
  h4_records <- lapply(seq_along(contrast_spec), function(i) list(
    split = contrast_spec[[i]]$split, label = contrast_spec[[i]]$label,
    res = c(contrast_res[[i]], list(method = "independence (disjoint subsamples)",
                                    n_units = NA_integer_)),
    family = "H4"))

  for (r in seq_len(nrow(diff_tab))) {
    message(sprintf("  %-36s diff = %8s  SE = %7s  z = %7s  p = %s",
                    diff_tab$Contrast[r], diff_tab$Difference[r],
                    diff_tab$SE[r], diff_tab$z[r], diff_tab$p[r]))
  }

  names(diff_tab) <- c("Split", "Contrast", "Difference", "SE", "$z$", "$p$-value")

  xtab_h5 <- xtable(diff_tab, label = "tab:het_difftests")
  align(xtab_h5) <- "lllcccc"

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
    "Each row tests equality of two subgroup ATTs from the paired heterogeneity table ",
    "(not re-estimated). ", est_label_h5, " on log(adaptation commitments), headline ",
    "specification. Differences = low- minus high-capacity (positive: larger effect where ",
    "capacity weaker). H4 rows (LDC/governance/income): SE $=\\sqrt{se_a^2+se_b^2}$ on ",
    "multiplier-bootstrap SEs (", BITERS, " reps, seed 1242), disjoint subsamples. ",
    "H3/mitigation rows: aligned analytical-IF SEs (not bootstrap), same recipient-years. ",
    "$z$ vs.\\ standard normal; joint stats vs.\\ $\\chi^2$; $p$ two-sided"
  )
  source_h5 <- paste0("OECD CRS (Rio adaptation markers); UNFCCC NAP Central; ",
                      "WGI; World Bank income classification")

  out_path_h5 <- file.path(dir_tabs_h5, "het_difference_tests.tex")
  write_tex_float(
    out_path_h5,
    paste0("Formal tests of subgroup differences in the NAP effect ",
           "(test of H4: is the effect larger where baseline capacity is weaker?)"),
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
# SECTION 10. §26 Minimum detectable differences for the heterogeneity
# family.
#
# The heterogeneity discussion opens with what the design
# could have detected rather than with what it failed to reject. For every
# contrast in Table~\ref{tab:het_difftests} we report
#     MDE = (z_{0.975} + z_{0.80}) * se(contrast) = 2.8016 * se,
# the smallest true difference a two-sided 5% test would reject with 80%
# probability. Nothing is estimated here; the standard errors are the ones
# already reported.
# ==============================================================================

message("\n=== Section 26: minimum detectable differences ===\n")

mde_records <- c(if (is.null(h4_records)) list() else h4_records,
                 if (length(extra_records) == 0L) list() else extra_records)

if (length(mde_records) == 0L) {
  message("  No contrasts available — MDE table not written.")
} else {
  mde_z <- qnorm(1 - 0.05 / 2) + qnorm(0.80)
  message(sprintf("  MDE multiplier (80%% power, 5%% two-sided) = %.4f", mde_z))

  mde_tab <- data.frame(
    Split    = vapply(mde_records, `[[`, character(1L), "split"),
    Contrast = vapply(mde_records, `[[`, character(1L), "label"),
    Estimate = vapply(mde_records, function(r)
      if (is.na(r$res$diff)) "---" else sprintf("%.4f", r$res$diff), character(1L)),
    SE       = vapply(mde_records, function(r)
      if (is.na(r$res$se)) "---" else sprintf("%.4f", r$res$se), character(1L)),
    MDE      = vapply(mde_records, function(r)
      if (is.na(r$res$se)) "---" else sprintf("%.4f", mde_from_se(r$res$se)),
      character(1L)),
    Powered  = vapply(mde_records, function(r) {
      if (is.na(r$res$se) || is.na(r$res$diff)) return("---")
      if (abs(r$res$diff) >= mde_from_se(r$res$se)) "Yes" else "No"
    }, character(1L)),
    Inference = vapply(mde_records, function(r) {
      if (identical(r$family, "H4")) "Independent subsamples" else
        "Aligned influence functions"
    }, character(1L)),
    stringsAsFactors = FALSE
  )
  for (r in seq_len(nrow(mde_tab)))
    message(sprintf("  %-46s diff = %8s  SE = %7s  MDE = %7s  powered: %s",
                    mde_tab$Contrast[r], mde_tab$Estimate[r], mde_tab$SE[r],
                    mde_tab$MDE[r], mde_tab$Powered[r]))

  names(mde_tab) <- c("Split", "Contrast", "Estimated difference", "SE",
                      "MDE (80\\% power)", "$|$Diff$| \\geq$ MDE", "Inference")

  xtab_mde <- xtable(mde_tab, label = "tab:het_mde")
  align(xtab_mde) <- "lllccccc"
  raw_mde <- capture.output(
    print(xtab_mde, include.rownames = FALSE, booktabs = TRUE,
          sanitize.text.function = identity, size = "\\small", floating = FALSE))

  notes_mde <- paste0(
    "MDE $= (z_{0.975}+z_{0.80}) \\times se = ", sprintf("%.4f", mde_z),
    " \\times se$, for the contrasts of Table~\\ref{tab:het_difftests} (not re-estimated). ",
    "``No'' means the null is uninformative about smaller differences, not equality. Last ",
    "column: independence (disjoint subsamples) or aligned unit-level IFs (same recipient-",
    "years). Outcome: log(adaptation commitments); donor-type rows use the donor-group ",
    "commitments"
  )

  out_path_mde <- file.path(dir_tabs_h5, "het_mde.tex")
  write_tex_float(
    out_path_mde,
    paste0("Minimum detectable differences for the heterogeneity and ",
           "falsification contrasts"),
    "tab:het_mde", raw_mde, notes_mde,
    paste0("OECD CRS (Rio adaptation markers); UNFCCC NAP Central; WGI; ",
           "World Bank income classification"))

  if (dir.exists(here("paper"))) {
    paper_path_mde <- here("paper", "Tables", "heterogeneity", "het_mde.tex")
    dir.create(dirname(paper_path_mde), recursive = TRUE, showWarnings = FALSE)
    file.copy(out_path_mde, paper_path_mde, overwrite = TRUE)
    message("Saved: ", paper_path_mde)
  } else message("paper/ not found -- exhibits are left in output/ only")
}

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
  message("  Columns missing from the panel (", paste(missing_zero_vars, collapse = ", "),
          ") — zero-share table not written.")
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
# different sample. Those losses are stated here. Counts are reported
# per specification: recipient-years and recipients lost to a missing outcome,
# to a missing control (WGI government effectiveness or log population), and in
# total.
# ==============================================================================

message("\n=== Section 28: listwise-deletion losses ===\n")

#' Count recipient-years and recipients lost to missing outcome/controls.
#'
#' @param panel data frame with country_id and the columns below
#' @param outcome_var name of the outcome column
#' @param label specification label for the table
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

lw_specs <- list(
  list(panel = did_panel_h1, y = "log_commits",
       lbl = "Main panel: log(adaptation commitments)"),
  list(panel = did_panel_h1, y = "log_commits_dac",   lbl = "Donor type: DAC bilateral"),
  list(panel = did_panel_h1, y = "log_commits_multi", lbl = "Donor type: multilateral"),
  list(panel = did_panel_h1, y = "log_commits_other", lbl = "Donor type: other donors")
)
if (exists("ldc_groups"))
  for (g in ldc_groups) lw_specs[[length(lw_specs) + 1L]] <- list(
    panel = g$panel %>% filter(!(cohort_year %in% thin_cohorts)),
    y = "log_commits", lbl = paste0("LDC split: ", g$label))
if (exists("gov_groups"))
  for (g in gov_groups) lw_specs[[length(lw_specs) + 1L]] <- list(
    panel = g$panel %>% filter(!(cohort_year %in% thin_cohorts)),
    y = "log_commits", lbl = paste0("Governance split: ", g$label))
if (exists("inc_groups"))
  for (g in inc_groups) lw_specs[[length(lw_specs) + 1L]] <- list(
    panel = g$panel %>% filter(!(cohort_year %in% thin_cohorts)),
    y = "log_commits", lbl = paste0("Income split: ", g$label))

lw_tab <- do.call(rbind, lapply(lw_specs, function(sp) {
  if (!sp$y %in% names(sp$panel)) return(NULL)
  listwise_losses(sp$panel, sp$y, sp$lbl)
}))

if (is.null(lw_tab)) {
  message("  No specifications available — listwise-loss table not written.")
} else {
  for (r in seq_len(nrow(lw_tab)))
    message(sprintf("  %-44s raw %7s -> used %7s  (outcome %s, controls %s)",
                    lw_tab$Specification[r], lw_tab$`Recipient-years (raw)`[r],
                    lw_tab$`Recipient-years (used)`[r],
                    lw_tab$`Dropped: outcome`[r], lw_tab$`Dropped: controls`[r]))

  xtab_lw <- xtable(lw_tab, label = "tab:listwise_losses")
  align(xtab_lw) <- "llcccccc"
  raw_lw <- capture.output(
    print(xtab_lw, include.rownames = FALSE, booktabs = TRUE,
          sanitize.text.function = identity, size = "\\small", floating = FALSE))

  notes_lw <- paste0(
    "Recipient-years/recipients lost to listwise deletion, by specification (outcome ",
    "vs.\\ controls: WGI government effectiveness, log population); \\texttt{did} drops ",
    "these cells silently, so ``used'' is the estimation sample (see text). Panels are ",
    "post cohort-$\\geq 5$ restriction, and post subgroup restriction for split rows"
  )

  write_tex_float(
    file.path(dir_tabs_h5, "listwise_losses.tex"),
    "Listwise-deletion losses by specification",
    "tab:listwise_losses", raw_lw, notes_lw,
    "OECD CRS (Rio adaptation markers); UNFCCC NAP Central; WGI; World Bank WDI")
}

message("\n=== 05_heterogeneity.R: complete ===\n")
