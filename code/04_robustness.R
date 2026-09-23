# ==============================================================================
# 04_robustness.R
# Robustness checks: retained-cohorts spec, HonestDiD, placebo, mitigation.
# Paper: Beaucoral, Goujon and Marchand (2026) — §5 robustness (thin-cohort sensitivity, HonestDiD,
#        placebo test, mitigation falsification)
#
# Inputs:
#   data/processed/simple_panel_wgi.csv
#   data/processed/mitigation_panel.csv   (raw recipient-year flows; built by 01)
#
# Outputs (all under output/):
#   figures/cohorts_retained/did_combined_cohort_wgi.png   (§19b robustness)
#   tables/cohorts_retained/att_combined_wide.tex          (§19c robustness)
#   figures/notyettreated/did_combined_cohort_wgi.png      (§19d not-yet-treated)
#   figures/notyettreated/did_notyettreated_es.png         (§19d event study)
#   tables/notyettreated/att_notyettreated_wide.tex        (§19d not-yet-treated)
#   tables/units_zeros/diagnostic_units_zeros.tex          (§19e units/zeros sensitivity)
#   tables/cohorts_dropped/honestdid_rm.tex                (§20 HonestDiD)
#   figures/placebo/did_placebo_es.png                     (§25)
#   tables/placebo/att_placebo.tex                         (§25)
#   figures/mitigation/did_mitigation_es.png               (§26)
#   tables/mitigation/att_mitigation.tex                   (§26)
#   tables/dcdh/att_dcdh.tex                                (§A dCDH estimator)
#   figures/dcdh/did_dcdh_es.png                            (§A dCDH event study)
#   tables/bacon/bacon_decomp.tex                           (§B Goodman-Bacon)
#   figures/bacon/bacon_scatter.png                         (§B Goodman-Bacon)
#
# NOTE: honestdid_combined.png is NOT written (not used in paper).
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
# $\bar{M}$           | Mbarvec_hd              | HonestDiD RM parameter
# ============================================================

# ARM-mac headless gotcha (verified): DIDmultiplegtDYN pulls in rgl, which hangs
# indefinitely on a headless run unless RGL is told to use the null device.  This
# MUST be the first executable line, before ANY library() call.
Sys.setenv(RGL_USE_NULL = TRUE)

# duplicated from 03/05 §1 — keep in sync
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
library(fixest)
library(countrycode)
# NOTE: MASS is NOT attached via library() -- MASS::select() would mask
# dplyr::select() used throughout this script. compute_pretrend_test() below
# calls MASS::ginv() by full namespace instead (singular-covariance fallback).
# §A/§B appendix estimators (de Chaisemartin & D'Haultfoeuille; Goodman-Bacon).
# polars MUST be loaded before DIDmultiplegtDYN, else did_multiplegt_dyn() errors
# with "objet 'pl' introuvable" (the package references the loaded `pl` object).
library(polars)
library(DIDmultiplegtDYN)
library(bacondecomp)

# HonestDiD (GitHub-only optional dependency; used in §20). Guarded with
# requireNamespace() so a missing install fails with an actionable message
# rather than a bare "there is no package called 'HonestDiD'" error. Moved
# here from inside §20 so that all packages are loaded at the top of the
# script, before any data loading or computation.
if (!requireNamespace("HonestDiD", quietly = TRUE)) {
  stop(
    "Package 'HonestDiD' is not installed.\n",
    "Run renv::restore() (renv.lock pins HonestDiD 0.2.8); ",
    "see README.md, Computational requirements."
  )
}
library(HonestDiD)

set.seed(20240601)  # global seed — local set.seed(1242) calls follow each estimator

# -----------------------------------------------------------------------
# One fit, one SE: single bootstrap-replication constant.
# Every multiplier-bootstrap call in this script uses BITERS. Before this
# revision the §25 placebo fit passed no `biters` at all and silently used
# did's default of 1000 while its table note claimed 999 replications.
# 999 is kept so that no published SE is re-randomised
# merely by changing the replication count.
# -----------------------------------------------------------------------
BITERS <- 999L

# -----------------------------------------------------------------------
# Saved headline fits written by 03_main_results.R. §19e and §20 read
# these instead of re-fitting, so the paper can never print two different
# standard errors for one specification.
# -----------------------------------------------------------------------
FITS_DIR <- here("output", "fits")

# canonical pre-trend wording — keep byte-identical across scripts
PRETREND_NOTE_AGG <- function(min_e, max_e, k) sprintf(
  "Pre-trend $\\chi^2$ / $p$: joint Wald test on the aggregated pre-treatment event-time coefficients ($%d \\leq e \\leq %d$; %d restrictions), using the influence-function covariance of the dynamic aggregation from the analytical (non-bootstrap) fit; a generalized inverse is used if the block is singular",
  min_e, max_e, k)
PRETREND_NOTE_DID <- "\\texttt{did} pre-test $p$: \\texttt{did}'s built-in Wald test over all pre-period $ATT(g,t)$ cells against each cohort's $g-1$ base year"
PRETREND_NOTE <- function(min_e, max_e, k) paste0(PRETREND_NOTE_AGG(min_e, max_e, k), ". ", PRETREND_NOTE_DID)

#' Read a saved headline fit written by 03_main_results.R.
#'
#' @param stem one of "adaptation", "share", "total", "nonadaptation",
#'   "disbursements"
#' @return the saved list (see 03_main_results.R §7b for its fields)
read_headline_fit <- function(stem) {
  path <- file.path(FITS_DIR, paste0("headline_", stem, "_dr_bs.rds"))
  if (!file.exists(path)) {
    stop("Saved headline fit not found: ", path, "\n",
         "  04_robustness.R consumes the fits written by 03_main_results.R ",
         "(one fit, one SE).\n",
         "  Fix: run `Rscript code/03_main_results.R` first ",
         "(run_all.R already orders 03 before 04).")
  }
  readRDS(path)
}

# Outcome -> saved-fit stem (kept in sync with 03_main_results.R fit_stem_map).
fit_stem_map <- c(
  log_commits           = "adaptation",
  share_adapt           = "share",
  lcommitments_all      = "total",
  lcommitments_nonadapt = "nonadaptation",
  ldisbursements        = "disbursements"
)

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
# Null-coalescing operator (used in HonestDiD helpers)
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
# Wraps a bare tabular block in the project's standard complete float:
#   \begin{table}[H]
#   \centering
#   \caption{<title>}
#   \label{<label>}
#   \adjustbox{max width=\textwidth}{%
#     <tabular>
#   }
#   \par\vspace{4pt}
#   \begin{minipage}{\linewidth}\footnotesize
#   Notes: <notes>.\par
#   Source: <source>.
#   \end{minipage}
#   \end{table}
# Called by every table-generating block in this script.
# ==============================================================================

write_tex_float <- function(out_path, caption_title, label,
                             tabular_lines, notes_text, source_text,
                             size = "\\small") {
  # Strip any existing \begin{table}/\end{table} wrapper that xtable emits
  # when floating = TRUE.  We rebuild it ourselves.
  inner <- tabular_lines
  # Remove leading/trailing table float lines if present
  is_table_open  <- grepl("^\\\\begin\\{table\\}", inner)
  is_table_close <- grepl("^\\\\end\\{table\\}", inner)
  if (any(is_table_open)) inner <- inner[!is_table_open]
  if (any(is_table_close)) inner <- inner[!is_table_close]
  # Remove centering / caption / label lines emitted by xtable
  inner <- inner[!grepl("^\\\\centering", inner)]
  inner <- inner[!grepl("^\\\\caption", inner)]
  inner <- inner[!grepl("^\\\\label", inner)]
  # Remove blank lines at top/bottom
  while (length(inner) > 0 && trimws(inner[1]) == "") inner <- inner[-1]
  while (length(inner) > 0 && trimws(inner[length(inner)]) == "")
    inner <- inner[-length(inner)]

  # Locate tabular block to wrap in adjustbox
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

  # Build the complete float
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

message("\n=== 04_robustness.R: loading panel ===\n")

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

# --- Combined wide table ---
# Seed rule (reproduces the published SEs): set.seed(1242) inside before outcomes loop
make_wide_table <- function(did_panel_in, retain_thin, outcomes,
                             dir_tabs, tex_label, caption_spec,
                             control_grp  = "nevertreated",
                             out_filename = "att_combined_wide.tex") {

  use_dr <- if (retain_thin) "reg" else "dr"

  # Control-group labels (parameterized so the same helper builds the
  # never-treated main table and the not-yet-treated robustness table).
  control_grp_label <- if (control_grp == "notyettreated")
    "Not-yet-treated" else "Never-treated"
  control_grp_note  <- if (control_grp == "notyettreated")
    "not-yet-treated control group" else "never-treated control group"

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
          control_group = control_grp,
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
          control_group = control_grp,
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
          control_group = control_grp,
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
      pt_wpval_did = if (is.na(pt$Wpval_did)) "---" else sprintf("%.3f", pt$Wpval_did)
    )
    message(sprintf("    ATT = %s  |  Pre-trend chi2(%d) = %.3f  p = %.3f  |  did Wpval = %s",
                    table_stats[[oc$var]]$att_fmt, pt$df, pt$stat, pt$pval,
                    table_stats[[oc$var]]$pt_wpval_did))
  }

  # Row label reflects the actual inference method: bootstrap for main spec
  # (retain_thin = FALSE), analytical for thin-cohort robustness spec.
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

  # DECISION 1: main spec uses multiplier-bootstrap; thin-cohort spec uses analytical.
  bootstrap_row_val <- if (!retain_thin) "Yes (multiplier, 999 reps)" else "No (analytical SE)"

  for (i in seq_along(outcomes)) {
    oc  <- outcomes[[i]]
    s   <- table_stats[[oc$var]]
    col <- if (is.null(s)) rep("---", length(row_labels)) else c(
      s$att_fmt, s$se_fmt, s$t_fmt,
      s$mean_pre_fmt, s$implied_fmt,
      s$n_obs, s$n_country,
      s$pt_stat, s$pt_pval, s$pt_wpval_did,
      "CS (2021)", control_grp_label, bootstrap_row_val
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
  se_label <- if (!retain_thin) "multiplier-bootstrap SE (999 reps, seed 1242; \\texttt{did} 2.5.0; CRS Apr.\\ 2026)" else "analytical (IF) SE; \\texttt{did} 2.5.0; CRS Apr.\\ 2026"

  notes_txt <- paste0(
    "CS\\,(2021) ", est_method_label,
    "; WGI gov.\\ effectiveness + log population; ",
    control_grp_note, "; ", se_label, ". ",
    wpval_reconciliation(pre_egt = first_pt_leads, df_did = first_pt_df_did,
                         n_clusters = first_n_country, wpval_did = col_wpval,
                         pval_wald = col_pwald, labels = col_labels_pt,
                         wpval_reason = first_reason, ginv_used = col_ginv),
    "Implied effects: naive back-transforms on the pre-treatment mean, not additive; ",
    "suppressed if insignificant at 5\\%. ",
    "* $p<0.10$, ** $p<0.05$, *** $p<0.01$"
  )
  source_txt <- "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"

  cap_title <- paste0(
    "Effect of NAP adoption on climate finance: simple ATT across outcomes (",
    caption_spec, ")"
  )

  out_path <- file.path(dir_tabs, out_filename)
  write_tex_float(out_path, cap_title, tex_label, raw_lines, notes_txt, source_txt)

  invisible(table_stats)
}

# ==============================================================================
# SECTION 4. Helper: make_cohort_plot()
# (duplicated from 03_main_results.R §9)
# ==============================================================================

make_cohort_plot <- function(results_group, outcomes, dir_figs,
                              spec_label, use_bstrap) {

  group_all <- bind_rows(results_group) %>%
    mutate(
      outcome = factor(outcome, levels = sapply(outcomes, `[[`, "label")),
      cohort  = as.integer(cohort)
    ) %>%
    filter(!is.na(ATT))

  if (nrow(group_all) == 0) {
    message("No cohort ATT data for: ", spec_label); return(invisible(NULL))
  }

  palette_vec <- setNames(sapply(outcomes, `[[`, "color"),
                           sapply(outcomes, `[[`, "label"))
  dodge_w <- 0.6

  p <- ggplot(group_all,
              aes(x = factor(cohort), y = ATT,
                  colour = outcome, shape = outcome, group = outcome)) +
    geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed",
               linewidth = 0.5) +
    geom_linerange(aes(ymin = Lower, ymax = Upper),
                   position = position_dodge(width = dodge_w),
                   linewidth = 0.6, alpha = 0.85) +
    geom_point(size = 2.5,
               position = position_dodge(width = dodge_w)) +
    scale_colour_manual(values = palette_vec) +
    scale_shape_manual(values  = c(16, 17, 15, 18, 8)) +
    # No title, subtitle, or caption — those go in LaTeX \caption{}
    labs(
      title    = NULL, subtitle = NULL, caption = NULL,
      x        = "Adoption cohort (year of NAP)",
      y        = "Cohort ATT estimate",
      colour   = NULL, shape = NULL
    ) +
    theme_minimal() +
    theme(
      text             = element_text(family = "serif", size = 12),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      legend.position  = "bottom",
      legend.text      = element_text(size = 10),
      panel.grid.minor = element_blank()
    )

  out_path <- file.path(dir_figs, "did_combined_cohort_wgi.png")
  ggsave(out_path, p, width = 12, height = 7, dpi = 300)
  message("Saved: ", out_path)
  invisible(p)
}

# ==============================================================================
# SECTION 5. §19 retained-cohorts robustness spec
# Seed rule (reproduces the published SEs): set.seed(1242) before outcomes loop.
# ==============================================================================

message("\n=== Section 19: retained-cohorts robustness ===\n")

dir.create(here("output", "figures", "cohorts_retained"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "tables",  "cohorts_retained"), recursive = TRUE, showWarnings = FALSE)

# Storage for retained-cohorts results
results_simple_r  <- list()
results_group_r   <- list()
results_dynamic_r <- list()

# Seed rule (reproduces the published SEs): set.seed(1242) before loop
set.seed(1242)

for (oc in outcomes) {
  message(sprintf("\n--- Retained spec: %s ---", oc$label))

  gt_r <- tryCatch(
    att_gt(
      yname         = oc$var,
      tname         = "year",
      idname        = "country_id",
      gname         = "cohort_year",
      xformla       = ~ ge_est + log_population,
      data          = did_panel_full,
      est_method    = "reg",
      bstrap        = FALSE,
      cband         = FALSE,
      control_group = "nevertreated",
      anticipation  = 0,
      base_period   = "universal",
      panel         = TRUE,
      allow_unbalanced_panel = TRUE
    ),
    error = function(e) {
      message("  att_gt failed: ", conditionMessage(e)); NULL
    }
  )
  if (is.null(gt_r)) next

  agg_s_r <- tryCatch(aggte(gt_r, type = "simple", na.rm = TRUE), error = function(e) NULL)
  agg_g_r <- tryCatch(aggte(gt_r, type = "group",  na.rm = TRUE), error = function(e) NULL)
  agg_d_r <- tryCatch(
    aggte(gt_r, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf),
    error = function(e) NULL
  )

  if (!is.null(agg_s_r)) {
    results_simple_r[[oc$var]] <- data.frame(
      outcome = oc$label,
      ATT     = round(agg_s_r$overall.att, 4),
      SE      = round(agg_s_r$overall.se,  4),
      t_stat  = round(agg_s_r$overall.att / agg_s_r$overall.se, 3)
    )
  }
  if (!is.null(agg_g_r)) {
    results_group_r[[oc$var]] <- data.frame(
      outcome = oc$label,
      cohort  = agg_g_r$egt,
      ATT     = round(agg_g_r$att.egt,  4),
      SE      = round(agg_g_r$se.egt,   4),
      Lower   = round(agg_g_r$att.egt - agg_g_r$crit.val.egt * agg_g_r$se.egt, 4),
      Upper   = round(agg_g_r$att.egt + agg_g_r$crit.val.egt * agg_g_r$se.egt, 4)
    )
  }
  if (!is.null(agg_d_r)) {
    results_dynamic_r[[oc$var]] <- data.frame(
      outcome    = oc$label,
      event_time = agg_d_r$egt,
      ATT        = round(agg_d_r$att.egt, 4),
      SE         = round(agg_d_r$se.egt,  4),
      Lower      = round(agg_d_r$att.egt - agg_d_r$crit.val.egt * agg_d_r$se.egt, 4),
      Upper      = round(agg_d_r$att.egt + agg_d_r$crit.val.egt * agg_d_r$se.egt, 4),
      color      = oc$color
    )
  }
}

# --- §19b cohort plot: cohorts_retained ---
make_cohort_plot(
  results_group = results_group_r,
  outcomes      = outcomes,
  dir_figs      = file.path(here("output", "figures"), "cohorts_retained"),
  spec_label    = "All cohorts retained (asymptotic SE)",
  use_bstrap    = FALSE
)

# --- §19c wide table: cohorts_retained ---
make_wide_table(
  did_panel_in  = did_panel_full,
  retain_thin   = TRUE,
  outcomes      = outcomes,
  dir_tabs      = file.path(here("output", "tables"), "cohorts_retained"),
  tex_label     = "tab:combined_wide_robust",
  caption_spec  = "all cohorts retained, asymptotic SE"
)

# ==============================================================================
# SECTION 5b. §19d Not-yet-treated control-group robustness
# Mirrors the MAIN specification (03_main_results.R) exactly — doubly-robust
# estimation, multiplier-bootstrap SE (999 reps), WGI gov. effectiveness +
# log population controls, universal base period, no anticipation, main panel
# (cohorts >= 5) — changing ONLY the comparison group from never-treated to
# not-yet-treated.  Under "notyettreated", a unit that adopts a NAP in a later
# year serves as a control for earlier-adopting cohorts up to its own adoption,
# and never-treated units remain in the control pool.  A stable ATT across the
# two control groups shows the estimate is not an artifact of the never-treated
# comparison set.
# Seed rule (reproduces the published SEs): set.seed(1242) before the outcomes loop and again before each
# bootstrap fit (matches make_wide_table / main-spec seeding for reproducibility).
# ==============================================================================

message("\n=== Section 19d: not-yet-treated control-group robustness ===\n")

dir_figs_nyt <- file.path(here("output", "figures"), "notyettreated")
dir_tabs_nyt <- file.path(here("output", "tables"),  "notyettreated")
dir.create(dir_figs_nyt, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tabs_nyt, recursive = TRUE, showWarnings = FALSE)

# Hold the estimation panel fixed at the MAIN-spec definition (cohorts_dropped):
# thin cohorts (< thin_threshold treated units) are removed, exactly as in
# 03_main_results.R run_did_estimation(retain_thin_cohorts = FALSE).  Only the
# control group differs from the headline result.
did_panel_nyt <- if (length(thin_cohorts) > 0)
  did_panel_full %>% filter(!(cohort_year %in% thin_cohorts)) else did_panel_full
message(sprintf("Not-yet-treated panel: dropped thin cohorts {%s}; N = %d rows, %d countries",
                paste(thin_cohorts, collapse = ", "),
                nrow(did_panel_nyt), length(unique(did_panel_nyt$country_id))))

results_group_nyt   <- list()
results_dynamic_nyt <- list()

# Seed rule (reproduces the published SEs): set.seed(1242) before loop
set.seed(1242)

for (oc in outcomes) {
  message(sprintf("\n--- Not-yet-treated spec: %s ---", oc$label))

  set.seed(1242)  # re-seed before each bootstrap fit (DECISION 1 reproducibility)
  gt_nyt <- tryCatch(
    att_gt(
      yname         = oc$var,
      tname         = "year",
      idname        = "country_id",
      gname         = "cohort_year",
      xformla       = ~ ge_est + log_population,
      data          = did_panel_nyt,
      est_method    = "dr",
      bstrap        = TRUE,
      biters        = BITERS,
      cband         = FALSE,
      control_group = "notyettreated",
      anticipation  = 0,
      base_period   = "universal",
      panel         = TRUE,
      allow_unbalanced_panel = TRUE
    ),
    error = function(e) { message("  att_gt failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(gt_nyt)) next

  agg_g_nyt <- tryCatch(aggte(gt_nyt, type = "group",  na.rm = TRUE), error = function(e) NULL)
  agg_d_nyt <- tryCatch(
    aggte(gt_nyt, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf),
    error = function(e) NULL
  )

  if (!is.null(agg_g_nyt)) {
    results_group_nyt[[oc$var]] <- data.frame(
      outcome = oc$label,
      cohort  = agg_g_nyt$egt,
      ATT     = round(agg_g_nyt$att.egt, 4),
      SE      = round(agg_g_nyt$se.egt,  4),
      Lower   = round(agg_g_nyt$att.egt - agg_g_nyt$crit.val.egt * agg_g_nyt$se.egt, 4),
      Upper   = round(agg_g_nyt$att.egt + agg_g_nyt$crit.val.egt * agg_g_nyt$se.egt, 4)
    )
  }
  if (!is.null(agg_d_nyt)) {
    results_dynamic_nyt[[oc$var]] <- data.frame(
      outcome    = oc$label,
      event_time = agg_d_nyt$egt,
      ATT        = round(agg_d_nyt$att.egt, 4),
      SE         = round(agg_d_nyt$se.egt,  4),
      Lower      = round(agg_d_nyt$att.egt - agg_d_nyt$crit.val.egt * agg_d_nyt$se.egt, 4),
      Upper      = round(agg_d_nyt$att.egt + agg_d_nyt$crit.val.egt * agg_d_nyt$se.egt, 4),
      color      = oc$color
    )
  }
}

# --- §19d-i cohort plot (all outcomes): notyettreated ---
make_cohort_plot(
  results_group = results_group_nyt,
  outcomes      = outcomes,
  dir_figs      = dir_figs_nyt,
  spec_label    = "Not-yet-treated controls (multiplier bootstrap)",
  use_bstrap    = TRUE
)

# --- §19d-ii headline event-study: log(adaptation commitments) ---
# No in-figure title/subtitle/caption — those live in LaTeX \caption{}.
es_nyt <- results_dynamic_nyt[["log_commits"]]
if (!is.null(es_nyt)) {
  p_es_nyt <- ggplot(es_nyt, aes(x = event_time, y = ATT)) +
    geom_hline(yintercept = 0,    colour = "grey50", linetype = "dashed") +
    geom_vline(xintercept = -0.5, colour = "grey30", linetype = "dotted") +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, fill = "#2E86C1") +
    geom_line(colour = "#2E86C1", linewidth = 0.8) +
    geom_point(colour = "#2E86C1", size = 2) +
    scale_x_continuous(breaks = sort(unique(es_nyt$event_time))) +
    labs(title = NULL, subtitle = NULL, caption = NULL,
         x = "Years relative to NAP adoption",
         y = "ATT on log(adaptation commitments)") +
    theme_minimal() +
    theme(text = element_text(family = "serif", size = 12),
          panel.grid.minor = element_blank())
  es_path <- file.path(dir_figs_nyt, "did_notyettreated_es.png")
  ggsave(es_path, p_es_nyt, width = 9, height = 5.5, dpi = 300)
  message("Saved: ", es_path)
}

# --- §19d-iii wide ATT table (all outcomes), mirroring main spec ---
make_wide_table(
  did_panel_in  = did_panel_nyt,
  retain_thin   = FALSE,
  outcomes      = outcomes,
  dir_tabs      = dir_tabs_nyt,
  tex_label     = "tab:combined_wide_notyettreated",
  caption_spec  = "not-yet-treated controls, multiplier-bootstrap SE",
  control_grp   = "notyettreated",
  out_filename  = "att_notyettreated_wide.tex"
)

# ==============================================================================
# SECTION 5c. §19e Outcome-units and zeros sensitivity
# One diagnostic table for the HEADLINE outcome only (log adaptation
# commitments).  It reports the simple ATT under five outcome-construction
# choices to show that the log1p transform does not manufacture the headline
# effect through its treatment of zeros / sub-million observations.
# All specs share the MAIN-spec estimation settings: never-treated control,
# WGI gov. effectiveness + log population controls, universal base period, no
# anticipation, and the main panel (thin cohorts < thin_threshold dropped).
# Five outcome definitions (commitments are in USD millions; raw USD = x 1e6):
#   1. log1p (USD millions)  -- log_commits as-is        (dr  + bootstrap 999)
#   2. log1p (raw USD)       -- log1p(commitments * 1e6) (dr  + bootstrap 999)
#   3. asinh (USD millions)  -- asinh(commitments)       (dr  + bootstrap 999)
#   4. asinh (raw USD)       -- asinh(commitments * 1e6) (dr  + bootstrap 999)
#   5. Pure log, zeros dropped (unit-invariant) -- filter commitments > 0,
#      outcome log(commitments)                          (reg + bootstrap 999)
#
# GOTCHA (do not "optimise" away): att_gt(est_method = "dr") SEGFAULTS inside
# fastglm's colMax_dense on the positives-only panel of spec 5.  A segfault
# kills the R process and CANNOT be caught by tryCatch, so spec 5 MUST use
# est_method = "reg".  Specs 1--4 keep "dr" (matching the headline result).
#
# Seed rule (reproduces the published SEs): set.seed(1242) once before the spec loop, and again
# immediately before each att_gt() call (matches make_wide_table / 5b seeding).
# ==============================================================================

message("\n=== Section 19e: outcome-units and zeros sensitivity ===\n")

dir_tabs_uz <- file.path(here("output", "tables"), "units_zeros")
dir.create(dir_tabs_uz, recursive = TRUE, showWarnings = FALSE)

# Estimation panel = MAIN-spec definition (thin cohorts dropped), identical to
# the not-yet-treated section's did_panel_nyt.  Guard for the no-thin case.
did_panel_main <- if (length(thin_cohorts) > 0)
  did_panel_full %>% filter(!(cohort_year %in% thin_cohorts)) else did_panel_full
message(sprintf("Units/zeros panel: dropped thin cohorts {%s}; N = %d rows, %d countries",
                paste(thin_cohorts, collapse = ", "),
                nrow(did_panel_main), length(unique(did_panel_main$country_id))))

# Spec definitions.  `transform` maps the panel to (data, yname) ready for
# att_gt(); `est_method` is "dr" for specs 1--4 and "reg" for spec 5 (the
# positives-only panel that segfaults under "dr").
uz_specs <- list(
  list(key   = "log1p_millions",
       label = "log1p (USD millions)",
       est_method = "dr",
       transform  = function(d) {
         d$uz_y <- log1p(d$commitments)
         list(data = d, yname = "uz_y")
       }),
  list(key   = "log1p_rawusd",
       label = "log1p (raw USD)",
       est_method = "dr",
       transform  = function(d) {
         d$uz_y <- log1p(d$commitments * 1e6)
         list(data = d, yname = "uz_y")
       }),
  list(key   = "asinh_millions",
       label = "asinh (USD millions)",
       est_method = "dr",
       transform  = function(d) {
         d$uz_y <- asinh(d$commitments)
         list(data = d, yname = "uz_y")
       }),
  list(key   = "asinh_rawusd",
       label = "asinh (raw USD)",
       est_method = "dr",
       transform  = function(d) {
         d$uz_y <- asinh(d$commitments * 1e6)
         list(data = d, yname = "uz_y")
       }),
  list(key   = "purelog_dropzeros",
       label = "Pure log, zeros dropped",
       est_method = "reg",   # GOTCHA: "dr" segfaults on positives-only panel
       transform  = function(d) {
         d <- d[!is.na(d$commitments) & d$commitments > 0, ]
         d$uz_y <- log(d$commitments)
         list(data = d, yname = "uz_y")
       })
)

uz_stats <- list()

# Seed rule (reproduces the published SEs): seed once before the spec loop.
set.seed(1242)

for (sp in uz_specs) {
  message(sprintf("\n--- Units/zeros spec: %s (est_method = %s) ---",
                  sp$label, sp$est_method))

  tr      <- sp$transform(did_panel_main)
  d_sp    <- tr$data
  yn      <- tr$yname

  # N = rows with a non-missing outcome actually used in this spec.
  n_obs_sp <- sum(!is.na(d_sp[[yn]]))

  # --------------------------------------------------------------------
  # The log1p(USD millions) row IS the headline specification:
  # log_commits = log1p(commitments), same panel, same estimator, same
  # controls, same control group, same seed. Re-fitting it here used to
  # produce a SECOND standard error for one specification (0.1173 in this
  # table vs 0.1245 in Table 2). The two fits were never different: aggte()
  # re-runs the multiplier bootstrap when the att_gt object was fitted with
  # bstrap = TRUE, so the aggregated SE depends on the RNG state at the
  # aggte() call, which differed between 03 (post-att_gt-mboot state) and
  # this script (fresh set.seed(1242) state). Seeding before att_gt() does
  # not pin the aggregated SE; sharing the fit does. This row is therefore
  # READ from output/fits/headline_adaptation_dr_bs.rds.
  # --------------------------------------------------------------------
  if (identical(sp$key, "log1p_millions")) {
    fit_head <- read_headline_fit("adaptation")
    fit_head_uz <- fit_head          # kept for the note built after the loop
    stopifnot(
      identical(fit_head$outcome, "log_commits"),
      identical(fit_head$spec$est_method, sp$est_method),
      identical(fit_head$spec$control_group, "nevertreated"),
      identical(as.integer(fit_head$spec$biters), as.integer(BITERS)),
      identical(as.integer(fit_head$n_obs), as.integer(n_obs_sp)),
      isTRUE(all.equal(d_sp[[yn]], d_sp$log_commits))
    )
    att_sp <- fit_head$att
    se_sp  <- fit_head$se
    t_sp   <- att_sp / se_sp
    stars_sp <- if (abs(t_sp) > 2.576) "***" else if (abs(t_sp) > 1.960) "**" else
                if (abs(t_sp) > 1.645) "*" else ""
    pt_sp <- fit_head$pretrend

    uz_stats[[sp$key]] <- list(
      label        = sp$label,
      att_fmt      = paste0(sprintf("%.4f", att_sp), stars_sp),
      se_fmt       = sprintf("(%.4f)", se_sp),
      t_fmt        = sprintf("%.3f", t_sp),
      pt_stat      = if (is.na(pt_sp$stat)) "---" else sprintf("%.3f", pt_sp$stat),
      pt_pval      = if (is.na(pt_sp$pval)) "---" else sprintf("%.3f", pt_sp$pval),
      pt_wpval_did = if (is.na(pt_sp$Wpval_did)) "---" else sprintf("%.3f", pt_sp$Wpval_did),
      n_obs        = format(n_obs_sp, big.mark = ",")
    )
    message(sprintf(paste0("    [read from output/fits/headline_adaptation_dr_bs.rds] ",
                           "ATT = %s  SE = %s  t = %s  |  Pre-trend chi2(%s) = %s  ",
                           "p = %s  |  did Wpval = %s  |  N = %s"),
                    uz_stats[[sp$key]]$att_fmt, uz_stats[[sp$key]]$se_fmt,
                    uz_stats[[sp$key]]$t_fmt, pt_sp$df,
                    uz_stats[[sp$key]]$pt_stat, uz_stats[[sp$key]]$pt_pval,
                    uz_stats[[sp$key]]$pt_wpval_did, uz_stats[[sp$key]]$n_obs))
    next
  }

  # Bootstrap fit — reported ATT and SE (multiplier bootstrap, BITERS reps).
  set.seed(1242)  # re-seed immediately before att_gt (reproducibility invariant)
  gt_sp <- tryCatch(
    att_gt(
      yname         = yn,
      tname         = "year",
      idname        = "country_id",
      gname         = "cohort_year",
      xformla       = ~ ge_est + log_population,
      data          = d_sp,
      est_method    = sp$est_method,
      bstrap        = TRUE,
      biters        = BITERS,
      cband         = FALSE,
      control_group = "nevertreated",
      anticipation  = 0,
      base_period   = "universal",
      panel         = TRUE,
      allow_unbalanced_panel = TRUE
    ),
    error = function(e) { message("  att_gt failed: ", conditionMessage(e)); NULL }
  )

  # Separate analytical fit (bstrap = FALSE) feeds the pre-trend Wald test only,
  # keeping the pre-trend test on an analytical influence-function covariance
  # (did returns no analytical variance matrix for a bootstrap fit).
  set.seed(1242)  # re-seed immediately before att_gt (reproducibility invariant)
  gt_sp_analytical <- tryCatch(
    att_gt(
      yname         = yn,
      tname         = "year",
      idname        = "country_id",
      gname         = "cohort_year",
      xformla       = ~ ge_est + log_population,
      data          = d_sp,
      est_method    = sp$est_method,
      bstrap        = FALSE,
      cband         = FALSE,
      control_group = "nevertreated",
      anticipation  = 0,
      base_period   = "universal",
      panel         = TRUE,
      allow_unbalanced_panel = TRUE
    ),
    error = function(e) NULL
  )

  agg_s_sp <- if (!is.null(gt_sp))
    tryCatch(aggte(gt_sp, type = "simple", na.rm = TRUE), error = function(e) NULL)
  else NULL

  agg_d_sp <- if (!is.null(gt_sp_analytical))
    tryCatch(aggte(gt_sp_analytical, type = "dynamic", na.rm = TRUE,
                   min_e = -5, max_e = Inf), error = function(e) NULL)
  else NULL

  att_sp <- if (!is.null(agg_s_sp)) agg_s_sp$overall.att else NA_real_
  se_sp  <- if (!is.null(agg_s_sp)) agg_s_sp$overall.se  else NA_real_
  t_sp   <- if (!is.na(att_sp) && !is.na(se_sp) && se_sp > 0) att_sp / se_sp else NA_real_
  stars_sp <- if (is.na(t_sp)) "" else
    if (abs(t_sp) > 2.576) "***" else
    if (abs(t_sp) > 1.960) "**"  else
    if (abs(t_sp) > 1.645) "*"   else ""

  pt_sp <- if (!is.null(agg_d_sp))
    tryCatch(compute_pretrend_test(agg_d_sp, gt_sp_analytical),
             error = function(e) list(stat = NA_real_, pval = NA_real_, df = 0L,
                                       W_did = NA_real_, Wpval_did = NA_real_,
                                       df_did = NA_integer_))
  else
    list(stat = NA_real_, pval = NA_real_, df = 0L,
         W_did = NA_real_, Wpval_did = NA_real_,
         df_did = NA_integer_)

  uz_stats[[sp$key]] <- list(
    label        = sp$label,
    att_fmt      = if (is.na(att_sp)) "---" else paste0(sprintf("%.4f", att_sp), stars_sp),
    se_fmt       = if (is.na(se_sp))  "---" else sprintf("(%.4f)", se_sp),
    t_fmt        = if (is.na(t_sp))   "---" else sprintf("%.3f", t_sp),
    pt_stat      = if (is.na(pt_sp$stat)) "---" else sprintf("%.3f", pt_sp$stat),
    pt_pval      = if (is.na(pt_sp$pval)) "---" else sprintf("%.3f", pt_sp$pval),
    pt_wpval_did = if (is.na(pt_sp$Wpval_did)) "---" else sprintf("%.3f", pt_sp$Wpval_did),
    n_obs        = format(n_obs_sp, big.mark = ",")
  )

  message(sprintf("    ATT = %s  SE = %s  t = %s  |  Pre-trend chi2(%s) = %s  p = %s  |  did Wpval = %s  |  N = %s",
                  uz_stats[[sp$key]]$att_fmt, uz_stats[[sp$key]]$se_fmt,
                  uz_stats[[sp$key]]$t_fmt, pt_sp$df,
                  uz_stats[[sp$key]]$pt_stat, uz_stats[[sp$key]]$pt_pval,
                  uz_stats[[sp$key]]$pt_wpval_did, uz_stats[[sp$key]]$n_obs))
}

# -----------------------------------------------------------------------
# SIXTH spec: PPML in LEVELS (Poisson) — the unit-invariant, zero-retaining
# benchmark (Chen & Roth remedy).  Estimated on `commitments` in USD millions
# via the Sun & Abraham (2021) interaction-weighted estimator with a Poisson
# link (fixest::fepois + sunab()), which is robust to staggered timing.
# Poisson in levels handles the 230 zeros natively and its coefficient is a
# proportional (semi-elasticity) effect, invariant to the units of the outcome.
# fepois drops perfectly-separated observations, so N may be < 2,014 — that is
# expected and correct; we report whatever fepois uses (model$nobs).
# Analytical clustered SE (clustered by country) — not bootstrap-based, so no
# seeding is required for this spec.
# -----------------------------------------------------------------------
message("\n--- Units/zeros spec: PPML, levels (Poisson) [Sun & Abraham 2021] ---")

# Sentinel cohort for never-treated units: must sit OUTSIDE the year range so
# sunab() treats them as the pure control group.
ppml_dat <- did_panel_main
ppml_dat$cohort_sa <- ifelse(ppml_dat$cohort_year == 0, 10000, ppml_dat$cohort_year)

ppml_model <- tryCatch(
  fepois(commitments ~ sunab(cohort_sa, year) + ge_est + log_population |
           country_id + year,
         data = ppml_dat, cluster = ~country_id),
  error = function(e) { message("  fepois failed: ", conditionMessage(e)); NULL }
)

# Overall ATT via Sun & Abraham aggregation.
ppml_att_row <- if (!is.null(ppml_model)) {
  tryCatch(summary(ppml_model, agg = "att")$coeftable["ATT", ],
           error = function(e) { message("  agg=att failed: ", conditionMessage(e)); NULL })
} else NULL

att_ppml <- if (!is.null(ppml_att_row)) ppml_att_row[["Estimate"]]   else NA_real_
se_ppml  <- if (!is.null(ppml_att_row)) ppml_att_row[["Std. Error"]] else NA_real_
t_ppml   <- if (!is.null(ppml_att_row)) ppml_att_row[[3L]]           else NA_real_
p_ppml   <- if (!is.null(ppml_att_row)) ppml_att_row[[4L]]           else NA_real_

stars_ppml <- {
  if (is.na(p_ppml)) "" else
  if (p_ppml < 0.01) "***" else
  if (p_ppml < 0.05) "**"  else
  if (p_ppml < 0.10) "*"   else ""
}

# N: observations actually used by fepois (after any separation drops).
n_obs_ppml <- if (!is.null(ppml_model)) as.integer(ppml_model$nobs) else NA_integer_

# Pre-trend: joint Wald test on PRE-treatment event-study coefficients
# (pre-period interaction terms carry negative event-time labels "year::-").
# Non-blocking: NA -> "---" cell if the test is not cleanly available.
#
# NOTE: sunab() exposes only the DISAGGREGATED cohort x period interaction
# terms, so fixest::wald(keep = "year::-") jointly tests ~50 individual
# pre-period cells (50 DoF) rather than the 4-DoF aggregated event-study ATTs
# reported for specs 1--5.  The resulting statistic is on a different scale and
# is NOT comparable to the CS\,(2021) pre-trend column, so we render this cell
# as "---" to avoid a misleading cross-row comparison.  (The test runs clean;
# we deliberately suppress the non-comparable value rather than display it.)
ppml_pt <- list(stat = NA_real_, pval = NA_real_, Wpval_did = NA_real_,
                df_did = NA_integer_, n_leads = 0L, leads = integer(0),
                ginv_used = FALSE, wpval_reason = NA_character_)

uz_stats[["ppml_levels"]] <- list(
  label        = "PPML, levels (Poisson)",
  att_fmt      = if (is.na(att_ppml)) "---" else paste0(sprintf("%.4f", att_ppml), stars_ppml),
  se_fmt       = if (is.na(se_ppml))  "---" else sprintf("(%.4f)", se_ppml),
  t_fmt        = if (is.na(t_ppml))   "---" else sprintf("%.3f", t_ppml),
  pt_stat      = if (is.na(ppml_pt$stat)) "---" else sprintf("%.3f", ppml_pt$stat),
  pt_pval      = if (is.na(ppml_pt$pval)) "---" else sprintf("%.3f", ppml_pt$pval),
  pt_wpval_did = if (is.na(ppml_pt$Wpval_did)) "---" else sprintf("%.3f", ppml_pt$Wpval_did),
  n_obs        = if (is.na(n_obs_ppml)) "---" else format(n_obs_ppml, big.mark = ",")
)

message(sprintf("    ATT = %s  SE = %s  t = %s  p = %s  |  Pre-trend chi2 = %s  p = %s  did Wpval = %s  |  N = %s",
                uz_stats[["ppml_levels"]]$att_fmt, uz_stats[["ppml_levels"]]$se_fmt,
                uz_stats[["ppml_levels"]]$t_fmt,
                if (is.na(p_ppml)) "---" else sprintf("%.4f", p_ppml),
                uz_stats[["ppml_levels"]]$pt_stat, uz_stats[["ppml_levels"]]$pt_pval,
                uz_stats[["ppml_levels"]]$pt_wpval_did,
                uz_stats[["ppml_levels"]]$n_obs))

# Degrees of freedom behind the two pre-trend columns, read from the stored
# headline fit so that the note quotes the same test the first row reports.
# The guard is not decorative: if the log1p-millions spec were ever removed or
# renamed, the note would silently quote whatever `fit_head_uz` happened to be.
stopifnot(exists("fit_head_uz"), !is.null(fit_head_uz$pretrend))
uz_pt_df       <- fit_head_uz$pretrend$df
uz_pt_df_did   <- fit_head_uz$pretrend$df_did
uz_n_country   <- as.integer(fit_head_uz$n_country)
uz_pt_leads    <- fit_head_uz$pretrend$leads
uz_wpval       <- fit_head_uz$pretrend$Wpval_did
uz_pval        <- fit_head_uz$pretrend$pval
uz_reason      <- fit_head_uz$pretrend$wpval_reason
uz_ginv        <- isTRUE(fit_head_uz$pretrend$ginv_used)

# --- §19e diagnostic table: tables/units_zeros/diagnostic_units_zeros.tex ---
uz_order <- c("log1p_millions", "log1p_rawusd", "asinh_millions",
              "asinh_rawusd", "purelog_dropzeros", "ppml_levels")

uz_tab <- data.frame(
  Specification            = vapply(uz_order, function(k) uz_stats[[k]]$label,   character(1L)),
  ATT                      = vapply(uz_order, function(k) uz_stats[[k]]$att_fmt, character(1L)),
  SE                       = vapply(uz_order, function(k) uz_stats[[k]]$se_fmt,  character(1L)),
  `$t$`                    = vapply(uz_order, function(k) uz_stats[[k]]$t_fmt,   character(1L)),
  `Pre-trend $\\chi^2$`    = vapply(uz_order, function(k) uz_stats[[k]]$pt_stat, character(1L)),
  `Pre-trend $p$`          = vapply(uz_order, function(k) uz_stats[[k]]$pt_pval, character(1L)),
  `\\texttt{did} pre-test $p$` = vapply(uz_order, function(k) uz_stats[[k]]$pt_wpval_did, character(1L)),
  N                        = vapply(uz_order, function(k) uz_stats[[k]]$n_obs,   character(1L)),
  check.names      = FALSE,
  stringsAsFactors = FALSE
)

xtab_uz <- xtable(uz_tab, label = "tab:units_zeros")
align(xtab_uz) <- "llccccccc"   # row-name col (dropped) + Specification (l) + 7 numeric cols (c)

raw_uz <- capture.output(
  print(xtab_uz, include.rownames = FALSE, booktabs = TRUE,
        sanitize.text.function = identity, size = "\\small",
        floating = FALSE)
)

notes_uz <- paste0(
  "Outcome: adaptation commitments. Specs 1--4: CS\\,(2021) DR, boot SE (",
  BITERS, " reps, seed 1242); 5: OR, boot SE, positive subsample. Row 1: headline ",
  "fit, not re-estimated. Row 6: Sun--Abraham ",
  "(2021) PPML, levels, country+year FE, clustered SE; retains zeros; no ",
  "comparable pre-trend statistic. Never-treated controls; WGI GE + log ",
  "population. Commitments in USD millions (raw USD $= \\times 10^6$). ",
  wpval_reconciliation(pre_egt = uz_pt_leads, df_did = uz_pt_df_did,
                       n_clusters = uz_n_country, wpval_did = uz_wpval,
                       pval_wald = uz_pval, wpval_reason = uz_reason,
                       ginv_used = uz_ginv),
  "* $p<0.10$, ** $p<0.05$, *** $p<0.01$"
)
source_uz <- "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"

write_tex_float(
  out_path      = file.path(dir_tabs_uz, "diagnostic_units_zeros.tex"),
  caption_title = "Sensitivity of the headline ATT to outcome units and the treatment of zeros",
  label         = "tab:units_zeros",
  tabular_lines = raw_uz,
  notes_text    = notes_uz,
  source_text   = source_uz
)

# ==============================================================================
# SECTION 6. §20 HonestDiD sensitivity analysis
#
# DESIGN:
#   The paper's headline numbers come from a doubly-robust fit with
#   multiplier-bootstrap standard errors, so the sensitivity analysis must
#   describe that same specification (not an analytical "reg" re-fit). §20 runs on the SAVED
#   headline fits (output/fits/headline_*_dr_bs.rds, written by
#   03_main_results.R), i.e. exactly the fits behind Table 2.
#
# HOW THE COVARIANCE IS BUILT (stated explicitly because HonestDiD needs a full
# matrix and `aggte()` returns only a vector of standard errors):
#   1. did's dynamic aggregation stores the unit-level influence-function
#      matrix at agg$inf.function$dynamic.inf.func.e (n x K, aligned with
#      agg$egt; the alignment is asserted).
#   2. did's bootstrap standard error is NOT sqrt(diag(cov(bres))/n): mboot()
#      reports an interquartile-range scale estimate,
#      se = bSigma * sqrt(n_clusters)/n. Its covariance matrix is discarded.
#   3. We re-run did::mboot() on that same influence-function matrix under an
#      explicit seed to recover the bootstrap covariance V (V = cov(bres),
#      bres already scaled by sqrt(n)), take R = cov2cor(V/n), and rescale:
#          Sigma = diag(se.egt) %*% R %*% diag(se.egt).
#      The DIAGONAL of Sigma is then exactly the event-study standard-error
#      vector the paper reports, and the off-diagonal correlations are the
#      multiplier-bootstrap ones. The gap between the covariance-based and the
#      published IQR-based standard errors is logged for every outcome.
#
# l_vec: the paper's headline number is did's SIMPLE ATT, which is the
#   group-size-weighted (NOT equal-weighted) average of the post-treatment
#   event-time effects. 03_main_results.R computes those weights and stores
#   them after verifying that sum(l_vec * att.egt_post) reproduces the reported
#   simple ATT; here we re-verify and use them, so the sensitivity analysis
#   targets the estimand Table 2 prints. Equal weights (0.25 each) are reported
#   in the console for comparison.
#
# THREE EXHIBITS:
#   tables/cohorts_dropped/honestdid_rm.tex         relative magnitudes, same
#                                                   layout as before
#   tables/cohorts_dropped/honestdid_prepriods.tex  breakdown Mbar for
#                                                   numPrePeriods in {4, 6, 8}
#   tables/cohorts_dropped/honestdid_sd.tex         smoothness restriction
#                                                   Delta^SD(M), M in
#                                                   {0, 0.01, 0.02, 0.05}
# Seed rule (reproduces the published SEs): set.seed(1242) immediately before every bootstrap draw.
# ==============================================================================

message("\n=== Section 20: HonestDiD (on the saved DR/bootstrap fits) ===\n")

Mbarvec_hd   <- seq(0, 2, by = 0.25)
# Outcomes for which the 6- and 8-lead ladders are computed (see the loop below).
PRIMARY_OUTCOMES <- c("log_commits", "share_adapt")
Mvec_sd      <- c(0, 0.01, 0.02, 0.05)
numPre_grid  <- c(4L, 6L, 8L)

#' Multiplier-bootstrap covariance of a dynamic aggregation.
#' (Byte-identical copy of the helper in 03_main_results.R §7b — duplicated
#'  per this project's no-sourced-utils convention.)
#'
#' @param agg_d dynamic AGGTEobj from aggte(type = "dynamic")
#' @param gt_obj the att_gt MP object the aggregation came from
#' @param seed integer seed set immediately before the bootstrap draw
#' @return list(sigma, keep, se_published, se_cov, max_rel_dev)
boot_vcov_dynamic <- function(agg_d, gt_obj, seed = 1242L) {
  IF <- agg_d$inf.function$dynamic.inf.func.e
  stopifnot(is.matrix(IF), ncol(IF) == length(agg_d$egt))
  keep <- which(!is.na(agg_d$se.egt) & agg_d$se.egt > 1e-10)
  n    <- nrow(IF)

  set.seed(seed)
  mb <- did::mboot(IF[, keep, drop = FALSE], DIDparams = gt_obj$DIDparams,
                   return_V = TRUE)
  stopifnot(is.matrix(mb$V), ncol(mb$V) == length(keep))

  V_theta <- mb$V / n                       # covariance of the estimator
  se_cov  <- sqrt(diag(V_theta))
  R       <- stats::cov2cor(V_theta)
  se_pub  <- agg_d$se.egt[keep]
  sigma   <- diag(se_pub, nrow = length(se_pub)) %*% R %*%
             diag(se_pub, nrow = length(se_pub))
  # Symmetrise against floating-point asymmetry (HonestDiD checks symmetry).
  sigma <- (sigma + t(sigma)) / 2

  list(sigma        = sigma,
       keep         = keep,
       se_published = se_pub,
       se_cov       = se_cov,
       max_rel_dev  = max(abs(se_cov - se_pub) / se_pub))
}

#' Smallest grid value at which a sensitivity CI first contains zero.
#'
#' @param lb numeric vector of robust lower bounds, ordered as `grid`
#' @param grid numeric vector of M or Mbar values
#' @return formatted LaTeX string
breakdown_value <- function(lb, grid) {
  if (all(is.na(lb))) return("---")
  if (all(lb > 0, na.rm = TRUE)) return(sprintf("$> %.2f$", max(grid)))
  idx <- which(lb <= 0)[1]
  if (is.na(idx)) sprintf("$> %.2f$", max(grid)) else {
    if (max(grid) < 0.1) sprintf("$%.3f$", grid[idx]) else sprintf("$%.2f$", grid[idx])
  }
}

#' Assemble HonestDiD inputs (betahat, Sigma, l_vec) from a saved headline fit.
#'
#' @param fit saved fit list from 03_main_results.R
#' @param numPre required number of pre-treatment event times (4, 6 or 8)
#' @return list(betahat, sigma, numPre, numPost, l_vec, egt, max_rel_dev) or NULL
hd_inputs_from_fit <- function(fit, numPre) {
  # numPre pre-treatment event times require min_e = -(numPre + 1): the
  # universal base period e = -1 has a zero standard error and is dropped.
  #
  # NO SECOND BOOTSTRAP DRAW FOR THE PRIMARY WINDOW. fit$agg_dynamic is the
  # stored dynamic aggregation (min_e = -5, i.e. exactly four pre-treatment
  # event times) whose standard errors the paper reports; re-running aggte()
  # here would draw a fresh multiplier bootstrap and produce a second SE for
  # one specification -- exactly what the saved fits exist to prevent. The primary window
  # therefore READS the stored aggregation. The 6- and 8-lead windows do not
  # exist in the stored object and must be re-aggregated; that is flagged in
  # the returned `reaggregated` field and disclosed in the table note.
  min_e_k <- -(numPre + 1)
  reagg   <- !identical(as.integer(numPre), 4L)
  if (!reagg) {
    agg_k <- fit$agg_dynamic
    if (is.null(agg_k)) {
      message("    stored dynamic aggregation missing — falling back to aggte().")
      reagg <- TRUE
    } else {
      stopifnot(identical(as.numeric(agg_k$min_e), -5))
    }
  }
  if (reagg) {
    set.seed(1242)
    agg_k <- tryCatch(
      aggte(fit$gt_boot, type = "dynamic", na.rm = TRUE, min_e = min_e_k, max_e = Inf),
      error = function(e) { message("    aggte(min_e = ", min_e_k, ") failed: ",
                                     conditionMessage(e)); NULL })
  }
  if (is.null(agg_k)) return(NULL)

  bvc  <- tryCatch(boot_vcov_dynamic(agg_k, fit$gt_boot, seed = 1242L),
                   error = function(e) { message("    boot vcov failed: ",
                                                  conditionMessage(e)); NULL })
  if (is.null(bvc)) return(NULL)

  keep    <- bvc$keep
  egt     <- agg_k$egt[keep]
  betahat <- agg_k$att.egt[keep]
  nPre    <- sum(egt < 0)
  nPost   <- sum(egt >= 0)
  if (nPre != numPre) {
    message(sprintf("    requested numPre = %d but only %d pre-treatment event ",
                    numPre, nPre), "times are estimable — skipping this cell.")
    return(NULL)
  }

  # l_vec: group-size weights that reproduce the simple ATT (verified in 03).
  l_vec <- fit$l_vec_simple
  if (is.null(l_vec) || length(l_vec) != nPost) {
    message("    l_vec unavailable or wrong length — falling back to equal weights.")
    l_vec <- rep(1 / nPost, nPost)
  } else {
    # Hard stop: the sensitivity analysis must target the estimand Table 2
    # prints. If the group-size weights stop reproducing the reported simple
    # ATT, the l_vec is wrong and the exhibit would silently describe a
    # different parameter.
    gap <- sum(l_vec * betahat[egt >= 0]) - fit$att
    stopifnot(abs(gap) < 1e-6)
  }

  list(betahat = betahat, sigma = bvc$sigma, numPre = nPre, numPost = nPost,
       l_vec = l_vec, egt = egt, max_rel_dev = bvc$max_rel_dev,
       reaggregated = reagg)
}

short_labels_hd <- c("log(Adapt.)", "Share (pp)", "log(Total)",
                     "log(Non-adapt.)", "log(Disb.)")

hd_all  <- list()   # primary RM results (numPre = 4)
hd_pre  <- list()   # RM breakdown by numPrePeriods
hd_sd   <- list()   # smoothness-restriction results

for (oc in outcomes) {
  message(sprintf("\n=== HonestDiD: %s ===", oc$label))
  stem <- fit_stem_map[[oc$var]]
  fit  <- read_headline_fit(stem)

  # Runtime: the 6- and 8-lead relative-magnitude ladders are expensive and are
  # only interpretatively interesting for the two primary outcomes (the three
  # aggregate outcomes already break down at Mbar = 0, so a longer pre-window
  # cannot change their verdict). Restricting them keeps §20 under ~5 minutes.
  ks <- if (oc$var %in% PRIMARY_OUTCOMES) numPre_grid else numPre_grid[1L]
  for (k in ks) {
    hin <- hd_inputs_from_fit(fit, k)
    if (is.null(hin)) next
    message(sprintf("  numPre = %d | numPost = %d | event times %s",
                    hin$numPre, hin$numPost,
                    paste(sprintf("%d", as.integer(hin$egt)), collapse = ",")))
    message(sprintf(paste0("  aggregation: %s | bootstrap vcov: max ",
                           "|se_cov - se_diag| / se_diag = %.4f"),
                    if (isTRUE(hin$reaggregated))
                      sprintf("re-aggregated at min_e = %d under seed 1242", -(k + 1L))
                    else "stored (the SEs Table 2 reports)",
                    hin$max_rel_dev))

    rm_sens <- tryCatch(
      createSensitivityResults_relativeMagnitudes(
        betahat = hin$betahat, sigma = hin$sigma,
        numPrePeriods = hin$numPre, numPostPeriods = hin$numPost,
        Mbarvec = Mbarvec_hd, l_vec = hin$l_vec, alpha = 0.05),
      error = function(e) { message("  RM sensitivity failed: ",
                                     conditionMessage(e)); NULL })

    if (!is.null(rm_sens)) {
      bd <- breakdown_value(rm_sens$lb, Mbarvec_hd[seq_len(nrow(rm_sens))])
      message(sprintf("  RM breakdown Mbar (numPre = %d) = %s", k, bd))
      hd_pre[[paste0(oc$var, "_", k)]] <- list(var = oc$var, numPre = k,
                                               rm = rm_sens, breakdown = bd)
      if (k == numPre_grid[1L]) {
        hd_all[[oc$var]] <- list(label = oc$label, rm = rm_sens,
                                 numPre = hin$numPre, numPost = hin$numPost,
                                 l_vec = hin$l_vec, breakdown = bd,
                                 max_rel_dev = hin$max_rel_dev,
                                 reaggregated = hin$reaggregated)
      }
      hd_pre[[paste0(oc$var, "_", k)]]$reaggregated <- hin$reaggregated
    }

    # Smoothness restriction Delta^SD(M) on the primary (numPre = 4) window only.
    if (k == numPre_grid[1L]) {
      sd_sens <- tryCatch(
        createSensitivityResults(
          betahat = hin$betahat, sigma = hin$sigma,
          numPrePeriods = hin$numPre, numPostPeriods = hin$numPost,
          Mvec = Mvec_sd, l_vec = hin$l_vec, alpha = 0.05),
        error = function(e) { message("  SD sensitivity failed: ",
                                       conditionMessage(e)); NULL })
      if (!is.null(sd_sens)) {
        hd_sd[[oc$var]] <- list(label = oc$label, sd = sd_sens,
                                breakdown = breakdown_value(sd_sens$lb, Mvec_sd),
                                method = unique(as.character(sd_sens$method)))
        message(sprintf("  SD breakdown M = %s (method: %s)",
                        hd_sd[[oc$var]]$breakdown,
                        paste(hd_sd[[oc$var]]$method, collapse = "/")))
      }
    }
  }
}

message("\n=== Section 20 estimation complete ===\n")

# --- §20b HonestDiD LaTeX tables ---------------------------------------------

fmt_ci_hd <- function(lb, ub) {
  if (is.na(lb) || is.na(ub)) return("---")
  sprintf("[%.3f, %.3f]", lb, ub)
}

# (i) Relative magnitudes — tables/cohorts_dropped/honestdid_rm.tex
#     Layout is unchanged (Mbar rows + breakdown row, one column per outcome)
#     so the manuscript compiles without edits.
make_honestdid_table <- function(hd_all, Mbarvec, outcomes, dir_tabs) {

  n_mbar   <- length(Mbarvec)
  body_mat <- matrix("---", nrow = n_mbar, ncol = length(outcomes))
  rownames(body_mat) <- sprintf("$\\bar{M} = %.2f$", Mbarvec)
  colnames(body_mat) <- short_labels_hd

  for (j in seq_along(outcomes)) {
    res <- hd_all[[outcomes[[j]]$var]]
    if (is.null(res$rm)) next
    for (i in seq_len(min(n_mbar, nrow(res$rm))))
      body_mat[i, j] <- fmt_ci_hd(res$rm$lb[i], res$rm$ub[i])
  }

  breakdown_row <- matrix(vapply(outcomes, function(oc) {
    res <- hd_all[[oc$var]]
    if (is.null(res$rm)) "---" else res$breakdown
  }, character(1L)), nrow = 1)
  rownames(breakdown_row) <- "Breakdown $\\bar{M}$"
  colnames(breakdown_row) <- short_labels_hd

  full_mat <- rbind(body_mat, breakdown_row)
  tab_df   <- as.data.frame(full_mat, stringsAsFactors = FALSE)
  tab_df   <- cbind(`$\\bar{M}$` = rownames(full_mat), tab_df, stringsAsFactors = FALSE)
  rownames(tab_df) <- NULL

  xt <- xtable(tab_df, label = "tab:honestdid")
  dir.create(dir_tabs, showWarnings = FALSE, recursive = TRUE)
  raw_lines <- capture.output(
    print(xt, include.rownames = FALSE, booktabs = TRUE,
          sanitize.text.function = identity, size = "\\small", floating = FALSE))

  l_vec_txt <- {
    lv <- hd_all[[outcomes[[1L]]$var]]$l_vec
    if (is.null(lv)) "equal weights" else
      paste0("(", paste(sprintf("%.4f", lv), collapse = ", "), ")")
  }
  # Largest covariance-vs-reported standard-error gap across outcomes, printed
  # in the note so the rescaling is quantified rather than merely asserted.
  mrd_vals <- vapply(outcomes, function(oc) {
    r <- hd_all[[oc$var]]
    if (is.null(r$max_rel_dev)) NA_real_ else r$max_rel_dev
  }, numeric(1L))
  mrd_txt <- if (all(is.na(mrd_vals))) "not available" else
    sprintf("%.1f\\%% (%s)", 100 * max(mrd_vals, na.rm = TRUE),
            short_labels_hd[which.max(replace(mrd_vals, is.na(mrd_vals), -Inf))])

  notes_txt <- paste0(
    "Sensitivity of \\citet{rambachan_roth_2023} on the headline spec (Table~",
    "\\ref{tab:combined_wide_main}): CS\\,(2021) DR, never-treated controls, ",
    "multiplier-bootstrap SE (", BITERS, " reps, seed 1242); event-study ",
    "covariance from the stored fit, not re-estimated. Target: group-size-weighted ",
    "average post-treatment effect ($e=0,\\dots,3$), weights ", l_vec_txt,
    ", reproducing the simple ATT. $\\Sigma$ rescaled to reported bootstrap SEs ",
    "(max gap ", mrd_txt, "). Four pre-treatment leads ($e=-5,\\dots,-2$). ",
    "``Breakdown $\\bar{M}$'' = smallest grid value where the robust CI contains zero"
  )
  source_txt <- "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"

  write_tex_float(
    file.path(dir_tabs, "honestdid_rm.tex"),
    caption_title = paste0("HonestDiD robust 95\\% CI under relative magnitude ",
                           "restriction (\\citealt{rambachan_roth_2023})"),
    label         = "tab:honestdid",
    tabular_lines = raw_lines,
    notes_text    = notes_txt,
    source_text   = source_txt)
}

make_honestdid_table(hd_all, Mbarvec_hd, outcomes,
                     here("output", "tables", "cohorts_dropped"))

# (ii) Breakdown Mbar by number of pre-treatment periods -----------------------
if (length(hd_pre) > 0) {
  pre_mat <- matrix("---", nrow = length(numPre_grid), ncol = length(outcomes))
  rownames(pre_mat) <- sprintf("%d pre-treatment periods", numPre_grid)
  colnames(pre_mat) <- short_labels_hd
  for (i in seq_along(numPre_grid)) {
    for (j in seq_along(outcomes)) {
      key <- paste0(outcomes[[j]]$var, "_", numPre_grid[i])
      if (!is.null(hd_pre[[key]])) pre_mat[i, j] <- hd_pre[[key]]$breakdown
    }
  }
  win_row <- matrix(sprintf("$e \\in [-%d, 3]$", numPre_grid + 1L), ncol = 1)
  tab_pre <- data.frame(
    `Pre-treatment window` = paste0(rownames(pre_mat), " ($e \\in [-",
                                    numPre_grid + 1L, ", -2]$)"),
    pre_mat, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(tab_pre) <- NULL

  xt_pre <- xtable(tab_pre, label = "tab:honestdid_prepriods")
  raw_pre <- capture.output(
    print(xt_pre, include.rownames = FALSE, booktabs = TRUE,
          sanitize.text.function = identity, size = "\\small", floating = FALSE))

  notes_pre <- paste0(
    "Each cell: breakdown value $\\bar{M}$, disciplined by 4, 6 or 8 pre-treatment ",
    "event times. Estimator/covariance as Table~\\ref{tab:honestdid}; only ",
    "\\texttt{min\\_e} changes ($-5$, $-7$, $-9$). 4-lead row: stored aggregation ",
    "behind Table~\\ref{tab:combined_wide_main}. 6-/8-lead rows: re-aggregated from ",
    "the same stored fit (seed 1242), two primary outcomes only; point estimates ",
    "unchanged, SEs a fresh bootstrap draw. Three aggregate outcomes already break ",
    "at $\\bar{M} = 0$ with four leads, so those cells are left blank"
  )

  write_tex_float(
    here("output", "tables", "cohorts_dropped", "honestdid_prepriods.tex"),
    caption_title = paste0("HonestDiD breakdown value $\\bar{M}$ by the number ",
                           "of pre-treatment periods used"),
    label         = "tab:honestdid_prepriods",
    tabular_lines = raw_pre,
    notes_text    = notes_pre,
    source_text   = "OECD CRS (Rio adaptation markers); UNFCCC NAP Central")
}

# (iii) Smoothness restriction Delta^SD(M) ------------------------------------
if (length(hd_sd) > 0) {
  sd_mat <- matrix("---", nrow = length(Mvec_sd), ncol = length(outcomes))
  rownames(sd_mat) <- sprintf("$M = %.2f$", Mvec_sd)
  colnames(sd_mat) <- short_labels_hd
  for (j in seq_along(outcomes)) {
    res <- hd_sd[[outcomes[[j]]$var]]
    if (is.null(res$sd)) next
    for (i in seq_len(min(length(Mvec_sd), nrow(res$sd))))
      sd_mat[i, j] <- fmt_ci_hd(res$sd$lb[i], res$sd$ub[i])
  }
  bd_sd <- matrix(vapply(outcomes, function(oc) {
    res <- hd_sd[[oc$var]]
    if (is.null(res$sd)) "---" else res$breakdown
  }, character(1L)), nrow = 1)
  rownames(bd_sd) <- "Breakdown $M$"
  colnames(bd_sd) <- short_labels_hd

  full_sd <- rbind(sd_mat, bd_sd)
  tab_sd  <- cbind(`$M$` = rownames(full_sd),
                   as.data.frame(full_sd, stringsAsFactors = FALSE),
                   stringsAsFactors = FALSE)
  rownames(tab_sd) <- NULL

  xt_sd  <- xtable(tab_sd, label = "tab:honestdid_sd")
  raw_sd <- capture.output(
    print(xt_sd, include.rownames = FALSE, booktabs = TRUE,
          sanitize.text.function = identity, size = "\\small", floating = FALSE))

  sd_methods <- unique(unlist(lapply(hd_sd, `[[`, "method")))
  notes_sd <- paste0(
    "Smoothness restriction $\\Delta^{SD}(M)$ of \\citet{rambachan_roth_2023}: the ",
    "differential trend's slope changes by at most $M$ between periods; $M = 0$ still ",
    "allows an exactly linear violation (not the conventional CI). Confidence sets: ",
    paste(sd_methods, collapse = "/"), " method, 5\\% level. Estimator, stored fit, ",
    "covariance and target weights as Table~\\ref{tab:honestdid}; four pre-treatment ",
    "leads ($e = -5,\\dots,-2$). ``Breakdown $M$'' = smallest grid value at which the ",
    "robust CI contains zero. $M$ is in outcome units per period, not comparable ",
    "across the log and share columns"
  )

  write_tex_float(
    here("output", "tables", "cohorts_dropped", "honestdid_sd.tex"),
    caption_title = paste0("HonestDiD robust 95\\% CI under the smoothness ",
                           "restriction $\\Delta^{SD}(M)$"),
    label         = "tab:honestdid_sd",
    tabular_lines = raw_sd,
    notes_text    = notes_sd,
    source_text   = "OECD CRS (Rio adaptation markers); UNFCCC NAP Central")
}

message("\n=== Section 20 complete ===\n")


# ==============================================================================
# SECTION 7. §25 Placebo test
# Seed rule (reproduces the published SEs): set.seed(1242) before gt_placebo call.
# ==============================================================================

message("\n=== Section 25: Placebo test (treatment shifted -2 years) ===\n")

dir_figs_placebo <- here("output", "figures", "placebo")
dir_tabs_placebo <- here("output", "tables",  "placebo")
dir.create(dir_figs_placebo, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tabs_placebo, recursive = TRUE, showWarnings = FALSE)

did_panel_placebo <- bind_rows(
  did_panel_full %>%
    filter(cohort_year > 0, year < cohort_year) %>%
    mutate(cohort_year = cohort_year - 2L),
  did_panel_full %>%
    filter(cohort_year == 0)
) %>%
  filter(cohort_year == 0 | cohort_year >= 2015)

message(sprintf("  Placebo panel: %d obs | %d countries | %d treated cohorts",
  nrow(did_panel_placebo),
  n_distinct(did_panel_placebo$country_id),
  n_distinct(did_panel_placebo$cohort_year[did_panel_placebo$cohort_year > 0])))

thin_placebo <- did_panel_placebo %>%
  filter(cohort_year > 0) %>%
  group_by(cohort_year) %>%
  summarise(n = n_distinct(country_id), .groups = "drop") %>%
  filter(n < 5L) %>%
  pull(cohort_year)

did_panel_placebo <- did_panel_placebo %>%
  filter(!(cohort_year %in% thin_placebo))

message(sprintf("  After dropping thin cohorts: %d obs | %d treated cohorts",
  nrow(did_panel_placebo),
  n_distinct(did_panel_placebo$cohort_year[did_panel_placebo$cohort_year > 0])))

n_treated_plac <- n_distinct(
  did_panel_placebo$country_id[did_panel_placebo$cohort_year > 0])
# Placebo uses outcome regression (reg): the doubly-robust estimator is not
# stable on this truncated panel. Inference is the multiplier bootstrap in every
# case (see the note immediately below) -- the old "bootstrap only when
# N >= 40" rule is gone, so the reported placebo fit is always the bootstrap
# one and the analytical twin exists solely to feed the pre-trend test.
em_plac <- "reg"
# The `n_treated >= 40` bootstrap switch is removed project-wide -- 40 is a
# convention, not an estimator requirement, and it disabled the bootstrap
# exactly where analytical influence-function SEs are most anti-conservative.
# The multiplier bootstrap is now always on, with an explicit replication count
# (previously this call passed no `biters` and silently used did's default of
# 1000 while the table note claimed 999).
bs_plac <- TRUE
message(sprintf("  N treated (placebo) = %d  ->  est_method = %s  bstrap = %s (%d reps)",
                n_treated_plac, em_plac, bs_plac, BITERS))

# Seed set immediately before the estimator (reproduces the published SE)
set.seed(1242)
gt_placebo <- tryCatch(
  att_gt(
    yname         = "log_commits",
    tname         = "year",
    idname        = "country_id",
    gname         = "cohort_year",
    xformla       = ~ ge_est + log_population,
    data          = did_panel_placebo,
    est_method    = em_plac,
    bstrap        = bs_plac,
    biters        = BITERS,
    cband         = FALSE,
    control_group = "nevertreated",
    anticipation  = 0,
    base_period   = "universal",
    panel         = TRUE,
    allow_unbalanced_panel = TRUE
  ),
  error = function(e) {
    message("  att_gt (", em_plac, ") failed: ", conditionMessage(e))
    NULL
  }
)

# Analytical twin (bstrap = FALSE, otherwise identical) so the pre-trend test
# below is never bootstrap-based, matching make_wide_table()'s convention
# (gt_tab_analytical). The reported placebo ATT and SE always come from the
# bootstrap fit gt_placebo above; this twin is never reported.
gt_placebo_analytical <- tryCatch(
  att_gt(
    yname         = "log_commits",
    tname         = "year",
    idname        = "country_id",
    gname         = "cohort_year",
    xformla       = ~ ge_est + log_population,
    data          = did_panel_placebo,
    est_method    = em_plac,
    bstrap        = FALSE,
    cband         = FALSE,
    control_group = "nevertreated",
    anticipation  = 0,
    base_period   = "universal",
    panel         = TRUE,
    allow_unbalanced_panel = TRUE
  ),
  error = function(e) {
    message("  att_gt (analytical twin) failed: ", conditionMessage(e))
    NULL
  }
)

if (!is.null(gt_placebo)) {

  agg_s_plac <- tryCatch(aggte(gt_placebo, type = "simple",  na.rm = TRUE),
                          error = function(e) NULL)
  agg_d_plac <- tryCatch(aggte(gt_placebo, type = "dynamic", na.rm = TRUE,
                                min_e = -4, max_e = Inf),
                          error = function(e) NULL)
  # Analytical dynamic aggregation feeding the pre-trend test only (figure /
  # CI ribbon above still use agg_d_plac from the primary, possibly-bootstrap
  # fit, unchanged).
  agg_d_plac_analytical <- if (!is.null(gt_placebo_analytical)) tryCatch(
    aggte(gt_placebo_analytical, type = "dynamic", na.rm = TRUE,
          min_e = -4, max_e = Inf),
    error = function(e) NULL
  ) else NULL

  if (!is.null(agg_s_plac)) {
    att_p   <- agg_s_plac$overall.att
    se_p    <- agg_s_plac$overall.se
    t_p     <- att_p / se_p
    stars_p <- if (abs(t_p) > 2.576) "***" else if (abs(t_p) > 1.960) "**" else
               if (abs(t_p) > 1.645) "*" else ""
    message(sprintf("  Placebo ATT = %.4f%s (SE = %.4f, t = %.3f)",
                    att_p, stars_p, se_p, t_p))


  }

  if (!is.null(agg_d_plac)) {
    cv_p <- agg_d_plac$crit.val.egt
    es_plac <- data.frame(
      event_time = agg_d_plac$egt,
      ATT        = agg_d_plac$att.egt,
      SE         = agg_d_plac$se.egt,
      Lower      = agg_d_plac$att.egt - cv_p * agg_d_plac$se.egt,
      Upper      = agg_d_plac$att.egt + cv_p * agg_d_plac$se.egt,
      stringsAsFactors = FALSE
    ) %>% filter(!is.na(SE), SE > 1e-10)

    p_placebo <- ggplot(es_plac, aes(x = event_time, y = ATT)) +
      geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed", linewidth = 0.4) +
      geom_vline(xintercept = -0.5, colour = "grey60", linetype = "dotted", linewidth = 0.4) +
      geom_linerange(aes(ymin = Lower, ymax = Upper),
                     linewidth = 0.6, alpha = 0.8, colour = "#2166ac") +
      geom_point(size = 2.5, colour = "#2166ac") +
      labs(
        title    = NULL, subtitle = NULL,  # titles go in the LaTeX caption
        x        = "Event time (years relative to placebo date)",
        y        = "ATT — log(adaptation commitments)"
      ) +
      theme_minimal() +
      theme(
        text             = element_text(family = "serif", size = 11),
        panel.grid.minor = element_blank()
      )
    ggsave(file.path(dir_figs_placebo, "did_placebo_es.png"),
           p_placebo, width = 10, height = 5, dpi = 300)
    message("Saved: ", file.path(dir_figs_placebo, "did_placebo_es.png"))
  }

  if (!is.null(agg_s_plac)) {

    pt_plac <- if (!is.null(agg_d_plac_analytical))
      tryCatch(compute_pretrend_test(agg_d_plac_analytical, gt_placebo_analytical),
               error = function(e) list(stat = NA_real_, pval = NA_real_, df = 0L,
                                         W_did = NA_real_, Wpval_did = NA_real_,
                                         df_did = NA_integer_))
    else
      list(stat = NA_real_, pval = NA_real_, df = 0L,
           W_did = NA_real_, Wpval_did = NA_real_,
           df_did = NA_integer_)

    rows_p     <- did_panel_placebo[!is.na(did_panel_placebo$log_commits), ]
    n_obs_p    <- nrow(rows_p)
    n_cty_p    <- length(unique(rows_p$country_id))

    pre_rows_p <- did_panel_placebo %>%
      filter(cohort_year > 0, year < cohort_year, !is.na(commitments))
    mean_pre_p <- mean(pre_rows_p$commitments, na.rm = TRUE)
    # Suppress the naive back-transform unless the ATT is significant
    # at 5% (|t| > 1.960) -- see make_wide_table() for the full rationale.
    is_sig5_p <- !is.na(t_p) && abs(t_p) > 1.960
    implied_p <- if (is_sig5_p) (exp(att_p) - 1) * mean_pre_p else NA_real_
    implied_fmt_p <- if (is.na(implied_p)) "---" else sprintf("%.1f", implied_p)

    row_labels_p <- c(
      "ATT", "SE", "$t$-statistic",
      "Mean (pre-treat., USD M)",
      "Implied effect (USD M, exp(ATT)$-$1 $\\times$ pre-treat.\\ mean)",
      "Observations", "Countries",
      "Pre-trend $\\chi^2$", "Pre-trend $p$", "\\texttt{did} pre-test $p$",
      "\\midrule Estimator", "Control group", "Bootstrap SE"
    )
    col_vals_p <- c(
      paste0(sprintf("%.4f", att_p), stars_p),
      sprintf("(%.4f)", se_p),
      sprintf("%.3f", t_p),
      sprintf("%.1f", mean_pre_p),
      implied_fmt_p,
      format(n_obs_p, big.mark = ","),
      as.character(n_cty_p),
      sprintf("%.3f", pt_plac$stat),
      sprintf("%.3f", pt_plac$pval),
      if (is.na(pt_plac$Wpval_did)) "---" else sprintf("%.3f", pt_plac$Wpval_did),
      "CS (2021)", "Never-treated",
      if (bs_plac) "Yes" else paste0("No (", em_plac, ")")
    )

    tab_plac <- data.frame(
      ` ` = row_labels_p,
      `Placebo (treatment shifted $-2$ years)` = col_vals_p,
      check.names = FALSE, stringsAsFactors = FALSE
    )

    xtab_plac <- xtable(tab_plac, label = "tab:placebo")
    align(xtab_plac) <- "llc"

    raw_plac <- capture.output(
      print(xtab_plac, include.rownames = FALSE, booktabs = TRUE,
            sanitize.text.function = identity, size = "\\small",
            floating = FALSE)
    )

    se_label_p <- if (bs_plac) paste0("multiplier-bootstrap SE (", BITERS, " reps), seed 1242") else "analytical SE"
    est_label_p <- if (em_plac == "dr") "DR" else "regression adjustment"

    notes_plac <- paste0(
      "CS\\,(2021) ", est_label_p,
      "; ", se_label_p,
      "; never-treated control; WGI GE + log population. ",
      wpval_reconciliation(pre_egt = pt_plac$leads, df_did = pt_plac$df_did,
                           n_clusters = n_cty_p, wpval_did = pt_plac$Wpval_did,
                           pval_wald = pt_plac$pval,
                           wpval_reason = pt_plac$wpval_reason,
                           ginv_used = pt_plac$ginv_used),
      "Implied effect: back-transform on pre-treatment mean; suppressed if ",
      "insignificant at 5\\%. ",
      "* $p<0.10$, ** $p<0.05$, *** $p<0.01$"
    )
    source_plac <- "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"

    out_path_plac <- file.path(dir_tabs_placebo, "att_placebo.tex")
    write_tex_float(
      out_path_plac,
      "Placebo test: treatment date assigned 2 years before actual NAP adoption",
      "tab:placebo",
      raw_plac,
      notes_plac,
      source_plac
    )
    # Persist this fit so that 13_cohort_anticipation.R can place the
    # published -2-year placebo at the foot of its -3/-4/-5 ladder without
    # re-estimating it (which would print a second SE for one specification).
    saveRDS(list(
      shift_years = 2L,
      att         = att_p,
      se          = se_p,
      n_obs       = n_obs_p,
      n_country   = n_cty_p,
      n_treated   = n_treated_plac,
      cohorts     = sort(unique(did_panel_placebo$cohort_year[
                          did_panel_placebo$cohort_year > 0])),
      spec        = list(estimator = "CS (2021)", est_method = em_plac,
                         control_group = "nevertreated", bstrap = bs_plac,
                         biters = BITERS, seed = 1242L,
                         xformla = "~ ge_est + log_population",
                         cohort_rule = paste0("shifted cohorts >= 2015 with at ",
                                              "least 5 treated recipients")),
      pretrend    = pt_plac,
      provenance  = list(script = "code/04_robustness.R",
                         created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    ), file.path(FITS_DIR, "placebo_m2_reg_bs.rds"))
    message("  [fit] Saved: ", file.path(FITS_DIR, "placebo_m2_reg_bs.rds"))
  }

} else {
  message("  Placebo estimation failed — skipping table and figure.")
}

message("\n=== Section 25 complete ===\n")

# ==============================================================================
# SECTION 8. §26 Mitigation falsification test
# Seed rule (reproduces the published SEs): set.seed(1242) before gt_mit call.
# ==============================================================================

message("\n=== Section 26: Mitigation falsification test ===\n")

dir_tabs_mit <- here("output", "tables",  "mitigation")
dir_figs_mit <- here("output", "figures", "mitigation")
dir.create(dir_tabs_mit, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figs_mit, recursive = TRUE, showWarnings = FALSE)

mit_panel_path <- here("data", "processed", "mitigation_panel.csv")
if (!file.exists(mit_panel_path)) {
  message("WARNING: mitigation_panel.csv not found — skipping Section 26.")
  message("  Expected path: ", mit_panel_path)
  message("  This file is produced by 01_prepare_data.R.")
} else {

  mitigation_ry <- as.data.frame(fread(mit_panel_path))

  required_mit_cols <- c("recipient_name", "year", "commits_mitigation")
  missing_mit <- setdiff(required_mit_cols, names(mitigation_ry))
  if (length(missing_mit) > 0) {
    message("  mitigation_panel.csv missing columns: ", paste(missing_mit, collapse = ", "))
    message("  Skipping Section 26.")
  } else {

  did_panel_mit_full <- did_panel_full %>%
    left_join(mitigation_ry, by = c("recipient_name", "year")) %>%
    mutate(
      commits_mitigation = if_else(is.na(commits_mitigation), 0, commits_mitigation),
      log_commits_mit    = log1p(commits_mitigation)
    )

  message(sprintf("  Merged: %d rows | %d countries | mean mitigation commits = %.1f M USD",
                  nrow(did_panel_mit_full),
                  n_distinct(did_panel_mit_full$country_id),
                  mean(did_panel_mit_full$commits_mitigation, na.rm = TRUE)))

    did_panel_mit_main <- did_panel_mit_full %>%
      filter(!(cohort_year %in% thin_cohorts))

    n_treated_mit <- n_distinct(
      did_panel_mit_main$country_id[did_panel_mit_main$cohort_year > 0])
    em_mit <- if (n_treated_mit >= 40L) "dr" else "reg"
    # Bootstrap always on (the >= 40 convention is removed project-wide).
    # With 40 treated units this evaluates to the same TRUE as before, so the
    # published mitigation SE is unchanged by the rule change itself.
    bs_mit <- TRUE
    message(sprintf("  N treated = %d  ->  est_method = '%s'  bstrap = %s",
                    n_treated_mit, em_mit, bs_mit))

    run_mit_att_gt <- function(est_method, bstrap) {
      att_gt(
        yname         = "log_commits_mit",
        tname         = "year",
        idname        = "country_id",
        gname         = "cohort_year",
        xformla       = ~ ge_est + log_population,
        data          = did_panel_mit_main,
        est_method    = est_method,
        bstrap        = bstrap,
        biters        = BITERS,  # One replication constant per script
        cband         = FALSE,
        control_group = "nevertreated",
        anticipation  = 0,
        base_period   = "universal",
        panel         = TRUE,
        allow_unbalanced_panel = TRUE
      )
    }

    # Seed set immediately before the estimator (reproduces the published SE)
    set.seed(1242)
    gt_mit <- tryCatch(
      run_mit_att_gt(em_mit, bs_mit),
      error = function(e) {
        message("  att_gt (", em_mit, ") failed: ", conditionMessage(e))
        if (em_mit == "dr") {
          message("  Retrying with est_method = 'reg' ...")
          tryCatch(run_mit_att_gt("reg", FALSE),
                   error = function(e2) {
                     message("  reg fallback also failed: ", conditionMessage(e2))
                     NULL
                   })
        } else NULL
      }
    )

    # Analytical twin (bstrap = FALSE, otherwise identical est_method) so the
    # pre-trend test below is never bootstrap-based, matching make_wide_table()'s
    # convention (gt_tab_analytical). The reported mitigation ATT and SE always
    # come from the bootstrap fit gt_mit above; this twin supplies the influence
    # function that 05_heterogeneity.R uses for the adaptation-minus-mitigation
    # contrast, and is never reported on its own.
    gt_mit_analytical <- tryCatch(
      run_mit_att_gt(em_mit, FALSE),
      error = function(e) {
        message("  att_gt (analytical twin) failed: ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(gt_mit)) {

      agg_s_mit <- tryCatch(aggte(gt_mit, type = "simple",  na.rm = TRUE),
                             error = function(e) NULL)
      agg_d_mit <- tryCatch(aggte(gt_mit, type = "dynamic", na.rm = TRUE,
                                   min_e = -5, max_e = Inf),
                             error = function(e) NULL)
      # Analytical dynamic aggregation feeding the pre-trend test only (figure
      # / CI ribbon above still use agg_d_mit from the primary, possibly-
      # bootstrap fit, unchanged).
      agg_d_mit_analytical <- if (!is.null(gt_mit_analytical)) tryCatch(
        aggte(gt_mit_analytical, type = "dynamic", na.rm = TRUE,
              min_e = -5, max_e = Inf),
        error = function(e) NULL
      ) else NULL

      if (!is.null(agg_s_mit)) {
        att_m   <- agg_s_mit$overall.att
        se_m    <- agg_s_mit$overall.se
        t_m     <- att_m / se_m
        stars_m <- if (abs(t_m) > 2.576) "***" else if (abs(t_m) > 1.960) "**" else
                   if (abs(t_m) > 1.645) "*" else ""
        message(sprintf("  Mitigation ATT = %.4f%s (SE = %.4f, t = %.3f)",
                        att_m, stars_m, se_m, t_m))

        # Persist the mitigation fit so that 05_heterogeneity.R can form
        # the adaptation-minus-mitigation contrast with ALIGNED unit-level
        # influence functions instead of re-estimating the falsification test
        # on its own (which would print a second SE for one specification).
        # The analytical twin supplies the influence function; the bootstrap fit
        # supplies the reported ATT/SE. `ids` records the unit order the
        # influence-function rows correspond to (did sorts panel units by
        # idname), so the alignment can be asserted rather than assumed.
        agg_s_mit_analytical <- if (!is.null(gt_mit_analytical))
          tryCatch(aggte(gt_mit_analytical, type = "simple", na.rm = TRUE),
                   error = function(e) NULL) else NULL
        saveRDS(list(
          outcome       = "log_commits_mit",
          outcome_label = "log(Mitigation commitments)",
          spec          = list(estimator = "CS (2021)", est_method = em_mit,
                               control_group = "nevertreated",
                               xformla = "~ ge_est + log_population",
                               anticipation = 0L, base_period = "universal",
                               cohort_rule = "adoption cohorts with >= 5 treated units",
                               bstrap = bs_mit, biters = BITERS, seed = 1242L),
          agg_simple       = agg_s_mit,
          agg_simple_analytic = agg_s_mit_analytical,
          att              = att_m,
          se               = se_m,
          ids              = sort(unique(did_panel_mit_main$country_id)),
          n_obs            = sum(!is.na(did_panel_mit_main$log_commits_mit)),
          provenance    = list(script = "code/04_robustness.R",
                               created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
        ), file.path(FITS_DIR, "mitigation_dr_bs.rds"))
        message("  [fit] Saved: ", file.path(FITS_DIR, "mitigation_dr_bs.rds"))
      } else {
        att_m <- NA_real_; se_m <- NA_real_; t_m <- NA_real_; stars_m <- ""
      }

      # Event-study figure
      if (!is.null(agg_d_mit)) {
        cv_m <- agg_d_mit$crit.val.egt
        es_mit <- data.frame(
          event_time = agg_d_mit$egt,
          ATT        = agg_d_mit$att.egt,
          SE         = agg_d_mit$se.egt,
          Lower      = agg_d_mit$att.egt - cv_m * agg_d_mit$se.egt,
          Upper      = agg_d_mit$att.egt + cv_m * agg_d_mit$se.egt,
          stringsAsFactors = FALSE
        ) %>% filter(!is.na(SE), SE > 1e-10)

        p_mit <- ggplot(es_mit, aes(x = event_time, y = ATT)) +
          geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed",
                     linewidth = 0.4) +
          geom_vline(xintercept = -0.5, colour = "grey60", linetype = "dotted",
                     linewidth = 0.4) +
          geom_linerange(aes(ymin = Lower, ymax = Upper),
                         linewidth = 0.6, alpha = 0.8, colour = "#8E44AD") +
          geom_point(size = 2.5, colour = "#8E44AD") +
          labs(
            title    = NULL, subtitle = NULL,  # titles go in the LaTeX caption
            x        = "Event time (years relative to NAP adoption)",
            y        = "ATT — log(mitigation commitments)"
          ) +
          theme_minimal() +
          theme(
            text             = element_text(family = "serif", size = 11),
            panel.grid.minor = element_blank()
          )
        ggsave(file.path(dir_figs_mit, "did_mitigation_es.png"),
               p_mit, width = 10, height = 5, dpi = 300)
        message("Saved: ", file.path(dir_figs_mit, "did_mitigation_es.png"))
      }

      # LaTeX table
      if (!is.null(agg_s_mit)) {

        pt_mit <- if (!is.null(agg_d_mit_analytical))
          tryCatch(compute_pretrend_test(agg_d_mit_analytical, gt_mit_analytical),
                   error = function(e) list(stat = NA_real_, pval = NA_real_, df = 0L,
                                             W_did = NA_real_, Wpval_did = NA_real_,
                                             df_did = NA_integer_))
        else
          list(stat = NA_real_, pval = NA_real_, df = 0L,
               W_did = NA_real_, Wpval_did = NA_real_,
               df_did = NA_integer_)

        rows_m     <- did_panel_mit_main[!is.na(did_panel_mit_main$log_commits_mit), ]
        n_obs_m    <- nrow(rows_m)
        n_cty_m    <- length(unique(rows_m$country_id))

        pre_rows_m <- did_panel_mit_main %>%
          filter(cohort_year > 0, year < cohort_year, !is.na(commits_mitigation))
        mean_pre_m <- mean(pre_rows_m$commits_mitigation, na.rm = TRUE)
        # Suppress the naive back-transform unless the ATT is
        # significant at 5% (|t| > 1.960) -- see make_wide_table() for the
        # full rationale.
        is_sig5_m <- !is.na(t_m) && abs(t_m) > 1.960
        implied_m <- if (is_sig5_m) (exp(att_m) - 1) * mean_pre_m else NA_real_
        implied_fmt_m <- if (is.na(implied_m)) "---" else sprintf("%.1f", implied_m)

        row_labels_m <- c(
          "ATT", "SE", "$t$-statistic",
          "Mean (pre-treat., USD M)",
          "Implied effect (USD M, exp(ATT)$-$1 $\\times$ pre-treat.\\ mean)",
          "Observations", "Countries",
          "Pre-trend $\\chi^2$", "Pre-trend $p$", "\\texttt{did} pre-test $p$",
          "\\midrule Estimator", "Control group", "Bootstrap SE"
        )
        col_vals_m <- c(
          paste0(sprintf("%.4f", att_m), stars_m),
          sprintf("(%.4f)", se_m),
          sprintf("%.3f", t_m),
          sprintf("%.1f", mean_pre_m),
          implied_fmt_m,
          format(n_obs_m, big.mark = ","),
          as.character(n_cty_m),
          sprintf("%.3f", pt_mit$stat),
          sprintf("%.3f", pt_mit$pval),
          if (is.na(pt_mit$Wpval_did)) "---" else sprintf("%.3f", pt_mit$Wpval_did),
          "CS (2021)",
          "Never-treated",
          if (bs_mit) "Yes (multiplier, 999 reps)" else paste0("No (", em_mit, ")")
        )

        tab_mit <- data.frame(
          ` ` = row_labels_m,
          `log(Mitigation commitments)` = col_vals_m,
          check.names = FALSE, stringsAsFactors = FALSE
        )

        xtab_mit <- xtable(tab_mit, label = "tab:mitigation")
        align(xtab_mit) <- "llc"

        raw_mit <- capture.output(
          print(xtab_mit, include.rownames = FALSE, booktabs = TRUE,
                sanitize.text.function = identity, size = "\\small",
                floating = FALSE)
        )

        est_label_m <- if (em_mit == "dr") "DR" else "regression adjustment"
        se_label_m  <- if (bs_mit) "multiplier-bootstrap SE (999 reps, seed 1242; \\texttt{did} 2.5.0; CRS Apr.\\ 2026)" else "analytical (IF) SE; \\texttt{did} 2.5.0; CRS Apr.\\ 2026"

        notes_mit <- paste0(
          "CS\\,(2021) ", est_label_m,
          "; ", se_label_m,
          "; never-treated control; cohorts $\\geq 5$ units; ",
          "controls: WGI gov.\\ effectiveness + log population. ",
          wpval_reconciliation(pre_egt = pt_mit$leads, df_did = pt_mit$df_did,
                               n_clusters = n_cty_m, wpval_did = pt_mit$Wpval_did,
                               pval_wald = pt_mit$pval,
                               wpval_reason = pt_mit$wpval_reason,
                               ginv_used = pt_mit$ginv_used),
          "Implied effect: naive back-transform on the pre-treatment mean; suppressed if ",
          "insignificant at 5\\%. ",
          "* $p<0.10$, ** $p<0.05$, *** $p<0.01$"
        )
        source_mit <- "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"

        out_path_mit <- file.path(dir_tabs_mit, "att_mitigation.tex")
        write_tex_float(
          out_path_mit,
          "Falsification test: effect of NAP adoption on log(mitigation commitments)",
          "tab:mitigation",
          raw_mit,
          notes_mit,
          source_mit
        )
      }

    } else {
      message("  Mitigation estimation failed — skipping table and figure.")
    }
  }
}

message("\n=== Section 26 complete ===\n")

# ==============================================================================
# SECTION A. §27 de Chaisemartin & D'Haultfoeuille (2024) estimator
# Alternative heterogeneity-robust DID estimator on the headline outcome
# (log adaptation commitments).  Per the project's locked decision, dCDH uses the
# SAME controls as the CS main spec (WGI gov. effectiveness + log population).
# Panel: thin-dropped main panel.  Treatment is the current-treatment dummy
# treated_post = 1{cohort > 0 & year >= cohort}; never-treated and not-yet-treated
# units serve as controls.  4 dynamic effects, 3 placebos.  This is an APPENDIX
# exhibit.
# Seed rule (reproduces the published SEs): set.seed(1242) immediately before the estimator call.
# ==============================================================================

message("\n=== Section A (27): de Chaisemartin & D'Haultfoeuille (2024) ===\n")

dir_tabs_dcdh <- file.path(here("output", "tables"),  "dcdh")
dir_figs_dcdh <- file.path(here("output", "figures"), "dcdh")
dir.create(dir_tabs_dcdh, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figs_dcdh, recursive = TRUE, showWarnings = FALSE)

# Thin-dropped main panel (same definition as §5b / §20) + current-treatment dummy.
did_panel_dcdh <- did_panel_full %>%
  filter(!(cohort_year %in% thin_cohorts)) %>%
  mutate(treated_post = as.integer(cohort_year > 0 & year >= cohort_year))

message(sprintf("dCDH panel: %d rows, %d countries, years %d-%d",
                nrow(did_panel_dcdh), length(unique(did_panel_dcdh$country_id)),
                min(did_panel_dcdh$year), max(did_panel_dcdh$year)))

set.seed(1242)  # Seed rule (reproduces the published SEs): seed immediately before the estimator
res_dcdh <- tryCatch(
  did_multiplegt_dyn(
    df        = did_panel_dcdh,
    outcome   = "log_commits",
    group     = "country_id",
    time      = "year",
    treatment = "treated_post",
    effects   = 4,
    placebo   = 3,
    controls  = c("ge_est", "log_population"),
    graph_off = TRUE
  ),
  error = function(e) { message("  did_multiplegt_dyn failed: ", conditionMessage(e)); NULL }
)

if (is.null(res_dcdh)) {
  message("  dCDH estimation returned NULL — skipping table and figure.")
} else {

  # --- Pull results.  Package column names: "Estimate","SE","LB CI","UB CI","N","Switchers".
  ate_row  <- res_dcdh$results$ATE          # single row "Av_tot_eff"
  eff_mat  <- res_dcdh$results$Effects      # rows Effect_1..Effect_4
  plac_mat <- res_dcdh$results$Placebos     # rows Placebo_1..Placebo_3

  # Helper: pull a named column from a results matrix/data.frame row by index.
  pull_col <- function(m, i, col) as.numeric(m[i, col])

  # Format a 95% CI as "[lb, ub]".
  fmt_ci <- function(lb, ub) sprintf("[%.3f, %.3f]", lb, ub)

  # --- Assemble the table rows: Av. total effect; Effect +1..+4; Placebo -1..-3.
  n_eff  <- nrow(eff_mat)
  n_plac <- nrow(plac_mat)

  est_vec <- c(
    pull_col(ate_row, 1L, "Estimate"),
    vapply(seq_len(n_eff),  function(i) pull_col(eff_mat,  i, "Estimate"), numeric(1L)),
    vapply(seq_len(n_plac), function(i) pull_col(plac_mat, i, "Estimate"), numeric(1L))
  )
  se_vec <- c(
    pull_col(ate_row, 1L, "SE"),
    vapply(seq_len(n_eff),  function(i) pull_col(eff_mat,  i, "SE"), numeric(1L)),
    vapply(seq_len(n_plac), function(i) pull_col(plac_mat, i, "SE"), numeric(1L))
  )
  lb_vec <- c(
    pull_col(ate_row, 1L, "LB CI"),
    vapply(seq_len(n_eff),  function(i) pull_col(eff_mat,  i, "LB CI"), numeric(1L)),
    vapply(seq_len(n_plac), function(i) pull_col(plac_mat, i, "LB CI"), numeric(1L))
  )
  ub_vec <- c(
    pull_col(ate_row, 1L, "UB CI"),
    vapply(seq_len(n_eff),  function(i) pull_col(eff_mat,  i, "UB CI"), numeric(1L)),
    vapply(seq_len(n_plac), function(i) pull_col(plac_mat, i, "UB CI"), numeric(1L))
  )
  n_vec <- c(
    pull_col(ate_row, 1L, "N"),
    vapply(seq_len(n_eff),  function(i) pull_col(eff_mat,  i, "N"), numeric(1L)),
    vapply(seq_len(n_plac), function(i) pull_col(plac_mat, i, "N"), numeric(1L))
  )

  row_labels_dcdh <- c(
    "Av. total effect",
    paste0("Effect $+", seq_len(n_eff), "$"),
    paste0("Placebo $-", seq_len(n_plac), "$")
  )

  dcdh_tab <- data.frame(
    ` `        = row_labels_dcdh,
    Estimate   = sprintf("%.3f", est_vec),
    SE         = sprintf("(%.3f)", se_vec),
    `95\\% CI` = mapply(fmt_ci, lb_vec, ub_vec),
    N          = format(round(n_vec), big.mark = ","),
    check.names = FALSE, stringsAsFactors = FALSE
  )

  # --- Joint placebo test returned by the package ------------------------
  # did_multiplegt_dyn() reports a joint test that all requested placebos are
  # zero. Field names have moved across package versions, so we search the
  # results list for the p-value and the statistic rather than hard-coding one
  # name, and report exactly what the package returns (no re-derivation).
  dcdh_res_names <- names(res_dcdh$results)
  message("  dCDH results fields: ", paste(dcdh_res_names, collapse = ", "))

  # Exactly one field must match, or we do not know what we are reporting:
  # silently taking the first of several matches is how a table ends up
  # labelling one statistic with another's name.
  pick_first <- function(nms, patterns) {
    hit <- nms[grepl(patterns, nms, ignore.case = TRUE)]
    if (length(hit) == 0L) return(NA_character_)
    if (length(hit) > 1L)
      stop("Ambiguous did_multiplegt_dyn result field: ",
           paste(hit, collapse = ", "), " all match /", patterns,
           "/. Update the field search in §A before reporting the joint ",
           "placebo test.")
    hit[1L]
  }
  nm_pjoint <- pick_first(dcdh_res_names, "p_?joint.*placebo|placebo.*p_?joint")
  nm_fjoint <- pick_first(dcdh_res_names, "^(F|chi|stat).*joint.*placebo|joint.*placebo.*(F|stat)")

  p_joint_plac <- if (!is.na(nm_pjoint))
    suppressWarnings(as.numeric(res_dcdh$results[[nm_pjoint]])[1L]) else NA_real_
  f_joint_plac <- if (!is.na(nm_fjoint))
    suppressWarnings(as.numeric(res_dcdh$results[[nm_fjoint]])[1L]) else NA_real_

  message(sprintf("  dCDH joint placebo test: field = %s | statistic = %s | p = %s",
                  ifelse(is.na(nm_pjoint), "not returned", nm_pjoint),
                  ifelse(is.na(f_joint_plac), "not returned", sprintf("%.4f", f_joint_plac)),
                  ifelse(is.na(p_joint_plac), "not returned", sprintf("%.4f", p_joint_plac))))

  dcdh_tab <- rbind(dcdh_tab, data.frame(
    ` `        = paste0("\\midrule Joint placebo test ($", n_plac, "$ placebos)"),
    Estimate   = if (is.na(f_joint_plac)) "---" else sprintf("%.3f", f_joint_plac),
    SE         = "---",
    `95\\% CI` = if (is.na(p_joint_plac)) "$p$ = ---" else sprintf("$p$ = %.3f", p_joint_plac),
    N          = "---",
    check.names = FALSE, stringsAsFactors = FALSE
  ))

  xtab_dcdh <- xtable(dcdh_tab, label = "tab:dcdh")
  align(xtab_dcdh) <- "llcccc"   # dropped row-name col + 5 body cols

  raw_dcdh <- capture.output(
    print(xtab_dcdh, include.rownames = FALSE, booktabs = TRUE,
          sanitize.text.function = identity, size = "\\small",
          floating = FALSE)
  )

  notes_dcdh <- paste0(
    "de Chaisemartin--D'Haultfoeuille (2024) estimator, current-treatment dummy; ",
    "never- and not-yet-treated controls; same controls as the CS main specification. ",
    "Last row: package's own joint test that all ", n_plac, " placebo estimates are ",
    "zero (statistic in Estimate, $p$-value in the CI column), accounting for the ",
    "covariance across placebo horizons; $p$-value taken from field \\texttt{",
    gsub("_", "\\\\_", ifelse(is.na(nm_pjoint), "not returned", nm_pjoint)),
    "} of \\texttt{did\\_multiplegt\\_dyn}'s results object (v",
    as.character(utils::packageVersion("DIDmultiplegtDYN")),
    "), which returns no accompanying statistic (Estimate column blank)"
  )
  source_dcdh <- "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"

  write_tex_float(
    out_path      = file.path(dir_tabs_dcdh, "att_dcdh.tex"),
    caption_title = "Robustness: de Chaisemartin and D'Haultfoeuille (2024) estimator",
    label         = "tab:dcdh",
    tabular_lines = raw_dcdh,
    notes_text    = notes_dcdh,
    source_text   = source_dcdh
  )

  message(sprintf("  dCDH Av. total effect = %.3f (SE %.3f, CI [%.3f, %.3f])",
                  est_vec[1L], se_vec[1L], lb_vec[1L], ub_vec[1L]))

  # --- §A event-study figure: built MANUALLY from Effects (+1..+4) and Placebos
  # (-1..-3), with a normalized point at event time 0 = 0.  Ribbon/linerange uses
  # the package 95% CI (LB CI / UB CI).  No in-figure title/subtitle/caption.
  es_dcdh <- data.frame(
    event_time = c(-rev(seq_len(n_plac)), 0L, seq_len(n_eff)),
    estimate   = c(rev(vapply(seq_len(n_plac), function(i) pull_col(plac_mat, i, "Estimate"), numeric(1L))),
                   0,
                   vapply(seq_len(n_eff), function(i) pull_col(eff_mat, i, "Estimate"), numeric(1L))),
    lb         = c(rev(vapply(seq_len(n_plac), function(i) pull_col(plac_mat, i, "LB CI"), numeric(1L))),
                   0,
                   vapply(seq_len(n_eff), function(i) pull_col(eff_mat, i, "LB CI"), numeric(1L))),
    ub         = c(rev(vapply(seq_len(n_plac), function(i) pull_col(plac_mat, i, "UB CI"), numeric(1L))),
                   0,
                   vapply(seq_len(n_eff), function(i) pull_col(eff_mat, i, "UB CI"), numeric(1L))),
    stringsAsFactors = FALSE
  )

  p_dcdh <- ggplot(es_dcdh, aes(x = event_time, y = estimate)) +
    geom_hline(yintercept = 0,    colour = "grey50", linetype = "dashed",
               linewidth = 0.5) +
    geom_vline(xintercept = -0.5, colour = "grey30", linetype = "dotted",
               linewidth = 0.5) +
    geom_linerange(aes(ymin = lb, ymax = ub),
                   linewidth = 0.6, alpha = 0.85, colour = "#2E86C1") +
    geom_line(colour = "#2E86C1", linewidth = 0.7) +
    geom_point(size = 2.5, colour = "#2E86C1") +
    scale_x_continuous(breaks = sort(unique(es_dcdh$event_time))) +
    labs(title = NULL, subtitle = NULL, caption = NULL,  # titles go in the LaTeX caption
         x = "Years relative to NAP adoption",
         y = "ATT on log(adaptation commitments)") +
    theme_minimal() +
    theme(text             = element_text(family = "serif", size = 12),
          panel.grid.minor = element_blank())

  ggsave(file.path(dir_figs_dcdh, "did_dcdh_es.png"),
         p_dcdh, width = 9, height = 5.5, dpi = 300)
  message("Saved: ", file.path(dir_figs_dcdh, "did_dcdh_es.png"))
}

message("\n=== Section A complete ===\n")

# ==============================================================================
# SECTION B. §28 Goodman-Bacon (2021) decomposition
# Decomposes the static two-way fixed-effects estimate of post_nap on log
# adaptation commitments into its 2x2 comparison-group weights, exposing the
# share carried by the problematic "Later vs Earlier Treated" forbidden
# comparison (already-treated units used as controls).  This is an APPENDIX
# exhibit.
#
# TWO traps (verified):
#  (a) simple_panel_wgi.csv already carries a column named `treated`; bacon()'s
#      internal merge collides with it and throws a misleading "Treatment not
#      weakly increasing with time".  We pass a CLEAN minimal frame containing
#      only (country_id, year, log_commits, post_nap) — NO other `treated`-named
#      column.
#  (b) bacon() needs a STRONGLY BALANCED panel — subset to countries observed in
#      all years first.
# Wrapped in tryCatch -> NULL/skip on failure.
# ==============================================================================

message("\n=== Section B (28): Goodman-Bacon decomposition ===\n")

dir_tabs_bacon <- file.path(here("output", "tables"),  "bacon")
dir_figs_bacon <- file.path(here("output", "figures"), "bacon")
dir.create(dir_tabs_bacon, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figs_bacon, recursive = TRUE, showWarnings = FALSE)

dpd <- did_panel_full %>% filter(!(cohort_year %in% thin_cohorts))

# (b) strongly-balanced subsample: keep only countries observed in all years.
ny      <- length(unique(dpd$year))
bal_ids <- (dpd %>% count(country_id) %>% filter(n == ny))$country_id

# (a) clean minimal frame — the treatment dummy is named post_nap (NOT `treated`).
clean_bacon <- dpd %>%
  filter(country_id %in% bal_ids) %>%
  transmute(country_id, year, log_commits,
            post_nap = as.integer(cohort_year > 0 & year >= cohort_year)) %>%
  filter(!is.na(log_commits)) %>%
  as.data.frame()

message(sprintf("Bacon panel: %d balanced countries x %d years = %d rows",
                length(bal_ids), ny, nrow(clean_bacon)))

bd <- tryCatch(
  bacon(log_commits ~ post_nap, data = clean_bacon,
        id_var = "country_id", time_var = "year"),
  error = function(e) { message("  bacon() failed: ", conditionMessage(e)); NULL }
)

if (is.null(bd)) {
  message("  Goodman-Bacon decomposition returned NULL — skipping table and figure.")
} else {

  # --- Aggregate by comparison type: total weight and weighted-mean estimate.
  # Compute the weighted mean safely as sum(estimate*weight)/sum(weight) within
  # each type (do NOT use weighted.mean with mismatched-length args).
  # NOTE: compute the weighted-mean estimate BEFORE collapsing `weight`, so that
  # sum(estimate * weight) uses the per-comparison weights (not the already-summed
  # scalar).  Reassigning `weight` first would silently corrupt the numerator.
  bacon_by_type <- bd %>%
    group_by(type) %>%
    summarise(
      estimate = sum(estimate * weight) / sum(weight),
      weight   = sum(weight),
      .groups  = "drop"
    )

  # Overall TWFE estimate = weighted sum across ALL 2x2 comparisons.
  overall_twfe <- sum(bd$estimate * bd$weight)

  message(sprintf("  Overall TWFE (Goodman-Bacon weighted) = %.4f", overall_twfe))
  for (i in seq_len(nrow(bacon_by_type))) {
    message(sprintf("    %-28s weight = %.4f  est = %.4f",
                    bacon_by_type$type[i],
                    bacon_by_type$weight[i], bacon_by_type$estimate[i]))
  }

  # --- §B table: one row per comparison type + final "Overall (TWFE)" row.
  # Map the package's type labels to readable comparison names; preserve a fixed
  # display order regardless of how bacon() orders them internally.
  type_order  <- c("Treated vs Untreated",
                   "Earlier vs Later Treated",
                   "Later vs Earlier Treated")
  bacon_by_type <- bacon_by_type %>%
    mutate(type = factor(type, levels = type_order)) %>%
    arrange(type)

  comp_labels <- as.character(bacon_by_type$type)
  weight_vals <- bacon_by_type$weight
  est_vals    <- bacon_by_type$estimate

  bacon_tab <- data.frame(
    Comparison        = c(comp_labels, "Overall (TWFE)"),
    Weight            = c(sprintf("%.4f", weight_vals), sprintf("%.4f", sum(weight_vals))),
    `Avg. estimate`   = c(sprintf("%.4f", est_vals),    sprintf("%.4f", overall_twfe)),
    check.names = FALSE, stringsAsFactors = FALSE
  )

  xtab_bacon <- xtable(bacon_tab, label = "tab:bacon")
  align(xtab_bacon) <- "llcc"   # dropped row-name col + 3 body cols

  raw_bacon <- capture.output(
    print(xtab_bacon, include.rownames = FALSE, booktabs = TRUE,
          sanitize.text.function = identity, size = "\\small",
          floating = FALSE)
  )

  notes_bacon <- paste0(
    "Decomposition of the static two-way fixed-effects estimate of \\texttt{post\\_nap} ",
    "(an indicator for years after NAP adoption) on log adaptation commitments into ",
    "$2\\times2$ comparison groups (Goodman-Bacon 2021); the weights sum to one. ",
    "The ``Later vs Earlier Treated'' group is the problematic forbidden comparison ",
    "that uses already-treated units as controls. ",
    "Estimated on the balanced-panel subsample"
  )
  source_bacon <- "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"

  write_tex_float(
    out_path      = file.path(dir_tabs_bacon, "bacon_decomp.tex"),
    caption_title = "Goodman-Bacon decomposition of the two-way fixed-effects estimator",
    label         = "tab:bacon",
    tabular_lines = raw_bacon,
    notes_text    = notes_bacon,
    source_text   = source_bacon
  )

  # --- §B scatter: weight (x) vs estimate (y), coloured/shaped by comparison type,
  # with a dashed horizontal line at the overall TWFE weighted estimate.
  # No in-figure title; serif font; legend at bottom.
  bd_plot <- bd %>%
    mutate(type = factor(type, levels = type_order))

  p_bacon <- ggplot(bd_plot, aes(x = weight, y = estimate,
                                 colour = type, shape = type)) +
    geom_hline(yintercept = overall_twfe, colour = "grey30",
               linetype = "dashed", linewidth = 0.5) +
    geom_point(size = 2.6, alpha = 0.85) +
    scale_colour_brewer(palette = "Set2") +
    labs(title = NULL, subtitle = NULL, caption = NULL,  # titles go in the LaTeX caption
         x = "Weight", y = "2x2 DID estimate",
         colour = NULL, shape = NULL) +
    theme_minimal() +
    theme(text             = element_text(family = "serif", size = 12),
          legend.position  = "bottom",
          panel.grid.minor = element_blank())

  ggsave(file.path(dir_figs_bacon, "bacon_scatter.png"),
         p_bacon, width = 9, height = 5.5, dpi = 300)
  message("Saved: ", file.path(dir_figs_bacon, "bacon_scatter.png"))
}

message("\n=== Section B complete ===\n")

# ==============================================================================
# SECTION C. Outlier robustness: excluding India
# India is the largest single recipient of adaptation-related commitments
# (never-treated control in the design), so it could dominate the estimated
# counterfactual. Re-run the MAIN specification (doubly-robust, multiplier-
# bootstrap SE, never-treated control, cohorts >= 5 treated units) on the
# panel excluding India.
# ==============================================================================

message("\n=== Section C: Outlier robustness — excluding India ===\n")

# Executable precondition: India must be never-treated (cohort_year == 0), so
# dropping it removes a control, not a treated unit (interpretation depends on it).
stopifnot(all(did_panel_full$cohort_year[
  did_panel_full$recipient_name == "India"] == 0))

# Main-spec panel: thin cohorts filtered OUT (make_wide_table does not filter
# internally — retain_thin only switches estimator/SE; mirrors §19d / 03 main).
did_panel_noindia <- did_panel_full %>%
  filter(!(cohort_year %in% thin_cohorts), recipient_name != "India")
stopifnot(n_distinct(did_panel_noindia$recipient_name) ==
            n_distinct(did_panel_full$recipient_name[
              !(did_panel_full$cohort_year %in% thin_cohorts)]) - 1L)

make_wide_table(
  did_panel_in  = did_panel_noindia,
  retain_thin   = FALSE,
  outcomes      = outcomes,
  dir_tabs      = file.path(here("output", "tables"), "outlier_india"),
  tex_label     = "tab:combined_wide_noindia",
  caption_spec  = "main specification excluding India (largest recipient)"
)

message("\n=== Section C complete ===\n")
message("\n=== 04_robustness.R: complete ===\n")
