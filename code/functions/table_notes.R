# ==============================================================================
# code/functions/table_notes.R
# The single source of the table-note wording shared across stages 03, 04, 05,
# 07, 09, 10, 11, 13 and 14. Table 2 (tab:combined_wide_main) states the
# estimator, covariates, control group, inference, pre-trend tests, pre-treatment mean, implied
# effect and significance stars in full; every other table refers to it with
# as_in_main() for the parts it shares verbatim and writes out only what
# differs (estimator, control group, SE type, pre-trend window, counts). The
# Source line is written in full on every table (SOURCE_MAIN).
# Base R only.
# ==============================================================================

MAIN_TAB_LABEL <- "tab:combined_wide_main"

# Pre-trend window of Table 2 (leads -5..-2, 4 restrictions). make_wide_table()
# checks Table 2 against it; pretrend_window_note() states any other window.
MAIN_PRETREND_WINDOW <- c(min_e = -5L, max_e = -2L, k = 4L)

PRETREND_NOTE_AGG <- function(min_e, max_e, k) sprintf(
  "Pre-trend $\\chi^2$ / $p$: joint Wald test on the aggregated pre-treatment event-time coefficients ($%d \\leq e \\leq %d$; %d restrictions), using the influence-function covariance of the dynamic aggregation from the analytical (non-bootstrap) fit; a generalized inverse is used if the block is singular",
  min_e, max_e, k)
PRETREND_NOTE_DID <- "\\texttt{did} pre-test $p$: \\texttt{did}'s built-in Wald test over all pre-period $ATT(g,t)$ cells against each cohort's $g-1$ base year"
STAR_NOTE   <- "* $p<0.10$, ** $p<0.05$, *** $p<0.01$"
SOURCE_MAIN <- "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"
# CRS data vintage, in the two forms the notes print (Table 2 spec line; stage 14
# source lines). Kept without parentheses so paste0("(", CRS_VINTAGE, ")") never nests.
CRS_VINTAGE      <- "April 2026 vintage"
CRS_VINTAGE_ABBR <- "Apr.\\ 2026"

#' Reference phrase to Table 2 for the note elements a table shares with it
#'
#' @param parts note elements, e.g. c("estimator", "covariates", "control group")
#' @param capitalize FALSE inside a sentence
#' @return "Estimator, controls and inference as in Table~\ref{tab:combined_wide_main}"
#'   (no trailing period)
as_in_main <- function(parts, capitalize = TRUE) {
  n <- length(parts)
  stopifnot(n >= 1L)
  s <- if (n == 1L) parts else
    paste0(paste(parts[-n], collapse = ", "), " and ", parts[n])
  if (capitalize) s <- paste0(toupper(substr(s, 1L, 1L)), substring(s, 2L))
  paste0(s, " as in Table~\\ref{", MAIN_TAB_LABEL, "}")
}

#' The pre-trend window, stated only when it differs from Table 2's
#'
#' @return "" for Table 2's window, else "Pre-trend window: ... . "
pretrend_window_note <- function(min_e, max_e, k) {
  if (identical(as.integer(c(min_e, max_e, k)), unname(MAIN_PRETREND_WINDOW))) return("")
  sprintf("Pre-trend window: $%d \\leq e \\leq %d$ (%d restrictions). ",
          as.integer(min_e), as.integer(max_e), as.integer(k))
}

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
#' @param full TRUE for Table 2, which defines both tests in full; FALSE (other
#'   tables, which refer to Table 2 for the definitions) states only the lead
#'   window when it differs from Table 2's, singular blocks, counts and verdicts
#' @return a character string for the table note
wpval_reconciliation <- function(pre_egt, df_did, n_clusters,
                                 wpval_did = NA_real_, pval_wald = NA_real_,
                                 labels = NULL, wpval_reason = NA_character_,
                                 ginv_used = FALSE, df = length(pre_egt),
                                 full = FALSE) {
  # Compact pre-trend disclosure for table notes. The averaging mechanism and the interpretation of each
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

  head_txt <- if (n_lead == 0L) {
    "No aggregated pre-trend test: no pre-treatment event time has a usable SE. "
  } else {
    # The lead window is disclosed as min/max of pre_egt.
    # Restrictions = the test's df (the rank actually inverted), not the lead
    # count; if the df differ across columns the lead count is quoted and the
    # singular columns are flagged just below.
    df_u <- unique(as.integer(df[!is.na(df)]))
    k    <- if (length(df_u) == 1L) df_u else n_lead
    ginv_txt <- if (isTRUE(any(ginv_used)))
      paste0("the block was singular",
             if (!is.null(labels) && length(labels) == length(ginv_used))
               paste0(" for ", paste(labels[ginv_used], collapse = ", ")) else "",
             " and a generalized inverse was used") else ""
    if (full) {
      paste0(PRETREND_NOTE_AGG(min(as.integer(pre_egt)), max(as.integer(pre_egt)), k),
             if (nzchar(ginv_txt)) paste0("; ", ginv_txt) else "", ". ")
    } else {
      paste0(pretrend_window_note(min(as.integer(pre_egt)), max(as.integer(pre_egt)), k),
             if (nzchar(ginv_txt)) paste0("Pre-trend test: ", ginv_txt, ". ") else "")
    }
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
    if (full) PRETREND_NOTE_DID else "\\texttt{did} pre-test $p$", " (",
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
