# ==============================================================================
# 03_main_results.R
# Main CS (2021) estimation — cohorts_dropped specification.
# Paper: Beaucoral, Goujon and Marchand (2026) — §5 (main
# results: tab:combined_wide_main, fig:did_combined_es, fig:did_combined_cohort)
#
# Inputs : data/processed/simple_panel_wgi.csv
# Outputs:
#   output/figures/cohorts_dropped/did_combined_es_wgi.png     (Figure 2, fig:did_combined_es)
#   output/figures/cohorts_dropped/did_combined_cohort_wgi.png (Figure C.3, fig:did_combined_cohort)
#   output/tables/cohorts_dropped/att_combined_wide.tex        (Table 2, tab:combined_wide_main)
#   output/tables/extensive_margin/att_extensive.tex           (tab:extensive_margin)
#   output/tables/nap_cohorts.tex                              (tab:nap_cohorts_app)
#   output/fits/headline_*_dr_bs.rds, extensive_margin_dr_bs.rds (read by later stages)
# Figures 2 and C.3 are drawn from the stored fits (Section 10b), never refitted.
# ==============================================================================

# ============================================================
# Paper-to-Code Naming Map
# ============================================================
# Paper Notation      | Code Name             | Description
# $Y_{it}$            | log_commits, etc.     | Five outcome variables
# $G_i$               | cohort_year           | Year of first NAP adoption (0 = never)
# $ATT(g,t)$          | gt_obj                | Group-time ATT from att_gt()
# $\hat\theta^{simp}$ | agg_s$overall.att     | Calendar-time simple ATT
# $\hat\theta^{dyn}$  | agg_d$att.egt         | Dynamic ATT by event time
# $X_{it}$            | ge_est, log_population| Controls: WGI GE + log pop.
# ============================================================

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
source(here("code", "functions", "pretrend_test.R"))  # compute_pretrend_test()
source(here("code", "functions", "two_line_head.R"))     # two_line_head()
source(here("code", "functions", "mde.R"))            # mde()
source(here("code", "functions", "make_country_id.R"))    # make_country_id()
source(here("code", "functions", "crs_positive.R"))       # crs_positive()
source(here("code", "functions", "make_wide_table.R"))    # make_wide_table()
source(here("code", "functions", "sup_t_crit.R"))         # sup_t_crit()

set.seed(20240601)  # global seed — local set.seed(1242) calls follow each estimator

# -----------------------------------------------------------------------
# One fit, one SE: single bootstrap-replication constant.
# Every multiplier-bootstrap call in this script uses BITERS. 999 is kept
# because sharing the saved fit across scripts already removes
# the 0.1245 / 0.1173 duplicate-SE problem, and raising the replication count
# would re-randomise every standard error in the paper for no inferential gain.
# -----------------------------------------------------------------------
BITERS <- 999L

# -----------------------------------------------------------------------
# Output directories
# -----------------------------------------------------------------------
dir.create(here("output", "figures", "cohorts_dropped"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "tables",  "cohorts_dropped"), recursive = TRUE, showWarnings = FALSE)
# Saved headline fits (one .rds per main outcome) consumed by later stages.
dir.create(here("output", "fits"), recursive = TRUE, showWarnings = FALSE)
# A failed run must never leave the previous run's fits for later stages to
# read: the fits this stage owns are deleted before they are re-estimated.
unlink(Sys.glob(here("output", "fits",
                     c("headline_*_dr_bs.rds", "extensive_margin_dr_bs.rds"))))
dir.create(here("output", "tables", "extensive_margin"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------
# Escape LaTeX-special characters in TABLE COLUMN HEADERS only.
# (Needed for xtable column names.)
# -----------------------------------------------------------------------
esc_header <- function(x) {
  x <- gsub("%", "\\\\%", x)
  x <- gsub("_", "\\\\_", x)
  x <- gsub("#", "\\\\#", x)
  x
}

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

# ==============================================================================
# Helper: write_tex_float()
# Wraps a bare tabular block in the project's standard complete float:
#   \begin{table}[H] ... \adjustbox ... \minipage{Notes/Source} \end{table}
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
# SECTION 1. Load and prepare the DiD panel
# ==============================================================================

message("\n=== 03_main_results.R: loading panel ===\n")

aggregated <- fread(here("data", "processed", "simple_panel_wgi.csv"))
aggregated  <- as.data.frame(aggregated)

# Convert to plain data.frame (did package does not handle data.table well)
did_panel <- aggregated

# Numeric unit ID required by att_gt()
did_panel <- did_panel %>%
  mutate(country_id = make_country_id(recipient_name))

# Build gname (cohort = year of first NAP adoption).
# att_gt() requires gname to be CONSTANT within each unit across all years.
# The panel window is fixed at 2009 (first year of Rio adaptation-marker
# reporting); 01_prepare_data.R filters and asserts it, and this re-checks it.
first_year <- 2009L
stopifnot(min(did_panel$year) == first_year)

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

# Merge share_adapt if not already present
if (!"share_adapt" %in% names(did_panel)) {
  share_lookup <- aggregated %>%
    select(recipient_name, year, share_adapt) %>%
    distinct()
  did_panel <- left_join(did_panel, share_lookup, by = c("recipient_name", "year"))
}

# Preserve full panel before thin-cohort filtering (used in Section 3 below)
did_panel_full <- did_panel

# ==============================================================================
# SECTION 2. Sample description
# ==============================================================================

message("\n=== FINAL DiD SAMPLE ===\n")
message(sprintf("Observations : %d", nrow(did_panel)))
message(sprintf("Countries    : %d", n_distinct(did_panel$recipient_name)))
message(sprintf("Years        : %d  (%d – %d)",
                n_distinct(did_panel$year),
                min(did_panel$year),
                max(did_panel$year)))

cohort_summary <- did_panel %>%
  distinct(recipient_name, cohort_year, always_treated) %>%
  mutate(group = case_when(
    always_treated   ~ "Always treated (used as never-treated control)",
    cohort_year == 0 ~ "Never treated",
    TRUE             ~ paste0("Cohort ", cohort_year)
  )) %>%
  count(group) %>%
  arrange(group)
message("\nCountries by treatment group:")
print(cohort_summary, row.names = FALSE)

message(sprintf("\nRemaining NAs — ge_est: %d | log_population: %d",
                sum(is.na(did_panel$ge_est)),
                sum(is.na(did_panel$log_population))))

region_summary <- did_panel %>%
  distinct(recipient_name, WB_region) %>%
  count(WB_region) %>%
  arrange(desc(n))
message("\nCountries by World Bank region:")
print(region_summary, row.names = FALSE)
message("\n=== END OF SAMPLE DESCRIPTION ===\n")

# Cohort sizes (diagnostic)
cohort_sizes <- did_panel %>%
  filter(cohort_year > 0) %>%
  distinct(recipient_name, cohort_year) %>%
  count(cohort_year, name = "n_treated") %>%
  arrange(cohort_year)
message("\n--- Cohort sizes (treated countries per adoption year) ---")
print(cohort_sizes, row.names = FALSE)

thin_threshold <- 5L
thin_cohorts   <- cohort_sizes$cohort_year[cohort_sizes$n_treated < thin_threshold]

# --- Appendix table: NAP adoption cohorts (tab:nap_cohorts_app) ---------------
# Documents the treated-cohort composition and the >= thin_threshold inclusion
# rule used by the main specification. Single source of truth: cohort_sizes +
# thin_threshold above. Exhibit: paper/Tables/nap_cohorts (\input).
n_adopters <- sum(cohort_sizes$n_treated)
n_never    <- did_panel %>% filter(cohort_year == 0) %>% distinct(recipient_name) %>% nrow()
coh_tab    <- cohort_sizes %>%
  mutate(in_main = ifelse(n_treated >= thin_threshold, "Yes", "No"))
yr_span <- function(y) paste(range(y), collapse = "--")  # e.g. "2021--2024"
nap_cohort_tabular <- c(
  "\\begin{tabular}{lcc}",
  "\\toprule",
  "Adoption year & Number of adopters & Included in main specification \\\\",
  "\\midrule",
  paste0(coh_tab$cohort_year, " & ", coh_tab$n_treated, " & ", coh_tab$in_main, " \\\\"),
  "\\midrule",
  paste0("Total adopters & ", n_adopters, " & \\\\"),
  "\\bottomrule",
  "\\end{tabular}"
)
write_tex_float(
  out_path      = here("output", "tables", "nap_cohorts.tex"),
  caption_title = "NAP adoption cohorts",
  label         = "tab:nap_cohorts_app",
  tabular_lines = nap_cohort_tabular,
  notes_text    = paste0(
    "The table reports the number of countries first submitting a National ",
    "Adaptation Plan in each adoption year, among the ", n_adopters,
    " adopters in the estimation sample; a further ", n_never,
    " never-adopting countries serve as controls. The main specification ",
    "retains adoption cohorts with at least ", thin_threshold,
    " treated units (", yr_span(coh_tab$cohort_year[coh_tab$in_main == "Yes"]),
    "); the retained-cohorts robustness specification ",
    "(Section~\\ref{sec:robust}) additionally includes the smaller ",
    yr_span(coh_tab$cohort_year[coh_tab$in_main == "No"]), " cohorts"),
  source_text   = "UNFCCC NAP Central tracking tool"
)

# ==============================================================================
# SECTION 4. Outcomes list
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
# SECTION 7b. Helpers — save the headline fit so that every downstream
# script reports ONE ATT and ONE SE per specification.
#
# Why: `aggte()` runs its own multiplier bootstrap when the att_gt object was
# fitted with bstrap = TRUE, so an aggregated SE depends on the RNG state at
# the aggte() call, not only at the att_gt() call. Two scripts that fit the
# same specification with the same seed can therefore report different SEs
# (for example, an analytical fit between att_gt() and aggte() leaves a
# different RNG state). Re-seeding before att_gt() does not pin the reported
# SE; fitting and aggregating once here, saving the result, and having every
# other script read the saved aggregation does.
# ==============================================================================

#' Event-time weights implied by did's `simple` aggregation.
#'
#' did's simple ATT is the group-size-weighted average of the post-treatment
#' ATT(g,t) cells, so it equals sum_e w_e * att.egt(e) with
#'   w_e proportional to the number of treated units in the cohorts that
#'   contribute a post-treatment cell at event time e.
#' The returned vector is the l_vec that makes a HonestDiD sensitivity analysis
#' target exactly the ATT printed in Table 2. The reconciliation against
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

#' Multiplier-bootstrap covariance of a dynamic aggregation.
#'
#' `aggte()` reports a multiplier-bootstrap standard error but discards the
#' bootstrap covariance matrix, and did's bootstrap SE is an interquartile-range
#' scale estimate (`mboot`: se = bSigma * sqrt(n_clusters) / n), not
#' sqrt(diag(cov(bres))/n). HonestDiD needs a full covariance. We therefore
#'   (i)  re-run did::mboot() on the SAME dynamic influence-function matrix the
#'        aggregation used, with an explicit seed, to obtain the bootstrap
#'        covariance `V` (V = cov(bres), bres already scaled by sqrt(n));
#'   (ii) take its correlation matrix R = cov2cor(V/n); and
#'   (iii) rescale it by the SEs that the paper actually reports:
#'        Sigma = diag(se.egt) %*% R %*% diag(se.egt).
#' The diagonal of Sigma is then EXACTLY the published event-study SE vector,
#' and the off-diagonal correlations are the bootstrap ones. The function also
#' returns the raw covariance-based SEs so the caller can report how far the
#' IQR-based published SEs sit from the covariance-based ones.
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

# Stable file stems for the saved fits (paper outcome -> file name).
fit_stem_map <- c(
  log_commits           = "adaptation",
  share_adapt           = "share",
  lcommitments_all      = "total",
  lcommitments_nonadapt = "nonadaptation",
  ldisbursements        = "disbursements"
)

# ==============================================================================
# SECTION 8. Helper: store_headline_fit()
# Table 2 is built by make_wide_table() (code/functions/make_wide_table.R),
# which calls this function once per outcome, after every reported quantity of
# that column has been computed, to store the fit behind the column.
# ==============================================================================

#' Store the fit behind one column of Table 2 in output/fits/
#'
#' Adds the bootstrap dynamic and group aggregations (Figure 2, Figure C.3,
#' HonestDiD), the event-time weights that reproduce the simple ATT, and the
#' bootstrap covariance of the dynamic aggregation. The next column re-seeds
#' with set.seed(1242) before its own bootstrap fit, so the draws made here do
#' not reach it.
#'
#' @param oc outcome entry of `outcomes`
#' @param col the column's fits and statistics, as passed by make_wide_table()
#' @return invisibly, the path of the stored fit (NULL for an outcome without a
#'   stem in fit_stem_map)
store_headline_fit <- function(oc, col) {
  stem <- fit_stem_map[oc$var]
  if (is.na(stem)) return(invisible(NULL))
  if (!is.finite(col$att) || !is.finite(col$se))
    stop("Headline simple aggregation failed for ", oc$var, ": no stored fit is written.")

  agg_dyn_boot <- aggte(col$gt_boot, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf)
  agg_grp_boot <- aggte(col$gt_boot, type = "group", na.rm = TRUE)

  # l_vec matching the ATT printed in Table 2 (group-size-weighted average of
  # the post-treatment event-time effects), reconciled against the simple ATT.
  l_vec_simple <- NULL
  l_vec_check  <- NA_real_
  keep_d   <- which(!is.na(agg_dyn_boot$se.egt) & agg_dyn_boot$se.egt > 1e-10)
  egt_post <- agg_dyn_boot$egt[keep_d][agg_dyn_boot$egt[keep_d] >= 0]
  if (length(egt_post) > 0) {
    l_vec_simple <- simple_att_event_weights(col$gt_boot, col$panel, egt_post)
    att_post     <- agg_dyn_boot$att.egt[match(egt_post, agg_dyn_boot$egt)]
    l_vec_check  <- sum(l_vec_simple * att_post) - col$att
    message(sprintf(
      "    [fit] l_vec = (%s) reproduces the simple ATT to %.2e",
      paste(sprintf("%.4f", l_vec_simple), collapse = ", "), abs(l_vec_check)))
  }

  # Bootstrap covariance of the dynamic aggregation (feeds HonestDiD in 04).
  bvc <- boot_vcov_dynamic(agg_dyn_boot, col$gt_boot, seed = 1242L)
  message(sprintf(paste0("    [fit] bootstrap vcov: max |se_cov - se_published| / ",
                         "se_published = %.4f"), bvc$max_rel_dev))

  fit_obj <- list(
    outcome        = oc$var,
    outcome_label  = oc$label,
    spec           = list(
      estimator      = "CS (2021)",
      est_method     = col$est_method,
      control_group  = "nevertreated",
      xformla        = "~ ge_est + log_population",
      anticipation   = 0L,
      base_period    = "universal",
      cohort_rule    = "adoption cohorts with >= 5 treated units",
      bstrap         = TRUE,
      biters         = BITERS,
      seed           = 1242L,
      min_e          = -5,
      max_e          = Inf
    ),
    gt_boot          = col$gt_boot,
    gt_analytical    = col$gt_analytical,
    agg_simple       = col$agg_simple,
    agg_simple_analytic = col$agg_simple_analytic,
    agg_dynamic      = agg_dyn_boot,
    agg_group        = agg_grp_boot,
    agg_dyn_analytic = col$agg_dyn_analytic,
    ids              = sort(unique(col$panel$country_id)),
    att              = col$att,
    se               = col$se,
    pretrend         = col$pretrend,
    l_vec_simple     = l_vec_simple,
    l_vec_egt        = egt_post,
    l_vec_recon_gap  = l_vec_check,
    boot_vcov        = bvc,
    n_obs            = col$n_obs,
    n_country        = col$n_country,
    provenance       = list(
      script       = "code/03_main_results.R",
      created      = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      r_version    = R.version.string,
      did_version  = as.character(utils::packageVersion("did"))
    )
  )
  fit_path <- here("output", "fits", paste0("headline_", stem, "_dr_bs.rds"))
  saveRDS(fit_obj, fit_path)
  message("    [fit] Saved: ", fit_path)
  invisible(fit_path)
}

# ==============================================================================
# SECTION 9. Helper: make_cohort_plot() — Figure C.3 (called in Section 10b)
# ==============================================================================

make_cohort_plot <- function(results_group, outcomes, dir_figs, spec_label) {

  group_all <- bind_rows(results_group) %>%
    mutate(
      outcome = factor(outcome, levels = vapply(outcomes, `[[`, character(1L), "label")),
      cohort  = as.integer(cohort)
    ) %>%
    filter(!is.na(ATT))

  if (nrow(group_all) == 0) stop("No cohort ATT data for: ", spec_label)

  palette_vec <- setNames(vapply(outcomes, `[[`, character(1L), "color"),
                           vapply(outcomes, `[[`, character(1L), "label"))
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
# SECTION 10. §19c combined wide table — cohorts_DROPPED half
# Seed rule (reproduces the published SEs): set.seed(1242) in make_wide_table before estimator loop.
# ==============================================================================

did_panel_tab_main <- did_panel_full %>% filter(!(cohort_year %in% thin_cohorts))
make_wide_table(
  did_panel_in  = did_panel_tab_main,
  retain_thin   = FALSE,
  outcomes      = outcomes,
  dir_tabs      = file.path(here("output", "tables"), "cohorts_dropped"),
  tex_label     = "tab:combined_wide_main",
  caption_spec  = "cohorts $\\geq 5$ units, DR + multiplier-bootstrap SE",
  after_column  = store_headline_fit
)

# ==============================================================================
# SECTION 10b. Figures 2 and C.3, drawn from the stored fits
# Both figures read the aggregations store_headline_fit() saved in output/fits/.
# Figure 2 uses the stored agg_dynamic draw, the same draw behind the
# event-time SEs quoted in the text and the HonestDiD inputs; Figure C.3 uses
# the stored agg_group draw. (Table 2's SE is the agg_simple draw.) There is
# no second aggte() call: aggte() re-runs the bootstrap, so a second call on
# the same fit reports different SEs. The intervals are simultaneous (sup-t)
# 95% bands: ATT +/- crit x stored SE, with crit from sup_t_crit() on the stored
# influence functions (what cband = TRUE would give, up to the RNG draw). One
# band family per outcome curve: uniform over that outcome's event times
# (Figure 2) or cohorts (Figure C.3), not jointly over the five outcomes.
# ==============================================================================

message("\n=== Section 10b: Figures 2 and C.3 from the stored fits ===\n")

dir_figs_main   <- here("output", "figures", "cohorts_dropped")
results_dynamic <- vector("list", length(outcomes))
results_group   <- vector("list", length(outcomes))
for (i in seq_along(outcomes)) {
  oc     <- outcomes[[i]]
  stored <- readRDS(here("output", "fits",
                         paste0("headline_", fit_stem_map[[oc$var]], "_dr_bs.rds")))
  d <- stored$agg_dynamic
  g <- stored$agg_group
  stopifnot(!is.null(d), !is.null(g))
  cv_d <- sup_t_crit(d$inf.function$dynamic.inf.func.e, d$se.egt, biters = BITERS)
  cv_g <- sup_t_crit(g$inf.function$selective.inf.func.g, g$se.egt, biters = BITERS)
  message(sprintf("  Fig 2 sup-t crit (%s): %.4f | Fig C.3 sup-t crit (%s): %.4f",
                  oc$label, cv_d, oc$label, cv_g))
  results_dynamic[[i]] <- data.frame(
    outcome = oc$label, event_time = d$egt, ATT = d$att.egt,
    Lower = d$att.egt - cv_d * d$se.egt,
    Upper = d$att.egt + cv_d * d$se.egt)
  results_group[[i]] <- data.frame(
    outcome = oc$label, cohort = g$egt, ATT = g$att.egt,
    Lower = g$att.egt - cv_g * g$se.egt,
    Upper = g$att.egt + cv_g * g$se.egt)
  # Cohort-level ATTs quoted in the text (Figure C.3 plots them).
  message(sprintf("  Cohort ATTs [%s]: %s", oc$label,
                  paste(sprintf("%d %.4f (%.4f)", as.integer(g$egt), g$att.egt, g$se.egt),
                        collapse = "; ")))
}

# --- Figure 2: combined event-study overlay ---
outcome_levels <- vapply(outcomes, `[[`, character(1L), "label")
palette_vec    <- setNames(vapply(outcomes, `[[`, character(1L), "color"), outcome_levels)
dynamic_all    <- bind_rows(results_dynamic) %>%
  mutate(outcome = factor(outcome, levels = outcome_levels))
dodge_w <- 0.4

p_combined <- ggplot(dynamic_all,
                     aes(x = event_time, y = ATT,
                         colour = outcome, shape = outcome, group = outcome)) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
  geom_vline(xintercept = -0.5, colour = "grey30", linetype = "dotted") +
  geom_linerange(aes(ymin = Lower, ymax = Upper),
                 position = position_dodge(width = dodge_w),
                 linewidth = 0.6, alpha = 0.8, na.rm = TRUE) +
  geom_point(size = 2.5, position = position_dodge(width = dodge_w)) +
  scale_x_continuous(breaks = seq(min(dynamic_all$event_time), max(dynamic_all$event_time))) +
  scale_colour_manual(values = palette_vec) +
  scale_shape_manual(values = c(16, 17, 15, 18, 8)) +
  # No title, subtitle, or caption — those go in LaTeX \caption{}
  labs(title = NULL, subtitle = NULL, caption = NULL,
       x = "Years relative to NAP adoption", y = "ATT estimate",
       colour = NULL, shape = NULL) +
  theme_minimal() +
  theme(text             = element_text(family = "serif", size = 12),
        legend.position  = "bottom",
        legend.text      = element_text(size = 10),
        panel.grid.minor = element_blank())
ggsave(file.path(dir_figs_main, "did_combined_es_wgi.png"),
       p_combined, width = 12, height = 7, dpi = 300)
message("Saved: ", file.path(dir_figs_main, "did_combined_es_wgi.png"))

# --- Figure C.3: ATT by adoption cohort ---
make_cohort_plot(
  results_group = results_group,
  outcomes      = outcomes,
  dir_figs      = dir_figs_main,
  spec_label    = "Cohorts >= 5 units (stored bootstrap fits)"
)

# ==============================================================================
# SECTION 11. Extensive margin: 1{adaptation commitments > 0}
# Does NAP adoption move the probability of receiving ANY
# adaptation-marked commitment, as opposed to the (log) amount conditional on
# receiving something. Identical main specification: CS (2021) doubly robust,
# never-treated controls, WGI gov. effectiveness + log population, universal
# base period, no anticipation, cohorts with >= 5 treated units, multiplier
# bootstrap with BITERS replications. Pooled (simple) and dynamic ATTs are
# reported in one table.
# Seed rule (reproduces the published SEs): set.seed(1242) immediately before each att_gt() call.
# ==============================================================================

message("\n=== Section 11: extensive-margin ATT ===\n")

did_panel_ext <- did_panel_tab_main %>%
  mutate(any_adapt = if_else(is.na(commitments), NA_integer_,
                             as.integer(crs_positive(commitments))))

n_ext_obs <- sum(!is.na(did_panel_ext$any_adapt))
message(sprintf("  Extensive-margin panel: %d rows with a non-missing indicator; ",
                n_ext_obs),
        sprintf("mean(any_adapt) = %.4f",
                mean(did_panel_ext$any_adapt, na.rm = TRUE)))

set.seed(1242)
gt_ext_analytical <- tryCatch(
  att_gt(yname = "any_adapt", tname = "year", idname = "country_id",
         gname = "cohort_year", xformla = ~ ge_est + log_population,
         data = did_panel_ext, est_method = "dr", bstrap = FALSE, cband = FALSE,
         control_group = "nevertreated", anticipation = 0,
         base_period = "universal", panel = TRUE, allow_unbalanced_panel = TRUE),
  error = function(e) { message("  att_gt (analytical) failed: ",
                                conditionMessage(e)); NULL })

set.seed(1242)
gt_ext <- tryCatch(
  att_gt(yname = "any_adapt", tname = "year", idname = "country_id",
         gname = "cohort_year", xformla = ~ ge_est + log_population,
         data = did_panel_ext, est_method = "dr", bstrap = TRUE, biters = BITERS,
         cband = FALSE, control_group = "nevertreated", anticipation = 0,
         base_period = "universal", panel = TRUE, allow_unbalanced_panel = TRUE),
  error = function(e) { message("  att_gt (bootstrap) failed: ",
                                conditionMessage(e)); NULL })

if (is.null(gt_ext) || is.null(gt_ext_analytical)) {
  stop("Extensive-margin estimation failed: no table written.")
} else {
  agg_s_ext <- aggte(gt_ext, type = "simple",  na.rm = TRUE)
  agg_d_ext <- aggte(gt_ext, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf)
  agg_d_ext_analytical <- aggte(gt_ext_analytical, type = "dynamic", na.rm = TRUE,
                                min_e = -5, max_e = Inf)

  att_ext <- agg_s_ext$overall.att
  se_ext  <- agg_s_ext$overall.se
  t_ext   <- att_ext / se_ext
  stars_ext <- if (abs(t_ext) > 2.576) "***" else if (abs(t_ext) > 1.960) "**" else
               if (abs(t_ext) > 1.645) "*" else ""
  pt_ext <- compute_pretrend_test(agg_d_ext_analytical, gt_ext_analytical)

  # MDE at 80% power, 5% two-sided.
  mde_ext <- mde(se_ext)
  message(sprintf("  Extensive-margin ATT = %.4f%s (SE = %.4f, t = %.3f) | MDE = %.4f",
                  att_ext, stars_ext, se_ext, t_ext, mde_ext))

  keep_ext <- which(!is.na(agg_d_ext$se.egt) & agg_d_ext$se.egt > 1e-10)
  dyn_rows <- vapply(keep_ext, function(i) sprintf(
    "\\quad $e = %d$ & %.4f & (%.4f) \\\\",
    as.integer(agg_d_ext$egt[i]), agg_d_ext$att.egt[i], agg_d_ext$se.egt[i]),
    character(1L))

  pre_mean_ext <- did_panel_ext %>%
    filter(cohort_year > 0, year < cohort_year, !is.na(any_adapt)) %>%
    summarise(m = mean(any_adapt)) %>% pull(m)

  n_cty_ext <- length(unique(did_panel_ext$country_id[!is.na(did_panel_ext$any_adapt)]))

  ext_tabular <- c(
    "\\begin{tabular}{lcc}",
    "\\toprule",
    " & ATT & SE \\\\",
    "\\midrule",
    "\\multicolumn{3}{l}{\\textit{Panel A: pooled}} \\\\",
    sprintf("Simple ATT & %.4f%s & (%.4f) \\\\", att_ext, stars_ext, se_ext),
    sprintf("$t$-statistic & %.3f & \\\\", t_ext),
    sprintf("Pre-treatment mean of $1\\{Y>0\\}$ (treated) & %.4f & \\\\", pre_mean_ext),
    "\\midrule",
    "\\multicolumn{3}{l}{\\textit{Panel B: dynamic (event time)}} \\\\",
    dyn_rows,
    "\\midrule",
    sprintf("Observations & %s & \\\\", format(n_ext_obs, big.mark = ",")),
    sprintf("Countries & %d & \\\\", n_cty_ext),
    sprintf("Pre-trend $\\chi^2$ & %.3f & \\\\", pt_ext$stat),
    sprintf("Pre-trend $p$ & %.3f & \\\\", pt_ext$pval),
    sprintf("\\texttt{did} pre-test $p$ & %s & \\\\",
            if (is.na(pt_ext$Wpval_did)) "---" else sprintf("%.3f", pt_ext$Wpval_did)),
    "\\bottomrule",
    "\\end{tabular}"
  )

  notes_ext <- paste0(
    "Outcome: $1\\{$adaptation commitments $>0\\}$. CS\\,(2021) DR; WGI GE + log ",
    "population; headline specification; boot SE (", BITERS,
    " reps, seed 1242). Panel B: dynamic aggregation, $e \\in [",
    paste(range(as.integer(agg_d_ext$egt[keep_ext])), collapse = ", "),
    "]$ ($e=-1$ omitted, base period). ",
    wpval_reconciliation(pre_egt = pt_ext$leads, df_did = pt_ext$df_did,
                         n_clusters = n_cty_ext, wpval_did = pt_ext$Wpval_did,
                         pval_wald = pt_ext$pval,
                         wpval_reason = pt_ext$wpval_reason,
                         ginv_used = pt_ext$ginv_used, df = pt_ext$df),
    "Outcome $=1$ for any positive marked amount, however small. ",
    "* $p<0.10$, ** $p<0.05$, *** $p<0.01$"
  )

  write_tex_float(
    out_path      = here("output", "tables", "extensive_margin", "att_extensive.tex"),
    caption_title = paste0("Extensive margin: effect of NAP adoption on the ",
                           "probability of receiving any adaptation-marked commitment"),
    label         = "tab:extensive_margin",
    tabular_lines = ext_tabular,
    notes_text    = notes_ext,
    source_text   = "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"
  )

  saveRDS(list(gt_boot = gt_ext, gt_analytical = gt_ext_analytical,
               agg_simple = agg_s_ext, agg_dynamic = agg_d_ext,
               att = att_ext, se = se_ext, pretrend = pt_ext),
          here("output", "fits", "extensive_margin_dr_bs.rds"))
}

message("\n=== 03_main_results.R: complete ===\n")
