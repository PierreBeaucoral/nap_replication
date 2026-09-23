# ==============================================================================
# 03_main_results.R
# Main CS (2021) estimation — cohorts_dropped specification + pre-trends figure.
# Paper: Beaucoral, Goujon and Marchand (2026) — §4 (pre-trends), §5 (main results, Table 2, Figures 2 and B.4)
#
# Inputs : data/processed/simple_panel_wgi.csv
# Outputs:
#   output/figures/pretrends_analysis.png                     (§4 pre-trends)
#   output/figures/cohorts_dropped/did_combined_es_wgi.png    (Figure 2)
#   output/figures/cohorts_dropped/did_combined_cohort_wgi.png(Figure B.4)
#   output/tables/cohorts_dropped/att_combined_wide.tex       (Table 2 main)
#   output/tables/nap_cohorts.tex                             (appendix cohorts)
#
# NOTE: per-outcome did_es_*, did_cohort_*, did_cohort_es_* figures
#   are built in memory but NOT written to disk (not used in paper).
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
# NOTE: MASS is NOT attached via library() -- MASS::select() would mask
# dplyr::select() used throughout this script. compute_pretrend_test() below
# calls MASS::ginv() by full namespace instead (singular-covariance fallback).

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

# Null-coalescing operator (used in compute_pretrend_test)
`%||%` <- function(a, b) if (!is.null(a)) a else b

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
  mutate(country_id = as.integer(factor(recipient_name)))

# Build gname (cohort = year of first NAP adoption).
# att_gt() requires gname to be CONSTANT within each unit across all years.
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
    " treated units (2021--2024); the retained-cohorts robustness specification ",
    "(Section~\\ref{sec:robust}) additionally includes the smaller 2015--2020 cohorts"),
  source_text   = "UNFCCC NAP Central tracking tool"
)

# ==============================================================================
# SECTION 3. Pre-trends figure — paper Figure pretrends_analysis
# Output: output/figures/pretrends_analysis.png  (final name — no rename shim)
# labs(title = NULL, subtitle = NULL): titles go in the LaTeX \caption{}
# ==============================================================================

message("\n=== Building pre-trends figure ===\n")

yearly_averages <- aggregated %>%
  group_by(year, treated) %>%
  summarise(
    avg_commitments  = mean(commitments,  na.rm = TRUE),
    avg_disbursements = mean(disbursements, na.rm = TRUE),
    .groups = "drop"
  )

p6 <- ggplot(yearly_averages, aes(x = year)) +
  geom_line(aes(y = avg_commitments,  color = "Commitments",  linetype = factor(treated)), linewidth = 1) +
  geom_line(aes(y = avg_disbursements, color = "Disbursements", linetype = factor(treated)), linewidth = 1) +
  geom_point(aes(y = avg_commitments,  color = "Commitments",  shape = factor(treated)), size = 2) +
  geom_point(aes(y = avg_disbursements, color = "Disbursements", shape = factor(treated)), size = 2) +
  scale_color_manual(values = c("Commitments" = "#2E86C1", "Disbursements" = "#E67E22")) +
  scale_linetype_manual(values = c("0" = "dashed", "1" = "solid"),
                        labels = c("Control", "Treated"), name = "Group") +
  scale_shape_manual(values = c("0" = 1, "1" = 16),
                     labels = c("Control", "Treated"), name = "Group") +
  # No title, subtitle, or caption — those go in LaTeX \caption{}
  labs(title = NULL, subtitle = NULL, caption = NULL,
       x = "Year", y = "USD (Millions)", color = "Type") +
  theme_minimal() +
  theme(text = element_text(family = "serif", size = 12),
        legend.position = "bottom")

ggsave(here("output", "figures", "pretrends_analysis.png"),
       p6, width = 12, height = 8, dpi = 300)
message("Saved: output/figures/pretrends_analysis.png")

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
# SECTION 5. Helper: run_did_estimation()
# Encapsulates the full CS(2021) multi-outcome loop for one cohort strategy.
# Seed rule (reproduces the published SEs): set.seed(1242) stays immediately before the att_gt loop.
# ==============================================================================

run_did_estimation <- function(did_panel_in, retain_thin_cohorts,
                               outcomes, thin_cohorts, thin_threshold = 5L) {

  cohort_tag <- if (retain_thin_cohorts) "cohorts_retained" else "cohorts_dropped"
  dir_figs   <- file.path(here("output", "figures"), cohort_tag)
  dir_tabs   <- file.path(here("output", "tables"),  cohort_tag)
  dir.create(dir_figs, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_tabs, recursive = TRUE, showWarnings = FALSE)
  message(sprintf("\n=== run_did_estimation: %s ===", cohort_tag))
  message(sprintf("    Figures -> %s  |  Tables -> %s", dir_figs, dir_tabs))

  did_panel <- did_panel_in
  if (!retain_thin_cohorts && length(thin_cohorts) > 0) {
    message(sprintf("Dropping cohorts with < %d treated units: %s",
                    thin_threshold, paste(thin_cohorts, collapse = ", ")))
    did_panel <- did_panel %>% filter(!(cohort_year %in% thin_cohorts))
  } else if (retain_thin_cohorts && length(thin_cohorts) > 0) {
    message(sprintf("Retaining thin cohorts (%s). Asymptotic SE (bstrap = FALSE).",
                    paste(thin_cohorts, collapse = ", ")))
  } else {
    message("No thin cohorts at threshold = ", thin_threshold)
  }

  use_bstrap <- if (retain_thin_cohorts) FALSE else TRUE
  use_dr     <- if (retain_thin_cohorts) "reg" else "dr"
  message(sprintf("Bootstrap: %s  |  est_method: %s", use_bstrap, use_dr))

  outcome_vars     <- sapply(outcomes, `[[`, "var")
  missing_outcomes <- setdiff(outcome_vars, names(did_panel))
  if (length(missing_outcomes) > 0) {
    stop("Outcome column(s) missing from did_panel: ",
         paste(missing_outcomes, collapse = ", "),
         "\nCheck that share_adapt and log-transformed columns were created and merged.")
  }
  message("Pre-flight OK — outcome columns present: ", paste(outcome_vars, collapse = ", "))

  results_simple  <- list()
  results_group   <- list()
  results_dynamic <- list()

  # Seed rule (reproduces the published SEs): seed immediately before estimator loop
  set.seed(1242)

  for (oc in outcomes) {

    message(sprintf("\n--- Estimating: %s ---", oc$label))

    # WHAT THIS FIT IS FOR. This pass exists ONLY to
    # draw Figures 2 and B.4, which show SIMULTANEOUS confidence bands
    # (cband = TRUE); the stored fits written by make_wide_table() below use
    # cband = FALSE and therefore cannot supply them. Every NUMBER the paper
    # reports comes from those stored fits, never from here. The replication
    # count is pinned to BITERS (999, as everywhere else), and the seed is set immediately before EACH
    # outcome's estimator call rather than once before the five-outcome loop,
    # so each panel of the figure is reproducible on its own. The point
    # estimates are identical to the stored fits BY CONSTRUCTION (the
    # multiplier bootstrap affects only standard errors), and §10b asserts
    # exactly that rather than trusting it; with the seed and replication count
    # aligned, the aggregated standard errors coincide too, so Figures 2 and B.4
    # and Table 2 report the same numbers.
    set.seed(1242)
    gt_obj <- tryCatch(
      att_gt(
        yname         = oc$var,
        tname         = "year",
        idname        = "country_id",
        gname         = "cohort_year",
        xformla       = ~ ge_est + log_population,
        data          = did_panel,
        est_method    = use_dr,
        bstrap        = use_bstrap,
        biters        = BITERS,
        cband         = use_bstrap,
        control_group = "nevertreated",
        anticipation  = 0,
        base_period   = "universal",
        panel         = TRUE,
        allow_unbalanced_panel = TRUE
      ),
      error = function(e) {
        message("  att_gt failed for ", oc$var, ": ", conditionMessage(e)); NULL
      }
    )
    if (is.null(gt_obj)) next

    # Simple ATT
    agg_s <- tryCatch(aggte(gt_obj, type = "simple", na.rm = TRUE), error = function(e) NULL)
    if (!is.null(agg_s)) {
      results_simple[[oc$var]] <- data.frame(
        outcome = oc$label,
        ATT     = round(agg_s$overall.att, 4),
        SE      = round(agg_s$overall.se,  4),
        t_stat  = round(agg_s$overall.att / agg_s$overall.se, 3)
      )
      message(sprintf("  Simple ATT = %.4f (SE = %.4f)", agg_s$overall.att, agg_s$overall.se))
    }

    # By-cohort ATT
    agg_g <- tryCatch(aggte(gt_obj, type = "group", na.rm = TRUE), error = function(e) NULL)
    if (!is.null(agg_g)) {
      results_group[[oc$var]] <- data.frame(
        outcome = oc$label,
        cohort  = agg_g$egt,
        ATT     = round(agg_g$att.egt,  4),
        SE      = round(agg_g$se.egt,   4),
        Lower   = round(agg_g$att.egt - agg_g$crit.val.egt * agg_g$se.egt, 4),
        Upper   = round(agg_g$att.egt + agg_g$crit.val.egt * agg_g$se.egt, 4)
      )
      # (Figure B.4 is drawn from results_group below.)
    }

    # Dynamic ATT (event study)
    agg_d <- tryCatch(
      aggte(gt_obj, type = "dynamic", na.rm = TRUE,
            min_e = -5,
            max_e = Inf),
      error = function(e) NULL
    )
    if (!is.null(agg_d)) {
      results_dynamic[[oc$var]] <- data.frame(
        outcome    = oc$label,
        event_time = agg_d$egt,
        ATT        = round(agg_d$att.egt, 4),
        SE         = round(agg_d$se.egt,  4),
        Lower      = round(agg_d$att.egt - agg_d$crit.val.egt * agg_d$se.egt, 4),
        Upper      = round(agg_d$att.egt + agg_d$crit.val.egt * agg_d$se.egt, 4),
        color      = oc$color
      )
      # Per-outcome event-study plot built in memory only — not used in paper
      p_dyn <- ggdid(agg_d) +
        geom_vline(xintercept = -0.5, linetype = "dashed", colour = "grey40") +
        labs(title    = paste("Event Study —", oc$label),
             subtitle = "Controls: WGI gov. effectiveness, log population",
             x        = "Years relative to NAP adoption", y = "ATT",
             caption  = "CS (2021); outcome regression; never-treated controls") +
        theme_minimal() +
        theme(text       = element_text(family = "serif", size = 11),
              plot.title = element_text(face = "bold", size = 13, hjust = 0.5))
    }

    # Cohort-specific event study (faceted)
    es_raw <- data.frame(
      cohort   = gt_obj$group,
      cal_time = gt_obj$t,
      att      = gt_obj$att,
      se       = gt_obj$se
    ) %>%
      filter(!is.na(se), se > 0) %>%
      mutate(
        event_time = cal_time - cohort,
        cohort_lab = paste0("Cohort ", cohort),
        ci_lo      = att - 1.96 * se,
        ci_hi      = att + 1.96 * se
      )
    # Per-outcome cohort event-study plot built in memory only — not used in paper
    p_cs_ev <- ggplot(es_raw, aes(x = event_time, y = att)) +
      geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
      geom_vline(xintercept = -0.5, colour = "grey30", linetype = "dotted") +
      geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, fill = oc$color) +
      geom_line(colour = oc$color, linewidth = 0.8) +
      geom_point(colour = oc$color, size = 2) +
      facet_wrap(~ cohort_lab, scales = "free_y") +
      labs(
        title    = paste("Cohort Event Studies —", oc$label),
        subtitle = "Controls: WGI gov. effectiveness, log population",
        x        = "Years relative to NAP adoption", y = "ATT",
        caption  = "CS (2021); outcome regression; never-treated controls; 95% CI (±1.96 SE)"
      ) +
      theme_minimal() +
      theme(text             = element_text(family = "serif", size = 11),
            plot.title       = element_text(face = "bold", size = 13, hjust = 0.5),
            strip.text       = element_text(face = "bold"),
            panel.grid.minor = element_blank())
  }

  # --- Combined event-study overlay (paper Figure 2, cohorts_dropped) ---
  if (length(results_dynamic) == 0) {
    message("No dynamic results — skipping combined plot.")
    return(invisible(list(simple  = results_simple,
                          group   = results_group,
                          dynamic = results_dynamic)))
  }

  dynamic_all <- bind_rows(results_dynamic) %>%
    mutate(outcome = factor(outcome, levels = sapply(outcomes, `[[`, "label")))

  palette_vec <- setNames(sapply(outcomes, `[[`, "color"),
                          sapply(outcomes, `[[`, "label"))
  dodge_w <- 0.4

  p_combined <- ggplot(dynamic_all,
                       aes(x = event_time, y = ATT,
                           colour = outcome, shape = outcome,
                           group  = outcome)) +
    geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
    geom_vline(xintercept = -0.5, colour = "grey30", linetype = "dotted") +
    geom_linerange(aes(ymin = Lower, ymax = Upper),
                   position = position_dodge(width = dodge_w),
                   linewidth = 0.6, alpha = 0.8) +
    geom_point(size = 2.5, position = position_dodge(width = dodge_w)) +
    scale_colour_manual(values = palette_vec) +
    scale_shape_manual(values = c(16, 17, 15, 18, 8)) +
    # No title, subtitle, or caption — those go in LaTeX \caption{}
    labs(
      title    = NULL, subtitle = NULL, caption = NULL,
      x        = "Years relative to NAP adoption",
      y        = "ATT estimate",
      colour   = NULL, shape = NULL
    ) +
    theme_minimal() +
    theme(
      text             = element_text(family = "serif", size = 12),
      legend.position  = "bottom",
      legend.text      = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  print(p_combined)
  ggsave(file.path(dir_figs, "did_combined_es_wgi.png"),
         p_combined, width = 12, height = 7, dpi = 300)

  # att_simple_all.tex and att_dynamic_all.tex are NOT written to disk —
  # they are not used by the manuscript. Objects remain in memory below.
  simple_all <- bind_rows(results_simple) %>%
    mutate(
      stars   = case_when(
        abs(t_stat) > 2.576 ~ "***",
        abs(t_stat) > 1.960 ~ "**",
        abs(t_stat) > 1.645 ~ "*",
        TRUE                ~ ""
      ),
      ATT_fmt = paste0(sprintf("%.4f", ATT), stars)
    )

  message(sprintf("\n--- Completed: %s | Figures -> %s | Tables -> %s ---",
                  cohort_tag, dir_figs, dir_tabs))

  invisible(list(simple  = results_simple,
                 group   = results_group,
                 dynamic = results_dynamic))
}

# ==============================================================================
# SECTION 6. Main estimation — cohorts_dropped (bootstrap SE)
# ==============================================================================

message("\n=== Main estimation: cohorts_dropped (bootstrap SE) ===\n")

res_main <- run_did_estimation(
  did_panel_in        = did_panel,
  retain_thin_cohorts = FALSE,
  outcomes            = outcomes,
  thin_cohorts        = thin_cohorts,
  thin_threshold      = thin_threshold
)

# ==============================================================================
# SECTION 7. Helper: compute_pretrend_test()
# Wald chi-sq on all pre-treatment event-time ATTs using influence-function
# covariance (analytical, not bootstrap).
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

# ==============================================================================
# SECTION 7b. Helpers — save the headline fit so that every downstream
# script reports ONE ATT and ONE SE per specification.
#
# WHY THIS EXISTS (diagnosed 2026-09-15; the discrepancy is not bootstrap
# Monte Carlo noise): Table 2 reported SE 0.1245 for log(adaptation
# commitments) while the §19e units/zeros diagnostic reported 0.1173 for the
# byte-identical specification. The two fits are identical -- same panel, same
# est_method, same seed immediately before att_gt(), same 999 replications --
# and they return the same ATT to 6 decimals. The SEs differ because
# `aggte()` RUNS ITS OWN multiplier bootstrap when the att_gt object was
# fitted with bstrap = TRUE, so the aggregated SE depends on the RNG state at
# the aggte() call, not on the RNG state at the att_gt() call:
#   03 make_wide_table : set.seed(1242) -> att_gt(bstrap) -> aggte(simple)
#                        [RNG state at aggte = post-att_gt-mboot]    -> 0.1245
#   04 §19e            : set.seed(1242) -> att_gt(bstrap)
#                        -> set.seed(1242) -> att_gt(analytical, consumes no RNG)
#                        -> aggte(simple) [RNG state = fresh 1242]   -> 0.1173
# "Re-seed immediately before each att_gt()" therefore does NOT pin the
# reported SE. The fix is structural, not a seeding tweak: fit and aggregate
# once here, saveRDS the result, and have every other script read the saved
# aggregation instead of re-running it.
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
# SECTION 8. Helper: make_wide_table()
# Produces att_combined_wide.tex for cohorts_dropped (§19c main half).
# Seed rule (reproduces the published SEs): set.seed(1242) immediately before the outcomes loop
# ==============================================================================

# Map from log outcome var name to its underlying raw column (USD millions)
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

make_wide_table <- function(did_panel_in, retain_thin, outcomes,
                             dir_tabs, tex_label, caption_spec) {

  use_dr <- if (retain_thin) "reg" else "dr"

  table_stats <- list()
  # Seed rule (reproduces the published SEs): seed immediately before outcomes loop
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
    message(sprintf("    Mean pre-treat (raw) = %s  |  Implied effect (USD M) = %s",
                    mean_pre_fmt, implied_fmt))

    # ------------------------------------------------------------------
    # Persist the fit behind this column of Table 2.
    # Everything below runs AFTER every reported quantity for this outcome
    # has been computed, so it cannot perturb the RNG stream that produced
    # them; the next loop iteration re-seeds with set.seed(1242) before its
    # own bootstrap fit, so it cannot perturb the next column either.
    # Saved only for the main never-treated / cohorts >= 5 specification.
    # ------------------------------------------------------------------
    if (!retain_thin && !is.na(fit_stem_map[oc$var])) {

      agg_dyn_boot <- tryCatch(
        aggte(gt_tab, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf),
        error = function(e) NULL
      )
      agg_grp_boot <- tryCatch(aggte(gt_tab, type = "group", na.rm = TRUE),
                               error = function(e) NULL)
      # Analytical simple aggregation: supplies the unit-level influence
      # function used by 05_heterogeneity.R for correlated-sample contrasts
      # `ids` records the unit order those rows correspond to.
      agg_s_analytical <- if (!is.null(gt_tab_analytical))
        tryCatch(aggte(gt_tab_analytical, type = "simple", na.rm = TRUE),
                 error = function(e) NULL) else NULL

      # l_vec matching the ATT printed in Table 2 (group-size-weighted average
      # of the post-treatment event-time effects), with an explicit
      # reconciliation against agg_s$overall.att.
      l_vec_simple <- NULL
      l_vec_check  <- NA_real_
      egt_post     <- NULL
      if (!is.null(agg_dyn_boot)) {
        keep_d   <- which(!is.na(agg_dyn_boot$se.egt) & agg_dyn_boot$se.egt > 1e-10)
        egt_post <- agg_dyn_boot$egt[keep_d][agg_dyn_boot$egt[keep_d] >= 0]
        if (length(egt_post) > 0) {
          l_vec_simple <- simple_att_event_weights(gt_tab, did_panel_in, egt_post)
          att_post     <- agg_dyn_boot$att.egt[match(egt_post, agg_dyn_boot$egt)]
          l_vec_check  <- sum(l_vec_simple * att_post) - att
          message(sprintf(
            "    [fit] l_vec = (%s) reproduces the simple ATT to %.2e",
            paste(sprintf("%.4f", l_vec_simple), collapse = ", "), abs(l_vec_check)))
        }
      }

      # Bootstrap covariance of the dynamic aggregation (feeds HonestDiD in 04).
      bvc <- if (!is.null(agg_dyn_boot))
        tryCatch(boot_vcov_dynamic(agg_dyn_boot, gt_tab, seed = 1242L),
                 error = function(e) { message("    [fit] boot vcov failed: ",
                                                conditionMessage(e)); NULL })
      else NULL
      if (!is.null(bvc))
        message(sprintf(paste0("    [fit] bootstrap vcov: max |se_cov - se_published| / ",
                               "se_published = %.4f"), bvc$max_rel_dev))

      fit_obj <- list(
        outcome        = oc$var,
        outcome_label  = oc$label,
        spec           = list(
          estimator      = "CS (2021)",
          est_method     = use_dr,
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
        gt_boot          = gt_tab,
        gt_analytical    = gt_tab_analytical,
        agg_simple       = agg_s,
        agg_simple_analytic = agg_s_analytical,
        agg_dynamic      = agg_dyn_boot,
        agg_group        = agg_grp_boot,
        agg_dyn_analytic = agg_dyn_analytical,
        ids              = sort(unique(did_panel_in$country_id)),
        att              = att,
        se               = se,
        pretrend         = pt,
        l_vec_simple     = l_vec_simple,
        l_vec_egt        = egt_post,
        l_vec_recon_gap  = l_vec_check,
        boot_vcov        = bvc,
        n_obs            = n_obs,
        n_country        = n_country,
        provenance       = list(
          script       = "code/03_main_results.R",
          created      = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          r_version    = R.version.string,
          did_version  = as.character(utils::packageVersion("did"))
        )
      )
      fit_path <- here("output", "fits",
                       paste0("headline_", fit_stem_map[[oc$var]], "_dr_bs.rds"))
      saveRDS(fit_obj, fit_path)
      message("    [fit] Saved: ", fit_path)
    }
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

  tab_wide <- data.frame(` ` = row_labels, check.names = FALSE,
                          stringsAsFactors = FALSE)

  # Row label reflects the actual inference method: bootstrap for main spec
  # (retain_thin = FALSE), analytical for thin-cohort robustness spec.
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
  se_label <- if (!retain_thin) "multiplier-bootstrap SE (999 reps, seed 1242; \\texttt{did} 2.5.0; CRS Apr.\\ 2026)" else "analytical (IF) SE; \\texttt{did} 2.5.0; CRS Apr.\\ 2026"

  notes_txt <- paste0(
    "CS\\,(2021) ", est_method_label,
    "; WGI gov.\\ effectiveness + log population; ",
    "never-treated control group; ", se_label, ". ",
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

  out_path <- file.path(dir_tabs, "att_combined_wide.tex")
  write_tex_float(out_path, cap_title, tex_label, raw_lines, notes_txt, source_txt)

  invisible(table_stats)
}

# ==============================================================================
# SECTION 9. §19b combined cohort plot — cohorts_DROPPED half
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

# cohorts_dropped cohort plot (main spec)
make_cohort_plot(
  results_group = res_main$group,
  outcomes      = outcomes,
  dir_figs      = file.path(here("output", "figures"), "cohorts_dropped"),
  spec_label    = "Cohorts ≥ 5 units (bootstrap SE)",
  use_bstrap    = TRUE
)

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
  caption_spec  = "cohorts $\\geq 5$ units, DR + multiplier-bootstrap SE"
)

# ==============================================================================
# SECTION 10b. Reconciliation: the figure fit against the stored fits
# The event-study and cohort figures (Figures 2 and B.4) are drawn from a second
# pass over the same specification, needed only because they display
# SIMULTANEOUS confidence bands (cband = TRUE) while the stored fits use
# cband = FALSE. Point estimates cannot differ -- the multiplier bootstrap
# affects only standard errors -- so this is asserted rather than assumed, and
# the standard-error gap between the two draws is printed so that the size of
# the simulation noise is on the record.
# ==============================================================================

message("\n=== Section 10b: figure fit vs stored fit reconciliation ===\n")

for (oc in outcomes) {
  stem <- fit_stem_map[[oc$var]]
  if (is.na(stem)) next
  fit_path <- here("output", "fits", paste0("headline_", stem, "_dr_bs.rds"))
  if (!file.exists(fit_path)) next
  stored <- readRDS(fit_path)
  fig    <- res_main$simple[[oc$var]]
  if (is.null(fig) || is.null(stored$att)) next
  # The figure pass rounds its stored ATT to 4 decimals, so compare at that
  # tolerance; the SE gap is reported, not asserted.
  stopifnot(abs(fig$ATT - round(stored$att, 4)) < 1e-9)
  message(sprintf(paste0("  %-38s ATT identical (%.4f) | SE figure draw %.4f ",
                         "vs stored %.4f (gap %+.4f)"),
                  oc$label, fig$ATT, fig$SE, stored$se, fig$SE - stored$se))
}

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
                             as.integer(commitments > 0)))

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

if (is.null(gt_ext)) {
  message("  Extensive-margin estimation failed — table not written.")
} else {
  agg_s_ext <- aggte(gt_ext, type = "simple",  na.rm = TRUE)
  agg_d_ext <- aggte(gt_ext, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf)
  agg_d_ext_analytical <- if (!is.null(gt_ext_analytical))
    tryCatch(aggte(gt_ext_analytical, type = "dynamic", na.rm = TRUE,
                   min_e = -5, max_e = Inf), error = function(e) NULL) else NULL

  att_ext <- agg_s_ext$overall.att
  se_ext  <- agg_s_ext$overall.se
  t_ext   <- att_ext / se_ext
  stars_ext <- if (abs(t_ext) > 2.576) "***" else if (abs(t_ext) > 1.960) "**" else
               if (abs(t_ext) > 1.645) "*" else ""
  pt_ext <- if (!is.null(agg_d_ext_analytical))
    compute_pretrend_test(agg_d_ext_analytical, gt_ext_analytical) else
    list(stat = NA_real_, pval = NA_real_, df = 0L,
         W_did = NA_real_, Wpval_did = NA_real_, df_did = NA_integer_)

  # MDE at 80% power, 5% two-sided (the constant is z_{0.975} + z_{0.80}).
  mde_ext <- (qnorm(0.975) + qnorm(0.80)) * se_ext
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
    " reps, seed 1242). Panel B: dynamic aggregation, $e \\in [-5, 3]$ ($e=-1$ ",
    "omitted, base period). ",
    wpval_reconciliation(pre_egt = pt_ext$leads, df_did = pt_ext$df_did,
                         n_clusters = n_cty_ext, wpval_did = pt_ext$Wpval_did,
                         pval_wald = pt_ext$pval,
                         wpval_reason = pt_ext$wpval_reason,
                         ginv_used = pt_ext$ginv_used),
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
