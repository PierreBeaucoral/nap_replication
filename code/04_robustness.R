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
#   tables/cohorts_retained/att_group_retained.tex         (§19 cohort ATTs, retained spec)
#   tables/cohorts_dropped/honestdid_rm.tex                (§20 HonestDiD)
#   tables/cohorts_dropped/honestdid_prepriods.tex         (§20 HonestDiD, 4/6/8 leads)
#   tables/cohorts_dropped/honestdid_sd.tex                (§20 HonestDiD, smoothness)
#   figures/placebo/did_placebo_es.png                     (§25)
#   tables/placebo/att_placebo.tex                         (§25)
#   figures/mitigation/did_mitigation_es.png               (§26)
#   tables/mitigation/att_mitigation.tex                   (§26)
#   tables/dcdh/att_dcdh.tex                                (§A dCDH estimator)
#   figures/dcdh/did_dcdh_es.png                            (§A dCDH event study)
#   tables/bacon/bacon_decomp.tex                           (§B Goodman-Bacon)
#   figures/bacon/bacon_scatter.png                         (§B Goodman-Bacon)
#   tables/outlier_india/att_combined_wide.tex             (§C excluding India)
#   tables/balanced_panel/att_combined_wide.tex            (§D balanced panel)
#   tables/panel_2010/att_combined_wide.tex                (§E panel from 2010)
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
source(here("code", "functions", "pretrend_test.R"))      # compute_pretrend_test()
source(here("code", "functions", "two_line_head.R"))     # two_line_head()
source(here("code", "functions", "read_headline_fit.R"))  # read_headline_fit()
source(here("code", "functions", "make_country_id.R"))    # make_country_id()
source(here("code", "functions", "crs_positive.R"))       # crs_positive()
source(here("code", "functions", "make_wide_table.R"))    # make_wide_table()
source(here("code", "functions", "sup_t_crit.R"))         # sup_t_crit()
# §A/§B appendix estimators (de Chaisemartin & D'Haultfoeuille; Goodman-Bacon).
# polars MUST be loaded before DIDmultiplegtDYN, else did_multiplegt_dyn() errors
# with "objet 'pl' introuvable" (the package references the loaded `pl` object).
library(polars)
library(DIDmultiplegtDYN)
library(bacondecomp)

# HonestDiD (required; used in §20; renv.lock pins 0.2.8). Guarded with
# requireNamespace() so a missing install fails with an actionable message
# rather than a bare "there is no package called 'HonestDiD'" error.
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
# Every multiplier-bootstrap call in this script passes BITERS explicitly
# (did's own default is 1000), so the replication count stated in the table
# notes is the one used.
# -----------------------------------------------------------------------
BITERS <- 999L

# -----------------------------------------------------------------------
# §19e and §20 read the headline fits written by 03_main_results.R
# (read_headline_fit()) instead of re-fitting, so the paper can never print two
# different standard errors for one specification. This stage writes its own
# placebo and mitigation fits to the same folder.
# -----------------------------------------------------------------------
FITS_DIR <- here("output", "fits")

# A failed run must not leave the previous run's results in place: the fits
# and exhibits this stage owns are deleted before anything is re-estimated, and
# an estimator failure below stops the stage.
unlink(file.path(FITS_DIR, c("mitigation_dr_bs.rds", "placebo_m2_reg_bs.rds")))
unlink(c(list.files(here("output", c("tables", "figures"),
                         rep(c("cohorts_retained", "notyettreated", "units_zeros",
                               "outlier_india", "balanced_panel", "panel_2010", "placebo",
                               "mitigation", "dcdh", "bacon"), each = 2L)),
                    pattern = "\\.(tex|png|pdf)$", full.names = TRUE),
         here("output", "tables", "cohorts_dropped",
              c("honestdid_rm.tex", "honestdid_sd.tex", "honestdid_prepriods.tex"))))

# canonical pre-trend wording — keep byte-identical across scripts
PRETREND_NOTE_AGG <- function(min_e, max_e, k) sprintf(
  "Pre-trend $\\chi^2$ / $p$: joint Wald test on the aggregated pre-treatment event-time coefficients ($%d \\leq e \\leq %d$; %d restrictions), using the influence-function covariance of the dynamic aggregation from the analytical (non-bootstrap) fit; a generalized inverse is used if the block is singular",
  min_e, max_e, k)
PRETREND_NOTE_DID <- "\\texttt{did} pre-test $p$: \\texttt{did}'s built-in Wald test over all pre-period $ATT(g,t)$ cells against each cohort's $g-1$ base year"
PRETREND_NOTE <- function(min_e, max_e, k) paste0(PRETREND_NOTE_AGG(min_e, max_e, k), ". ", PRETREND_NOTE_DID)

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
# make_wide_table() is sourced from code/functions/make_wide_table.R and
# compute_pretrend_test() from code/functions/pretrend_test.R.
# ==============================================================================

# ==============================================================================
# SECTION 4. Helper: make_cohort_plot()
# (duplicated from 03_main_results.R §9)
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
# SECTION 5. §19 retained-cohorts robustness spec
# Seed rule (reproduces the published SEs): set.seed(1242) before outcomes loop.
# ==============================================================================

message("\n=== Section 19: retained-cohorts robustness ===\n")

dir.create(here("output", "figures", "cohorts_retained"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "tables",  "cohorts_retained"), recursive = TRUE, showWarnings = FALSE)

# Storage for retained-cohorts results, one slot per outcome
outcome_vars      <- vapply(outcomes, `[[`, character(1L), "var")
results_simple_r  <- setNames(vector("list", length(outcomes)), outcome_vars)
results_group_r   <- setNames(vector("list", length(outcomes)), outcome_vars)
results_dynamic_r <- setNames(vector("list", length(outcomes)), outcome_vars)
agg_group_r       <- setNames(vector("list", length(outcomes)), outcome_vars)

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
  if (is.null(gt_r)) stop("Retained-cohorts att_gt() failed for ", oc$var)

  agg_s_r <- tryCatch(aggte(gt_r, type = "simple", na.rm = TRUE), error = function(e) NULL)
  agg_g_r <- tryCatch(aggte(gt_r, type = "group",  na.rm = TRUE), error = function(e) NULL)
  agg_d_r <- tryCatch(
    aggte(gt_r, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf),
    error = function(e) NULL
  )
  if (is.null(agg_s_r) || is.null(agg_g_r) || is.null(agg_d_r))
    stop("Retained-cohorts aggte() failed for ", oc$var)

  if (!is.null(agg_s_r)) {
    results_simple_r[[oc$var]] <- data.frame(
      outcome = oc$label,
      ATT     = round(agg_s_r$overall.att, 4),
      SE      = round(agg_s_r$overall.se,  4),
      t_stat  = round(agg_s_r$overall.att / agg_s_r$overall.se, 3)
    )
  }
  if (!is.null(agg_g_r)) {
    agg_group_r[[oc$var]] <- agg_g_r
    # Figure C.10 band: simultaneous (sup-t) 95% over this outcome's cohorts.
    cv_g_r <- sup_t_crit(agg_g_r$inf.function$selective.inf.func.g, agg_g_r$se.egt,
                         biters = BITERS)
    message(sprintf("  Retained cohort fig sup-t crit (%s): %.4f", oc$label, cv_g_r))
    # Cohort-level ATTs quoted in the text (Figure C.10 plots them).
    message(sprintf("  Cohort ATTs [%s, retained]: %s", oc$label,
                    paste(sprintf("%d %.4f (%.4f)", as.integer(agg_g_r$egt),
                                  agg_g_r$att.egt, agg_g_r$se.egt), collapse = "; ")))
    results_group_r[[oc$var]] <- data.frame(
      outcome = oc$label,
      cohort  = agg_g_r$egt,
      ATT     = round(agg_g_r$att.egt,  4),
      SE      = round(agg_g_r$se.egt,   4),
      Lower   = round(agg_g_r$att.egt - cv_g_r * agg_g_r$se.egt, 4),
      Upper   = round(agg_g_r$att.egt + cv_g_r * agg_g_r$se.egt, 4)
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
  spec_label    = "All cohorts retained (asymptotic SE)"
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

# --- §19f cohort ATTs and pooled early / late adopters, retained spec --------
# The group aggregations above come from the retained-cohorts fits behind
# Table tab:combined_wide_robust (outcome regression, analytical SEs, all
# cohorts, never-treated controls). A pooled row averages the ATTs of a range
# of cohorts with weights proportional to each cohort's number of recipients,
# held fixed. Its influence function is the same weighted sum of the cohort
# influence functions (inf.function$selective.inf.func.g), and its SE uses the
# formula did applies to this analytical fit, sqrt(mean(IF^2) / n), which is
# checked below against did's own cohort SEs.
grp_outcomes <- c(log_commits           = "log(Adaptation commitments)",
                  lcommitments_all      = "log(Total commitments)",
                  lcommitments_nonadapt = "log(Non-adaptation commitments)")
grp_pools <- list(list(label = "2015--2019 adopters", years = 2015:2019),
                  list(label = "2021--2024 adopters", years = 2021:2024))
grp_cohorts <- cohort_sizes$cohort_year
stopifnot(all(unlist(lapply(grp_pools, `[[`, "years")) %in% grp_cohorts))

#' Fixed-weight pooled ATT over a set of cohorts from a group aggregation
#'
#' @param agg_g AGGTEobj from aggte(type = "group") of an analytical fit
#' @param years cohorts to pool
#' @return list(att, se, p, n_rec)
pool_cohort_atts <- function(agg_g, years) {
  IF_g <- agg_g$inf.function$selective.inf.func.g
  stopifnot(is.matrix(IF_g), ncol(IF_g) == length(agg_g$egt))
  se_chk <- sqrt(colMeans(IF_g^2) / nrow(IF_g))
  stopifnot(max(abs(se_chk - agg_g$se.egt) / agg_g$se.egt) < 1e-10)
  j <- match(years, agg_g$egt)
  stopifnot(!anyNA(j), !anyNA(agg_g$att.egt[j]))
  n_rec <- cohort_sizes$n_treated[match(years, cohort_sizes$cohort_year)]
  w     <- n_rec / sum(n_rec)
  IF_p  <- as.numeric(IF_g[, j, drop = FALSE] %*% w)
  att   <- sum(w * agg_g$att.egt[j])
  se    <- sqrt(mean(IF_p^2) / length(IF_p))
  list(att = att, se = se, p = 2 * pnorm(-abs(att / se)), n_rec = sum(n_rec))
}

stars_of_p <- function(p) if (is.na(p)) "" else if (p < 0.01) "***" else
  if (p < 0.05) "**" else if (p < 0.10) "*" else ""
grp_cell <- function(att, se) {
  p <- 2 * pnorm(-abs(att / se))
  c(paste0(sprintf("%.4f", att), stars_of_p(p)), sprintf("(%.4f)", se), sprintf("%.3f", p))
}

stopifnot(all(!vapply(agg_group_r[names(grp_outcomes)], is.null, logical(1L))))
grp_rows <- vapply(seq_along(grp_cohorts), function(i) {
  g <- grp_cohorts[i]
  cells <- unlist(lapply(names(grp_outcomes), function(v) {
    a <- agg_group_r[[v]]
    j <- match(g, a$egt)
    stopifnot(!is.na(j))
    message(sprintf("  Retained cohort ATT [%s] %d (%d recipients): %.4f (SE %.4f, p %.3f)",
                    v, as.integer(g), cohort_sizes$n_treated[i], a$att.egt[j], a$se.egt[j],
                    2 * pnorm(-abs(a$att.egt[j] / a$se.egt[j]))))
    grp_cell(a$att.egt[j], a$se.egt[j])
  }))
  paste0(g, " & ", cohort_sizes$n_treated[i], " & ", paste(cells, collapse = " & "), " \\\\")
}, character(1L))
pool_rows <- vapply(grp_pools, function(pl) {
  res <- lapply(names(grp_outcomes), function(v) pool_cohort_atts(agg_group_r[[v]], pl$years))
  for (k in seq_along(res))
    message(sprintf("  Retained pooled ATT [%s] %s (%d recipients): %.4f (SE %.4f, p %.4f)",
                    names(grp_outcomes)[k], pl$label, res[[k]]$n_rec,
                    res[[k]]$att, res[[k]]$se, res[[k]]$p))
  cells <- unlist(lapply(res, function(r) grp_cell(r$att, r$se)))
  paste0(pl$label, " & ", res[[1L]]$n_rec, " & ", paste(cells, collapse = " & "), " \\\\")
}, character(1L))

n_col_grp <- 2L + 3L * length(grp_outcomes)
grp_tabular <- c(
  paste0("\\begin{tabular}{l", strrep("c", n_col_grp - 1L), "}"),
  "\\toprule",
  paste0(" & & ", paste(sprintf("\\multicolumn{3}{c}{%s}", grp_outcomes), collapse = " & "),
         " \\\\"),
  paste(sprintf("\\cmidrule(lr){%d-%d}", 3L + 3L * (seq_along(grp_outcomes) - 1L),
                5L + 3L * (seq_along(grp_outcomes) - 1L)), collapse = " "),
  paste0("Adoption cohort & Recipients & ",
         paste(rep("ATT & SE & $p$", length(grp_outcomes)), collapse = " & "), " \\\\"),
  "\\midrule",
  grp_rows,
  "\\midrule",
  paste0("\\multicolumn{", n_col_grp, "}{l}{\\textit{Pooled cohorts (weights proportional to recipients)}} \\\\"),
  pool_rows,
  "\\bottomrule",
  "\\end{tabular}"
)
write_tex_float(
  out_path      = file.path(here("output", "tables"), "cohorts_retained", "att_group_retained.tex"),
  caption_title = "ATT by adoption cohort, all cohorts retained",
  label         = "tab:att_group_retained",
  tabular_lines = grp_tabular,
  notes_text    = paste0(
    "Specification of Table~\\ref{tab:combined_wide_robust}: CS\\,(2021) outcome regression, ",
    "all adoption cohorts retained (", format(nrow(did_panel_full), big.mark = ","),
    " recipient-years, ",
    n_distinct(did_panel_full$recipient_name), " recipients), never-treated control group, ",
    "WGI gov.\\ effectiveness + log population, universal base period. Each cohort row is ",
    "that cohort's average post-adoption ATT; ``Recipients'' counts its treated recipients. ",
    "Pooled rows average the cohort ATTs with weights proportional to their recipients, held ",
    "fixed (no term for estimated weights). SEs: analytical influence-function SEs clustered by recipient; a pooled row's SE ",
    "combines the cohorts' influence functions with the same weights. $p$: two-sided, normal. ",
    "Analytical SEs on cohorts of ",
    paste(range(cohort_sizes$n_treated[cohort_sizes$cohort_year %in% thin_cohorts]),
          collapse = "--"),
    " treated recipients are unreliable, so individual thin-cohort $p$-values and stars ",
    "should be read with caution; the pooled rows are the estimates of interest. ",
    "* $p<0.10$, ** $p<0.05$, *** $p<0.01$"),
  source_text   = "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"
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
# the headline panel of 03_main_results.R (Section 10).  Only the control group
# differs from the headline result.
did_panel_nyt <- if (length(thin_cohorts) > 0)
  did_panel_full %>% filter(!(cohort_year %in% thin_cohorts)) else did_panel_full
message(sprintf("Not-yet-treated panel: dropped thin cohorts {%s}; N = %d rows, %d countries",
                paste(thin_cohorts, collapse = ", "),
                nrow(did_panel_nyt), length(unique(did_panel_nyt$country_id))))

results_group_nyt   <- setNames(vector("list", length(outcomes)), outcome_vars)
results_dynamic_nyt <- setNames(vector("list", length(outcomes)), outcome_vars)

# Seed rule (reproduces the published SEs): set.seed(1242) before loop
set.seed(1242)

for (oc in outcomes) {
  message(sprintf("\n--- Not-yet-treated spec: %s ---", oc$label))

  set.seed(1242)  # re-seed before each bootstrap fit (reproduces the published SEs)
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
  if (is.null(gt_nyt)) stop("Not-yet-treated att_gt() failed for ", oc$var)

  agg_g_nyt <- tryCatch(aggte(gt_nyt, type = "group",  na.rm = TRUE), error = function(e) NULL)
  agg_d_nyt <- tryCatch(
    aggte(gt_nyt, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf),
    error = function(e) NULL
  )
  if (is.null(agg_g_nyt) || is.null(agg_d_nyt))
    stop("Not-yet-treated aggte() failed for ", oc$var)
  # Figure bands: simultaneous (sup-t) 95%, one band family per outcome curve
  # (uniform over its cohorts / its event times), on the SEs of these draws.
  cv_g_nyt <- sup_t_crit(agg_g_nyt$inf.function$selective.inf.func.g, agg_g_nyt$se.egt,
                         biters = BITERS)
  cv_d_nyt <- sup_t_crit(agg_d_nyt$inf.function$dynamic.inf.func.e, agg_d_nyt$se.egt,
                         biters = BITERS)
  message(sprintf("  Not-yet-treated sup-t crit (%s): cohort fig %.4f | event-study %.4f",
                  oc$label, cv_g_nyt, cv_d_nyt))

  if (!is.null(agg_g_nyt)) {
    results_group_nyt[[oc$var]] <- data.frame(
      outcome = oc$label,
      cohort  = agg_g_nyt$egt,
      ATT     = round(agg_g_nyt$att.egt, 4),
      SE      = round(agg_g_nyt$se.egt,  4),
      Lower   = round(agg_g_nyt$att.egt - cv_g_nyt * agg_g_nyt$se.egt, 4),
      Upper   = round(agg_g_nyt$att.egt + cv_g_nyt * agg_g_nyt$se.egt, 4)
    )
  }
  if (!is.null(agg_d_nyt)) {
    results_dynamic_nyt[[oc$var]] <- data.frame(
      outcome    = oc$label,
      event_time = agg_d_nyt$egt,
      ATT        = round(agg_d_nyt$att.egt, 4),
      SE         = round(agg_d_nyt$se.egt,  4),
      Lower      = round(agg_d_nyt$att.egt - cv_d_nyt * agg_d_nyt$se.egt, 4),
      Upper      = round(agg_d_nyt$att.egt + cv_d_nyt * agg_d_nyt$se.egt, 4),
      color      = oc$color
    )
  }
}

# --- §19d-i cohort plot (all outcomes): notyettreated ---
make_cohort_plot(
  results_group = results_group_nyt,
  outcomes      = outcomes,
  dir_figs      = dir_figs_nyt,
  spec_label    = "Not-yet-treated controls (multiplier bootstrap)"
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
# Text numbers (units/zeros section): zeros and how many fall in 2009.
is_zero_uz <- !is.na(did_panel_main$commitments) & !crs_positive(did_panel_main$commitments)
message(sprintf(paste0("Exact zeros in adaptation commitments (main sample): %d of %d ",
                       "observations (%.1f%%), %d of them in 2009"),
                sum(is_zero_uz), sum(!is.na(did_panel_main$commitments)),
                100 * mean(is_zero_uz[!is.na(did_panel_main$commitments)]),
                sum(is_zero_uz & did_panel_main$year == 2009L)))

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
         d <- d[!is.na(d$commitments) & crs_positive(d$commitments), ]
         d$uz_y <- log(d$commitments)
         list(data = d, yname = "uz_y")
       })
)

uz_stats <- setNames(vector("list", length(uz_specs) + 1L),
                     c(vapply(uz_specs, `[[`, character(1L), "key"), "ppml_levels"))

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

  # Separate analytical fit (bstrap = FALSE) feeds the pre-trend Wald test only:
  # its dynamic aggregation reports the influence-function SEs the test's
  # covariance reproduces (a bootstrap aggregation reports bootstrap SEs).
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
    error = function(e) { message("  att_gt (analytical) failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(gt_sp) || is.null(gt_sp_analytical))
    stop("Units/zeros estimation failed: ", sp$label)

  agg_s_sp <- aggte(gt_sp, type = "simple", na.rm = TRUE)
  agg_d_sp <- aggte(gt_sp_analytical, type = "dynamic", na.rm = TRUE,
                    min_e = -5, max_e = Inf)

  att_sp <- if (!is.null(agg_s_sp)) agg_s_sp$overall.att else NA_real_
  se_sp  <- if (!is.null(agg_s_sp)) agg_s_sp$overall.se  else NA_real_
  t_sp   <- if (!is.na(att_sp) && !is.na(se_sp) && se_sp > 0) att_sp / se_sp else NA_real_
  stars_sp <- if (is.na(t_sp)) "" else
    if (abs(t_sp) > 2.576) "***" else
    if (abs(t_sp) > 1.960) "**"  else
    if (abs(t_sp) > 1.645) "*"   else ""

  pt_sp <- compute_pretrend_test(agg_d_sp, gt_sp_analytical)

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
# Poisson in levels keeps the exact zeros (counted in the log above) and its
# coefficient is a proportional (semi-elasticity) effect, invariant to the units
# of the outcome. fepois drops perfectly-separated observations, so its N can be
# below the panel's; we report whatever fepois uses (model$nobs).
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

if (is.null(ppml_att_row)) stop("Units/zeros PPML estimation failed.")
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
  BITERS, " reps, seed 1242); 5: OR, boot SE, positive subsample. Row 1: the headline ",
  "estimate (Table~\\ref{tab:combined_wide_main}). Row 6: Sun--Abraham ",
  "(2021) PPML, levels, country+year FE, clustered SE; retains zeros; no ",
  "comparable pre-trend statistic. Never-treated controls; WGI GE + log ",
  "population. Commitments in USD millions (raw USD $= \\times 10^6$). ",
  wpval_reconciliation(pre_egt = uz_pt_leads, df_did = uz_pt_df_did,
                       n_clusters = uz_n_country, wpval_did = uz_wpval,
                       pval_wald = uz_pval, wpval_reason = uz_reason,
                       ginv_used = uz_ginv, df = uz_pt_df),
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
#' (Copy of the helper in 03_main_results.R §7b; keep the two in sync.)
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
#' A CI contains zero when lb <= 0 <= ub; testing lb alone would call an
#' interval lying wholly below zero a breakdown.
#' @param lb,ub numeric vectors of robust lower/upper bounds, ordered as `grid`
#' @param grid numeric vector of M or Mbar values
#' @return formatted LaTeX string
breakdown_value <- function(lb, ub, grid) {
  if (all(is.na(lb) | is.na(ub))) return("---")
  idx <- which(lb <= 0 & ub >= 0)[1]
  if (is.na(idx)) return(sprintf("$> %.2f$", max(grid)))
  if (max(grid) < 0.1) sprintf("$%.3f$", grid[idx]) else sprintf("$%.2f$", grid[idx])
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
    agg_k <- fit$agg_dynamic  # read_headline_fit() guarantees it is present
    stopifnot(identical(as.numeric(agg_k$min_e), -5))
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
  # Hard stop: the sensitivity analysis must target the estimand Table 2
  # prints. If the group-size weights do not reproduce the reported simple
  # ATT, the exhibit would silently describe a different parameter.
  l_vec <- fit$l_vec_simple
  if (length(l_vec) != nPost)
    stop("l_vec_simple has ", length(l_vec), " entries for ", nPost, " post-treatment event times")
  gap <- sum(l_vec * betahat[egt >= 0]) - fit$att
  stopifnot(abs(gap) < 1e-6)

  list(betahat = betahat, sigma = bvc$sigma, numPre = nPre, numPost = nPost,
       l_vec = l_vec, egt = egt, max_rel_dev = bvc$max_rel_dev,
       reaggregated = reagg)
}

short_labels_hd <- c("log(Adapt.)", "Share (pp)", "log(Total)",
                     "log(Non-adapt.)", "log(Disb.)")

# One slot per outcome (and per outcome x lead count for hd_pre).
hd_all  <- setNames(vector("list", length(outcomes)), outcome_vars)  # primary RM results (numPre = 4)
hd_pre_keys <- unlist(lapply(outcomes, function(oc)
  paste0(oc$var, "_", if (oc$var %in% PRIMARY_OUTCOMES) numPre_grid else numPre_grid[1L])))
hd_pre  <- setNames(vector("list", length(hd_pre_keys)), hd_pre_keys)  # RM breakdown by numPrePeriods
hd_sd   <- setNames(vector("list", length(outcomes)), outcome_vars)  # smoothness-restriction results

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
    if (is.null(hin)) stop("HonestDiD inputs unavailable for ", oc$var, ", numPre = ", k)
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
    if (is.null(rm_sens)) stop("HonestDiD RM sensitivity failed for ", oc$var, ", numPre = ", k)

    if (!is.null(rm_sens)) {
      bd <- breakdown_value(rm_sens$lb, rm_sens$ub, Mbarvec_hd[seq_len(nrow(rm_sens))])
      message(sprintf("  RM breakdown Mbar (numPre = %d) = %s", k, bd))
      hd_pre[[paste0(oc$var, "_", k)]] <- list(var = oc$var, numPre = k,
                                               rm = rm_sens, breakdown = bd)
      if (k == numPre_grid[1L]) {
        hd_all[[oc$var]] <- list(label = oc$label, rm = rm_sens, egt = hin$egt,
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
      if (is.null(sd_sens)) stop("HonestDiD SD sensitivity failed for ", oc$var)
      if (!is.null(sd_sens)) {
        hd_sd[[oc$var]] <- list(label = oc$label, sd = sd_sens,
                                breakdown = breakdown_value(sd_sens$lb, sd_sens$ub, Mvec_sd),
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

  # Event-time windows of the fit behind the table (first outcome), not typed.
  egt_hd  <- hd_all[[outcomes[[1L]]$var]]$egt
  post_hd <- egt_hd[egt_hd >= 0]
  pre_hd  <- egt_hd[egt_hd < 0]
  notes_txt <- paste0(
    "Sensitivity of \\citet{rambachan_roth_2023} on the headline spec (Table~",
    "\\ref{tab:combined_wide_main}): CS\\,(2021) DR, never-treated controls, ",
    "multiplier-bootstrap SE (", BITERS, " reps, seed 1242); event-study ",
    "estimates and covariance of that headline estimate. Target: group-size-weighted ",
    "average post-treatment effect ($e=", min(post_hd), ",\\dots,", max(post_hd),
    "$), weights ", l_vec_txt,
    ", reproducing the simple ATT. $\\Sigma$ rescaled to reported bootstrap SEs ",
    "(max gap ", mrd_txt, "). Pre-treatment window: ", length(pre_hd), " leads ($e=",
    min(pre_hd), ",\\dots,", max(pre_hd), "$). ",
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
if (!all(vapply(hd_pre, is.null, logical(1L)))) {
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

  # Lead counts, windows and outcome lists come from the objects, not typed.
  np_vars  <- setdiff(outcome_vars, PRIMARY_OUTCOMES)
  lab_of   <- function(v) {
    l <- short_labels_hd[match(v, outcome_vars)]
    if (length(l) < 2L) l else paste(paste(head(l, -1L), collapse = ", "), "and", tail(l, 1L))
  }
  # TRUE if the robust CI at the smallest Mbar (0) already contains zero.
  np_break0 <- vapply(np_vars, function(v) {
    r <- hd_pre[[paste0(v, "_", numPre_grid[1L])]]$rm
    r$lb[1L] <= 0 && r$ub[1L] >= 0
  }, logical(1L))
  stopifnot(Mbarvec_hd[1L] == 0)
  notes_pre <- paste0(
    "Each cell: breakdown value $\\bar{M}$, disciplined by ",
    paste(head(numPre_grid, -1L), collapse = ", "), " or ", tail(numPre_grid, 1L),
    " pre-treatment event times. Estimator/covariance as Table~\\ref{tab:honestdid}; ",
    "only the first event time of the dynamic aggregation changes (",
    paste(sprintf("$%d$", -(numPre_grid + 1L)), collapse = ", "), "). ",
    numPre_grid[1L], "-lead row: the event-study estimates behind ",
    "Table~\\ref{tab:combined_wide_main}. ",
    paste(numPre_grid[-1L], collapse = "-/"), "-lead rows: the same estimates aggregated ",
    "over a longer pre-treatment window (seed 1242), for ", lab_of(PRIMARY_OUTCOMES),
    " only; point estimates unchanged, SEs a new bootstrap draw. ",
    if (all(np_break0))
      paste0(lab_of(np_vars), " already break at $\\bar{M} = 0$ with ", numPre_grid[1L],
             " leads, so their longer-window cells are left blank")
    else
      paste0("Longer windows are not computed for ", lab_of(np_vars),
             " (cells left blank)")
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
if (!all(vapply(hd_sd, is.null, logical(1L)))) {
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
    paste(sd_methods, collapse = "/"), " method, 5\\% level. Estimator, estimates, ",
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
# case (see the note immediately below). The reported placebo fit is the
# bootstrap one; the analytical twin exists solely to feed the pre-trend test.
em_plac <- "reg"
# The multiplier bootstrap is used whatever the number of treated units (a
# 40-unit cut-off is a convention, not an estimator requirement, and thin
# cells are where analytical influence-function SEs are most
# anti-conservative), with the replication count passed explicitly (BITERS).
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
  if (is.null(agg_s_plac) || is.null(agg_d_plac) || is.null(gt_placebo_analytical))
    stop("Placebo estimation failed: no table or figure written.")
  agg_d_plac_analytical <- aggte(gt_placebo_analytical, type = "dynamic", na.rm = TRUE,
                                 min_e = -4, max_e = Inf)

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
    # Simultaneous (sup-t) 95% band over the plotted event times.
    cv_p <- sup_t_crit(agg_d_plac$inf.function$dynamic.inf.func.e, agg_d_plac$se.egt,
                       biters = BITERS)
    message(sprintf("  Placebo ES sup-t crit (log adaptation commitments): %.4f", cv_p))
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

    pt_plac <- compute_pretrend_test(agg_d_plac_analytical, gt_placebo_analytical)

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
                           ginv_used = pt_plac$ginv_used, df = pt_plac$df),
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
  stop("Placebo estimation failed: no table or figure written.")
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
  stop("mitigation_panel.csv not found at ", mit_panel_path,
       " (it is produced by 01_prepare_data.R).")
} else {

  mitigation_ry <- as.data.frame(fread(mit_panel_path))

  required_mit_cols <- c("recipient_name", "year", "commits_mitigation")
  missing_mit <- setdiff(required_mit_cols, names(mitigation_ry))
  if (length(missing_mit) > 0) {
    stop("mitigation_panel.csv missing columns: ", paste(missing_mit, collapse = ", "))
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
    # Headline estimator, fixed (no treated-count switch): the falsification
    # test must use the same estimator as the headline, so a failed DR fit
    # stops the stage (no fallback to another estimator).
    em_mit <- "dr"
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
        NULL
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
      if (is.null(agg_s_mit) || is.null(agg_d_mit) || is.null(gt_mit_analytical))
        stop("Mitigation estimation failed: no table or figure written.")
      agg_d_mit_analytical <- aggte(gt_mit_analytical, type = "dynamic", na.rm = TRUE,
                                    min_e = -5, max_e = Inf)

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
        agg_s_mit_analytical <- aggte(gt_mit_analytical, type = "simple", na.rm = TRUE)
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
        # Simultaneous (sup-t) 95% band over the plotted event times.
        cv_m <- sup_t_crit(agg_d_mit$inf.function$dynamic.inf.func.e, agg_d_mit$se.egt,
                           biters = BITERS)
        message(sprintf("  Mitigation ES sup-t crit (log mitigation commitments): %.4f", cv_m))
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

        pt_mit <- compute_pretrend_test(agg_d_mit_analytical, gt_mit_analytical)

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
                               ginv_used = pt_mit$ginv_used, df = pt_mit$df),
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
      stop("Mitigation estimation failed: no table or figure written.")
    }
  }
}

message("\n=== Section 26 complete ===\n")

# ==============================================================================
# SECTION A. §27 de Chaisemartin & D'Haultfoeuille (2026, REStat) estimator
# Alternative heterogeneity-robust DID estimator on the headline outcome
# (log adaptation commitments).  dCDH uses the SAME controls as the CS main
# spec (WGI gov. effectiveness + log population), so the two are comparable.
# Panel: thin-dropped main panel.  Treatment is the current-treatment dummy
# treated_post = 1{cohort > 0 & year >= cohort}; never-treated and not-yet-treated
# units serve as controls.  4 dynamic effects, 3 placebos.  This is an APPENDIX
# exhibit.
# Seed rule (reproduces the published SEs): set.seed(1242) immediately before the estimator call.
# ==============================================================================

message("\n=== Section A (27): de Chaisemartin & D'Haultfoeuille (2026) ===\n")

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
  stop("dCDH estimation failed: no table or figure written.")
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
    "de Chaisemartin--D'Haultf{\\oe}uille (2026) estimator, current-treatment dummy; ",
    "never- and not-yet-treated controls; same controls as the CS main specification. ",
    "Effect $+\\ell$ is event time $e = \\ell - 1$ of the Callaway--Sant'Anna event studies ",
    "($\\ell = 1$ is the adoption year, $e = 0$), and Placebo $-\\ell$ is $e = -1 - \\ell$. ",
    "Last row: package's own joint test that all ", n_plac, " placebo estimates are ",
    "zero ($p$-value in the CI column), accounting for the ",
    "covariance across placebo horizons. The package (DIDmultiplegtDYN ",
    as.character(utils::packageVersion("DIDmultiplegtDYN")),
    ") reports only the $p$-value of this test, so the Estimate column is blank",
    if (is.na(nm_pjoint)) " (not returned by this run)" else ""
  )
  source_dcdh <- "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"

  write_tex_float(
    out_path      = file.path(dir_tabs_dcdh, "att_dcdh.tex"),
    caption_title = "Robustness: de Chaisemartin and D'Haultf{\\oe}uille (2026) estimator",
    label         = "tab:dcdh",
    tabular_lines = raw_dcdh,
    notes_text    = notes_dcdh,
    source_text   = source_dcdh
  )

  message(sprintf("  dCDH Av. total effect = %.3f (SE %.3f, CI [%.3f, %.3f])",
                  est_vec[1L], se_vec[1L], lb_vec[1L], ub_vec[1L]))

  # --- §A event-study figure: built MANUALLY from Effects (+1..+4) and Placebos
  # (-1..-3), plotted on Figure 2's event-time axis e = t - F_g. In
  # did_multiplegt_dyn(), Effect_l compares F_g - 1 + l with F_g - 1 (Effect_1 is
  # the adoption period F_g) and Placebo_k compares F_g - 1 - k with F_g - 1
  # (checked on simulated data with known effects). So Effect_l sits at
  # e = l - 1, the normalized reference period F_g - 1 at e = -1, and Placebo_k
  # at e = -1 - k. Ribbon/linerange uses the package 95% CI (LB CI / UB CI).
  # No in-figure title/subtitle/caption.
  es_dcdh <- data.frame(
    event_time = c(-1L - rev(seq_len(n_plac)), -1L, seq_len(n_eff) - 1L),
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
# A failure stops the stage (no table or figure is written).
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
  stop("Goodman-Bacon decomposition failed: no table or figure written.")
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
    "Decomposition of the static two-way fixed-effects estimate of a post-adoption ",
    "indicator (1 in the years after NAP adoption) on log adaptation commitments into ",
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

# ==============================================================================
# SECTION D. Balanced panel
# The estimation panel is unbalanced (South Sudan enters in 2011; recipients
# that left the DAC List have no CRS record afterwards), so did's
# allow_unbalanced_panel path estimates it as repeated cross-sections, with
# covariates at their current-year values. Restricting to recipients observed
# in every year makes did use its panel estimator, with covariates at g-1.
# ==============================================================================

message("\n=== Section D: balanced panel ===\n")

n_years_panel <- n_distinct(did_panel_full$year)
did_panel_balanced <- did_panel_full %>%
  filter(!(cohort_year %in% thin_cohorts)) %>%
  group_by(recipient_name) %>%
  filter(n() == n_years_panel) %>%
  ungroup()
message("Recipients dropped (not observed every year): ",
        paste(setdiff(unique(did_panel_full$recipient_name[!(did_panel_full$cohort_year %in% thin_cohorts)]),
                      unique(did_panel_balanced$recipient_name)), collapse = ", "))
stopifnot(nrow(did_panel_balanced) == n_distinct(did_panel_balanced$recipient_name) * n_years_panel)

make_wide_table(
  did_panel_in  = did_panel_balanced,
  retain_thin   = FALSE,
  outcomes      = outcomes,
  dir_tabs      = file.path(here("output", "tables"), "balanced_panel"),
  tex_label     = "tab:combined_wide_balanced",
  caption_spec  = "main specification, balanced panel, covariates at $g-1$"
)

# ==============================================================================
# SECTION E. Panel starting in 2010
# 2009 precedes routine reporting of the adaptation marker: almost every
# recipient records zero adaptation finance that year. It enters no
# post-treatment ATT (each cohort is compared with its own g-1 >= 2020), only
# pre-period cells; this block re-estimates the main table without it.
# ==============================================================================

message("\n=== Section E: panel 2010-2024 ===\n")

did_panel_2010 <- did_panel_full %>%
  filter(!(cohort_year %in% thin_cohorts), year >= 2010)
# Text number: recipients with any adaptation-marked commitment in 2009.
panel_2009 <- did_panel_full %>% filter(year == 2009L)
message(sprintf("2009: %d of %d panel recipients record any adaptation-marked commitment",
                sum(crs_positive(panel_2009$commitments), na.rm = TRUE),
                n_distinct(panel_2009$recipient_name)))

make_wide_table(
  did_panel_in  = did_panel_2010,
  retain_thin   = FALSE,
  outcomes      = outcomes,
  dir_tabs      = file.path(here("output", "tables"), "panel_2010"),
  tex_label     = "tab:combined_wide_2010",
  caption_spec  = "main specification, panel 2010--2024"
)
message("\n=== 04_robustness.R: complete ===\n")
