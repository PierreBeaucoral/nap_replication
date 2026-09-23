# ==============================================================================
# 10_model_tests.R
# Two testable implications of the CES allocation model (Section 3 of the manuscript).
#
# TEST 1 (Lemma 2, l.352-360): envelope elasticity d ln B_i*/d ln phi_i =
#   eta*(1-omega_i), where omega_i is recipient i's pre-treatment share of the
#   donor's total envelope. Under the gross-substitutes maintained assumption
#   used throughout Sec. 3 for the envelope channel (sigma_b > 1, eta > 0),
#   the envelope (total-commitments) ATT should be larger for small recipients
#   (low omega_i) than for large recipients (high omega_i).
#
# TEST 2 (Proposition 1 footnote, l.427): sign(d phi_i/d alpha_i) =
#   sign(alpha_i - 1/2). Because s_i* = s_i*(alpha_i) is a strictly increasing
#   bijection with s_i*(alpha_i=1/2) = 1/2 for ANY sigma_c > 0 (Lemma 1 gives
#   s_i* = q_i/(1+q_i), q_i = (alpha_i/(1-alpha_i))^sigma_c, and q_i = 1 when
#   alpha_i = 1/2 regardless of sigma_c), alpha_i >= 1/2 iff s_i* >= 1/2 EXACTLY.
#   If the empirical proxy for s_i* (pre-treatment within-country adaptation
#   share) is far below 1/2 for every recipient, the model predicts a weakly
#   NEGATIVE envelope response to NAP adoption (again under eta > 0), i.e. the
#   ATT on total commitments should be <= 0.
#
# Inputs : data/processed/simple_panel_wgi.csv
# Outputs:
#   output/tables/model_tests/lemma2_size.tex     (tab:lemma2_size)
#   output/tables/model_tests/alpha_half.tex      (tab:alpha_half)
#   output/figures/model_tests/fig_lemma2_size.png
#   -> copied to paper/Tables/model_tests/ and paper/Figures/model_tests/
#
# This script is standalone (not called by run_all.R) and self-contained given
# that 01_prepare_data.R has produced data/processed/simple_panel_wgi.csv.
# ==============================================================================

# ============================================================
# Paper-to-Code Naming Map
# ============================================================
# Paper Notation        | Code Name            | Description
# $B_i^*$ / $S_i$       | lcommitments_all      | log(total commitments) — envelope
# $E_i$                 | log_commits           | log(adaptation commitments)
# $O_i$                 | lcommitments_nonadapt | log(non-adaptation commitments)
# Share$_i^E$           | share_adapt           | recipient's % share of global
#                        |                       | adaptation pool
# $\omega_i = B_i^*/B$  | omega_i               | pre-treatment (2009-2012 mean)
#                        |                       | share of commitments_all across
#                        |                       | the 144-panel countries
# $s_i^*(\alpha_i)$     | s_i                   | within-country adaptation share,
#                        |                       | commitments / commitments_all
# $\alpha_i$             | (not observed)       | proxied via sign(s_i - 0.5), exact
#                        |                       | equivalence per Lemma 1 (see header)
# $\eta=\rho_b/(1-\rho_b)$| (not estimated)      | sign unidentified in this reduced-
#                        |                       | form design; eta > 0 (gross
#                        |                       | substitutes) is the maintained
#                        |                       | assumption inherited from
#                        |                       | Proposition 1's own envelope test
# $G_i$                  | cohort_year          | NAP adoption cohort (0 = never)
# $ATT(g,t)$             | gt_obj / gt_g        | Group-time ATT from att_gt()
# $\Delta=\theta_a-\theta_b$| diff              | Subgroup ATT difference
# ============================================================

library(data.table)
library(dplyr)
library(ggplot2)
library(here)
library(did)

# Version guard: SEs on unbalanced panels changed in did 2.5.0 (renv.lock pins it).
if (utils::packageVersion("did") < "2.5.0") {
  stop(sprintf(paste0(
    "did %s is installed but this script's results require did >= 2.5.0.\n",
    "  Fix: run renv::restore() from the project root (renv.lock pins did 2.5.0)."),
    utils::packageVersion("did")))
}

Sys.setenv(RGL_USE_NULL = "TRUE")  # headless-safe; project convention, no rgl-dependent pkg used here

set.seed(20240601)  # global seed — local set.seed(1242) calls follow each bootstrap fit

# -----------------------------------------------------------------------
# Output directories
# -----------------------------------------------------------------------
dir_tabs <- here("output", "tables",  "model_tests")
dir_figs <- here("output", "figures", "model_tests")
dir.create(dir_tabs, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figs, recursive = TRUE, showWarnings = FALSE)

has_paper  <- dir.exists(here("paper"))  # FALSE in the stand-alone replication package
paper_tabs <- here("paper", "Tables",  "model_tests")
paper_figs <- here("paper", "Figures", "model_tests")
if (has_paper) {
  dir.create(paper_tabs, recursive = TRUE, showWarnings = FALSE)
  dir.create(paper_figs, recursive = TRUE, showWarnings = FALSE)
}

# ==============================================================================
# SECTION 1. Helpers (copied/adapted from code/03_main_results.R and
# code/05_heterogeneity.R per project convention — kept local, no shared
# functions/ directory in this repo).
# ==============================================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

# canonical pre-trend wording — keep byte-identical across scripts
PRETREND_NOTE_AGG <- function(min_e, max_e, k) sprintf(
  "Pre-trend $\\chi^2$ / $p$: joint Wald test on the aggregated pre-treatment event-time coefficients ($%d \\leq e \\leq %d$; %d restrictions), using the influence-function covariance of the dynamic aggregation from the analytical (non-bootstrap) fit; a generalized inverse is used if the block is singular",
  min_e, max_e, k)
PRETREND_NOTE_DID <- "\\texttt{did} pre-test $p$: \\texttt{did}'s built-in Wald test over all pre-period $ATT(g,t)$ cells against each cohort's $g-1$ base year"
PRETREND_NOTE <- function(min_e, max_e, k) paste0(PRETREND_NOTE_AGG(min_e, max_e, k), ". ", PRETREND_NOTE_DID)

# --- write_tex_float(): identical to 03/05's helper -------------------------
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

# --- compute_pretrend_test(): identical to 03's revised version -------------
# Joint Wald test on pre-treatment event-time ATTs, using the true dynamic-
# aggregation influence-function covariance (crossprod(IF)/n^2), read from
# agg_d$inf.function$dynamic.inf.func.e. did's own gt_obj$W / gt_obj$Wpval
# (group-time pre-test) is captured alongside for comparison.
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

# --- att_difference_test(): identical to 05's helper ------------------------
#' Wald test that two subgroup ATTs are equal.
#' Disjoint-sample independence: Var(theta_a - theta_b) = Var(theta_a) + Var(theta_b).
att_difference_test <- function(stats_a, stats_b) {
  stopifnot(is.list(stats_a), is.list(stats_b))
  att_a <- stats_a$att
  att_b <- stats_b$att
  se_a  <- stats_a$se
  se_b  <- stats_b$se
  if (any(is.na(c(att_a, att_b, se_a, se_b))))
    return(list(diff = NA_real_, se = NA_real_, z = NA_real_, pval = NA_real_))

  diff    <- att_a - att_b
  se_diff <- sqrt(se_a^2 + se_b^2)          # independence across disjoint samples
  z       <- if (se_diff > 1e-12) diff / se_diff else NA_real_
  pval    <- if (is.na(z)) NA_real_ else 2 * pnorm(-abs(z))
  list(diff = diff, se = se_diff, z = z, pval = pval)
}

# --- fit_subgroup_att(): mirrors 05's make_het_wide_table() convention ------
# est_method = "dr" if the subgroup has >= 40 treated units, else "reg";
# bstrap = FALSE always (analytical / influence-function SE), not-yet-treated
# controls. Matches code/05_heterogeneity.R l.596-640 exactly.
fit_subgroup_att <- function(panel_g, outcome_var, thin_cohorts_vec,
                              xformla = ~ ge_est + log_population,
                              control_group = "notyettreated") {
  stopifnot(is.data.frame(panel_g), outcome_var %in% names(panel_g))
  panel_g <- panel_g[!(panel_g$cohort_year %in% thin_cohorts_vec), ]

  n_treated_g <- n_distinct(panel_g$country_id[panel_g$cohort_year > 0])
  n_cohorts_g <- n_distinct(panel_g$cohort_year[panel_g$cohort_year > 0])
  if (n_cohorts_g < 2L) {
    message("  fit_subgroup_att: fewer than 2 treated cohorts — skipping")
    return(list(att = NA_real_, se = NA_real_, n_treated = n_treated_g,
                n_obs = NA_integer_, n_country = NA_integer_,
                em = NA_character_, agg_d = NULL, pretrend = NULL))
  }

  em_g <- if (n_treated_g >= 40L) "dr" else "reg"
  message(sprintf("  fit_subgroup_att[%s]: N treated = %d -> est_method = %s",
                  outcome_var, n_treated_g, em_g))

  gt_g <- tryCatch(
    att_gt(
      yname = outcome_var, tname = "year", idname = "country_id",
      gname = "cohort_year", xformla = xformla, data = panel_g,
      est_method = em_g, bstrap = FALSE, cband = FALSE,
      control_group = control_group, anticipation = 0,
      base_period = "universal", panel = TRUE, allow_unbalanced_panel = TRUE
    ),
    error = function(e) { message("  att_gt failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(gt_g)) {
    return(list(att = NA_real_, se = NA_real_, n_treated = n_treated_g,
                n_obs = NA_integer_, n_country = NA_integer_,
                em = em_g, agg_d = NULL, pretrend = NULL))
  }

  agg_s <- tryCatch(aggte(gt_g, type = "simple",  na.rm = TRUE), error = function(e) NULL)
  agg_d <- tryCatch(aggte(gt_g, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf),
                    error = function(e) NULL)

  att <- if (!is.null(agg_s)) agg_s$overall.att else NA_real_
  se  <- if (!is.null(agg_s)) agg_s$overall.se  else NA_real_
  pt  <- if (!is.null(agg_d)) tryCatch(compute_pretrend_test(agg_d), error = function(e) NULL) else NULL

  rows_g <- panel_g[!is.na(panel_g[[outcome_var]]), ]

  list(
    att = att, se = se, n_treated = n_treated_g,
    n_obs = nrow(rows_g), n_country = n_distinct(rows_g$country_id),
    em = em_g, agg_d = agg_d, pretrend = pt
  )
}

# --- fit_main_spec(): main-spec CS fit, matching code/03_main_results.R ----
# bootstrap fit for the ATT (bstrap=TRUE, biters=999L) + a parallel analytical
# (bstrap=FALSE) fit for the pre-trend Wald test. est_method="dr",
# control_group="nevertreated", base_period="universal" (cohorts_dropped spec).
fit_main_spec <- function(panel_in, outcome_var,
                          xformla = ~ ge_est + log_population) {
  stopifnot(is.data.frame(panel_in), outcome_var %in% names(panel_in))

  set.seed(1242)  # Seed rule (reproduces the published SEs): immediately before each bootstrap fit
  gt_boot <- tryCatch(
    att_gt(
      yname = outcome_var, tname = "year", idname = "country_id",
      gname = "cohort_year", xformla = xformla, data = panel_in,
      est_method = "dr", bstrap = TRUE, biters = 999L, cband = FALSE,
      control_group = "nevertreated", anticipation = 0,
      base_period = "universal", panel = TRUE, allow_unbalanced_panel = TRUE
    ),
    error = function(e) { message("  att_gt (bootstrap) failed: ", conditionMessage(e)); NULL }
  )
  agg_s <- if (!is.null(gt_boot))
    tryCatch(aggte(gt_boot, type = "simple", na.rm = TRUE), error = function(e) NULL) else NULL

  set.seed(1242)  # Seed rule (reproduces the published SEs): immediately before each bootstrap fit
  gt_analytical <- tryCatch(
    att_gt(
      yname = outcome_var, tname = "year", idname = "country_id",
      gname = "cohort_year", xformla = xformla, data = panel_in,
      est_method = "dr", bstrap = FALSE, cband = FALSE,
      control_group = "nevertreated", anticipation = 0,
      base_period = "universal", panel = TRUE, allow_unbalanced_panel = TRUE
    ),
    error = function(e) { message("  att_gt (analytical) failed: ", conditionMessage(e)); NULL }
  )
  agg_d <- if (!is.null(gt_analytical))
    tryCatch(aggte(gt_analytical, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf),
             error = function(e) NULL) else NULL
  pt <- if (!is.null(agg_d)) tryCatch(compute_pretrend_test(agg_d, gt_obj = gt_analytical),
                                       error = function(e) NULL) else NULL

  rows_oc <- panel_in[!is.na(panel_in[[outcome_var]]), ]

  list(
    att = if (!is.null(agg_s)) agg_s$overall.att else NA_real_,
    se  = if (!is.null(agg_s)) agg_s$overall.se  else NA_real_,
    n_obs = nrow(rows_oc), n_country = n_distinct(rows_oc$country_id),
    pretrend = pt
  )
}

# ==============================================================================
# SECTION 2. Load and prepare the DiD panel
# (verbatim from code/03_main_results.R Section 1, l.132-180)
# ==============================================================================

message("\n=== 10_model_tests.R: loading panel ===\n")

aggregated <- fread(here("data", "processed", "simple_panel_wgi.csv"))
aggregated  <- as.data.frame(aggregated)

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

stopifnot(
  n_distinct(did_panel_full$recipient_name) >= readRDS(here("output", "fits", "headline_adaptation_dr_bs.rds"))$n_country,
  sum(did_panel_full$always_treated) == 0L
)

cohort_sizes <- did_panel_full %>%
  filter(cohort_year > 0) %>%
  distinct(recipient_name, cohort_year) %>%
  count(cohort_year, name = "n_treated") %>%
  arrange(cohort_year)

thin_threshold <- 5L
thin_cohorts   <- cohort_sizes$cohort_year[cohort_sizes$n_treated < thin_threshold]
message("Thin cohorts (< ", thin_threshold, " treated units, dropped): ",
        paste(thin_cohorts, collapse = ", "))

# cohorts_dropped estimation sample (matches Table 2 in the paper: 126 countries)
did_panel_main <- did_panel_full %>% filter(!(cohort_year %in% thin_cohorts))
n_est_sample   <- n_distinct(did_panel_main$recipient_name)
message(sprintf("Estimation sample (cohorts_dropped): %d countries", n_est_sample))
stopifnot(n_est_sample == readRDS(here("output", "fits", "headline_adaptation_dr_bs.rds"))$n_country)

# ==============================================================================
# SECTION 3. Recipient-size measure omega_i and composition share s_i
# ==============================================================================

message("\n=== Section 3: omega_i and s_i construction ===\n")

PRE_YEARS <- 2009:2012

# --- omega_i: pre-treatment share of total commitments across the 144-panel
#     countries (proxy for B_i*/B, recipient i's share of the donor's envelope).
#     Alternative size proxy: pre-treatment share of the global adaptation pool.
size_measures <- did_panel_full %>%
  filter(year %in% PRE_YEARS) %>%
  group_by(recipient_name, country_id) %>%
  summarise(
    mean_commitments_all = mean(commitments_all, na.rm = TRUE),
    mean_share_adapt     = mean(share_adapt,     na.rm = TRUE),  # already in pp, sums to 100/yr
    .groups = "drop"
  ) %>%
  mutate(
    omega_i     = mean_commitments_all / sum(mean_commitments_all, na.rm = TRUE),
    omega_i_alt = mean_share_adapt / 100  # alt size proxy: share of global adaptation pool
  )

stopifnot(!anyNA(size_measures$omega_i), all(is.finite(size_measures$omega_i)))
message(sprintf(
  "omega_i (share of total commitments, N=%d): mean=%.5f median=%.5f sd=%.5f min=%.5f max=%.5f",
  nrow(size_measures), mean(size_measures$omega_i), median(size_measures$omega_i),
  sd(size_measures$omega_i), min(size_measures$omega_i), max(size_measures$omega_i)))
message(sprintf(
  "omega_i_alt (share of global adaptation pool, N=%d): mean=%.5f median=%.5f sd=%.5f min=%.5f max=%.5f",
  nrow(size_measures), mean(size_measures$omega_i_alt), median(size_measures$omega_i_alt),
  sd(size_measures$omega_i_alt), min(size_measures$omega_i_alt), max(size_measures$omega_i_alt)))
message(sprintf("corr(omega_i, omega_i_alt) = %.3f",
                cor(size_measures$omega_i, size_measures$omega_i_alt)))

# --- s_i: within-country adaptation share, ratio of 4-year means -----------
s_i_pre <- did_panel_full %>%
  filter(year %in% PRE_YEARS) %>%
  group_by(recipient_name) %>%
  summarise(
    mean_commitments     = mean(commitments,     na.rm = TRUE),
    mean_commitments_all = mean(commitments_all, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(s_i = mean_commitments / mean_commitments_all)

stopifnot(all(is.finite(s_i_pre$s_i)))
n_si_pre <- nrow(s_i_pre)

s_i_g1 <- did_panel_full %>%
  filter(cohort_year > 0) %>%
  distinct(recipient_name, cohort_year) %>%
  mutate(g1_year = cohort_year - 1L) %>%
  left_join(
    did_panel_full %>% select(recipient_name, year, commitments, commitments_all),
    by = c("recipient_name", "g1_year" = "year")
  ) %>%
  mutate(s_i_g1 = commitments / commitments_all)

n_si_g1_valid <- sum(!is.na(s_i_g1$s_i_g1) & is.finite(s_i_g1$s_i_g1))

s_i_dist <- function(x) {
  x <- x[is.finite(x)]
  list(n = length(x), mean = mean(x), median = median(x),
       p90 = as.numeric(quantile(x, 0.90)), p99 = as.numeric(quantile(x, 0.99)),
       max = max(x), n_ge_half = sum(x >= 0.5))
}
dist_pre <- s_i_dist(s_i_pre$s_i)
dist_g1  <- s_i_dist(s_i_g1$s_i_g1)

message(sprintf(
  "s_i (2009-2012 mean, N=%d): mean=%.4f median=%.4f p90=%.4f p99=%.4f max=%.4f N(s_i>=0.5)=%d",
  dist_pre$n, dist_pre$mean, dist_pre$median, dist_pre$p90, dist_pre$p99,
  dist_pre$max, dist_pre$n_ge_half))
message(sprintf(
  "s_i at g-1 (adopters, N=%d valid): mean=%.4f median=%.4f p90=%.4f p99=%.4f max=%.4f N(s_i>=0.5)=%d",
  dist_g1$n, dist_g1$mean, dist_g1$median, dist_g1$p90, dist_g1$p99,
  dist_g1$max, dist_g1$n_ge_half))

# ==============================================================================
# SECTION 4. TEST 1 — Lemma 2: envelope response decreasing in recipient size
# ==============================================================================

message("\n=== Section 4: Test 1 (Lemma 2, recipient size) ===\n")

# Median/tercile computed over the 126-country estimation sample, as specified.
est_sample_countries <- did_panel_main %>% distinct(recipient_name)
omega_est <- size_measures %>%
  inner_join(est_sample_countries, by = "recipient_name")
stopifnot(nrow(omega_est) == n_est_sample, !anyNA(omega_est$omega_i))

omega_median <- median(omega_est$omega_i)
omega_terciles <- quantile(omega_est$omega_i, probs = c(1 / 3, 2 / 3))
message(sprintf("omega_i median (126-country sample) = %.6f", omega_median))
message(sprintf("omega_i tercile cuts = %.6f, %.6f", omega_terciles[1], omega_terciles[2]))

omega_est <- omega_est %>%
  mutate(
    size_half  = if_else(omega_i <= omega_median, "Small", "Large"),
    size_terc  = case_when(
      omega_i <= omega_terciles[1] ~ "T1",
      omega_i <= omega_terciles[2] ~ "T2",
      TRUE                         ~ "T3"
    )
  )

half_counts <- omega_est %>% count(size_half)
terc_counts <- omega_est %>% count(size_terc)
message("Median-split country counts: ",
        paste(sprintf("%s=%d", half_counts$size_half, half_counts$n), collapse = ", "))
message("Tercile country counts: ",
        paste(sprintf("%s=%d", terc_counts$size_terc, terc_counts$n), collapse = ", "))

small_countries <- omega_est$recipient_name[omega_est$size_half == "Small"]
large_countries <- omega_est$recipient_name[omega_est$size_half == "Large"]
panel_small <- did_panel_full %>% filter(recipient_name %in% small_countries)
panel_large <- did_panel_full %>% filter(recipient_name %in% large_countries)

t1_outcomes <- c(
  lcommitments_all = "log(Total commitments)",
  log_commits       = "log(Adaptation commitments)",
  share_adapt       = "Adaptation share (\\% of global)"
)

set.seed(1242)  # Seed rule (reproduces the published SEs): before the subgroup-fit loop (05 convention)
t1_fits <- lapply(names(t1_outcomes), function(oc) {
  list(
    outcome = oc,
    small   = fit_subgroup_att(panel_small, oc, thin_cohorts),
    large   = fit_subgroup_att(panel_large, oc, thin_cohorts)
  )
})
names(t1_fits) <- names(t1_outcomes)

t1_diff <- lapply(t1_fits, function(f) att_difference_test(f$small, f$large))

for (oc in names(t1_outcomes)) {
  f <- t1_fits[[oc]]; d <- t1_diff[[oc]]
  message(sprintf(
    "  [%s] Small ATT=%.4f (SE=%.4f, N_tr=%d) | Large ATT=%.4f (SE=%.4f, N_tr=%d) | Diff=%.4f SE=%.4f z=%.3f p=%.4f",
    oc, f$small$att, f$small$se, f$small$n_treated,
    f$large$att, f$large$se, f$large$n_treated,
    d$diff, d$se, d$z, d$pval))
}

# --- Continuous version: ATT on lcommitments_all by tercile of omega_i -----
terc_countries <- split(omega_est$recipient_name, omega_est$size_terc)
panel_terc <- lapply(terc_countries, function(cc) did_panel_full %>% filter(recipient_name %in% cc))

set.seed(1242)  # Seed rule (reproduces the published SEs): before the subgroup-fit loop
t1_terc_fits <- lapply(names(panel_terc), function(tg) {
  fit_subgroup_att(panel_terc[[tg]], "lcommitments_all", thin_cohorts)
})
names(t1_terc_fits) <- names(panel_terc)
t1_terc_fits <- t1_terc_fits[c("T1", "T2", "T3")]  # enforce order

for (tg in names(t1_terc_fits)) {
  f <- t1_terc_fits[[tg]]
  message(sprintf("  [tercile %s] ATT (log total commitments)=%.4f SE=%.4f N_tr=%d",
                  tg, f$att, f$se, f$n_treated))
}

# ==============================================================================
# SECTION 5. TEST 1 — table tab:lemma2_size
# ==============================================================================

fmt4 <- function(x) if (is.na(x)) "---" else sprintf("%.4f", x)
fmt3 <- function(x) if (is.na(x)) "---" else sprintf("%.3f", x)
fmt_p <- function(x) {
  if (is.na(x)) return("---")
  if (x < 0.001) return("$<$0.001")
  sprintf("%.3f", x)
}
paren4 <- function(x) if (is.na(x)) "" else sprintf("(%.4f)", x)

panel_a_rows <- character(0)
for (oc in names(t1_outcomes)) {
  f <- t1_fits[[oc]]; d <- t1_diff[[oc]]
  panel_a_rows <- c(panel_a_rows,
    paste0(t1_outcomes[[oc]], " & ", fmt4(f$small$att), " & ", fmt4(f$large$att),
           " & ", fmt4(d$diff), " & ", fmt4(d$se), " & ", fmt3(d$z), " & ", fmt_p(d$pval), " \\\\"),
    paste0(" & ", paren4(f$small$se), " & ", paren4(f$large$se), " & & & & \\\\")
  )
}

panel_a_n <- paste0(
  "\\quad $N$ treated (small/large) & ", t1_fits[[1]]$small$n_treated, " & ",
  t1_fits[[1]]$large$n_treated, " & & & & \\\\"
)
panel_a_ncty <- paste0(
  "\\quad $N$ countries (small/large) & ", t1_fits[[1]]$small$n_country, " & ",
  t1_fits[[1]]$large$n_country, " & & & & \\\\"
)

panel_b_hdr <- paste0(" & T1 (smallest) & T2 (middle) & T3 (largest) & & & \\\\")
panel_b_att <- paste0("ATT & ", fmt4(t1_terc_fits$T1$att), " & ", fmt4(t1_terc_fits$T2$att),
                      " & ", fmt4(t1_terc_fits$T3$att), " & & & \\\\")
panel_b_se  <- paste0(" & ", paren4(t1_terc_fits$T1$se), " & ", paren4(t1_terc_fits$T2$se),
                      " & ", paren4(t1_terc_fits$T3$se), " & & & \\\\")
panel_b_n   <- paste0("\\quad $N$ treated & ", t1_terc_fits$T1$n_treated, " & ",
                      t1_terc_fits$T2$n_treated, " & ", t1_terc_fits$T3$n_treated,
                      " & & & \\\\")

lemma2_tabular <- c(
  "\\begin{tabular}{lcccccc}",
  "\\toprule",
  " & Small $\\omega_i$ & Large $\\omega_i$ & Difference & SE & $z$ & $p$ \\\\",
  "\\midrule",
  "\\multicolumn{7}{l}{\\textit{Panel A: ATT by recipient-size half (median split of $\\omega_i$)}} \\\\",
  "\\midrule",
  panel_a_rows,
  "\\midrule",
  panel_a_n,
  panel_a_ncty,
  "\\\\[0.5em]",
  "\\multicolumn{7}{l}{\\textit{Panel B: ATT on log(Total commitments) by tercile of $\\omega_i$}} \\\\",
  "\\midrule",
  panel_b_hdr,
  panel_b_att,
  panel_b_se,
  panel_b_n,
  "\\bottomrule",
  "\\end{tabular}"
)

notes_lemma2 <- paste0(
  "$\\omega_i$ = pre-treatment (", min(PRE_YEARS), "--", max(PRE_YEARS),
  ") mean share of commitments\\_all across the ", nrow(size_measures),
  "-country panel; median/terciles from the ", n_est_sample,
  "-country estimation sample (see text for the Lemma~\\ref{lem:envelope_elast} prediction). Panel~A: CS\\,(2021) ",
  "OR, analytical (IF) SE, not-yet-treated controls, cohorts $<$ ", thin_threshold,
  " dropped; difference row SE $=\\sqrt{se_a^2+se_b^2}$ (disjoint, country-clustered), ",
  "$z$ vs.\\ standard normal, two-sided $p$. Panel~B: same specification by tercile, ",
  "log(total commitments) only"
)
source_lemma2 <- "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"

write_tex_float(
  out_path      = file.path(dir_tabs, "lemma2_size.tex"),
  caption_title = "Envelope response by recipient size (test of Lemma 2)",
  label         = "tab:lemma2_size",
  tabular_lines = lemma2_tabular,
  notes_text    = notes_lemma2,
  source_text   = source_lemma2
)

# ==============================================================================
# SECTION 6. TEST 1 — figure fig_lemma2_size.png (ATT by tercile, 95% CI)
# ==============================================================================

crit_95 <- qnorm(0.975)
terc_df <- data.frame(
  tercile = factor(c("T1\n(smallest)", "T2\n(middle)", "T3\n(largest)"),
                   levels = c("T1\n(smallest)", "T2\n(middle)", "T3\n(largest)")),
  att = c(t1_terc_fits$T1$att, t1_terc_fits$T2$att, t1_terc_fits$T3$att),
  se  = c(t1_terc_fits$T1$se,  t1_terc_fits$T2$se,  t1_terc_fits$T3$se)
) %>%
  mutate(lower = att - crit_95 * se, upper = att + crit_95 * se)

p_lemma2 <- ggplot(terc_df, aes(x = tercile, y = att)) +
  geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed", linewidth = 0.4) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.7, linewidth = 0.8,
                  colour = "#2E86C1") +
  labs(title = NULL, subtitle = NULL, caption = NULL,
       x = "Tercile of pre-treatment recipient size (omega_i, share of total commitments)",
       y = "ATT \u2014 log(Total commitments)") +
  theme_minimal(base_family = "serif", base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(dir_figs, "fig_lemma2_size.png"), p_lemma2,
       width = 7, height = 5, dpi = 300)
message("Saved: ", file.path(dir_figs, "fig_lemma2_size.png"))

# ==============================================================================
# SECTION 7. TEST 2 — sign(alpha_i - 1/2) via s_i, and the main-spec sign test
# ==============================================================================

message("\n=== Section 7: Test 2 (Proposition 1 footnote, sign(alpha-1/2)) ===\n")

fit_env  <- fit_main_spec(did_panel_main, "lcommitments_all")
fit_nona <- fit_main_spec(did_panel_main, "lcommitments_nonadapt")

message(sprintf("Main-spec ATT log(Total commitments)      = %.4f (SE = %.4f, N=%d, countries=%d)",
                fit_env$att, fit_env$se, fit_env$n_obs, fit_env$n_country))
message(sprintf("Main-spec ATT log(Non-adaptation commits.) = %.4f (SE = %.4f, N=%d, countries=%d)",
                fit_nona$att, fit_nona$se, fit_nona$n_obs, fit_nona$n_country))

# Sanity check against the stored 03 fit for log(total commitments) (never a
# literal: the panel revision of 2026-09-15 moved every published value).
ref_total <- readRDS(here("output", "fits", "headline_total_dr_bs.rds"))
if (is.finite(fit_env$att) && abs(fit_env$att - ref_total$att) > 0.01) {
  warning(sprintf(
    "Reproduced ATT on lcommitments_all (%.4f) deviates from the stored 03 fit (%.4f) by more than 0.01.",
    fit_env$att, ref_total$att))
}

t_env  <- if (is.finite(fit_env$se)  && fit_env$se  > 0) fit_env$att  / fit_env$se  else NA_real_
t_nona <- if (is.finite(fit_nona$se) && fit_nona$se > 0) fit_nona$att / fit_nona$se else NA_real_

# One-sided test of the model's prediction (H0: ATT <= 0, i.e. model-consistent;
# H1: ATT > 0, i.e. envelope response violates the model's sign prediction).
p_onesided_env <- if (is.na(t_env)) NA_real_ else 1 - pnorm(t_env)

message(sprintf("t (log total commitments) = %.3f | one-sided p (H1: ATT>0) = %.4f",
                t_env, p_onesided_env))
message(sprintf("t (log non-adaptation commitments) = %.3f | two-sided p = %.4f",
                t_nona, if (is.na(t_nona)) NA_real_ else 2 * pnorm(-abs(t_nona))))

verdict_env <- if (is.na(p_onesided_env)) {
  "inconclusive (fit failed)"
} else if (fit_env$att <= 0) {
  "consistent with the model's weakly negative prediction"
} else if (p_onesided_env >= 0.10) {
  "positive but not significant at conventional levels; not a clear rejection of the model's prediction"
} else if (p_onesided_env >= 0.05) {
  "positive and marginally significant (10% but not 5% level, one-sided); a soft rejection of the model's weakly negative prediction"
} else {
  "positive and significant at the 5% level (one-sided); rejects the model's weakly negative prediction"
}
message("Verdict (Test 2, envelope sign): ", verdict_env)

# ==============================================================================
# SECTION 8. TEST 2 — table tab:alpha_half
# ==============================================================================

pt_env  <- fit_env$pretrend  %||% list(stat = NA_real_, pval = NA_real_, df = 0L)
pt_nona <- fit_nona$pretrend %||% list(stat = NA_real_, pval = NA_real_, df = 0L)

p_twosided_nona <- if (is.na(t_nona)) NA_real_ else 2 * pnorm(-abs(t_nona))

alpha_rows <- c(
  paste0("Mean $s_i$ (", min(PRE_YEARS), "--", max(PRE_YEARS), " mean, $N=", dist_pre$n, "$) & ",
        sprintf("%.4f", dist_pre$mean), " \\\\"),
  paste0("Median $s_i$ & ", sprintf("%.4f", dist_pre$median), " \\\\"),
  paste0("$p_{90}$ $s_i$ & ", sprintf("%.4f", dist_pre$p90), " \\\\"),
  paste0("$p_{99}$ $s_i$ & ", sprintf("%.4f", dist_pre$p99), " \\\\"),
  paste0("Max $s_i$ & ", sprintf("%.4f", dist_pre$max), " \\\\"),
  paste0("$N$ countries with $s_i \\geq 0.5$ & ", dist_pre$n_ge_half, " \\\\"),
  "\\midrule",
  paste0("Mean $s_i$ at $g-1$ (adopters, $N=", dist_g1$n, "$) & ",
        sprintf("%.4f", dist_g1$mean), " \\\\"),
  paste0("Median $s_i$ at $g-1$ & ", sprintf("%.4f", dist_g1$median), " \\\\"),
  paste0("Max $s_i$ at $g-1$ & ", sprintf("%.4f", dist_g1$max), " \\\\"),
  paste0("$N$ adopters with $s_i(g-1) \\geq 0.5$ & ", dist_g1$n_ge_half, " \\\\"),
  "\\midrule",
  "\\multicolumn{2}{l}{\\textit{ATT: log(Total commitments) --- main spec}} \\\\",
  paste0("ATT & ", sprintf("%.4f", fit_env$att), " \\\\"),
  paste0("SE & ", sprintf("(%.4f)", fit_env$se), " \\\\"),
  paste0("$t$-statistic & ", fmt3(t_env), " \\\\"),
  paste0("One-sided $p$ ($H_1$: ATT $>0$) & ", fmt_p(p_onesided_env), " \\\\"),
  paste0("Pre-trend $\\chi^2$ & ", fmt3(pt_env$stat), " \\\\"),
  paste0("Pre-trend $p$ & ", fmt_p(pt_env$pval), " \\\\"),
  "\\midrule",
  "\\multicolumn{2}{l}{\\textit{ATT: log(Non-adaptation commitments) --- main spec (crowd-out check)}} \\\\",
  paste0("ATT & ", sprintf("%.4f", fit_nona$att), " \\\\"),
  paste0("SE & ", sprintf("(%.4f)", fit_nona$se), " \\\\"),
  paste0("$t$-statistic & ", fmt3(t_nona), " \\\\"),
  paste0("Two-sided $p$ & ", fmt_p(p_twosided_nona), " \\\\"),
  paste0("Pre-trend $\\chi^2$ & ", fmt3(pt_nona$stat), " \\\\"),
  paste0("Pre-trend $p$ & ", fmt_p(pt_nona$pval), " \\\\"),
  "\\midrule",
  paste0("Observations / Countries & ", format(fit_env$n_obs, big.mark = ","),
        " / ", fit_env$n_country, " \\\\")
)

alpha_tabular <- c(
  "\\begin{tabular}{lc}",
  "\\toprule",
  "Quantity & Value \\\\",
  "\\midrule",
  alpha_rows,
  "\\bottomrule",
  "\\end{tabular}"
)

notes_alpha <- paste0(
  "$s_i = $ \\texttt{commitments}/\\texttt{commitments\\_all}, proxy for $\\alpha_i$ via ",
  "$s_i^*=q_i/(1+q_i)$ (Lemma~\\ref{lem:monotone_s}); $\\alpha_i \\geq 0.5 \\iff s_i^* \\geq 0.5$. ",
  "Pre-treatment $s_i$ is far below 0.5 for every recipient (max ",
  sprintf("%.3f", max(dist_pre$max, dist_g1$max)),
  "), so under gross substitutes the model predicts ATT $\\leq 0$ (sign chain: see text, ",
  "Lemma~\\ref{lem:envelope_elast}). One-sided test: $H_0$: ATT $\\leq 0$ vs.\\ $H_1$: ATT $>0$. ",
  "CS\\,(2021) DR, never-treated controls, headline specification; SE: multiplier-bootstrap ",
  "(999 reps, seed 1242; \\texttt{did} 2.5.0; CRS Apr.\\ 2026) for the ATT, analytical fit for ",
  "the pre-trend Wald. ",
  PRETREND_NOTE(-5L, -2L, pt_env$df)
)
source_alpha <- "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"

write_tex_float(
  out_path      = file.path(dir_tabs, "alpha_half.tex"),
  caption_title = "Recipient adaptation share and the sign of the envelope response (test of the $\\mathrm{sign}(\\alpha_i-1/2)$ condition)",
  label         = "tab:alpha_half",
  tabular_lines = alpha_tabular,
  notes_text    = notes_alpha,
  source_text   = source_alpha
)

# ==============================================================================
# SECTION 9. Copy exhibits into paper/Tables/model_tests and
#            paper/Figures/model_tests (run_all.R is not modified — this
#            standalone script performs its own copy step).
# ==============================================================================

message("\n=== Section 9: copying exhibits into paper/ ===\n")

copy_dir <- function(src, dst) {
  files <- list.files(src, full.names = FALSE)
  for (f in files) {
    file.copy(file.path(src, f), file.path(dst, f), overwrite = TRUE)
  }
  length(files)
}

if (has_paper) {
  n_tab_copied <- copy_dir(dir_tabs, paper_tabs)
  n_fig_copied <- copy_dir(dir_figs, paper_figs)
  message(sprintf("Copied %d table(s) to %s", n_tab_copied, paper_tabs))
  message(sprintf("Copied %d figure(s) to %s", n_fig_copied, paper_figs))
} else message("paper/ not found -- exhibits are left in output/ only")

message("\n=== 10_model_tests.R complete ===\n")
