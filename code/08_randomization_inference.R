# ==============================================================================
# 08_randomization_inference.R
# Fisher randomization inference (sharp null of no effect for any unit) for the
# actual main CS(2021) specification, recomputed on the current main-spec panel
# (superseding the June-2026 old-revision-battery RI exercise, which used
# est_method = "reg" and an older spec).
#
# Main specification (must match 03_main_results.R Section 6 exactly):
#   att_gt(yname = <outcome>, tname = "year", idname = "country_id",
#          gname = "cohort_year", xformla = ~ ge_est + log_population,
#          est_method = "dr", control_group = "nevertreated", anticipation = 0,
#          base_period = "universal", panel = TRUE, allow_unbalanced_panel = TRUE)
#   Pooled ATT: aggte(type = "simple", na.rm = TRUE)$overall.att
#
# Design: Fisher sharp null of no effect for any unit, implemented by permuting
#   the treatment assignment three ways, each with N_DRAWS = 2000 replications:
#     A. Timing permutation  — reshuffle cohort_year among the 40 treated
#        countries (cohort sizes preserved by construction).
#     B. Assignment permutation — draw 40 of the 126 countries at random and
#        assign them the observed (permuted) cohort-year vector.
#     C. Stratified assignment permutation — as B, but the 40 draws are made
#        within World Bank region strata, matching the observed number of
#        treated countries per region.
#   For every draw the thin-cohort rule (< 5 treated -> dropped) is re-applied
#   identically to the main specification. If a draw yields < 2 remaining
#   treated cohorts it is recorded as skipped (not a valid draw).
#
# Estimation per draw: est_method = "dr", bstrap = FALSE, cband = FALSE (point
#   estimate only — aggte() needs the influence function for SEs, so bootstrap
#   SEs are not computed for permutation draws, only for the one-time observed
#   ATT reproduction). If "dr" errors OR the forked worker crashes outright
#   (documented gotcha: att_gt(est_method = "dr") can segfault on reduced /
#   unbalanced panels via fastglm::colMax_dense — a crash that a plain
#   tryCatch() cannot catch because it kills the process), the draw is retried
#   with est_method = "reg" in a second pass. Using parallel::mclapply()
#   (fork-per-task, mc.preschedule = FALSE) means a segfaulting fork only
#   loses that one task — the parent session and all other cores are safe —
#   which is why mclapply is used here instead of future_lapply().
#
# Inputs : data/processed/simple_panel_wgi.csv
# Outputs:
#   output/tables/randomization/ri_pvalues.tex   (tab:ri_pvalues)
#   output/tables/randomization/ri_draws.csv     (raw permutation draws; also
#                                                the cache read on later runs)
#   output/figures/randomization/fig_ri_distributions.png
#   (copied to paper/Tables/randomization/ and paper/Figures/randomization/)
# ==============================================================================

# ============================================================
# Paper-to-Code Naming Map
# ============================================================
# Paper Notation           | Code Name             | Description
# $Y_{it}$                 | log_commits, share_adapt | Two RI outcomes
# $G_i$                    | cohort_year            | NAP adoption cohort (0=never)
# $ATT(g,t)$                | gt_obj (inside fit_att_simple) | Group-time ATT
# $\hat\theta^{simp}_{obs}$ | att_obs               | Observed pooled ATT (bstrap SE)
# $\hat\theta^{simp}_{(b)}$ | att_perm[b]           | Permutation-draw pooled ATT
# $p^{RI}_{two}$            | p_two                 | 2-sided RI p-value
# $p^{RI}_{one}$            | p_one                 | 1-sided (>=) RI p-value
# $X_{it}$                  | ge_est, log_population| Controls: WGI GE + log pop.
# ============================================================

# ARM-mac headless gotcha (verified in 04): must be set before any library()
# call that might pull in rgl transitively.
Sys.setenv(RGL_USE_NULL = TRUE)

library(data.table)
library(dplyr)
library(ggplot2)
library(xtable)
library(here)
library(did)
library(parallel)

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

# ------------------------------------------------------------------------
# REPRODUCIBILITY INVARIANT (RNG + forking) -- documented, not changed
# ------------------------------------------------------------------------
# All permutation draws (Design A/B/C, Section 4) are produced by calling
# set.seed(1242) once per design, immediately before a *sequential* for-loop
# of sample() calls (generate_perm_timing/_assignment/_stratified). The
# resulting draws (the perm_list objects, and therefore everything saved in
# ri_draws.csv) are fully determined by the default RNG kind BEFORE any
# forking happens.
#
# Estimation (Section 5, run_design()) then forks these ALREADY-MATERIALIZED
# draws out to worker processes via parallel::mclapply(). Each worker calls
# fit_att_simple() -> att_gt(bstrap = FALSE) on the data.table it is handed;
# att_gt() draws no random numbers when bstrap = FALSE, so there is no RNG
# stream running inside the forked workers to desynchronize. This is what
# makes the permutation draws reproducible across machines/core counts: the
# RNG only ever advances sequentially in the parent process, before
# mclapply() is invoked, and the workers are read-only with respect to
# random state.
#
# We deliberately do NOT call RNGkind("L'Ecuyer-CMRG") here. Doing so would
# change the default generator's actual sample() stream (L'Ecuyer-CMRG is a
# different algorithm from the Mersenne-Twister default), which would
# silently change every draw produced by the sequential set.seed(1242) +
# sample() calls above -- a fresh run would then no longer reproduce the
# already-saved output/tables/randomization/ri_draws.csv. RNGkind("L'Ecuyer-
# CMRG") + mc.reset.stream() is the correct fix ONLY when random draws are
# made *inside* forked workers, which is not this script's design (workers
# are deterministic given the panel they are handed). Left undone
# intentionally so a future editor does not "fix" this by adding it.
# ------------------------------------------------------------------------

set.seed(20240601)  # global seed — local set.seed(1242) calls follow, per README

SCRIPT_T0 <- Sys.time()

# -----------------------------------------------------------------------
# Output directories
# -----------------------------------------------------------------------
dir.create(here("output", "tables",  "randomization"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "figures", "randomization"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------
# Helpers duplicated from 03_main_results.R §1 (project convention: each
# stage script is self-contained — see 04/05 headers "duplicated from 03/04").
# -----------------------------------------------------------------------
esc_header <- function(x) {
  x <- gsub("%", "\\\\%", x)
  x <- gsub("_", "\\\\_", x)
  x <- gsub("#", "\\\\#", x)
  x
}

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
# SECTION 1. Load and prepare the DiD panel — replicates 03_main_results.R §1
# exactly, so the estimation sample is identical to the main specification
# (126 countries, 40 treated, N = 2,014 after the thin-cohort rule).
# ==============================================================================

message("\n=== 08_randomization_inference.R: loading panel ===\n")

aggregated <- fread(here("data", "processed", "simple_panel_wgi.csv"))
aggregated <- as.data.frame(aggregated)

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
    nap_year_c = if_else(is.infinite(nap_year_c), NA_real_, nap_year_c),
    cohort_year = case_when(
      is.na(nap_year_c)                              ~ 0,
      nap_year_c < first_year                         ~ 0,
      nap_year_c > max(did_panel$year, na.rm = TRUE)  ~ 0,
      TRUE                                             ~ as.numeric(nap_year_c)
    )
  )

did_panel <- did_panel %>%
  select(-any_of("cohort_year")) %>%
  left_join(country_gname %>% select(recipient_name, cohort_year),
            by = "recipient_name") %>%
  mutate(log_population = log(population))

if (!"share_adapt" %in% names(did_panel)) {
  stop("share_adapt column missing from simple_panel_wgi.csv — check 01_prepare_data.R output.")
}

thin_threshold <- 5L

cohort_sizes <- did_panel %>%
  filter(cohort_year > 0) %>%
  distinct(recipient_name, cohort_year) %>%
  count(cohort_year, name = "n_treated")

thin_cohorts <- cohort_sizes$cohort_year[cohort_sizes$n_treated < thin_threshold]

did_panel_main <- did_panel %>%
  filter(!(cohort_year %in% thin_cohorts))

# Keep only the columns needed downstream (lighter copies inside the
# permutation loop; country_id/year identify rows, ge_est/log_population are
# the controls, log_commits/share_adapt are the two RI outcomes, WB_region
# is used for the stratified-permutation design).
base_panel <- as.data.table(
  did_panel_main %>%
    select(country_id, recipient_name, year, cohort_year, ge_est,
           log_population, log_commits, share_adapt, WB_region)
)

n_countries <- uniqueN(base_panel$country_id)
n_treated   <- uniqueN(base_panel[cohort_year > 0, country_id])
n_obs       <- nrow(base_panel)

message(sprintf("Sample check — countries: %d | treated: %d | N: %d",
                n_countries, n_treated, n_obs))

# Reference sample and ATTs come from the stored 03 fits (output/fits/), not
# from literals, so a panel revision cannot silently invalidate this stage.
ref_head  <- readRDS(here("output", "fits", "headline_adaptation_dr_bs.rds"))
ref_share <- readRDS(here("output", "fits", "headline_share_dr_bs.rds"))
stopifnot(
  "Country count differs from the stored 03 headline fit" = n_countries == ref_head$n_country,
  "N differs from the stored 03 headline fit"             = n_obs == ref_head$n_obs
)
message(sprintf("Treated countries in the main-spec sample: %d", n_treated))

# ------------------------------------------------------------------------
# Cache file I/O helpers (the permutation draws are stored as CSV)
# ------------------------------------------------------------------------
# Doubles are stored twice: a readable decimal column and an exact hexadecimal
# column (sprintf("%a")) that reads back bit-for-bit. Decimal text does not
# always round-trip exactly through R's parser, and the permutation p-values
# compare each draw with the observed ATT exactly, so the reader uses the
# hexadecimal columns.
RI_NUM_COLS <- c("att_perm", "att_obs", "se_obs", "target_att", "generated_at")

#' Write the randomization-inference draws to one long CSV
#'
#' One row per design x outcome x draw; scalar metadata are repeated on every
#' row so the file is self-describing.
#' @param ri list with the structure built in Section 7
#' @param path output .csv path
write_ri_draws_csv <- function(ri, path) {
  rows <- lapply(names(ri$designs), function(dk) {
    d <- ri$designs[[dk]]
    do.call(rbind, lapply(names(d$outcomes), function(oc) {
      r <- d$outcomes[[oc]]
      data.frame(design_key = dk, design = d$design, n_skipped_thin = d$n_skipped_thin,
                 outcome = oc, draw = seq_along(r$att_perm),
                 att_perm = r$att_perm, fallback_used = r$fallback_used,
                 att_obs = unname(ri$att_obs[oc]), se_obs = unname(ri$se_obs[oc]),
                 target_att = unname(ri$target_att[oc]),
                 n_draws = ri$n_draws, thin_threshold = ri$thin_threshold,
                 seed_global = ri$seed_global, seed_permutation = ri$seed_permutation,
                 n_countries = ri$n_countries, n_treated = ri$n_treated, n_obs = ri$n_obs,
                 generated_at = as.numeric(ri$generated_at),
                 stringsAsFactors = FALSE)
    }))
  })
  out <- do.call(rbind, rows)
  for (v in RI_NUM_COLS) {
    hex <- ifelse(is.na(out[[v]]), NA_character_, sprintf("%a", out[[v]]))
    out[[paste0(v, "_hex")]] <- hex
    out[[v]] <- ifelse(is.na(out[[v]]), NA_character_, sprintf("%.17g", out[[v]]))
  }
  data.table::fwrite(out, path, na = "NA")
  invisible(path)
}

#' Read the long randomization-inference CSV back into the Section 7 list
#'
#' @param path .csv written by write_ri_draws_csv()
#' @return list identical to the object write_ri_draws_csv() was given
read_ri_draws_csv <- function(path) {
  x <- as.data.frame(data.table::fread(path, na.strings = "NA", encoding = "UTF-8",
                                       colClasses = "character"))
  for (v in RI_NUM_COLS) x[[v]] <- as.numeric(x[[paste0(v, "_hex")]])
  for (v in c("n_skipped_thin", "draw", "n_draws", "thin_threshold", "seed_global",
              "seed_permutation", "n_countries", "n_treated", "n_obs")) {
    x[[v]] <- as.integer(x[[v]])
  }
  x$fallback_used <- as.logical(x$fallback_used)
  outcome_vars <- unique(x$outcome)
  first_by_oc <- function(col) {
    setNames(x[[col]][match(outcome_vars, x$outcome)], outcome_vars)
  }
  designs <- lapply(setNames(nm = unique(x$design_key)), function(dk) {
    xd <- x[x$design_key == dk, ]
    outs <- lapply(setNames(nm = outcome_vars), function(oc) {
      xo <- xd[xd$outcome == oc, ]
      xo <- xo[order(xo$draw), ]
      list(att_perm = xo$att_perm, fallback_used = xo$fallback_used,
           n_fallback = sum(xo$fallback_used), n_valid = sum(!is.na(xo$att_perm)))
    })
    list(design = xd$design[1L], n_draws = sum(xd$outcome == outcome_vars[1L]),
         n_skipped_thin = xd$n_skipped_thin[1L], outcomes = outs)
  })
  list(seed_global = x$seed_global[1L], seed_permutation = x$seed_permutation[1L],
       n_draws = x$n_draws[1L], thin_threshold = x$thin_threshold[1L],
       outcome_vars = outcome_vars, att_obs = first_by_oc("att_obs"),
       se_obs = first_by_oc("se_obs"), target_att = first_by_oc("target_att"),
       n_countries = x$n_countries[1L], n_treated = x$n_treated[1L], n_obs = x$n_obs[1L],
       designs = designs, generated_at = .POSIXct(x$generated_at[1L]))
}

# ==============================================================================
# CACHE-AWARE FAST PATH -- table-only regeneration
# The permutation-fitting step below (Sections 2, 4-7) is the ~28-minute cost
# of this script (2,000 draws x 3 designs x 2 outcomes x 2 att_gt() passes).
# If output/tables/randomization/ri_draws.csv already exists, skip straight
# to Section 8 (p-values) / Section 9 (table) / Section 10 (figure) using the
# saved draws -- do NOT re-run att_gt() for the observed ATTs or for any
# permutation draw. Delete ri_draws.csv to force a full recompute (~29 min).
# ==============================================================================

ri_draws_cache_path <- here("output", "tables", "randomization", "ri_draws.csv")
USE_CACHED_DRAWS     <- file.exists(ri_draws_cache_path)

if (USE_CACHED_DRAWS) {

  message("\n=== Cache hit: ", ri_draws_cache_path, " exists -- loading saved permutation ",
          "draws and skipping the att_gt() reproduction + permutation-fitting steps ===\n")

  ri_draws       <- read_ri_draws_csv(ri_draws_cache_path)
  att_obs        <- ri_draws$att_obs
  se_obs         <- ri_draws$se_obs
  N_DRAWS        <- ri_draws$n_draws
  thin_threshold <- ri_draws$thin_threshold
  outcome_vars   <- ri_draws$outcome_vars
  all_designs    <- ri_draws$designs
  N_CORES        <- max(1L, min(8L, parallel::detectCores() - 1L))  # unused on this path; kept for the Section 12 wall-time message

  message(sprintf("  Loaded: %d draws/design, seed_permutation=%d, generated_at=%s",
                  N_DRAWS, ri_draws$seed_permutation, format(ri_draws$generated_at)))

} else {

# ==============================================================================
# SECTION 2. Reproduce the observed ATTs (main specification, bootstrap SE)
# Seed rule (reproduces the published SEs): set.seed(1242) immediately before the
# bootstrap call.
# ==============================================================================

message("\n=== Reproducing observed ATTs (main specification) ===\n")

outcome_vars <- c("log_commits", "share_adapt")
TARGET_ATT   <- c(log_commits = ref_head$att, share_adapt = ref_share$att)
ATT_TOL      <- 0.001

att_obs <- setNames(vector("numeric", length(outcome_vars)), outcome_vars)
se_obs  <- setNames(vector("numeric", length(outcome_vars)), outcome_vars)

for (oc in outcome_vars) {
  set.seed(1242)
  gt_obs <- att_gt(
    yname         = oc,
    tname         = "year",
    idname        = "country_id",
    gname         = "cohort_year",
    xformla       = ~ ge_est + log_population,
    data          = as.data.frame(base_panel),
    est_method    = "dr",
    bstrap        = TRUE,
    biters        = 999L,
    cband         = FALSE,
    control_group = "nevertreated",
    anticipation  = 0,
    base_period   = "universal",
    panel         = TRUE,
    allow_unbalanced_panel = TRUE
  )
  agg_obs <- aggte(gt_obs, type = "simple", na.rm = TRUE)
  att_obs[oc] <- agg_obs$overall.att
  se_obs[oc]  <- agg_obs$overall.se
  message(sprintf("  %-12s ATT = %.4f (SE = %.4f)  [target %.4f]",
                  oc, att_obs[oc], se_obs[oc], TARGET_ATT[oc]))
}

stopifnot(
  "log_commits observed ATT does not match the stored 03 fit" =
    abs(att_obs["log_commits"] - TARGET_ATT["log_commits"]) < ATT_TOL,
  "share_adapt observed ATT does not match the stored 03 fit" =
    abs(att_obs["share_adapt"] - TARGET_ATT["share_adapt"]) < ATT_TOL
)
message("Observed ATTs match target values within tolerance ", ATT_TOL, ".")

# ==============================================================================
# SECTION 3. Permutation-fitting primitives
# ==============================================================================

#' Rebuild the panel under a permuted treatment assignment and re-apply the
#' thin-cohort rule identically to the main specification.
#'
#' @param perm_dt data.table with columns country_id, cohort_year — the new
#'   treated-country assignment for this draw (untreated countries are
#'   implicit: any country_id in base_panel not present in perm_dt is coded
#'   cohort_year = 0, i.e. never-treated).
#' @param base_panel_in data.table, the fixed main-spec panel (cohort_year
#'   column will be replaced).
#' @param thin_threshold_in integer, minimum treated units per cohort to keep.
#' @return named list: panel_dt (data.table or NULL if invalid), valid
#'   (logical), n_cohorts (integer, post-filter treated-cohort count).
rebuild_panel <- function(perm_dt, base_panel_in, thin_threshold_in) {
  stopifnot(is.data.table(perm_dt), all(c("country_id", "cohort_year") %in% names(perm_dt)))

  dt <- copy(base_panel_in)
  dt[, cohort_year := NULL]
  dt <- merge(dt, perm_dt, by = "country_id", all.x = TRUE, sort = FALSE)
  dt[is.na(cohort_year), cohort_year := 0]

  cs   <- dt[cohort_year > 0, .(n_treated = uniqueN(country_id)), by = cohort_year]
  thin <- cs[n_treated < thin_threshold_in, cohort_year]
  if (length(thin) > 0L) dt <- dt[!(cohort_year %in% thin)]

  n_cohorts_post <- uniqueN(dt[cohort_year > 0, cohort_year])
  is_valid       <- n_cohorts_post >= 2L

  list(panel_dt = if (is_valid) dt else NULL,
       valid    = is_valid,
       n_cohorts = n_cohorts_post)
}

#' Fit att_gt() + aggte(type = "simple") for one outcome on one permuted
#' panel; point estimate only (bstrap = FALSE, cband = FALSE, per spec — the
#' permutation exercise needs only the point estimate to build the RI
#' reference distribution, not its standard error).
#'
#' @param panel_dt data.table, permuted+filtered panel from rebuild_panel().
#' @param outcome_var character, "log_commits" or "share_adapt".
#' @param est_method_in character, "dr" or "reg".
#' @return numeric pooled ATT, or NA_real_ on any estimation failure.
fit_att_simple <- function(panel_dt, outcome_var, est_method_in) {
  gt <- tryCatch(
    suppressWarnings(suppressMessages(
      att_gt(
        yname         = outcome_var,
        tname         = "year",
        idname        = "country_id",
        gname         = "cohort_year",
        xformla       = ~ ge_est + log_population,
        data          = as.data.frame(panel_dt),
        est_method    = est_method_in,
        bstrap        = FALSE,
        cband         = FALSE,
        control_group = "nevertreated",
        anticipation  = 0,
        base_period   = "universal",
        panel         = TRUE,
        allow_unbalanced_panel = TRUE
      )
    )),
    error = function(e) NULL
  )
  if (is.null(gt)) return(NA_real_)

  agg <- tryCatch(
    suppressWarnings(suppressMessages(aggte(gt, type = "simple", na.rm = TRUE))),
    error = function(e) NULL
  )
  if (is.null(agg) || is.null(agg$overall.att)) return(NA_real_)
  as.numeric(agg$overall.att)
}

# ==============================================================================
# SECTION 4. Permutation-draw generators
# Seed rule (reproduces the published SEs): set.seed(1242) immediately before each design's draw loop,
# so each design is independently reproducible.
# ==============================================================================

N_DRAWS <- 2000L

treated_ids        <- base_panel[cohort_year > 0, unique(country_id)]
observed_cohort_dt  <- unique(base_panel[country_id %in% treated_ids,
                                          .(country_id, cohort_year)])
observed_cohort_yrs <- observed_cohort_dt$cohort_year
all_ids             <- unique(base_panel$country_id)

stopifnot(length(treated_ids) == 40L, nrow(observed_cohort_dt) == 40L,
          length(all_ids) == ref_head$n_country)

#' Design A — timing permutation within the 40 treated countries. Cohort
#' sizes are preserved by construction (same multiset of adoption years
#' reassigned to a shuffled set of the same 40 countries).
generate_perm_timing <- function(observed_cohort_dt_in, n_draws) {
  set.seed(1242)
  perm_list <- vector("list", n_draws)
  for (b in seq_len(n_draws)) {
    shuffled_years <- sample(observed_cohort_dt_in$cohort_year)
    perm_list[[b]] <- data.table(
      country_id  = observed_cohort_dt_in$country_id,
      cohort_year = shuffled_years
    )
  }
  perm_list
}

#' Design B — assignment permutation. Draw n_treated countries at random
#' from the full pool of all_ids and assign them the observed (permuted)
#' cohort-year vector; all other countries are implicitly never-treated.
generate_perm_assignment <- function(all_ids_in, observed_cohort_yrs_in,
                                      n_treated_in, n_draws) {
  set.seed(1242)
  perm_list <- vector("list", n_draws)
  for (b in seq_len(n_draws)) {
    drawn_ids      <- sample(all_ids_in, n_treated_in)
    shuffled_years <- sample(observed_cohort_yrs_in)
    perm_list[[b]] <- data.table(country_id = drawn_ids, cohort_year = shuffled_years)
  }
  perm_list
}

#' Design C — stratified assignment permutation. As Design B, but the drawn
#' treated countries are sampled within World Bank region strata, matching
#' the observed number of treated countries per region (missing WB_region
#' values are pooled into an "Unknown" stratum).
generate_perm_stratified <- function(base_panel_in, treated_ids_in,
                                      observed_cohort_yrs_in, n_draws) {
  region_map <- unique(base_panel_in[, .(country_id, WB_region)])
  region_map[, region_grp := fifelse(is.na(WB_region), "Unknown", WB_region)]

  treated_region_counts <- region_map[country_id %in% treated_ids_in,
                                       .(n_needed = .N), by = region_grp]
  stopifnot(sum(treated_region_counts$n_needed) == length(treated_ids_in))

  # Precompute the eligible-country pool per stratum once (outside the loop).
  region_pools <- split(region_map$country_id, region_map$region_grp)
  for (r in treated_region_counts$region_grp) {
    stopifnot(length(region_pools[[r]]) >=
                treated_region_counts[region_grp == r, n_needed])
  }

  set.seed(1242)
  perm_list <- vector("list", n_draws)
  for (b in seq_len(n_draws)) {
    drawn_ids <- unlist(lapply(seq_len(nrow(treated_region_counts)), function(r) {
      reg      <- treated_region_counts$region_grp[r]
      n_needed <- treated_region_counts$n_needed[r]
      sample(region_pools[[reg]], n_needed)
    }), use.names = FALSE)
    shuffled_years <- sample(observed_cohort_yrs_in)
    perm_list[[b]] <- data.table(country_id = drawn_ids, cohort_year = shuffled_years)
  }
  perm_list
}

# ==============================================================================
# SECTION 5. Two-pass, crash-safe permutation fitting
# Pass 1: est_method = "dr" for every valid draw, forked one-task-per-core
#   (mc.preschedule = FALSE) so a segfault in one fork only loses that draw.
# Pass 2: est_method = "reg" fallback, only for draws where pass 1 failed
#   (ordinary error OR crashed fork — both surface as NULL/NA from mclapply).
# ==============================================================================

N_CORES <- max(1L, min(8L, parallel::detectCores() - 1L))
message(sprintf("\nUsing %d cores for permutation fitting (mclapply, fork-per-task).\n",
                N_CORES))

#' Run one full RI design: generate + rebuild once, fit both outcomes.
#'
#' @return named list: design, n_draws, n_skipped_thin, per-outcome results
#'   (att_perm numeric vector length n_draws, fallback_used logical vector,
#'   n_fallback integer, n_valid integer).
run_design <- function(design_name, perm_list, base_panel_in, thin_threshold_in,
                        outcome_vars_in, n_cores) {

  t0 <- Sys.time()
  message(sprintf("=== Design: %s (%d draws) ===", design_name, length(perm_list)))

  rebuilt_list <- lapply(perm_list, rebuild_panel,
                          base_panel_in = base_panel_in,
                          thin_threshold_in = thin_threshold_in)
  valid_flags  <- vapply(rebuilt_list, `[[`, logical(1L), "valid")
  valid_idx    <- which(valid_flags)
  n_skipped    <- sum(!valid_flags)
  message(sprintf("  %d/%d draws valid (>= 2 treated cohorts post thin-filter); %d skipped.",
                  length(valid_idx), length(perm_list), n_skipped))

  outcome_results <- vector("list", length(outcome_vars_in))
  names(outcome_results) <- outcome_vars_in

  for (oc in outcome_vars_in) {
    att_perm      <- rep(NA_real_, length(perm_list))
    fallback_used <- rep(FALSE, length(perm_list))

    # --- Pass 1: est_method = "dr" ---
    pass1 <- parallel::mclapply(
      valid_idx,
      function(i) fit_att_simple(rebuilt_list[[i]]$panel_dt, oc, "dr"),
      mc.cores = n_cores, mc.preschedule = FALSE
    )
    pass1_failed_local <- vapply(pass1, function(x) is.null(x) || is.na(x), logical(1L))
    for (k in seq_along(valid_idx)) {
      if (!pass1_failed_local[k]) att_perm[valid_idx[k]] <- pass1[[k]]
    }
    fail_idx <- valid_idx[pass1_failed_local]
    fallback_used[fail_idx] <- TRUE

    # --- Pass 2: est_method = "reg" fallback ---
    if (length(fail_idx) > 0L) {
      pass2 <- parallel::mclapply(
        fail_idx,
        function(i) fit_att_simple(rebuilt_list[[i]]$panel_dt, oc, "reg"),
        mc.cores = n_cores, mc.preschedule = FALSE
      )
      for (k in seq_along(fail_idx)) {
        val <- pass2[[k]]
        if (!(is.null(val) || is.na(val))) att_perm[fail_idx[k]] <- val
      }
    }

    n_fallback <- length(fail_idx)
    n_valid    <- sum(!is.na(att_perm))
    message(sprintf("  [%s] valid ATTs: %d | reg fallbacks: %d | still-NA after fallback: %d",
                    oc, n_valid, n_fallback, length(fail_idx) - sum(!is.na(att_perm[fail_idx]))))

    outcome_results[[oc]] <- list(
      att_perm      = att_perm,
      fallback_used = fallback_used,
      n_fallback    = n_fallback,
      n_valid       = n_valid
    )
  }

  message(sprintf("  Design %s done in %.1f sec.\n",
                  design_name, as.numeric(Sys.time() - t0, units = "secs")))

  list(design = design_name, n_draws = length(perm_list),
       n_skipped_thin = n_skipped, outcomes = outcome_results)
}

# ==============================================================================
# SECTION 6. Generate the three designs and run the permutation fits
# ==============================================================================

message("\n=== Generating permutation draws ===\n")

perm_timing     <- generate_perm_timing(observed_cohort_dt, N_DRAWS)
perm_assignment <- generate_perm_assignment(all_ids, observed_cohort_yrs, n_treated, N_DRAWS)
perm_stratified <- generate_perm_stratified(base_panel, treated_ids, observed_cohort_yrs, N_DRAWS)

message("\n=== Running permutation fits (this is the slow step) ===\n")

res_timing     <- run_design("Timing permutation (within adopters)",
                              perm_timing, base_panel, thin_threshold, outcome_vars, N_CORES)
res_assignment <- run_design("Assignment permutation (across countries)",
                              perm_assignment, base_panel, thin_threshold, outcome_vars, N_CORES)
res_stratified <- run_design("Stratified assignment permutation (by WB region)",
                              perm_stratified, base_panel, thin_threshold, outcome_vars, N_CORES)

all_designs <- list(timing = res_timing, assignment = res_assignment,
                    stratified = res_stratified)

# ==============================================================================
# SECTION 7. Save raw draws — no re-run needed for downstream inspection.
# ==============================================================================

ri_draws <- list(
  seed_global      = 20240601L,
  seed_permutation = 1242L,
  n_draws          = N_DRAWS,
  thin_threshold   = thin_threshold,
  outcome_vars     = outcome_vars,
  att_obs          = att_obs,
  se_obs           = se_obs,
  target_att       = TARGET_ATT,
  n_countries      = n_countries,
  n_treated        = n_treated,
  n_obs            = n_obs,
  designs          = all_designs,
  generated_at     = Sys.time()
)
write_ri_draws_csv(ri_draws, ri_draws_cache_path)
message("Saved: ", ri_draws_cache_path)

}  # end of `else` -- full recompute branch (USE_CACHED_DRAWS == FALSE)

# ==============================================================================
# SECTION 8. Summary statistics — RI p-values, N valid, N fallback, SD, 95% range
# ==============================================================================

message("\n=== Computing RI p-values ===\n")

design_labels <- c(timing = "Timing (within adopters)",
                   assignment = "Assignment (across countries)",
                   stratified = "Assignment, stratified (WB region)")

summary_rows <- vector("list", length(all_designs) * length(outcome_vars))
row_i <- 1L
for (dname in names(all_designs)) {
  d <- all_designs[[dname]]
  for (oc in outcome_vars) {
    r          <- d$outcomes[[oc]]
    valid_atts <- r$att_perm[!is.na(r$att_perm)]
    obs        <- att_obs[oc]

    p_two <- if (length(valid_atts) > 0L) mean(abs(valid_atts) >= abs(obs)) else NA_real_
    p_one <- if (length(valid_atts) > 0L) mean(valid_atts >= obs)          else NA_real_
    sd_p  <- if (length(valid_atts) > 1L) sd(valid_atts)                   else NA_real_
    q_lo  <- if (length(valid_atts) > 0L) as.numeric(quantile(valid_atts, 0.025)) else NA_real_
    q_hi  <- if (length(valid_atts) > 0L) as.numeric(quantile(valid_atts, 0.975)) else NA_real_

    summary_rows[[row_i]] <- data.frame(
      outcome       = oc,
      design        = unname(design_labels[dname]),
      att_obs       = obs,
      p_two         = p_two,
      p_one         = p_one,
      n_valid       = r$n_valid,
      n_fallback    = r$n_fallback,
      n_skipped     = d$n_skipped_thin,
      sd_perm       = sd_p,
      q_lo          = q_lo,
      q_hi          = q_hi,
      stringsAsFactors = FALSE
    )
    message(sprintf("  %-12s | %-34s | ATT=%.4f | p(2s)=%.3f | p(1s,>=)=%.3f | N=%d | fallback=%d",
                    oc, design_labels[dname], obs, p_two, p_one, r$n_valid, r$n_fallback))
    row_i <- row_i + 1L
  }
}
summary_df <- do.call(rbind, summary_rows)

# ==============================================================================
# SECTION 9. Table tab:ri_pvalues — panel per outcome, row per design
# ==============================================================================

outcome_labels <- c(log_commits = "log(Adaptation commitments)",
                    share_adapt = "Adaptation share (\\% of global)")
# Plain-text variant for ggplot facet strips (LaTeX escapes above are only
# valid inside the .tex table, not inside a plotted string).
outcome_labels_fig <- c(log_commits = "log(Adaptation commitments)",
                        share_adapt = "Adaptation share (% of global)")

fmt_row <- function(row) {
  paste0(
    row$design, " & ",
    sprintf("%.4f", row$att_obs), " & ",
    sprintf("%.3f", row$p_two), " & ",
    sprintf("%.3f", row$p_one), " & ",
    row$n_valid, " & ",
    row$n_fallback, " & ",
    sprintf("%.4f", row$sd_perm), " & ",
    sprintf("[%.4f, %.4f]", row$q_lo, row$q_hi),
    " \\\\"
  )
}

tab_lines <- c(
  "\\begin{tabular}{lccccccc}",
  "\\toprule",
  "Design & ATT (obs.) & $p$ (two-sided) & $p$ (one-sided, $\\geq$) & $N$ valid & $N$ fallback & SD (perm.) & 95\\% range (perm.) \\\\",
  "\\midrule"
)
for (oc in outcome_vars) {
  tab_lines <- c(tab_lines,
    sprintf("\\multicolumn{8}{l}{\\textit{Panel %s: %s}} \\\\",
            if (oc == "log_commits") "A" else "B", outcome_labels[oc]),
    "\\midrule")
  oc_rows <- summary_df[summary_df$outcome == oc, ]
  for (i in seq_len(nrow(oc_rows))) {
    tab_lines <- c(tab_lines, fmt_row(oc_rows[i, ]))
  }
  if (oc != outcome_vars[length(outcome_vars)]) {
    tab_lines <- c(tab_lines, "\\\\[0.25em]")
  }
}
tab_lines <- c(tab_lines, "\\bottomrule", "\\end{tabular}")

n_skip_note <- unique(summary_df$n_skipped)
notes_txt <- paste0(
  "Fisher RI (sharp null), headline specification, ", N_DRAWS, " draws/design. ",
  "\\textit{Timing} permutes adoption year among ", n_treated, " treated countries; ",
  "\\textit{Assignment} redraws which ", n_treated, " of ", n_countries, " are treated; ",
  "\\textit{Stratified} draws within World Bank region strata. $N$ skipped: draws dropped by the ",
  "$<$", thin_threshold, " thin-cohort rule (", paste(n_skip_note, collapse = "/"),
  "; never binds under permutation). $N$ fallback: \"dr\" failed, \"reg\" used. $p$-value: two-sided $=$ share of $|ATT_{perm}| \\geq |ATT_{obs}|$; ",
  "one-sided $=$ share of $ATT_{perm} \\geq ATT_{obs}$, over valid draws, observed excluded. ",
  "SD/95\\% range describe the permutation distribution"
)
source_txt <- "OECD CRS (Rio adaptation markers); UNFCCC NAP Central; World Bank WGI/WDI"

write_tex_float(
  out_path      = here("output", "tables", "randomization", "ri_pvalues.tex"),
  caption_title = "Randomization inference: permutation $p$-values for the main specification",
  label         = "tab:ri_pvalues",
  tabular_lines = tab_lines,
  notes_text    = notes_txt,
  source_text   = source_txt
)

# ==============================================================================
# SECTION 10. Figure fig_ri_distributions.png — permutation histograms
# ==============================================================================

message("\n=== Building RI distribution figure ===\n")

hist_rows <- vector("list", length(all_designs) * length(outcome_vars))
row_i <- 1L
for (dname in names(all_designs)) {
  d <- all_designs[[dname]]
  for (oc in outcome_vars) {
    valid_atts <- d$outcomes[[oc]]$att_perm
    valid_atts <- valid_atts[!is.na(valid_atts)]
    if (length(valid_atts) == 0L) next
    hist_rows[[row_i]] <- data.frame(
      att     = valid_atts,
      design  = unname(design_labels[dname]),
      outcome = unname(outcome_labels_fig[oc]),
      stringsAsFactors = FALSE
    )
    row_i <- row_i + 1L
  }
}
hist_df <- do.call(rbind, hist_rows[!vapply(hist_rows, is.null, logical(1L))])
hist_df$design  <- factor(hist_df$design, levels = unname(design_labels))
hist_df$outcome <- factor(hist_df$outcome, levels = unname(outcome_labels_fig))

obs_df <- data.frame(
  outcome = factor(unname(outcome_labels_fig[outcome_vars]), levels = unname(outcome_labels_fig)),
  att_obs = as.numeric(att_obs[outcome_vars])
)

p_ri <- ggplot(hist_df, aes(x = att)) +
  geom_histogram(bins = 40, fill = "#2E86C1", colour = "white", alpha = 0.85) +
  geom_vline(data = obs_df, aes(xintercept = att_obs),
             colour = "#C0392B", linewidth = 0.9, linetype = "dashed") +
  facet_grid(design ~ outcome, scales = "free") +
  # No title, subtitle, or caption inside the plot
  labs(title = NULL, subtitle = NULL, caption = NULL,
       x = "Permutation-draw pooled ATT", y = "Count") +
  theme_minimal() +
  theme(text = element_text(family = "serif", size = 11),
        strip.text = element_text(face = "bold", size = 9),
        panel.grid.minor = element_blank())

ggsave(here("output", "figures", "randomization", "fig_ri_distributions.png"),
       p_ri, width = 11, height = 9, dpi = 300)
message("Saved: ", here("output", "figures", "randomization", "fig_ri_distributions.png"))

# ==============================================================================
# SECTION 11. Copy this script's own outputs into paper/ — self-contained
# (mirrors run_all.R's copy_tree(), scoped only to output/*/randomization).
# ==============================================================================

copy_randomization_outputs <- function(sub, out_root, paper_root) {
  src <- here("output", out_root, sub)
  dst <- here("paper", paper_root, sub)
  if (!dir.exists(src)) return(invisible(0L))
  dir.create(dst, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(src, full.names = FALSE)
  for (f in files) {
    file.copy(file.path(src, f), file.path(dst, f), overwrite = TRUE)
  }
  length(files)
}

if (dir.exists(here("paper"))) {  # no-op in the stand-alone replication package
  n_tab_copied <- copy_randomization_outputs("randomization", "tables",  "Tables")
  n_fig_copied <- copy_randomization_outputs("randomization", "figures", "Figures")
  message(sprintf("Copied %d table file(s) -> paper/Tables/randomization/", n_tab_copied))
  message(sprintf("Copied %d figure file(s) -> paper/Figures/randomization/", n_fig_copied))
} else message("paper/ not found -- exhibits are left in output/ only")

# ==============================================================================
# SECTION 12. Wall-time summary
# ==============================================================================

total_secs <- as.numeric(Sys.time() - SCRIPT_T0, units = "secs")
message(sprintf("\n=== 08_randomization_inference.R: complete in %.1f min (%d cores) ===\n",
                total_secs / 60, N_CORES))
