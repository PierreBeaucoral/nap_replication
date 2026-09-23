# ==============================================================================
# 09_remarking_decomposition.R
# APPENDIX ANALYSIS -- addresses the relabeling concern: is the
# post-NAP rise in adaptation-marked commitments new money, or re-marking of
# continuing activities / NAP-support-and-readiness pipelines?
#
# Stage 09 of run_all.R. It is the only analysis script that reads CRS
# activity-level microdata (via its cached extract); every
# other script consumes data/processed/simple_panel_wgi.csv (built by
# code/01_prepare_data.R). Its outputs reconcile exactly against that panel's
# `commitments` column (Section 3 below).
#
# Inputs :
#   data/raw/CRS/CRS <year> data.txt   (2009-2024, read ONE YEAR AT A TIME)
#   data/processed/simple_panel_wgi.csv (built by code/01_prepare_data.R)
#
# Outputs:
#   data/processed/crs_adaptation_activities/crs_adaptation_activities_<year>.csv.gz
#                                                   (cached activity extract;
#     ALL CRS activities 2009-2024 for the 144 simple_panel_wgi recipients,
#     regional/multi-country flows excluded -- NOT filtered to
#     ClimateAdaptation-marked rows only, because the margin decomposition
#     needs each activity's marking status in years it was UNMARKED too.
#     One gzip-compressed CSV per year, 2009-2024. Delete the folder to
#     force a raw-data reread.)
#   output/tables/remarking/att_remarking_margins.tex      (tab:remarking_margins)
#   output/tables/remarking/att_remarking_counts.tex       (tab:remarking_counts)
#   output/tables/remarking/att_remarking_exclusions.tex   (tab:remarking_exclusions)
#   output/tables/remarking/tab_remarking_sectors.tex      (tab:remarking_sectors)
#   output/tables/remarking/remarking_flag_shares_by_year.tex   (descriptive)
#   output/tables/remarking/remarking_flag_shares_by_group.tex  (descriptive)
#   output/tables/remarking/id_linkage_diagnostic.csv      (diagnostic)
#   output/tables/remarking/reconciliation_check.csv       (diagnostic)
#   output/tables/remarking/remarking_results.rds          (all result objects)
#   output/figures/remarking/fig_remarking_sectors.png     (Fig: sector composition)
#   -- all .tex and .png outputs above are additionally copied to
#      paper/Tables/remarking/ and paper/Figures/remarking/ (this script does
#      its own copy since it is not wired into run_all.R).
#
# ============================================================
# Paper-to-Code Naming Map
# ============================================================
# Concept                              | Code name(s)
# Activity identifier (composite key)  | activity_id
# New activity (first year observed)   | commit_new / log_commit_new
# Continuing, already marked at first  |
#   observation ("old money")          | commit_cont_marked / log_commit_cont_marked
# Re-marked (first obs. unmarked)      | commit_remarked / log_commit_remarked
# Count of marked activities           | n_marked / log_n_marked
# Mean commitment per marked activity  | mean_commit_marked / log_mean_commit_marked
# Marker penetration (pp of activities)| share_marked_pp
# Exclusion outcomes (title/purpose/   |
#   fund/all/multilateral-ex-fund)     | commit_excl_* / log_commit_excl_* ,
#                                      | commit_multi_excl_fund / log_commit_multi_excl_fund
# Headline (all marked, reference)     | log_commits (already in simple_panel_wgi.csv)
# Group-time ATT                       | gt_obj  (did::att_gt())
# Correct joint pre-trend Wald         | compute_pretrend_wald_correct()
# did's own built-in pre-test p-value  | did_wpval  (att_gt()$Wpval)
# Cohort (year of first NAP adoption)  | cohort_year (0 = never-treated)
# Main-sample adopters (treated)       | cohort_year %in% 2021:2024 (n = 40)
# ============================================================

# ==============================================================================
# SECTION 0. PACKAGES, SEED, VERSION GUARD, DIRECTORIES, LOGGING
# ==============================================================================

library(data.table)
library(MASS)   # ginv() fallback for a near-singular pre-trend covariance; loaded
                # BEFORE dplyr so dplyr::select() (not MASS::select()) wins the mask
library(dplyr)
library(tidyr)
library(ggplot2)
library(xtable)
library(here)
library(did)

if (utils::packageVersion("did") < "2.5.0") {
  stop(sprintf(paste0(
    "did %s is installed but this script's SEs / pre-trend tests require ",
    "did >= 2.5.0 (renv.lock pins it). Run renv::restore() or update did."),
    utils::packageVersion("did")))
}

set.seed(20240601)  # global seed; local set.seed(1242) immediately before each att_gt()

t_script_start <- Sys.time()

dir.create(here("output", "tables",  "remarking"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "figures", "remarking"), recursive = TRUE, showWarnings = FALSE)
has_paper <- dir.exists(here("paper"))  # FALSE in the stand-alone replication package
if (has_paper) {
  dir.create(here("paper",  "Tables",  "remarking"), recursive = TRUE, showWarnings = FALSE)
  dir.create(here("paper",  "Figures", "remarking"), recursive = TRUE, showWarnings = FALSE)
}
dir.create(here("data",   "processed"),            recursive = TRUE, showWarnings = FALSE)

# Run log (timings, row counts, reconciliation result); not a paper exhibit.
log_dir  <- here("output", "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, "09_remarking_decomposition_log.txt")

log_msg <- function(...) {
  txt <- sprintf(...)
  message(txt)
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", txt, "\n",
      file = log_file, append = TRUE)
}

log_msg("=== 09_remarking_decomposition.R: START ===")

# ------------------------------------------------------------------------
# Helpers duplicated from 01_prepare_data.R / 03_main_results.R (additive
# script; not editing the source scripts, so the small shared helpers are
# copied here verbatim, matching the project's own "additive-only" pattern
# used for §3c in 01_prepare_data.R).
# ------------------------------------------------------------------------

esc_header <- function(x) {
  x <- gsub("%",  "\\\\%",  x)
  x <- gsub("_",  "\\\\_",  x)
  x <- gsub("#",  "\\\\#",  x)
  x
}

write_tex_float <- function(out_path, caption_title, label,
                             tabular_lines, notes_text, source_text,
                             size = "\\small") {
  inner <- tabular_lines
  is_table_open  <- grepl("^\\\\begin\\{table\\}", inner)
  is_table_close <- grepl("^\\\\end\\{table\\}", inner)
  if (any(is_table_open))  inner <- inner[!is_table_open]
  if (any(is_table_close)) inner <- inner[!is_table_close]
  inner <- inner[!grepl("^\\\\centering", inner)]
  inner <- inner[!grepl("^\\\\caption", inner)]
  inner <- inner[!grepl("^\\\\label", inner)]
  while (length(inner) > 0L && trimws(inner[1]) == "") inner <- inner[-1]
  while (length(inner) > 0L && trimws(inner[length(inner)]) == "")
    inner <- inner[-length(inner)]

  tab_start <- which(grepl("^\\\\begin\\{tabular", inner))[1]
  tab_end   <- which(grepl("^\\\\end\\{tabular",   inner))[1]

  if (!is.na(tab_start) && !is.na(tab_end)) {
    inner <- c(
      if (tab_start > 1L) inner[seq_len(tab_start - 1L)] else character(0),
      "\\adjustbox{max width=\\textwidth}{%",
      inner[tab_start:tab_end],
      "}",
      if (tab_end < length(inner)) inner[seq(tab_end + 1L, length(inner))] else character(0)
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

build_wide_tex_table <- function(row_labels, col_labels, col_data, tex_label) {
  tab <- data.frame(` ` = row_labels, check.names = FALSE, stringsAsFactors = FALSE)
  for (i in seq_along(col_labels)) {
    tab[[esc_header(col_labels[i])]] <- col_data[[i]]
  }
  xtab <- xtable(tab, label = tex_label)
  align(xtab) <- paste0("ll", paste(rep("c", length(col_labels)), collapse = ""))
  capture.output(
    print(xtab, include.rownames = FALSE, booktabs = TRUE,
          sanitize.text.function = identity, size = "\\small", floating = FALSE)
  )
}

crs_path <- function(y) {
  for (nm in c(paste0("CRS ", y, " Data.txt"), paste0("CRS ", y, " data.txt"))) {
    p <- here("data", "raw", "CRS", nm)
    if (file.exists(p)) return(p)
  }
  stop("CRS file for year ", y, " not found in data/raw/CRS/.")
}

# Regional/multi-country flow exclusion list -- IDENTICAL to
# 01_prepare_data.R §3c/§4 `regionalflows_dry` / `regionalflows` (duplicated,
# not reordered, so this script's reconciliation is exact).
# Mirrors the single `regionalflows` constant in 01_prepare_data.R. Code 860
# (Federated States of Micronesia) is a country, not a regional aggregate, and
# was removed from both lists on 2026-09-15; 1034 is the genuine "Micronesia,
# regional" code and stays excluded.
regionalflows_dry <- c(88, 89, 189, 237, 289, 298, 389, 489, 498, 589, 619,
                       679, 689, 789, 798, 889, 1027:1035, 9998)

# Donor-type code lists -- IDENTICAL to 01_prepare_data.R §3 `donor_list`
# (needed only for the "multilateral cell excluding fund-flagged" cut,
# objective 3(v)).
donor_list <- list(
  DAC_members = c(
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 18, 20, 21, 22, 40, 50, 61, 68,
    69, 75, 76, 301, 302, 701, 742, 801, 820, 918, 30, 45, 55, 62, 70, 72,
    77, 82, 83, 84, 87, 130, 133, 358, 543, 546, 552, 561, 566, 576, 611,
    613, 732, 764, 765
  ),
  Multilateral_donors = c(
    104, 807, 811, 812, 901, 902, 903, 905, 906, 907, 909, 913, 914, 915,
    921, 923, 926, 928, 932, 940, 944, 948, 951, 952, 953, 954, 956, 958,
    959, 960, 963, 964, 966, 967, 971, 974, 976, 978, 979, 980, 981, 982,
    983, 988, 990, 992, 997, 1011, 1012, 1013, 1014, 1015, 1016, 1017,
    1018, 1019, 1020, 1023, 1024, 1025, 1037, 1038, 1039, 1058,
    1311, 1312, 1403
  ),
  Private_donors = c(1601:1639)
)

# ==============================================================================
# SECTION 1. LOAD PANEL, RECONSTRUCT cohort_year / thin cohorts / treated
# (self-contained duplicate of code/03_main_results.R lines ~132-230, so this
# script never sources 03 and never risks re-running its heavy exhibits)
# ==============================================================================

log_msg("Section 1: loading simple_panel_wgi.csv and reconstructing cohorts")

panel_raw <- fread(here("data", "processed", "simple_panel_wgi.csv"))
panel_raw <- as.data.frame(panel_raw)

did_panel <- panel_raw %>%
  mutate(country_id = as.integer(factor(recipient_name)))

first_year <- min(did_panel$year)
last_year  <- max(did_panel$year)

country_gname <- did_panel %>%
  group_by(recipient_name) %>%
  summarise(
    nap_year_c = suppressWarnings(min(nap_year, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    nap_year_c = if_else(is.infinite(nap_year_c), NA_real_, nap_year_c),
    cohort_year = case_when(
      is.na(nap_year_c)                    ~ 0,
      nap_year_c < first_year              ~ 0,
      nap_year_c > last_year               ~ 0,
      TRUE                                 ~ as.numeric(nap_year_c)
    )
  )

did_panel <- did_panel %>%
  select(-any_of("cohort_year")) %>%
  left_join(country_gname %>% select(recipient_name, cohort_year), by = "recipient_name") %>%
  mutate(log_population = log(population))

cohort_sizes <- did_panel %>%
  filter(cohort_year > 0) %>%
  distinct(recipient_name, cohort_year) %>%
  count(cohort_year, name = "n_treated") %>%
  arrange(cohort_year)

thin_threshold <- 5L
thin_cohorts   <- cohort_sizes$cohort_year[cohort_sizes$n_treated < thin_threshold]

did_panel_126 <- did_panel %>% filter(!(cohort_year %in% thin_cohorts))

n_countries_144 <- n_distinct(did_panel$recipient_name)
n_countries_126 <- n_distinct(did_panel_126$recipient_name)

main_adopters <- did_panel_126 %>%
  filter(cohort_year >= 2021) %>%
  distinct(recipient_name) %>%
  pull(recipient_name)

never_treated_126 <- did_panel_126 %>%
  filter(cohort_year == 0) %>%
  distinct(recipient_name) %>%
  pull(recipient_name)

stopifnot(
  # Reference counts come from the stored 03 headline fit (output/fits/), not
  # literals (the 2026-09-15 panel revision restored recipient code 860).
  "estimation-sample size differs from the stored 03 fit" =
    n_countries_126 == readRDS(here("output", "fits", "headline_adaptation_dr_bs.rds"))$n_country,
  "full panel must contain the estimation sample" = n_countries_144 >= n_countries_126,
  "main-sample adopters must be non-empty" = length(main_adopters) > 0L
)

log_msg("Panel: %d countries (144-panel), %d in 126-country estimation sample, %d main adopters (cohort>=2021), %d never-treated",
        n_countries_144, n_countries_126, length(main_adopters), length(never_treated_126))

panel_recipients <- unique(panel_raw$recipient_name)  # 144 names, exact CRS RecipientName strings

# ==============================================================================
# SECTION 2. RAW CRS ACTIVITY-LEVEL EXTRACTION (cache-aware, one year at a time)
# ==============================================================================

cache_dir <- here("data", "processed", "crs_adaptation_activities")

crs_cols <- c("Year", "DonorCode", "DonorName", "RecipientCode", "RecipientName",
              "CrsID", "ProjectNumber", "ProjectTitle", "PurposeCode", "PurposeName",
              "ClimateAdaptation", "ClimateMitigation", "USD_Commitment_Defl", "FlowCode")

# Column classes of the activity extract as produced by the raw read below.
# The cache is read back with exactly these classes (and "NA" as the only
# missing-value string) so the cached object is identical to a fresh raw read.
activity_col_classes <- c(
  Year = "integer", DonorCode = "integer", DonorName = "character",
  RecipientCode = "integer", RecipientName = "character", CrsID = "character",
  ProjectNumber = "character", ProjectTitle = "character", PurposeCode = "integer",
  PurposeName = "character", ClimateAdaptation = "integer", ClimateMitigation = "integer",
  USD_Commitment_Defl = "numeric", FlowCode = "integer")

activity_cache_file <- function(y) {
  file.path(cache_dir, sprintf("crs_adaptation_activities_%d.csv.gz", y))
}

#' Read one year of the cached activity extract (gzip-compressed CSV)
#'
#' Decompressed with base R's gzfile() so that data.table does not need the
#' optional R.utils package to read .gz files.
#' @param y integer year
#' @return data.table with the columns and classes in activity_col_classes
read_activity_cache_year <- function(y) {
  con <- gzfile(activity_cache_file(y), encoding = "UTF-8")
  on.exit(close(con))
  fread(text = readLines(con, warn = FALSE), na.strings = "NA",
        colClasses = activity_col_classes, encoding = "UTF-8")
}

#' Read one CRS year, restricted to the columns and recipients we need
#'
#' @param y integer year
#' @param recipients character vector of RecipientName values to keep
#' @param regional_codes integer vector of RecipientCode values to exclude
#' @return data.table, one row per activity-year record
read_crs_year_activities <- function(y, recipients, regional_codes) {
  t0 <- Sys.time()
  dt <- fread(crs_path(y), select = crs_cols, encoding = "UTF-8", showProgress = FALSE)
  n_raw <- nrow(dt)
  dt <- dt[!RecipientCode %in% regional_codes]
  dt <- dt[RecipientName %in% recipients]
  dt[, CrsID         := trimws(CrsID)]
  dt[, ProjectNumber := trimws(ProjectNumber)]
  dt[CrsID == "",         CrsID := NA_character_]
  dt[ProjectNumber == "", ProjectNumber := NA_character_]
  log_msg("  Year %d: read %d raw rows -> kept %d after recipient filter (%.1fs)",
          y, n_raw, nrow(dt), as.numeric(Sys.time() - t0, units = "secs"))
  dt[]
}

years_needed <- 2009:2024

# The cache is valid only if it was built for exactly the current recipient set
# (a panel revision must trigger a rebuild, not a silent stale hit).
cache_valid <- FALSE
if (all(file.exists(vapply(years_needed, activity_cache_file, character(1L))))) {
  activities <- rbindlist(lapply(years_needed, read_activity_cache_year))
  cache_valid <- setequal(unique(activities$RecipientName), panel_recipients)
  if (!cache_valid) {
    log_msg("Section 2: cache found but its recipient set (%d) differs from the current panel (%d) -- rebuilding",
            uniqueN(activities$RecipientName), length(panel_recipients))
    rm(activities); gc(verbose = FALSE)
  }
}
if (cache_valid) {
  log_msg("Section 2: cache hit -- loaded %s (recipient set verified)", cache_dir)
} else {
  log_msg("Section 2: cache miss -- reading %d raw CRS years one at a time",
          length(years_needed))
  t_read_start <- Sys.time()
  chunks <- vector("list", length(years_needed))
  for (i in seq_along(years_needed)) {
    chunks[[i]] <- read_crs_year_activities(years_needed[i], panel_recipients, regionalflows_dry)
    gc(verbose = FALSE)  # never hold more than one raw year in memory
  }
  activities <- rbindlist(chunks, use.names = TRUE, fill = TRUE)
  rm(chunks); gc(verbose = FALSE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  for (y in years_needed) {
    fwrite(activities[activities$Year == y], activity_cache_file(y), na = "NA")
  }
  log_msg("Section 2: raw read complete in %.1f min; %d rows cached to %s",
          as.numeric(Sys.time() - t_read_start, units = "mins"),
          nrow(activities), cache_dir)
}

log_msg("Activity extract: %d rows, %d years, %d recipients",
        nrow(activities), uniqueN(activities$Year), uniqueN(activities$RecipientName))

# ==============================================================================
# SECTION 3. RECONCILIATION -- activity-level sums must equal panel$commitments
# ==============================================================================

log_msg("Section 3: reconciliation check against simple_panel_wgi.csv$commitments")

recon_activity <- activities[
  ClimateAdaptation %in% c(1L, 2L),
  .(my_commitments = sum(USD_Commitment_Defl, na.rm = TRUE)),
  by = .(recipient_name = RecipientName, year = Year)
]

recon_panel <- as.data.table(panel_raw)[, .(recipient_name, year, commitments)]
recon_merged <- merge(recon_panel, recon_activity, by = c("recipient_name", "year"), all.x = TRUE)
recon_merged[is.na(my_commitments), my_commitments := 0]
recon_merged[, diff := abs(my_commitments - commitments)]

max_diff <- max(recon_merged$diff)
n_mismatch <- sum(recon_merged$diff >= 1e-6)

fwrite(recon_merged[order(-diff)], here("output", "tables", "remarking", "reconciliation_check.csv"))

log_msg("Reconciliation: max |diff| = %.10f USD m across %d recipient-years (%d exceed 1e-6)",
        max_diff, nrow(recon_merged), n_mismatch)

stopifnot(
  "Reconciliation FAILED: activity-level sums do not match simple_panel_wgi.csv$commitments to 1e-6 USD m" =
    max_diff < 1e-6
)

log_msg("Reconciliation PASSED (tolerance 1e-6 USD m).")

# ==============================================================================
# SECTION 4. IDENTIFIER LINKAGE DIAGNOSTIC -- CrsID vs ProjectNumber
# ==============================================================================

log_msg("Section 4: identifier linkage diagnostic (2015-2024)")

id_link_check <- function(dt, idcol) {
  sub <- dt[Year %in% 2015:2024]
  ids <- sub[[idcol]]
  ok  <- !is.na(ids)
  na_share <- 1 - mean(ok)
  sub_ok <- sub[ok]
  by_id <- sub_ok[, .(n_years = uniqueN(Year)), by = c(idcol)]
  list(
    n_activities_total   = nrow(sub),
    na_share              = na_share,
    n_distinct_ids        = nrow(by_id),
    share_multi_year      = mean(by_id$n_years >= 2L)
  )
}

res_crsid <- id_link_check(activities, "CrsID")
res_pn    <- id_link_check(activities, "ProjectNumber")

id_diag <- data.frame(
  identifier          = c("CrsID", "ProjectNumber"),
  n_activities_2015_24 = c(res_crsid$n_activities_total, res_pn$n_activities_total),
  na_share             = round(c(res_crsid$na_share, res_pn$na_share), 4),
  n_distinct_ids       = c(res_crsid$n_distinct_ids, res_pn$n_distinct_ids),
  share_multi_year     = round(c(res_crsid$share_multi_year, res_pn$share_multi_year), 4)
)
fwrite(id_diag, here("output", "tables", "remarking", "id_linkage_diagnostic.csv"))
log_msg("Identifier diagnostic:\n%s", paste(capture.output(print(id_diag, row.names = FALSE)), collapse = "\n"))

chosen_id <- if (res_crsid$share_multi_year >= res_pn$share_multi_year) "CrsID" else "ProjectNumber"
log_msg("Chosen linkage identifier: %s (higher share of ids observed in >=2 distinct years, 2015-2024)",
        chosen_id)

# Composite activity_id = DonorCode_RecipientCode_<chosen id>. Rows with a
# missing chosen id cannot be linked across years and are conservatively
# assigned a unique singleton id (always classified as "new" if marked, since
# continuity cannot be verified) -- documented limitation.
activities[, id_raw := get(chosen_id)]
activities[, row_seq := .I]
activities[, is_singleton := is.na(id_raw)]
activities[, activity_id := ifelse(
  is_singleton,
  paste0("SINGLETON_", row_seq),
  paste(DonorCode, RecipientCode, id_raw, sep = "_")
)]
n_singleton <- sum(activities$is_singleton)
log_msg("%d / %d activity rows (%.2f%%) lack a resolvable %s and are treated as singleton (always-new) activities",
        n_singleton, nrow(activities), 100 * n_singleton / nrow(activities), chosen_id)

# ==============================================================================
# SECTION 5. ROW-LEVEL FLAGS (NAP-support / readiness / fund exclusion, §3)
# ==============================================================================

log_msg("Section 5: constructing NAP-support / readiness / fund flags")

# Title regex (co-occurrence rule). The original single regex included the
# bare token "planning", which matches family planning / urban planning /
# land-use planning titles that have nothing to do with NAP support. Fix:
# flag directly on the NAP/readiness core phrasing; flag the generic word
# "planning" only when it co-occurs with an adaptation/climate/NAP/resilience
# term in the SAME title (title_pattern_planning_only tracks the marginal
# contribution of this second branch alone, for the flag-share diagnostics).
title_pattern_core     <- "national adaptation plan|adaptation plan|\\bnap\\b|nap-|readiness"
title_pattern_planning <- "planning"
title_pattern_context  <- "adaptation|climate|resilien|nap"

title_lc                 <- tolower(activities$ProjectTitle)
flag_title_core_hit      <- grepl(title_pattern_core, title_lc)
flag_title_planning_hit  <- grepl(title_pattern_planning, title_lc) & grepl(title_pattern_context, title_lc)
activities[, flag_title               := flag_title_core_hit | flag_title_planning_hit]
activities[, flag_title_planning_only := flag_title_planning_hit & !flag_title_core_hit]

# Pre-fix (broad) title regex, retained ONLY to build the before/after
# regex-sensitivity comparison in Section 9b -- not used for any headline
# exhibit.
title_pattern_old <- "national adaptation plan|nap |readiness|nap-|adaptation plan|planning"
activities[, flag_title_old := grepl(title_pattern_old, title_lc)]

# Purpose flag narrowed to 41010 (environmental policy & administrative
# management) only. 15110 (public sector policy & admin. management) and
# 43010 (multisector aid) are too generic on their own to indicate
# NAP-support/readiness activity -- kept as separate diagnostic flags so the
# flag-share tables can report the commitment share attributable to each
# code alone, and the exclusion table can show what the narrower rule now
# retains that the old rule used to exclude.
activities[, flag_purpose       := PurposeCode == 41010L]
activities[, flag_purpose_15110 := PurposeCode == 15110L]
activities[, flag_purpose_43010 := PurposeCode == 43010L]
activities[, flag_purpose_old   := PurposeCode %in% c(41010L, 15110L, 43010L)]

fund_donors <- c("Green Climate Fund", "Adaptation Fund", "Global Environment Facility",
                 "Least Developed Countries Fund")
activities[, flag_fund := DonorName %in% fund_donors]

activities[, donor_type := fcase(
  DonorCode %in% donor_list$DAC_members,        "Bilateral_members",
  DonorCode %in% donor_list$Multilateral_donors, "Multilateral_donors",
  default = "Other"
)]

# ==============================================================================
# SECTION 6. ACTIVITY-YEAR AGGREGATES (commit_sum, marked status, flags)
# ==============================================================================

log_msg("Section 6: collapsing to activity-year level")

# commit_sum must total ONLY the marked-row commitment amount for an
# (activity_id, Year), not USD_Commitment_Defl across ALL rows for that
# activity-year: unmarked CRS line-items (amendments/co-financing rows without
# the Rio adaptation marker) would leak non-adaptation value into the
# margin/exclusion outcomes. 1,852 / 2,998,080 activity-years mix marked and
# unmarked rows (pooled excess ~USD 4,475m if summed indiscriminately; the
# Section 7 reconciliation below would fail). Restricting commit_sum -- and
# the NAP-support/fund flags -- to is_marked_row rows closes this exactly.
activities[, is_marked_row := ClimateAdaptation %in% c(1L, 2L)]

act_year <- activities[, .(
  commit_sum   = sum(USD_Commitment_Defl[is_marked_row], na.rm = TRUE),
  marked       = as.integer(any(is_marked_row)),
  flag_title   = as.integer(any(flag_title[is_marked_row])),
  flag_purpose = as.integer(any(flag_purpose[is_marked_row])),
  flag_fund    = as.integer(any(flag_fund[is_marked_row])),
  flag_title_planning_only = as.integer(any(flag_title_planning_only[is_marked_row])),
  flag_purpose_15110       = as.integer(any(flag_purpose_15110[is_marked_row])),
  flag_purpose_43010       = as.integer(any(flag_purpose_43010[is_marked_row])),
  flag_title_old           = as.integer(any(flag_title_old[is_marked_row])),
  flag_purpose_old         = as.integer(any(flag_purpose_old[is_marked_row])),
  donor_type_multi = as.integer(any(donor_type[is_marked_row] == "Multilateral_donors" &
                                    !flag_fund[is_marked_row]))
), by = .(activity_id, RecipientName, Year)]

first_obs <- act_year[, .(first_year = min(Year)), by = activity_id]
act_year  <- merge(act_year, first_obs, by = "activity_id")

first_status <- act_year[Year == first_year, .(activity_id, first_marked = marked)]
act_year <- merge(act_year, first_status, by = "activity_id")

act_year[, margin := fcase(
  marked == 0L,               NA_character_,
  Year == first_year,         "new",
  first_marked == 1L,         "continuing_marked",
  first_marked == 0L,         "remarked",
  default = NA_character_
)]

# Four-bucket margin (Section 7b): splits "new" into "new_linked" (a
# resolvable cross-year identifier exists) and "unlinked" (no resolvable id
# -- SINGLETON_* activity_id, Section 4 -- so the activity is classified new
# by construction, since continuity cannot be verified). A singleton
# activity_id is, by construction, observed in exactly one Year, so
# Year == first_year always holds for singletons when marked; unlinked can
# therefore never overlap continuing_marked or remarked.
act_year[, is_singleton_activity := grepl("^SINGLETON_", activity_id)]
act_year[, margin4 := fcase(
  marked == 0L,                    NA_character_,
  is_singleton_activity == TRUE,   "unlinked",
  Year == first_year,              "new_linked",
  first_marked == 1L,              "continuing_marked",
  first_marked == 0L,              "remarked",
  default = NA_character_
)]

log_msg("Activity-year table: %d rows, %d distinct activities", nrow(act_year), uniqueN(act_year$activity_id))

# ==============================================================================
# SECTION 6b. SINGLETON-ID SHARE DIAGNOSTIC (feeds tab:remarking_margins and
# tab:remarking_margins_unlinked notes)
# ==============================================================================

log_msg("Section 6b: singleton-id share among marked rows (overall and post-adoption treated cells)")

marked_rows_dt <- activities[is_marked_row == TRUE]

singleton_share_rows_overall   <- mean(marked_rows_dt$is_singleton)
singleton_share_commit_overall <- sum(marked_rows_dt$USD_Commitment_Defl[marked_rows_dt$is_singleton], na.rm = TRUE) /
                                   sum(marked_rows_dt$USD_Commitment_Defl, na.rm = TRUE)

cohort_lookup_act <- as.data.table(did_panel_126 %>% distinct(recipient_name, cohort_year))
marked_rows_dt <- merge(marked_rows_dt, cohort_lookup_act,
                         by.x = "RecipientName", by.y = "recipient_name", all.x = TRUE)
marked_rows_dt[, is_post_treated_cell := !is.na(cohort_year) & cohort_year >= 2021 & Year >= cohort_year]

post_treated_rows <- marked_rows_dt[is_post_treated_cell == TRUE]
n_post_treated_rows <- nrow(post_treated_rows)
singleton_share_rows_post   <- if (n_post_treated_rows > 0L) mean(post_treated_rows$is_singleton) else NA_real_
singleton_share_commit_post <- if (n_post_treated_rows > 0L) {
  sum(post_treated_rows$USD_Commitment_Defl[post_treated_rows$is_singleton], na.rm = TRUE) /
    sum(post_treated_rows$USD_Commitment_Defl, na.rm = TRUE)
} else NA_real_

log_msg("Singleton share -- overall marked rows: %.2f%% of rows, %.2f%% of commitment value (n rows = %d)",
        100 * singleton_share_rows_overall, 100 * singleton_share_commit_overall, nrow(marked_rows_dt))
log_msg("Singleton share -- post-adoption treated recipient-years: %.2f%% of rows, %.2f%% of commitment value (n rows = %d)",
        100 * singleton_share_rows_post, 100 * singleton_share_commit_post, n_post_treated_rows)
# Exported so the appendix text quoting these shares is traceable.
fwrite(data.table(
  sample      = c("all marked rows 2009-2024", "post-adoption treated recipient-years"),
  share_rows  = c(singleton_share_rows_overall, singleton_share_rows_post),
  share_value = c(singleton_share_commit_overall, singleton_share_commit_post),
  n_rows      = c(nrow(marked_rows_dt), n_post_treated_rows)
), here("output", "tables", "remarking", "singleton_shares.csv"))

# ==============================================================================
# SECTION 7. OBJECTIVE 1 -- MARGIN DECOMPOSITION (new / continuing / remarked)
# ==============================================================================

log_msg("Section 7: margin decomposition by recipient-year")

margins_long <- act_year[marked == 1L & !is.na(margin),
                          .(value = sum(commit_sum, na.rm = TRUE)),
                          by = .(recipient_name = RecipientName, year = Year, margin)]

margins_wide <- dcast(margins_long, recipient_name + year ~ margin,
                      value.var = "value", fill = 0)
for (nm in c("new", "continuing_marked", "remarked")) {
  if (!nm %in% names(margins_wide)) margins_wide[[nm]] <- 0
}
setnames(margins_wide,
         old = c("new", "continuing_marked", "remarked"),
         new = c("commit_new", "commit_cont_marked", "commit_remarked"))

did_panel_126 <- did_panel_126 %>%
  select(-any_of(c("commit_new", "commit_cont_marked", "commit_remarked"))) %>%
  left_join(as.data.frame(margins_wide), by = c("recipient_name", "year")) %>%
  mutate(
    commit_new          = coalesce(commit_new, 0),
    commit_cont_marked  = coalesce(commit_cont_marked, 0),
    commit_remarked     = coalesce(commit_remarked, 0),
    log_commit_new         = log1p(commit_new),
    log_commit_cont_marked = log1p(commit_cont_marked),
    log_commit_remarked    = log1p(commit_remarked)
  )

# Internal consistency check: the three margins must sum to total marked
# commitments (== log_commits' underlying `commitments` column).
margin_sum_check <- did_panel_126 %>%
  mutate(margin_total = commit_new + commit_cont_marked + commit_remarked,
         diff = abs(margin_total - commitments)) %>%
  summarise(max_diff = max(diff, na.rm = TRUE)) %>%
  pull(max_diff)
log_msg("Margin components sum to headline commitments: max |diff| = %.6f USD m", margin_sum_check)
stopifnot("Margin decomposition does not sum to headline commitments" = margin_sum_check < 1e-6)

# ==============================================================================
# SECTION 7b. FOUR-BUCKET MARGIN DIAGNOSTIC -- splits "new" into "new_linked"
# and "unlinked" (singleton) so the reader can see how much of the headline
# "new" ATT is attributable to activities whose continuity cannot be
# verified at all (Section 6b). "Continuing, already marked" and
# "re-marked" are unaffected (a singleton is, by construction, always
# classified as new -- see margin4 above), so only the "new" split is new.
# ==============================================================================

log_msg("Section 7b: four-bucket margin decomposition (new_linked / unlinked / continuing / re-marked)")

margins4_long <- act_year[marked == 1L & !is.na(margin4),
                           .(value = sum(commit_sum, na.rm = TRUE)),
                           by = .(recipient_name = RecipientName, year = Year, margin4)]
margins4_wide <- dcast(margins4_long, recipient_name + year ~ margin4,
                       value.var = "value", fill = 0)
for (nm in c("new_linked", "unlinked", "continuing_marked", "remarked")) {
  if (!nm %in% names(margins4_wide)) margins4_wide[[nm]] <- 0
}
setnames(margins4_wide, old = c("new_linked", "unlinked"),
         new = c("commit_new_linked", "commit_unlinked"))

did_panel_126 <- did_panel_126 %>%
  select(-any_of(c("commit_new_linked", "commit_unlinked"))) %>%
  left_join(as.data.frame(margins4_wide[, .(recipient_name, year, commit_new_linked, commit_unlinked)]),
            by = c("recipient_name", "year")) %>%
  mutate(
    commit_new_linked     = coalesce(commit_new_linked, 0),
    commit_unlinked       = coalesce(commit_unlinked, 0),
    log_commit_new_linked = log1p(commit_new_linked),
    log_commit_unlinked   = log1p(commit_unlinked)
  )

# Internal consistency check: new_linked + unlinked must reconstruct the
# 3-bucket "new" from Section 7 exactly.
split_check <- did_panel_126 %>%
  mutate(diff = abs((commit_new_linked + commit_unlinked) - commit_new)) %>%
  summarise(max_diff = max(diff, na.rm = TRUE)) %>%
  pull(max_diff)
log_msg("New-linked + unlinked reconstructs the 3-bucket 'new': max |diff| = %.6f USD m", split_check)
stopifnot("4-bucket split of 'new' does not reconcile with the 3-bucket margin" = split_check < 1e-6)

# ==============================================================================
# SECTION 8. OBJECTIVE 2 -- COUNTS VS VALUES
# ==============================================================================

log_msg("Section 8: counts-vs-values outcomes by recipient-year")

counts_marked <- act_year[marked == 1L, .(
  n_marked            = uniqueN(activity_id),
  total_marked_commit = sum(commit_sum, na.rm = TRUE)
), by = .(recipient_name = RecipientName, year = Year)]

counts_total <- act_year[, .(
  n_total = uniqueN(activity_id)
), by = .(recipient_name = RecipientName, year = Year)]

counts_all <- merge(counts_total, counts_marked, by = c("recipient_name", "year"), all.x = TRUE)
counts_all[is.na(n_marked), n_marked := 0L]
counts_all[is.na(total_marked_commit), total_marked_commit := 0]
counts_all[, mean_commit_marked := ifelse(n_marked > 0L, total_marked_commit / n_marked, 0)]
counts_all[, share_marked_pp := ifelse(n_total > 0L, 100 * n_marked / n_total, NA_real_)]

did_panel_126 <- did_panel_126 %>%
  select(-any_of(c("n_marked", "mean_commit_marked", "share_marked_pp"))) %>%
  left_join(
    as.data.frame(counts_all[, .(recipient_name, year, n_marked, mean_commit_marked, share_marked_pp)]),
    by = c("recipient_name", "year")
  ) %>%
  mutate(
    n_marked            = coalesce(n_marked, 0),
    mean_commit_marked  = coalesce(mean_commit_marked, 0),
    share_marked_pp     = coalesce(share_marked_pp, 0),
    log_n_marked            = log1p(n_marked),
    log_mean_commit_marked  = log1p(mean_commit_marked)
  )

# ==============================================================================
# SECTION 9. OBJECTIVE 3 -- NAP-SUPPORT / READINESS / FUND EXCLUSION
# ==============================================================================

log_msg("Section 9: NAP-support / readiness / fund exclusion outcomes and shares")

marked_act <- act_year[marked == 1L]

excl_ry <- marked_act[, .(
  commit_total          = sum(commit_sum, na.rm = TRUE),
  commit_flag_title     = sum(commit_sum[flag_title == 1L],   na.rm = TRUE),
  commit_flag_purpose   = sum(commit_sum[flag_purpose == 1L], na.rm = TRUE),
  commit_flag_fund      = sum(commit_sum[flag_fund == 1L],    na.rm = TRUE),
  commit_flag_title_planning_only = sum(commit_sum[flag_title_planning_only == 1L], na.rm = TRUE),
  commit_flag_purpose_15110       = sum(commit_sum[flag_purpose_15110 == 1L],       na.rm = TRUE),
  commit_flag_purpose_43010       = sum(commit_sum[flag_purpose_43010 == 1L],       na.rm = TRUE),
  commit_excl_title     = sum(commit_sum[flag_title == 0L],   na.rm = TRUE),
  commit_excl_purpose   = sum(commit_sum[flag_purpose == 0L], na.rm = TRUE),
  commit_excl_fund      = sum(commit_sum[flag_fund == 0L],    na.rm = TRUE),
  commit_excl_all       = sum(commit_sum[flag_title == 0L & flag_purpose == 0L & flag_fund == 0L], na.rm = TRUE),
  commit_multi_excl_fund = sum(commit_sum[donor_type_multi == 1L], na.rm = TRUE),
  # Pre-fix (old, broad) title/purpose flags -- diagnostic only, feeds the
  # regex-sensitivity comparison in Section 9b.
  commit_excl_title_old   = sum(commit_sum[flag_title_old == 0L], na.rm = TRUE),
  commit_excl_purpose_old = sum(commit_sum[flag_purpose_old == 0L], na.rm = TRUE),
  commit_excl_all_old     = sum(commit_sum[flag_title_old == 0L & flag_purpose_old == 0L & flag_fund == 0L], na.rm = TRUE)
), by = .(recipient_name = RecipientName, year = Year)]

did_panel_126 <- did_panel_126 %>%
  select(-any_of(c("commit_excl_title", "commit_excl_purpose", "commit_excl_fund",
                    "commit_excl_all", "commit_multi_excl_fund",
                    "commit_excl_title_old", "commit_excl_purpose_old", "commit_excl_all_old"))) %>%
  left_join(
    as.data.frame(excl_ry[, .(recipient_name, year, commit_excl_title, commit_excl_purpose,
                              commit_excl_fund, commit_excl_all, commit_multi_excl_fund,
                              commit_excl_title_old, commit_excl_purpose_old, commit_excl_all_old)]),
    by = c("recipient_name", "year")
  ) %>%
  mutate(
    commit_excl_title      = coalesce(commit_excl_title, 0),
    commit_excl_purpose    = coalesce(commit_excl_purpose, 0),
    commit_excl_fund       = coalesce(commit_excl_fund, 0),
    commit_excl_all        = coalesce(commit_excl_all, 0),
    commit_multi_excl_fund = coalesce(commit_multi_excl_fund, 0),
    commit_excl_title_old   = coalesce(commit_excl_title_old, 0),
    commit_excl_purpose_old = coalesce(commit_excl_purpose_old, 0),
    commit_excl_all_old     = coalesce(commit_excl_all_old, 0),
    log_commit_excl_title      = log1p(commit_excl_title),
    log_commit_excl_purpose    = log1p(commit_excl_purpose),
    log_commit_excl_fund       = log1p(commit_excl_fund),
    log_commit_excl_all        = log1p(commit_excl_all),
    log_commit_multi_excl_fund = log1p(commit_multi_excl_fund),
    log_commit_excl_title_old   = log1p(commit_excl_title_old),
    log_commit_excl_purpose_old = log1p(commit_excl_purpose_old),
    log_commit_excl_all_old     = log1p(commit_excl_all_old)
  )

# --- Descriptive flag-coverage shares: by year ------------------------------
flag_shares_year <- excl_ry[, .(
  commit_total        = sum(commit_total),
  commit_flag_title   = sum(commit_flag_title),
  commit_flag_purpose = sum(commit_flag_purpose),
  commit_flag_fund    = sum(commit_flag_fund),
  commit_flag_title_planning_only = sum(commit_flag_title_planning_only),
  commit_flag_purpose_15110       = sum(commit_flag_purpose_15110),
  commit_flag_purpose_43010       = sum(commit_flag_purpose_43010)
), by = year][order(year)]
flag_shares_year[, `:=`(
  share_title_pct   = 100 * commit_flag_title   / commit_total,
  share_purpose_pct = 100 * commit_flag_purpose / commit_total,
  share_fund_pct    = 100 * commit_flag_fund    / commit_total,
  share_title_planning_only_pct = 100 * commit_flag_title_planning_only / commit_total,
  share_purpose_15110_pct       = 100 * commit_flag_purpose_15110       / commit_total,
  share_purpose_43010_pct       = 100 * commit_flag_purpose_43010       / commit_total
)]

# --- Descriptive flag-coverage shares: by treated / never-treated ----------
group_lookup <- did_panel_126 %>%
  distinct(recipient_name, cohort_year) %>%
  mutate(group = case_when(
    cohort_year >= 2021 ~ "Treated (main sample, cohort >= 2021)",
    cohort_year == 0    ~ "Never-treated",
    TRUE                 ~ "Thin-cohort adopter (excluded from main sample)"
  ))

flag_shares_group <- merge(excl_ry, as.data.table(group_lookup), by = "recipient_name")
flag_shares_group <- flag_shares_group[, .(
  commit_total        = sum(commit_total),
  commit_flag_title   = sum(commit_flag_title),
  commit_flag_purpose = sum(commit_flag_purpose),
  commit_flag_fund    = sum(commit_flag_fund),
  commit_flag_title_planning_only = sum(commit_flag_title_planning_only),
  commit_flag_purpose_15110       = sum(commit_flag_purpose_15110),
  commit_flag_purpose_43010       = sum(commit_flag_purpose_43010)
), by = group]
flag_shares_group[, `:=`(
  share_title_pct   = 100 * commit_flag_title   / commit_total,
  share_purpose_pct = 100 * commit_flag_purpose / commit_total,
  share_fund_pct    = 100 * commit_flag_fund    / commit_total,
  share_title_planning_only_pct = 100 * commit_flag_title_planning_only / commit_total,
  share_purpose_15110_pct       = 100 * commit_flag_purpose_15110       / commit_total,
  share_purpose_43010_pct       = 100 * commit_flag_purpose_43010       / commit_total
)]

log_msg("Flag-coverage shares by group:\n%s",
        paste(capture.output(print(flag_shares_group)), collapse = "\n"))

# Write descriptive tex tables (supplementary; not paper-mandatory table names)
flag_row_labels <- c(
  "Total marked commitments (USD M)", "Title-flagged (\\%)",
  "Purpose-flagged (\\%)", "Fund-flagged (\\%)",
  "  of which: planning-only title (\\%)",
  "  of which: purpose 15110 alone (\\%)",
  "  of which: purpose 43010 alone (\\%)"
)

flag_year_tex <- build_wide_tex_table(
  row_labels = flag_row_labels,
  col_labels = as.character(flag_shares_year$year),
  col_data   = lapply(seq_len(nrow(flag_shares_year)), function(i) c(
    sprintf("%.1f", flag_shares_year$commit_total[i]),
    sprintf("%.1f", flag_shares_year$share_title_pct[i]),
    sprintf("%.1f", flag_shares_year$share_purpose_pct[i]),
    sprintf("%.1f", flag_shares_year$share_fund_pct[i]),
    sprintf("%.1f", flag_shares_year$share_title_planning_only_pct[i]),
    sprintf("%.1f", flag_shares_year$share_purpose_15110_pct[i]),
    sprintf("%.1f", flag_shares_year$share_purpose_43010_pct[i])
  )),
  tex_label = "tab:remarking_flag_shares_by_year"
)
write_tex_float(
  out_path      = here("output", "tables", "remarking", "remarking_flag_shares_by_year.tex"),
  caption_title = "Share of adaptation-marked commitments flagged as NAP-support, readiness, or climate-fund activity, by year",
  label         = "tab:remarking_flag_shares_by_year",
  tabular_lines = flag_year_tex,
  notes_text    = paste0(
    "Title-flag: title matches core NAP/adaptation-plan/readiness phrasing, or generic ",
    "``planning'' co-occurring with an adaptation/climate/NAP/resilience term. Purpose-flag: ",
    "CRS 41010 only. Fund-flag: GCF, Adaptation Fund, or GEF (LDC Fund under GEF). ``Of which'' ",
    "rows are diagnostic: ``planning-only title'' = matched only via the generic co-occurrence; ",
    "``purpose 15110/43010 alone'' = shares under the narrowed rule's excluded codes ",
    "(Table~\\ref{tab:remarking_exclusions_regex_comparison} has the ATT sensitivity). ",
    "Categories not mutually exclusive"
  ),
  source_text = "OECD CRS activity-level microdata, Rio adaptation marker 1 or 2"
)

flag_group_tex <- build_wide_tex_table(
  row_labels = flag_row_labels,
  col_labels = flag_shares_group$group,
  col_data   = lapply(seq_len(nrow(flag_shares_group)), function(i) c(
    sprintf("%.1f", flag_shares_group$commit_total[i]),
    sprintf("%.1f", flag_shares_group$share_title_pct[i]),
    sprintf("%.1f", flag_shares_group$share_purpose_pct[i]),
    sprintf("%.1f", flag_shares_group$share_fund_pct[i]),
    sprintf("%.1f", flag_shares_group$share_title_planning_only_pct[i]),
    sprintf("%.1f", flag_shares_group$share_purpose_15110_pct[i]),
    sprintf("%.1f", flag_shares_group$share_purpose_43010_pct[i])
  )),
  tex_label = "tab:remarking_flag_shares_by_group"
)
write_tex_float(
  out_path      = here("output", "tables", "remarking", "remarking_flag_shares_by_group.tex"),
  caption_title = "Share of adaptation-marked commitments flagged as NAP-support, readiness, or climate-fund activity, by treatment group",
  label         = "tab:remarking_flag_shares_by_group",
  tabular_lines = flag_group_tex,
  notes_text    = paste0(
    "Groups follow the main specification's cohort definition: treated = 40 main-sample adopters with adoption year ",
    ">= 2021; thin-cohort adopters (2015-2020, < 5 units) are excluded from the main estimation sample but shown ",
    "for completeness. See Table~\\ref{tab:remarking_flag_shares_by_year} for the title/purpose flag definitions ",
    "and the diagnostic ``of which'' rows"
  ),
  source_text   = "OECD CRS activity-level microdata; UNFCCC NAP Central"
)

# ==============================================================================
# SECTION 10. CS ESTIMATION HELPERS
# ==============================================================================

log_msg("Section 10: CS estimation helpers")

# canonical pre-trend wording — keep byte-identical across scripts
PRETREND_NOTE_AGG <- function(min_e, max_e, k) sprintf(
  "Pre-trend $\\chi^2$ / $p$: joint Wald test on the aggregated pre-treatment event-time coefficients ($%d \\leq e \\leq %d$; %d restrictions), using the influence-function covariance of the dynamic aggregation from the analytical (non-bootstrap) fit; a generalized inverse is used if the block is singular",
  min_e, max_e, k)
PRETREND_NOTE_DID <- "\\texttt{did} pre-test $p$: \\texttt{did}'s built-in Wald test over all pre-period $ATT(g,t)$ cells against each cohort's $g-1$ base year"
PRETREND_NOTE <- function(min_e, max_e, k) paste0(PRETREND_NOTE_AGG(min_e, max_e, k), ". ", PRETREND_NOTE_DID)

#' Correct joint Wald test for pre-treatment event-time ATTs
#'
#' Uses the dynamic-aggregation influence function stored by did >= 2.5.0 at
#' agg_d$inf.function$dynamic.inf.func.e (Sigma = crossprod(IF)/n^2), the
#' same covariance that 03_main_results.R's compute_pretrend_test() uses since
#' its 2026-09-14 revision (it previously used a diagonal approximation).
#'
#' @param agg_d an AGGTEobj from aggte(..., type = "dynamic", min_e = -5)
#' @return list(stat, pval, df)
compute_pretrend_wald_correct <- function(agg_d) {
  keep <- which(!is.na(agg_d$se.egt) & agg_d$se.egt > 1e-10)
  pre_pos <- which(agg_d$egt[keep] < 0L)
  if (length(pre_pos) == 0L) return(list(stat = NA_real_, pval = NA_real_, df = 0L))

  IF <- agg_d$inf.function$dynamic.inf.func.e
  if (is.null(IF)) return(list(stat = NA_real_, pval = NA_real_, df = 0L))

  n <- nrow(IF)
  sigma_full <- crossprod(IF) / n^2
  beta      <- agg_d$att.egt[keep][pre_pos]
  sigma_pre <- sigma_full[keep, keep][pre_pos, pre_pos, drop = FALSE]

  W <- tryCatch(
    as.numeric(t(beta) %*% solve(sigma_pre) %*% beta),
    error = function(e) as.numeric(t(beta) %*% MASS::ginv(sigma_pre) %*% beta)
  )
  df <- length(pre_pos)
  list(stat = round(W, 3), pval = round(pchisq(W, df = df, lower.tail = FALSE), 3), df = df)
}

#' Estimate one outcome with the main-spec CS(2021) settings (matches
#' 03_main_results.R make_wide_table(): analytical fit (bstrap = FALSE) feeds
#' the pre-trend Wald test only; bootstrap fit (bstrap = TRUE, biters = 999)
#' feeds the reported ATT/SE. seed = 1242 immediately before each att_gt().
#'
#' @param panel data.frame; the 126-country estimation sample
#' @param yvar character; outcome column name
#' @return named list of estimation results
estimate_cs_outcome <- function(panel, yvar, seed = 1242L) {
  stopifnot(is.data.frame(panel), yvar %in% names(panel))

  set.seed(seed)
  gt_analytical <- tryCatch(
    att_gt(
      yname = yvar, tname = "year", idname = "country_id", gname = "cohort_year",
      xformla = ~ ge_est + log_population, data = panel,
      est_method = "dr", bstrap = FALSE, cband = FALSE,
      control_group = "nevertreated", anticipation = 0L,
      base_period = "universal", panel = TRUE, allow_unbalanced_panel = TRUE
    ),
    error = function(e) {
      message("  att_gt (analytical) failed for ", yvar, ": ", conditionMessage(e)); NULL
    }
  )

  set.seed(seed)
  gt_boot <- tryCatch(
    att_gt(
      yname = yvar, tname = "year", idname = "country_id", gname = "cohort_year",
      xformla = ~ ge_est + log_population, data = panel,
      est_method = "dr", bstrap = TRUE, biters = 999L, cband = FALSE,
      control_group = "nevertreated", anticipation = 0L,
      base_period = "universal", panel = TRUE, allow_unbalanced_panel = TRUE
    ),
    error = function(e) {
      message("  att_gt (bootstrap) failed for ", yvar, ": ", conditionMessage(e)); NULL
    }
  )

  if (is.null(gt_boot)) {
    return(list(att = NA_real_, se = NA_real_, t = NA_real_,
                n_obs = 0L, n_country = 0L,
                pretrend_stat = NA_real_, pretrend_pval = NA_real_, pretrend_df = 0L,
                did_wpval = NA_real_))
  }

  agg_s <- tryCatch(aggte(gt_boot, type = "simple", na.rm = TRUE), error = function(e) NULL)
  att <- if (!is.null(agg_s)) agg_s$overall.att else NA_real_
  se  <- if (!is.null(agg_s)) agg_s$overall.se  else NA_real_
  tst <- if (!is.na(att) && !is.na(se) && se > 0) att / se else NA_real_

  pt <- list(stat = NA_real_, pval = NA_real_, df = 0L)
  did_wpval <- NA_real_
  if (!is.null(gt_analytical)) {
    did_wpval <- round(gt_analytical$Wpval, 3)
    agg_d <- tryCatch(
      aggte(gt_analytical, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = Inf),
      error = function(e) NULL
    )
    if (!is.null(agg_d)) pt <- compute_pretrend_wald_correct(agg_d)
  }

  rows_oc <- panel[!is.na(panel[[yvar]]), ]

  list(
    att = att, se = se, t = tst,
    n_obs = nrow(rows_oc), n_country = length(unique(rows_oc$country_id)),
    pretrend_stat = pt$stat, pretrend_pval = pt$pval, pretrend_df = pt$df,
    did_wpval = did_wpval
  )
}

fmt_stars <- function(t_v) {
  if (is.na(t_v)) return("")
  if (abs(t_v) > 2.576) "***" else if (abs(t_v) > 1.960) "**" else if (abs(t_v) > 1.645) "*" else ""
}

pretreat_mean <- function(panel, raw_var) {
  rows <- panel %>% filter(cohort_year > 0, year < cohort_year, !is.na(.data[[raw_var]]))
  if (nrow(rows) == 0L) return(NA_real_)
  mean(rows[[raw_var]], na.rm = TRUE)
}

did_panel_126_df <- as.data.frame(did_panel_126)

# ==============================================================================
# SECTION 11. RUN ESTIMATIONS
# ==============================================================================

t_est_start <- Sys.time()

## --- Objective 1: margin decomposition -------------------------------------
log_msg("Section 11a: margin decomposition estimation")

margin_specs <- list(
  list(var = "log_commits",            raw = "commitments",         label = "Headline (all marked)"),
  list(var = "log_commit_new",         raw = "commit_new",          label = "New activities"),
  list(var = "log_commit_cont_marked", raw = "commit_cont_marked",  label = "Continuing, already marked"),
  list(var = "log_commit_remarked",    raw = "commit_remarked",     label = "Re-marked (previously unmarked)")
)

margin_results <- vector("list", length(margin_specs))
pretreat_headline <- pretreat_mean(did_panel_126_df, "commitments")
for (i in seq_along(margin_specs)) {
  sp <- margin_specs[[i]]
  log_msg("  Estimating: %s (%s)", sp$label, sp$var)
  res <- estimate_cs_outcome(did_panel_126_df, sp$var)
  res$label <- sp$label
  res$mean_pre <- pretreat_mean(did_panel_126_df, sp$raw)
  res$share_pre <- if (!is.na(res$mean_pre) && !is.na(pretreat_headline) && pretreat_headline > 0) {
    100 * res$mean_pre / pretreat_headline
  } else NA_real_
  margin_results[[i]] <- res
  log_msg("    ATT=%.4f SE=%.4f t=%.3f pretrend p=%.3f (df=%d) did-Wpval=%.3f mean_pre=%.2f share_pre=%.1f%%",
          res$att, res$se, res$t, res$pretrend_pval, res$pretrend_df, res$did_wpval,
          res$mean_pre, res$share_pre)
}
names(margin_results) <- sapply(margin_specs, `[[`, "label")

## --- Objective 1b: four-bucket margin diagnostic (unlinked/singleton) ------
log_msg("Section 11a-ii: four-bucket margin diagnostic estimation (new-linked vs.\\ unlinked)")

margin4_specs <- list(
  list(var = "log_commit_new_linked", raw = "commit_new_linked", label = "New (linked)"),
  list(var = "log_commit_unlinked",   raw = "commit_unlinked",   label = "Unlinked (singleton, always new)")
)

margin4_results <- vector("list", length(margin4_specs))
for (i in seq_along(margin4_specs)) {
  sp <- margin4_specs[[i]]
  log_msg("  Estimating: %s (%s)", sp$label, sp$var)
  res <- estimate_cs_outcome(did_panel_126_df, sp$var)
  res$label <- sp$label
  res$mean_pre <- pretreat_mean(did_panel_126_df, sp$raw)
  res$share_pre <- if (!is.na(res$mean_pre) && !is.na(pretreat_headline) && pretreat_headline > 0) {
    100 * res$mean_pre / pretreat_headline
  } else NA_real_
  margin4_results[[i]] <- res
  log_msg("    ATT=%.4f SE=%.4f t=%.3f pretrend p=%.3f (df=%d) did-Wpval=%.3f mean_pre=%.2f share_pre=%.1f%%",
          res$att, res$se, res$t, res$pretrend_pval, res$pretrend_df, res$did_wpval,
          res$mean_pre, res$share_pre)
}
names(margin4_results) <- sapply(margin4_specs, `[[`, "label")

## --- Objective 2: counts vs values ------------------------------------------
log_msg("Section 11b: counts-vs-values estimation")

counts_specs <- list(
  list(var = "log_n_marked",            raw = "n_marked",            label = "N marked activities"),
  list(var = "log_mean_commit_marked",  raw = "mean_commit_marked",  label = "Mean commitment / marked activity"),
  list(var = "share_marked_pp",         raw = "share_marked_pp",     label = "Share of activities marked (pp)")
)

counts_results <- vector("list", length(counts_specs))
for (i in seq_along(counts_specs)) {
  sp <- counts_specs[[i]]
  log_msg("  Estimating: %s (%s)", sp$label, sp$var)
  res <- estimate_cs_outcome(did_panel_126_df, sp$var)
  res$label <- sp$label
  res$mean_pre <- pretreat_mean(did_panel_126_df, sp$raw)
  counts_results[[i]] <- res
  log_msg("    ATT=%.4f SE=%.4f t=%.3f pretrend p=%.3f (df=%d) did-Wpval=%.3f mean_pre=%.3f",
          res$att, res$se, res$t, res$pretrend_pval, res$pretrend_df, res$did_wpval, res$mean_pre)
}
names(counts_results) <- sapply(counts_specs, `[[`, "label")

## --- Objective 3: NAP-support / readiness / fund exclusion ------------------
log_msg("Section 11c: NAP-support/readiness/fund exclusion estimation")

exclusion_specs <- list(
  list(var = "log_commits",                  raw = "commitments",            label = "Headline (all marked)"),
  list(var = "log_commit_excl_title",        raw = "commit_excl_title",      label = "Excl. title-flagged"),
  list(var = "log_commit_excl_purpose",      raw = "commit_excl_purpose",    label = "Excl. purpose-flagged"),
  list(var = "log_commit_excl_fund",         raw = "commit_excl_fund",       label = "Excl. fund-flagged"),
  list(var = "log_commit_excl_all",          raw = "commit_excl_all",        label = "Excl. all three"),
  list(var = "log_commit_multi_excl_fund",   raw = "commit_multi_excl_fund", label = "Multilateral cell, excl. fund-flagged")
)

exclusion_results <- vector("list", length(exclusion_specs))
for (i in seq_along(exclusion_specs)) {
  sp <- exclusion_specs[[i]]
  log_msg("  Estimating: %s (%s)", sp$label, sp$var)
  res <- estimate_cs_outcome(did_panel_126_df, sp$var)
  res$label <- sp$label
  res$mean_pre <- pretreat_mean(did_panel_126_df, sp$raw)
  exclusion_results[[i]] <- res
  log_msg("    ATT=%.4f SE=%.4f t=%.3f pretrend p=%.3f (df=%d) did-Wpval=%.3f mean_pre=%.2f",
          res$att, res$se, res$t, res$pretrend_pval, res$pretrend_df, res$did_wpval, res$mean_pre)
}
names(exclusion_results) <- sapply(exclusion_specs, `[[`, "label")

## --- Objective 3b: regex-sensitivity diagnostic -- old (pre-fix) vs.\ new --
## title/purpose flags for the three affected exclusion outcomes only. "Excl.
## fund-flagged" and "Multilateral cell, excl. fund-flagged" do not depend on
## the title/purpose flags, so they are identical before and after by
## construction and are not re-estimated here.
log_msg("Section 11c-ii: regex-sensitivity diagnostic (old vs.\\ new title/purpose flags)")

exclusion_old_specs <- list(
  list(var = "log_commit_excl_title_old",   raw = "commit_excl_title_old",   label = "Excl. title-flagged"),
  list(var = "log_commit_excl_purpose_old", raw = "commit_excl_purpose_old", label = "Excl. purpose-flagged"),
  list(var = "log_commit_excl_all_old",     raw = "commit_excl_all_old",     label = "Excl. all three")
)

exclusion_old_results <- vector("list", length(exclusion_old_specs))
for (i in seq_along(exclusion_old_specs)) {
  sp <- exclusion_old_specs[[i]]
  log_msg("  Estimating [old regex]: %s (%s)", sp$label, sp$var)
  res <- estimate_cs_outcome(did_panel_126_df, sp$var)
  res$label <- sp$label
  res$mean_pre <- pretreat_mean(did_panel_126_df, sp$raw)
  exclusion_old_results[[i]] <- res
  log_msg("    [old] ATT=%.4f SE=%.4f t=%.3f pretrend p=%.3f (df=%d) did-Wpval=%.3f mean_pre=%.2f",
          res$att, res$se, res$t, res$pretrend_pval, res$pretrend_df, res$did_wpval, res$mean_pre)
}
names(exclusion_old_results) <- sapply(exclusion_old_specs, `[[`, "label")

for (lbl in names(exclusion_old_results)) {
  old_r <- exclusion_old_results[[lbl]]
  new_r <- exclusion_results[[lbl]]
  log_msg("  Before/after [%s]: old ATT=%.4f (SE=%.4f) -> new ATT=%.4f (SE=%.4f), delta ATT=%+.4f",
          lbl, old_r$att, old_r$se, new_r$att, new_r$se, new_r$att - old_r$att)
}

log_msg("Section 11: all estimations complete in %.1f min",
        as.numeric(Sys.time() - t_est_start, units = "mins"))

# ==============================================================================
# SECTION 12. BUILD TABLES 1-3
# ==============================================================================

log_msg("Section 12: writing tables 1-3")

row_labels_est <- c(
  "ATT", "SE", "$t$-statistic", "Mean (pre-treat., USD M)",
  "Share of pre-treat.\\ marked commitments (\\%)",
  "Observations", "Countries",
  "Pre-trend $\\chi^2$", "Pre-trend $p$",
  "\\texttt{did} pre-test $p$"
)

fmt_col_margin <- function(res) c(
  paste0(sprintf("%.4f", res$att), fmt_stars(res$t)),
  sprintf("(%.4f)", res$se),
  sprintf("%.3f", res$t),
  if (is.na(res$mean_pre)) "---" else sprintf("%.2f", res$mean_pre),
  if (is.na(res$share_pre)) "---" else sprintf("%.1f", res$share_pre),
  format(res$n_obs, big.mark = ","),
  as.character(res$n_country),
  sprintf("%.3f", res$pretrend_stat),
  sprintf("%.3f", res$pretrend_pval),
  sprintf("%.3f", res$did_wpval)
)

margins_tex <- build_wide_tex_table(
  row_labels = row_labels_est,
  col_labels = sapply(margin_specs, `[[`, "label"),
  col_data   = lapply(margin_results, fmt_col_margin),
  tex_label  = "tab:remarking_margins"
)
write_tex_float(
  out_path      = here("output", "tables", "remarking", "att_remarking_margins.tex"),
  caption_title = "Margin decomposition of adaptation-marked commitments: new vs.\\ continuing vs.\\ re-marked activities",
  label         = "tab:remarking_margins",
  tabular_lines = margins_tex,
  notes_text    = paste0(
    "CS\\,(2021) DR, never-treated controls, headline specification; cohorts $<5$ treated ",
    "dropped (", n_countries_126, "-country sample). SE: multiplier-bootstrap (999 reps, ",
    "seed 1242; \\texttt{did} 2.5.0; CRS Apr.\\ 2026). ",
    PRETREND_NOTE(-5L, -2L, margin_results[[1L]]$pretrend_df), ". ",
    "Category definitions and singleton shares: see Section~\\ref{sec:remarking} and ",
    "Appendix~\\ref{app:remarking}. ``New'' ATT is an upper bound on genuinely new activity. ",
    "* $p<0.10$, ** $p<0.05$, *** $p<0.01$"
  ),
  source_text = "OECD CRS activity-level microdata (Rio adaptation marker 1 or 2); UNFCCC NAP Central"
)

## --- Table 1b (diagnostic): four-bucket margin, unlinked isolated ----------
combined4_labels <- c(
  "Headline (all marked)", "New (linked)", "Unlinked (singleton, always new)",
  "Continuing, already marked", "Re-marked (previously unmarked)"
)
combined4_results <- list(
  margin_results[["Headline (all marked)"]],
  margin4_results[["New (linked)"]],
  margin4_results[["Unlinked (singleton, always new)"]],
  margin_results[["Continuing, already marked"]],
  margin_results[["Re-marked (previously unmarked)"]]
)

margins4_tex <- build_wide_tex_table(
  row_labels = row_labels_est,
  col_labels = combined4_labels,
  col_data   = lapply(combined4_results, fmt_col_margin),
  tex_label  = "tab:remarking_margins_unlinked"
)
write_tex_float(
  out_path      = here("output", "tables", "remarking", "att_remarking_margins_unlinked.tex"),
  caption_title = "Margin decomposition with unlinked (singleton) activities isolated as a fourth bucket",
  label         = "tab:remarking_margins_unlinked",
  tabular_lines = margins4_tex,
  notes_text    = paste0(
    "Same sample, controls, estimator as Table~\\ref{tab:remarking_margins}. Splits its ",
    "``new'' bucket into ``new (linked)'' (resolvable cross-year identifier, ", chosen_id,
    ") and ``unlinked'' (singleton, new by construction). ``Continuing''/``re-marked'' are ",
    "identical to Table~\\ref{tab:remarking_margins}. Singleton shares: see ",
    "Appendix~\\ref{app:remarking}. * $p<0.10$, ** $p<0.05$, *** $p<0.01$"
  ),
  source_text = "OECD CRS activity-level microdata (Rio adaptation marker 1 or 2); UNFCCC NAP Central"
)

fmt_col_counts <- function(res) c(
  paste0(sprintf("%.4f", res$att), fmt_stars(res$t)),
  sprintf("(%.4f)", res$se),
  sprintf("%.3f", res$t),
  if (is.na(res$mean_pre)) "---" else sprintf("%.2f", res$mean_pre),
  format(res$n_obs, big.mark = ","),
  as.character(res$n_country),
  sprintf("%.3f", res$pretrend_stat),
  sprintf("%.3f", res$pretrend_pval),
  sprintf("%.3f", res$did_wpval)
)
row_labels_counts <- c(
  "ATT", "SE", "$t$-statistic", "Mean (pre-treat.)",
  "Observations", "Countries",
  "Pre-trend $\\chi^2$", "Pre-trend $p$",
  "\\texttt{did} pre-test $p$"
)

counts_tex <- build_wide_tex_table(
  row_labels = row_labels_counts,
  col_labels = sapply(counts_specs, `[[`, "label"),
  col_data   = lapply(counts_results, fmt_col_counts),
  tex_label  = "tab:remarking_counts"
)
write_tex_float(
  out_path      = here("output", "tables", "remarking", "att_remarking_counts.tex"),
  caption_title = "Counts vs.\\ values: extensive and intensive margins of adaptation marking",
  label         = "tab:remarking_counts",
  tabular_lines = counts_tex,
  notes_text    = paste0(
    "CS\\,(2021) DR, same specification as Table~\\ref{tab:remarking_margins}. SE: ",
    "multiplier-bootstrap (999 reps, seed 1242; \\texttt{did} 2.5.0; CRS Apr.\\ 2026). ",
    "``N marked activities'' and ``mean commitment/marked activity'' are log(1+$x$); ",
    "``share of activities marked'' is the \\% of the recipient's CRS activities (all ",
    "purposes) carrying an adaptation marker, untransformed (pp). ",
    "* $p<0.10$, ** $p<0.05$, *** $p<0.01$"
  ),
  source_text = "OECD CRS activity-level microdata (all purposes for the denominator; Rio adaptation marker 1 or 2 for the numerator)"
)

exclusions_tex <- build_wide_tex_table(
  row_labels = row_labels_est[row_labels_est != "Share of pre-treat.\\ marked commitments (\\%)"],
  col_labels = sapply(exclusion_specs, `[[`, "label"),
  col_data   = lapply(exclusion_results, function(res) c(
    paste0(sprintf("%.4f", res$att), fmt_stars(res$t)),
    sprintf("(%.4f)", res$se),
    sprintf("%.3f", res$t),
    if (is.na(res$mean_pre)) "---" else sprintf("%.2f", res$mean_pre),
    format(res$n_obs, big.mark = ","),
    as.character(res$n_country),
    sprintf("%.3f", res$pretrend_stat),
    sprintf("%.3f", res$pretrend_pval),
    sprintf("%.3f", res$did_wpval)
  )),
  tex_label = "tab:remarking_exclusions"
)
write_tex_float(
  out_path      = here("output", "tables", "remarking", "att_remarking_exclusions.tex"),
  caption_title = "Re-estimating the headline effect excluding NAP-support, readiness, and climate-fund activity",
  label         = "tab:remarking_exclusions",
  tabular_lines = exclusions_tex,
  notes_text    = paste0(
    "CS\\,(2021) DR, same specification as Table~\\ref{tab:remarking_margins}. SE: ",
    "multiplier-bootstrap (999 reps, seed 1242; \\texttt{did} 2.5.0; CRS Apr.\\ 2026). ",
    "Title/purpose/fund-flag definitions: see notes to ",
    "Table~\\ref{tab:remarking_flag_shares_by_year}. ``Multilateral cell, excl.\\ ",
    "fund-flagged'' recomputes log\\_commits\\_multi dropping fund-flagged activities. ",
    "Flag coverage: Tables~\\ref{tab:remarking_flag_shares_by_year}--",
    "\\ref{tab:remarking_flag_shares_by_group}. * $p<0.10$, ** $p<0.05$, *** $p<0.01$"
  ),
  source_text = "OECD CRS activity-level microdata (Rio adaptation marker 1 or 2); UNFCCC NAP Central"
)

## --- Table 3b (diagnostic): old vs.\ new title/purpose regex, before/after -
exclusions_compare_labels <- c(
  "Excl. title-flagged [old regex]", "Excl. title-flagged [new regex]",
  "Excl. purpose-flagged [old codes]", "Excl. purpose-flagged [new codes]",
  "Excl. all three [old]", "Excl. all three [new]"
)
exclusions_compare_results <- list(
  exclusion_old_results[["Excl. title-flagged"]],   exclusion_results[["Excl. title-flagged"]],
  exclusion_old_results[["Excl. purpose-flagged"]], exclusion_results[["Excl. purpose-flagged"]],
  exclusion_old_results[["Excl. all three"]],       exclusion_results[["Excl. all three"]]
)

exclusions_compare_tex <- build_wide_tex_table(
  row_labels = row_labels_est[row_labels_est != "Share of pre-treat.\\ marked commitments (\\%)"],
  col_labels = exclusions_compare_labels,
  col_data   = lapply(exclusions_compare_results, function(res) c(
    paste0(sprintf("%.4f", res$att), fmt_stars(res$t)),
    sprintf("(%.4f)", res$se),
    sprintf("%.3f", res$t),
    if (is.na(res$mean_pre)) "---" else sprintf("%.2f", res$mean_pre),
    format(res$n_obs, big.mark = ","),
    as.character(res$n_country),
    sprintf("%.3f", res$pretrend_stat),
    sprintf("%.3f", res$pretrend_pval),
    sprintf("%.3f", res$did_wpval)
  )),
  tex_label = "tab:remarking_exclusions_regex_comparison"
)
write_tex_float(
  out_path      = here("output", "tables", "remarking", "att_remarking_exclusions_regex_comparison.tex"),
  caption_title = "Sensitivity of the exclusion ATTs to the title/purpose flag regex fix",
  label         = "tab:remarking_exclusions_regex_comparison",
  tabular_lines = exclusions_compare_tex,
  notes_text    = paste0(
    "CS\\,(2021) DR, same specification/sample as Table~\\ref{tab:remarking_exclusions}. ",
    "SE: multiplier-bootstrap (999 reps, seed 1242; \\texttt{did} 2.5.0; CRS Apr.\\ 2026). ",
    "``[old]'' = pre-fix broad flags: title matches ``national adaptation plan''/``NAP''/",
    "``readiness''/``NAP-''/``adaptation plan''/``planning'' (bare); purpose ",
    "$\\in\\{41010,15110,43010\\}$. ``[new]'' reproduces Table~\\ref{tab:remarking_exclusions}. ",
    "Fund-flagged cells omitted (unchanged by this fix). ",
    "* $p<0.10$, ** $p<0.05$, *** $p<0.01$"
  ),
  source_text = "OECD CRS activity-level microdata (Rio adaptation marker 1 or 2); UNFCCC NAP Central"
)

# ==============================================================================
# SECTION 13. OBJECTIVE 4 -- SECTOR COMPOSITION (treated vs never-treated)
# ==============================================================================

log_msg("Section 13: sector composition around adoption")

sector_lookup <- c(
  "11" = "Education", "12" = "Health", "13" = "Population \\& reproductive health",
  "14" = "Water supply \\& sanitation", "15" = "Government \\& civil society",
  "16" = "Other social infrastructure", "21" = "Transport \\& storage",
  "22" = "Communications", "23" = "Energy", "24" = "Banking \\& financial services",
  "25" = "Business \\& other services", "31" = "Agriculture, forestry, fishing",
  "32" = "Industry, mining, construction", "33" = "Trade policy \\& regulations",
  "41" = "General environment protection", "43" = "Other multisector",
  "51" = "General budget support", "52" = "Food aid / security",
  "53" = "Other commodity assistance", "60" = "Action relating to debt",
  "72" = "Emergency response", "73" = "Reconstruction relief \\& rehabilitation",
  "74" = "Disaster prevention \\& preparedness", "91" = "Administrative costs of donors",
  "92" = "Support to NGOs", "93" = "Refugees in donor countries",
  "99" = "Unallocated / unspecified"
)

activities[, sector2 := substr(formatC(PurposeCode, width = 3, flag = "0", format = "d"), 1, 2)]
activities[, sector_label := ifelse(sector2 %in% names(sector_lookup),
                                    sector_lookup[sector2], paste0("Sector ", sector2))]

adapt_rows <- activities[ClimateAdaptation %in% c(1L, 2L)]

cohort_map <- did_panel_126 %>% distinct(recipient_name, cohort_year)
treated_cohorts <- cohort_map %>% filter(cohort_year >= 2021)

# Per-cohort calendar windows: pre = [c-3, c-1], post = [c, c+2]
window_tbl <- treated_cohorts %>%
  mutate(
    pre_start  = cohort_year - 3, pre_end  = cohort_year - 1,
    post_start = cohort_year,     post_end = cohort_year + 2
  )

pre_years_union  <- sort(unique(unlist(Map(seq, window_tbl$pre_start,  window_tbl$pre_end))))
post_years_union <- sort(unique(unlist(Map(seq, window_tbl$post_start, window_tbl$post_end))))
post_years_available <- post_years_union[post_years_union <= last_year]
n_post_truncated <- sum(post_years_union > last_year)
if (n_post_truncated > 0L) {
  log_msg("Sector composition: %d post-window calendar year(s) beyond %d are unavailable (right-censoring for 2023/2024 cohorts)",
          n_post_truncated, last_year)
}

treated_dt <- as.data.table(window_tbl)[, .(recipient_name, pre_start, pre_end, post_start, post_end)]

treated_sector_pre <- adapt_rows[treated_dt, on = .(RecipientName = recipient_name),
                                  nomatch = 0][Year >= pre_start & Year <= pre_end]
treated_sector_post <- adapt_rows[treated_dt, on = .(RecipientName = recipient_name),
                                   nomatch = 0][Year >= post_start & Year <= post_end]

treated_pre_comp  <- treated_sector_pre[,  .(commit = sum(USD_Commitment_Defl, na.rm = TRUE)), by = sector_label]
treated_post_comp <- treated_sector_post[, .(commit = sum(USD_Commitment_Defl, na.rm = TRUE)), by = sector_label]

never_rows <- adapt_rows[RecipientName %in% never_treated_126]
never_pre_comp  <- never_rows[Year %in% pre_years_union,  .(commit = sum(USD_Commitment_Defl, na.rm = TRUE)), by = sector_label]
never_post_comp <- never_rows[Year %in% post_years_union, .(commit = sum(USD_Commitment_Defl, na.rm = TRUE)), by = sector_label]

make_share <- function(dt) {
  dt <- copy(dt)
  dt[, share := 100 * commit / sum(commit)]
  dt
}
treated_pre_s  <- make_share(treated_pre_comp)
treated_post_s <- make_share(treated_post_comp)
never_pre_s    <- make_share(never_pre_comp)
never_post_s   <- make_share(never_post_comp)

all_sectors <- unique(c(treated_pre_s$sector_label, treated_post_s$sector_label,
                        never_pre_s$sector_label, never_post_s$sector_label))

wide_shares <- data.table(sector_label = all_sectors)
wide_shares <- merge(wide_shares, treated_pre_s[,  .(sector_label, treated_pre  = share)], by = "sector_label", all.x = TRUE)
wide_shares <- merge(wide_shares, treated_post_s[, .(sector_label, treated_post = share)], by = "sector_label", all.x = TRUE)
wide_shares <- merge(wide_shares, never_pre_s[,    .(sector_label, never_pre    = share)], by = "sector_label", all.x = TRUE)
wide_shares <- merge(wide_shares, never_post_s[,   .(sector_label, never_post   = share)], by = "sector_label", all.x = TRUE)
for (cl in c("treated_pre", "treated_post", "never_pre", "never_post")) {
  wide_shares[is.na(get(cl)), (cl) := 0]
}
wide_shares[, total_rank := treated_pre + treated_post + never_pre + never_post]
setorder(wide_shares, -total_rank)

top8 <- wide_shares[seq_len(min(8L, nrow(wide_shares)))]
other_row <- wide_shares[-seq_len(min(8L, nrow(wide_shares))),
                          .(sector_label = "Other",
                            treated_pre = sum(treated_pre), treated_post = sum(treated_post),
                            never_pre = sum(never_pre), never_post = sum(never_post),
                            total_rank = sum(total_rank))]
sector_table <- rbindlist(list(top8, other_row))

# Dissimilarity index D = 0.5 * sum(|share_post - share_pre|), computed on the
# FULL (pre-collapse) sector distribution, not the top-8+Other collapse.
dissim <- function(pre_dt, post_dt) {
  m <- merge(pre_dt[, .(sector_label, pre = share)], post_dt[, .(sector_label, post = share)],
             by = "sector_label", all = TRUE)
  m[is.na(pre), pre := 0]; m[is.na(post), post := 0]
  0.5 * sum(abs(m$post - m$pre))
}
D_treated <- dissim(treated_pre_s, treated_post_s)
D_never   <- dissim(never_pre_s, never_post_s)

log_msg("Dissimilarity index: treated pre-vs-post D=%.2f, never-treated pre-vs-post D=%.2f",
        D_treated, D_never)

sector_row_labels <- c(sector_table$sector_label,
                       "Dissimilarity index $D$ (pre vs.\\ post, this column pair)")
sector_col_data <- list(
  c(sprintf("%.1f", sector_table$treated_pre),  sprintf("%.1f", D_treated)),
  c(sprintf("%.1f", sector_table$treated_post), ""),
  c(sprintf("%.1f", sector_table$never_pre),    sprintf("%.1f", D_never)),
  c(sprintf("%.1f", sector_table$never_post),   "")
)
sector_tex <- build_wide_tex_table(
  row_labels = sector_row_labels,
  col_labels = c("Treated, pre", "Treated, post", "Never-treated, pre", "Never-treated, post"),
  col_data   = sector_col_data,
  tex_label  = "tab:remarking_sectors"
)
write_tex_float(
  out_path      = here("output", "tables", "remarking", "tab_remarking_sectors.tex"),
  caption_title = "Sector composition of adaptation-marked commitments, before vs.\\ after NAP adoption",
  label         = "tab:remarking_sectors",
  tabular_lines = sector_tex,
  notes_text    = paste0(
    "Sector = leading 2 digits of the CRS purpose code. Treated: 40 main-sample adopters ",
    "(cohort $\\geq$ 2021); pre = 3 years before adoption, post = adoption year + 2 ",
    "(right-censored for 2023--24 cohorts). Never-treated: pre/post pool commitments over ",
    "the union of treated cohorts' pre-/post-window years (can overlap). Shares are \\% of ",
    "column total; top 8 sectors shown, remainder in ``Other''. ",
    "$D = 0.5\\sum_i|\\text{share}_{i,\\text{post}}-\\text{share}_{i,\\text{pre}}|$, full ",
    "(uncollapsed) distribution"
  ),
  source_text = "OECD CRS activity-level microdata (Rio adaptation marker 1 or 2)"
)

## --- Figure: stacked bar of sector composition ------------------------------
plot_df <- sector_table %>%
  select(sector_label, treated_pre, treated_post, never_pre, never_post) %>%
  pivot_longer(cols = -sector_label, names_to = "cell", values_to = "share") %>%
  mutate(
    group_lab = case_when(
      cell == "treated_pre"  ~ "Treated -- pre",
      cell == "treated_post" ~ "Treated -- post",
      cell == "never_pre"    ~ "Never-treated -- pre",
      cell == "never_post"   ~ "Never-treated -- post"
    ),
    group_lab = factor(group_lab, levels = c("Treated -- pre", "Treated -- post",
                                             "Never-treated -- pre", "Never-treated -- post")),
    sector_label = gsub("\\\\&", "&", sector_label),
    sector_label = factor(sector_label, levels = rev(unique(sector_label)))
  )

p_sectors <- ggplot(plot_df, aes(x = group_lab, y = share, fill = sector_label)) +
  geom_col(position = "stack", colour = "white", linewidth = 0.15) +
  scale_fill_viridis_d(option = "D") +
  # No title/subtitle/caption inside the plot
  labs(title = NULL, subtitle = NULL, caption = NULL,
       x = NULL, y = "Share of adaptation-marked commitments (%)", fill = NULL) +
  theme_minimal(base_family = "serif", base_size = 12) +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 9),
    axis.text.x      = element_text(angle = 20, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 3))

ggsave(here("output", "figures", "remarking", "fig_remarking_sectors.png"),
       p_sectors, width = 12, height = 7, dpi = 300)
log_msg("Saved: output/figures/remarking/fig_remarking_sectors.png")

# ==============================================================================
# SECTION 14. SAVE DIAGNOSTIC/RESULT OBJECTS AND COPY EXHIBITS TO paper/
# ==============================================================================

log_msg("Section 14: saving diagnostic objects and copying exhibits to paper/")

saveRDS(
  list(
    id_diagnostic          = id_diag,
    reconciliation         = recon_merged,
    reconciliation_max_diff = max_diff,
    margin_results         = margin_results,
    margin4_results        = margin4_results,
    counts_results         = counts_results,
    exclusion_results      = exclusion_results,
    exclusion_old_results  = exclusion_old_results,
    flag_shares_year       = flag_shares_year,
    flag_shares_group      = flag_shares_group,
    sector_table           = sector_table,
    D_treated              = D_treated,
    D_never                = D_never,
    chosen_id              = chosen_id,
    singleton_share_rows_overall   = singleton_share_rows_overall,
    singleton_share_commit_overall = singleton_share_commit_overall,
    singleton_share_rows_post      = singleton_share_rows_post,
    singleton_share_commit_post    = singleton_share_commit_post
  ),
  here("output", "tables", "remarking", "remarking_results.rds")
)

if (has_paper) {
  tex_files <- list.files(here("output", "tables", "remarking"), pattern = "\\.tex$", full.names = TRUE)
  file.copy(tex_files, here("paper", "Tables", "remarking"), overwrite = TRUE)

  png_files <- list.files(here("output", "figures", "remarking"), pattern = "\\.png$", full.names = TRUE)
  file.copy(png_files, here("paper", "Figures", "remarking"), overwrite = TRUE)

  log_msg("Copied %d table(s) to paper/Tables/remarking/ and %d figure(s) to paper/Figures/remarking/",
          length(tex_files), length(png_files))
} else log_msg("paper/ not found -- exhibits are left in output/ only")

log_msg("=== 09_remarking_decomposition.R: COMPLETE in %.1f min ===",
        as.numeric(Sys.time() - t_script_start, units = "mins"))
