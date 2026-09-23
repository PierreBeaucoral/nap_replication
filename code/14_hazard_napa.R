##############################################################################
# 14_hazard_napa.R
# Author : Pierre Beaucoral
# Date   : 2026-09-15
# Design:  (1) lagged emergency-response AID is a bad control, so hazard-
#          control specs use EM-DAT only; (2) EM-DAT is an event-level
#          extract, aggregated to country-year here (see §2).
# Purpose: hazard-realisation controls (EM-DAT) + two NAP-timing
#          orthogonality checks (EM-DAT hazard; CRS emergency-response aid),
#          plus the prior-NAPA split (which doubles as a placebo) and the
#          NAPA-adoption falsification check.
#
# Section 1 (panel construction) is copied VERBATIM from
# 03_main_results.R Section 1, so the estimation sample here is identical to the headline spec's
# sample. 03_main_results.R itself is NOT modified or sourced by this script.
# Section 6's aligned-influence-function difference test is copied VERBATIM
# from 05_heterogeneity.R:428-475 (id-matched version); 05 itself is NOT
# modified or sourced by this script.
#
# Inputs (reads only — never touches raw CRS):
#   data/processed/simple_panel_wgi.csv           (03's estimation panel)
#   data/processed/emergency_response_panel.csv   (01 §23; CRS aid proxy,
#                                                   timing check only)
#   data/raw/emdat/emdat.csv                       (EM-DAT event-level
#                                                   extract, v2026-09-11)
#   data/raw/napa/napa_list_unfccc.csv             (NAPA list verified by the authors)
#
# Outputs:
#   data/processed/napa_list.csv
#   data/processed/emdat_panel.csv
#   output/tables/hazard/att_hazard_controls.tex          (EM-DAT only;
#                                                           not written if
#                                                           EM-DAT absent)
#   output/tables/hazard/nap_timing_vs_humanitarian_aid.tex
#   output/tables/hazard/nap_timing_vs_emdat_hazard.tex   (not written if
#                                                           EM-DAT absent)
#   output/tables/napa/att_prior_napa_split.tex
#   output/tables/napa/att_napa_falsification.tex
##############################################################################

##############################################################################
# §0. PACKAGES
##############################################################################

library(data.table)
library(dplyr)
library(tibble)
library(broom)
library(here)
library(did)
library(fixest)
library(countrycode)
library(stringr)

# Version guard identical to 03_main_results.R: SEs on unbalanced panels
# changed in did 2.5.0 (renv.lock pins it).
if (utils::packageVersion("did") < "2.5.0") {
  stop(sprintf(paste0(
    "did %s is installed but this script's results require did >= 2.5.0.\n",
    "  Fix: run renv::restore() from the project root (renv.lock pins did 2.5.0)."),
    utils::packageVersion("did")))
}

# Single set.seed() at top (global). A local set.seed(1242) precedes
# EVERY att_gt()/bootstrap call below, matching the seed rule documented
# in 03/05 ("global seed -- local set.seed(1242) calls follow each estimator").
set.seed(20240601)

# Kept short (no embedded parenthetical) so paste0("... (", CRS_VINTAGE, ")") call sites never produce nested parens.
CRS_VINTAGE  <- "April 2026, DATA\\_AVAILABILITY.md"
DID_VERSION  <- as.character(utils::packageVersion("did"))
EMDAT_VINTAGE <- "EM-DAT v2026-09-11 (CRED/UCLouvain public custom request, downloaded 2026-09-15)"

##############################################################################
# §0b. OUTPUT DIRECTORIES
##############################################################################

dir.create(here("output", "tables", "hazard"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "tables", "napa"),   recursive = TRUE, showWarnings = FALSE)
dir.create(here("data", "processed"),          recursive = TRUE, showWarnings = FALSE)

# Superseded filenames from earlier iterations of this script (renamed when
# the tables were relabelled).
# Removed UNCONDITIONALLY at the top of every run, regardless of which
# branches execute below, so a stale file never survives a rerun.
stale_output_files <- c(
  here("output", "tables", "hazard", "nap_timing_orthogonality.tex"),
  here("output", "tables", "napa",   "att_napa_placebo.tex")
)
for (sf in stale_output_files) {
  if (file.exists(sf)) {
    file.remove(sf)
    message("Deleted stale output file (superseded name): ", sf)
  }
}

##############################################################################
# §0c. HELPERS
# write_tex_float() and esc_tex() are duplicated verbatim from
# 01_prepare_data.R / 03_main_results.R -- each pipeline script is
# self-contained (not sourced) per the project's existing convention.
# esc_tex() is needed here: hand-labelled
# text and recipient_name never contain LaTeX-special characters in this
# panel, but the DYNAMICALLY CAPTURED did/fixest warning and error texts
# embedded in several table notes can --
# e.g. \texttt{fixest::feglm()}'s error message quotes the R call verbatim,
# including underscored variable names -- so those specific strings are
# escaped before insertion; everything else is left as-is.
##############################################################################

esc_tex <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([&%$#_{}])", "\\\\\\1", x)
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
  message("Wrote: ", out_path)
  invisible(out_path)
}

stars_of <- function(p) {
  case_when(
    is.na(p)   ~ "",
    p < 0.01   ~ "***",
    p < 0.05   ~ "**",
    p < 0.10   ~ "*",
    TRUE       ~ ""
  )
}
fmt4 <- function(x) ifelse(is.na(x), "--", sprintf("%.4f", x))
fmt3 <- function(x) ifelse(is.na(x), "--", sprintf("%.3f", x))
fmt1 <- function(x) ifelse(is.na(x), "--", sprintf("%.1f", x))

# Runs att_gt() with BOTH error and warning capture (with a tryCatch on
# error only, package warnings like "very few observations" would be printed
# to the console and then lost). Warnings are collected, not suppressed from the
# log (still print via message()), and returned so callers can fold a count
# into the table note.
run_att_gt_captured <- function(...) {
  warns <- character(0)
  gt_obj <- withCallingHandlers(
    tryCatch(did::att_gt(...), error = function(e) {
      message("  att_gt failed: ", conditionMessage(e)); NULL
    }),
    warning = function(w) {
      # Sanctioned exception to the project's `<<-` rule: the only way to
      # accumulate conditions out of a withCallingHandlers() handler; the
      # target is a local of the enclosing function, not a global.
      warns <<- c(warns, conditionMessage(w))
      message("  att_gt warning: ", conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(gt_obj = gt_obj, warnings = warns)
}

##############################################################################
# §0d. Aligned unit-level influence-function difference test.
#
# Copied VERBATIM from 05_heterogeneity.R:428-475 (the id-matched version;
# 05 itself is NOT modified or sourced by this script). That function's own
# header notes it is ported (logic and scaling convention) from
# 07_principal_and_share.R::compute_att_difference. The scaling is sqrt(mean((IFa - IFb)^2) / n), the UNCENTRED
# mean-square convention did itself uses (verified in 07 against
# overall.se), not sqrt(var(IFa - IFb)) / n and not sqrt(var(IFa - IFb)).
# 05 additionally added id-matching: when the two fits' influence-function vectors don't
# align in length, restrict to the intersection of the two fits' unit ids
# (sorted-ascending country_id, matching did's own internal unit ordering)
# rather than falling straight to the cruder independence approximation.
#
# @param att_a,att_b point estimates
# @param IFa,IFb unit-level influence-function vectors from
#   aggte(type = "simple")$inf.function$simple.att of the ANALYTICAL fits
# @param se_a,se_b reported standard errors (used only for the fallback)
# @param ids_a,ids_b unit identifiers the influence-function rows refer to
#   (sort(unique(data$country_id)) of the data each fit was estimated on)
# @return list(diff, se, z, pval, method, n_units)
##############################################################################

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

##############################################################################
# SECTION 1. Panel construction -- copied verbatim from
# 03_main_results.R Section 1 (03 itself is not modified or sourced).
##############################################################################

message("\n=== 14_hazard_napa.R: loading panel (copy of 03 Section 1) ===\n")

aggregated <- fread(here("data", "processed", "simple_panel_wgi.csv"))
aggregated  <- as.data.frame(aggregated)

did_panel <- aggregated
did_panel <- did_panel %>%
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

# Cohort sizes / thin-cohort rule -- identical to 03 (thin_threshold = 5L).
cohort_sizes <- did_panel %>%
  filter(cohort_year > 0) %>%
  distinct(recipient_name, cohort_year) %>%
  count(cohort_year, name = "n_treated") %>%
  arrange(cohort_year)

thin_threshold <- 5L
thin_cohorts   <- cohort_sizes$cohort_year[cohort_sizes$n_treated < thin_threshold]

did_panel_main <- did_panel_full %>% filter(!(cohort_year %in% thin_cohorts))

n_recipients_panel <- n_distinct(did_panel_full$recipient_name)

message(sprintf(
  "Main panel (cohorts >= %d treated units): %d obs, %d countries, %d treated cohorts (%s).",
  thin_threshold, nrow(did_panel_main), n_distinct(did_panel_main$recipient_name),
  n_distinct(did_panel_main$cohort_year[did_panel_main$cohort_year > 0]),
  paste(sort(unique(did_panel_main$cohort_year[did_panel_main$cohort_year > 0])), collapse = ", ")
))
message(sprintf("Full panel (all cohorts, for orthogonality/placebo checks): %d countries.",
                n_recipients_panel))

##############################################################################
# SECTION 2. Hazard data preparation
#
# (i) EM-DAT (data/raw/emdat/emdat.csv) -- the ONLY regressor used in any
#     hazard-CONTROL specification (§3). Design decision: an aid-flow
#     proxy from the same CRS/donors/recipients as the outcome is a bad
#     control and must not enter xformla or a residualisation regression.
#     EM-DAT is event-level (one row per disaster, EM-DAT v2026-09-11,
#     10,896 natural-disaster records, 2000-2026, columns DisNo./ISO/Country/
#     Region/Disaster Group ("Natural" for all rows)/Disaster Subgroup
#     (Hydrological/Meteorological/Climatological/Geophysical/Biological/
#     Extra-terrestrial)/Disaster Type/Disaster Subtype/Start Year/Start
#     Month/End Year/Total Deaths/No. Injured/No. Affected/No. Homeless/
#     Total Affected/Total Damage ('000 US$)/Total Damage, Adjusted
#     ('000 US$)/Magnitude/Magnitude Scale). This section aggregates it to a
#     recipient-year panel (2008-2024, zero-filled) for TWO hazard sets:
#       ALL natural            -- every Disaster Group = Natural event
#       CLIMATE-RELATED        -- Hydrological + Meteorological +
#                                  Climatological subgroups only (excludes
#                                  Geophysical, Biological, Extra-terrestrial
#                                  -- earthquakes, epidemics, etc. are not
#                                  plausibly related to the adaptation-finance
#                                  channel this paper studies)
#     assigned by Start Year, with n_events, log(1+Total Affected), and
#     log(1+Total Damage, Adjusted) for each set. Written to
#     data/processed/emdat_panel.csv.
# (ii) CRS emergency-response-AID proxy (01 §23) -- used ONLY by §4's
#     separate "NAP timing vs. humanitarian aid" orthogonality check, never
#     as a hazard-control covariate (see the design decision above).
##############################################################################

message("\n=== SECTION 2: Hazard data preparation ===\n")

recipients_est <- unique(did_panel_full$recipient_iso)
recipients_est <- recipients_est[!is.na(recipients_est)]

emdat_path      <- here("data", "raw", "emdat", "emdat.csv")
emdat_available <- file.exists(emdat_path)
emdat_panel     <- NULL

# EM-DAT's terms of use do not allow redistribution, so the replication package
# does not include it. Without it the two EM-DAT exhibits below cannot be
# regenerated; the rest of this script (humanitarian-aid timing check, NAPA
# tables) runs normally. Any existing copies of the two exhibits are left
# untouched, so they are the authors' versions, not regenerated ones.
if (!emdat_available) {
  message(paste0(
    "\n*** EM-DAT NOT FOUND (data/raw/emdat/emdat.csv) ***\n",
    "  The following exhibits CANNOT be regenerated in this run:\n",
    "    output/tables/hazard/att_hazard_controls.tex        (hazard-control specifications)\n",
    "    output/tables/hazard/nap_timing_vs_emdat_hazard.tex (NAP timing vs. lagged EM-DAT hazard)\n",
    "  Any copies already in output/tables/hazard/ are the authors' versions and are left as is.\n",
    "  To regenerate them, request the EM-DAT extract described in data/raw/emdat/README.md\n",
    "  (https://public.emdat.be, free registration, non-commercial use) and save it as\n",
    "  data/raw/emdat/emdat.csv, then rerun this stage.\n"))
}

if (emdat_available) {
  message("EM-DAT extract found: data/raw/emdat/emdat.csv (", EMDAT_VINTAGE, ").")
  emdat_raw <- fread(emdat_path, encoding = "UTF-8")
  req_cols <- c("ISO", "Disaster Group", "Disaster Subgroup", "Start Year",
               "Total Affected", "Total Damage, Adjusted ('000 US$)")
  if (!all(req_cols %in% names(emdat_raw))) {
    stop("data/raw/emdat/emdat.csv is missing expected column(s): ",
        paste(setdiff(req_cols, names(emdat_raw)), collapse = ", "),
        ". See the \u00a72 header comment for the expected EM-DAT layout.")
  }
  data.table::setnames(
    emdat_raw,
    old = req_cols,
    new = c("iso3", "disaster_group", "disaster_subgroup", "start_year",
           "total_affected", "total_damage_adj")
  )
  emdat_events <- emdat_raw[disaster_group == "Natural",
                            .(iso3, disaster_subgroup, start_year,
                              total_affected = suppressWarnings(as.numeric(total_affected)),
                              total_damage_adj = suppressWarnings(as.numeric(total_damage_adj)))]

  climate_subgroups <- c("Hydrological", "Meteorological", "Climatological")

  # --- Join diagnostic: recipients with
  # NO matching entity anywhere in EM-DAT's own ISO universe (not just no
  # events in our window -- EM-DAT never tracks them at all, e.g. Kosovo's
  # events are booked under Serbia). For these, "0 events" would silently
  # misrepresent a join failure as a confirmed absence of disasters, so their
  # ENTIRE hazard series (n_events, affected, damage, and everything derived
  # from them) is set to NA, not 0.
  emdat_iso_universe <- unique(emdat_raw$iso3)
  recipients_no_emdat_entity <- setdiff(recipients_est, emdat_iso_universe)
  message(sprintf(
    paste0("EM-DAT join diagnostic: %d / %d panel recipients have NO matching EM-DAT entity ",
          "(%s) -- their hazard series is set to NA, not 0."),
    length(recipients_no_emdat_entity), length(recipients_est),
    if (length(recipients_no_emdat_entity) == 0) "none" else
      paste(recipients_no_emdat_entity, collapse = ", ")
  ))

  # --- Missing vs. zero: a
  # recipient-year with n_events > 0 but where NONE of those events reported
  # a Total Affected / Total Damage value must resolve to NA on that
  # variable, not 0 (sum(..., na.rm = TRUE) over an all-NA vector silently
  # returns 0, which is indistinguishable from "confirmed zero affected/
  # damage" without tracking n_events_with_value separately).
  build_emdat_set <- function(events, subgroup_filter = NULL) {
    ev <- events
    if (!is.null(subgroup_filter)) ev <- ev[disaster_subgroup %in% subgroup_filter]
    agg <- ev[, .(n_events           = .N,
                  n_events_w_affected = sum(!is.na(total_affected)),
                  n_events_w_damage   = sum(!is.na(total_damage_adj)),
                  total_affected      = sum(total_affected,   na.rm = TRUE),
                  total_damage_adj    = sum(total_damage_adj, na.rm = TRUE)),
              by = .(iso3, start_year)]
    data.table::setnames(agg, c("iso3", "start_year"), c("recipient_iso", "year"))
    grid <- expand.grid(recipient_iso = recipients_est, year = 2008:2024,
                        stringsAsFactors = FALSE)
    out <- merge(grid, agg, by = c("recipient_iso", "year"), all.x = TRUE)
    out$n_events            <- coalesce(out$n_events, 0L)
    out$n_events_w_affected <- coalesce(out$n_events_w_affected, 0L)
    out$n_events_w_damage   <- coalesce(out$n_events_w_damage, 0L)
    # Missing/zero split: NA only when events occurred but none reported that
    # particular value; a true absence of events is a genuine 0.
    out$total_affected <- dplyr::case_when(
      out$n_events == 0             ~ 0,
      out$n_events_w_affected == 0  ~ NA_real_,
      TRUE                          ~ coalesce(out$total_affected, 0)
    )
    out$total_damage_adj <- dplyr::case_when(
      out$n_events == 0           ~ 0,
      out$n_events_w_damage == 0  ~ NA_real_,
      TRUE                        ~ coalesce(out$total_damage_adj, 0)
    )
    # No-EM-DAT-entity recipients: entire series NA, not 0 (see join diagnostic).
    no_entity_rows <- out$recipient_iso %in% recipients_no_emdat_entity
    out$n_events[no_entity_rows]            <- NA_integer_
    out$n_events_w_affected[no_entity_rows] <- NA_integer_
    out$n_events_w_damage[no_entity_rows]   <- NA_integer_
    out$total_affected[no_entity_rows]      <- NA_real_
    out$total_damage_adj[no_entity_rows]    <- NA_real_
    as.data.frame(out)
  }

  emdat_all <- build_emdat_set(emdat_events, NULL) %>%
    rename(n_events_all = n_events, affected_all = total_affected, damage_all = total_damage_adj) %>%
    select(-n_events_w_affected, -n_events_w_damage)
  emdat_climate <- build_emdat_set(emdat_events, climate_subgroups) %>%
    rename(n_events_climate = n_events, affected_climate = total_affected,
          damage_climate = total_damage_adj,
          n_events_climate_w_affected = n_events_w_affected,
          n_events_climate_w_damage   = n_events_w_damage)

  emdat_panel <- emdat_all %>%
    left_join(emdat_climate, by = c("recipient_iso", "year")) %>%
    mutate(
      log_affected_all     = log1p(affected_all),
      log_damage_all       = log1p(damage_all),
      log_affected_climate = log1p(affected_climate),
      log_damage_climate   = log1p(damage_climate)
    )
  write.csv(emdat_panel, here("data", "processed", "emdat_panel.csv"), row.names = FALSE)

  # Coverage: n_events itself is clean (0 IS a genuine zero, except for the
  # no-EM-DAT-entity recipients, which are NA and excluded via na.rm below).
  cov_all     <- mean(emdat_panel$n_events_all > 0, na.rm = TRUE)
  cov_climate <- mean(emdat_panel$n_events_climate > 0, na.rm = TRUE)
  # Missing-value audit, CLIMATE set: among recipient-years with >= 1 event,
  # what share have NO event reporting affected / damage (i.e. would have
  # been silently coded 0 under the old sum(na.rm=TRUE) logic)?
  event_years_climate <- emdat_panel %>% filter(n_events_climate > 0)
  share_missing_affected_climate <- mean(event_years_climate$n_events_climate_w_affected == 0)
  share_missing_damage_climate   <- mean(event_years_climate$n_events_climate_w_damage == 0)
  message(sprintf(
    paste0("Wrote: data/processed/emdat_panel.csv (%d recipient-year rows, %d recipients; %d ",
          "recipients NA throughout for lack of an EM-DAT entity). Coverage (share of non-NA ",
          "recipient-years with >= 1 event): ALL natural = %.1f%%; CLIMATE-RELATED = %.1f%%. ",
          "Missing-value audit (CLIMATE-RELATED event-years, n = %d): %.1f%% have NO event ",
          "reporting Total Affected (coded NA, not 0); %.1f%% have NO event reporting Total ",
          "Damage (coded NA, not 0)."),
    nrow(emdat_panel), length(recipients_est), length(recipients_no_emdat_entity),
    100 * cov_all, 100 * cov_climate, nrow(event_years_climate),
    100 * share_missing_affected_climate, 100 * share_missing_damage_climate
  ))

  # Lags at g-1 and g-2, plus a 3-year backward moving average ending at g-1,
  # built without an extra rolling-window package dependency.
  emdat_panel <- emdat_panel %>%
    arrange(recipient_iso, year) %>%
    group_by(recipient_iso) %>%
    mutate(
      log_affected_climate_lag1 = dplyr::lag(log_affected_climate, 1L),
      log_affected_climate_lag2 = dplyr::lag(log_affected_climate, 2L),
      log_damage_climate_lag1   = dplyr::lag(log_damage_climate, 1L),
      n_events_climate_lag1     = dplyr::lag(n_events_climate, 1L),
      n_events_climate_lag2     = dplyr::lag(n_events_climate, 2L),
      log_affected_all_lag1     = dplyr::lag(log_affected_all, 1L),
      # 3-year backward MA of the RAW (pre-log) affected count over
      # (t-3, t-2, t-1), then log1p'd -- "the hazard at g-1", smoothed.
      affected_climate_ma3      = (dplyr::lag(affected_climate, 1L) +
                                    dplyr::lag(affected_climate, 2L) +
                                    dplyr::lag(affected_climate, 3L)) / 3,
      log_affected_climate_ma3  = log1p(affected_climate_ma3)
    ) %>%
    ungroup()
} else {
  message("data/raw/emdat/emdat.csv NOT found -- \u00a73 (hazard-control specifications) ",
          "will be SKIPPED and no att_hazard_controls.tex will be written (message only, ",
          "no placeholder table -- see \u00a73 below). The EM-DAT branch of \u00a74's ",
          "orthogonality check is skipped for the same reason.")
}

# --- CRS emergency-response-AID proxy (timing check only, NOT a hazard
#     control -- see the design-decision note in the header). --------------
aid_proxy_raw <- fread(here("data", "processed", "emergency_response_panel.csv")) %>%
  rename(aid_value = emergency_commitments)

aid_proxy_panel <- aid_proxy_raw %>%
  arrange(recipient_iso, year) %>%
  group_by(recipient_iso) %>%
  mutate(
    aid_lag1 = dplyr::lag(aid_value, 1L),
    aid_lag2 = dplyr::lag(aid_value, 2L)
  ) %>%
  ungroup() %>%
  mutate(
    log_aid_lag1 = log1p(pmax(aid_lag1, 0, na.rm = FALSE)),
    log_aid_lag2 = log1p(pmax(aid_lag2, 0, na.rm = FALSE))
  ) %>%
  select(recipient_iso, year, aid_value, log_aid_lag1, log_aid_lag2)

did_panel_full <- did_panel_full %>% left_join(aid_proxy_panel, by = c("recipient_iso", "year"))
did_panel_main <- did_panel_main %>% left_join(aid_proxy_panel, by = c("recipient_iso", "year"))
if (!is.null(emdat_panel)) {
  emdat_join_cols <- emdat_panel %>%
    select(recipient_iso, year, log_affected_climate_lag1, log_affected_climate_lag2,
           log_damage_climate_lag1, n_events_climate_lag1, n_events_climate_lag2,
           log_affected_all_lag1, log_affected_climate_ma3, log_affected_climate,
           log_damage_climate)
  did_panel_full <- did_panel_full %>% left_join(emdat_join_cols, by = c("recipient_iso", "year"))
  did_panel_main <- did_panel_main %>% left_join(emdat_join_cols, by = c("recipient_iso", "year"))
}

##############################################################################
# SECTION 3. Hazard-realisation controls (EM-DAT only)
#
# att_gt()'s doubly-robust estimator DOES accept genuinely time-varying
# per-unit-per-period covariates in xformla: the
# regression-adjustment/propensity-score step uses each unit's covariate
# value for the specific (g,t) comparison being formed. Row (b)'s
# residualised-outcome spec is offered as an ADDITIONAL, complementary
# check -- not a workaround for a package limitation -- because it purges
# the hazard-outcome relationship from the dependent variable entirely
# (using only never-treated variation to estimate the slope), which is
# informative on its own even though xformla-conditioning already works.
#
# Headline hazard set: CLIMATE-RELATED (Hydrological+Meteorological+
# Climatological). ALL-natural is included as one robustness row.
# Regressors: log(1+Total Affected) at g-1 (primary, row a); log(1+Total
# Damage, Adjusted) at g-1 (second column, row a-damage); a 3-year moving
# average of log(1+Total Affected) ending at g-1 (row c); a residualised-
# outcome version of row (a) (row b). Every row's difference vs. the
# baseline (no hazard control) is tested with the aligned-influence-function
# helper (§6). Row (a) additionally gets a propensity-score overlap
# diagnostic and a trimmed re-estimate.
##############################################################################

message("\n=== SECTION 3: Hazard-realisation controls (EM-DAT only) ===\n")

if (!emdat_available) {

  message("SKIPPED (EM-DAT not found): output/tables/hazard/att_hazard_controls.tex ",
          "NOT regenerated -- see data/raw/emdat/README.md.")

} else {

  hazard_headline_spec <- function(yname, xformla, data) {
    set.seed(1242)  # Seed rule (reproduces the published SEs): seed immediately before the estimator call
    cap <- run_att_gt_captured(
      yname = yname, tname = "year", idname = "country_id", gname = "cohort_year",
      xformla = xformla, data = data, est_method = "dr", bstrap = TRUE, biters = 999L,
      cband = FALSE, control_group = "nevertreated", anticipation = 0,
      base_period = "universal", panel = TRUE, allow_unbalanced_panel = TRUE
    )
    if (is.null(cap$gt_obj)) return(list(gt_obj = NULL, agg_s = NULL, n_obs = nrow(data),
                                         n_treated = NA_integer_, warnings = cap$warnings,
                                         ids = integer(0)))
    agg_s <- tryCatch(aggte(cap$gt_obj, type = "simple", na.rm = TRUE), error = function(e) NULL)
    # Treated units that actually contribute to the fit: units in cohorts retained by att_gt(), counted from the
    # estimator's own data, not from the input frame.
    # Contributing = member of a cohort with at least one non-NA POST-treatment
    # ATT(g,t) in this fit (cohorts whose post cells are all NA add nothing to
    # aggte(type = "simple", na.rm = TRUE)). The input-frame count is kept
    # separately as n_treated_input and both are disclosed in the table note.
    gt <- cap$gt_obj
    ok_groups <- unique(gt$group[gt$t >= gt$group & !is.na(gt$att)])
    n_treated_fit <- tryCatch({
      dd <- gt$DIDparams$data
      n_distinct(dd$country_id[dd$cohort_year %in% ok_groups])
    }, error = function(e) NA_integer_)
    n_treated_input <- n_distinct(data$country_id[data$cohort_year > 0])
    n_treated_source <- if (is.na(n_treated_fit)) "input-frame fallback" else "fitted object"
    if (is.na(n_treated_fit)) n_treated_fit <- n_treated_input
    list(gt_obj = cap$gt_obj, agg_s = agg_s, n_obs = nrow(data),
        n_treated = n_treated_fit, n_treated_input = n_treated_input,
        n_treated_source = n_treated_source, n_ok_groups = length(ok_groups),
        warnings = cap$warnings, ids = sort(unique(data$country_id)))
  }

  did_panel_hazard_ok <- did_panel_main %>% filter(!is.na(log_affected_climate_lag1))
  n_dropped_hazard_na <- nrow(did_panel_main) - nrow(did_panel_hazard_ok)
  message(sprintf(
    "Rows dropped for missing EM-DAT one-year lag (affected): %d / %d (%.1f%%).",
    n_dropped_hazard_na, nrow(did_panel_main), 100 * n_dropped_hazard_na / nrow(did_panel_main)
  ))

  # Damage row: after the §2 missing/
  # zero fix, log_damage_climate_lag1 is genuinely NA (not log1p(0) = 0) for
  # recipient-years where events occurred but none reported a damage value.
  # This additionally restricts the damage row's sample -- reported here and
  # in the table note so the row's smaller N is not silently absorbed into
  # the shared "0 lacking a valid lag" line the other rows use.
  did_panel_damage_ok <- did_panel_hazard_ok %>% filter(!is.na(log_damage_climate_lag1))
  n_dropped_damage_na <- nrow(did_panel_hazard_ok) - nrow(did_panel_damage_ok)
  message(sprintf(
    paste0("Rows additionally dropped for missing damage value (event occurred, none ",
          "reported damage): %d / %d (%.1f%%)."),
    n_dropped_damage_na, nrow(did_panel_hazard_ok),
    100 * n_dropped_damage_na / nrow(did_panel_hazard_ok)
  ))

  # Headline reference row: the same
  # specification on the FULL main panel (Table 2 sample), so the table
  # carries its own reference point and the sample-restriction effect of the
  # EM-DAT merge is visible inside the exhibit.
  fit_headline_ref <- hazard_headline_spec("log_commits", ~ ge_est + log_population,
                                           did_panel_main)
  n_units_main   <- n_distinct(did_panel_main$country_id)
  n_units_hazard <- n_distinct(did_panel_hazard_ok$country_id)
  n_units_lost_hazard <- n_units_main - n_units_hazard
  fit_baseline <- hazard_headline_spec("log_commits", ~ ge_est + log_population,
                                       did_panel_hazard_ok)
  fit_a_affected <- hazard_headline_spec(
    "log_commits", ~ ge_est + log_population + log_affected_climate_lag1, did_panel_hazard_ok)
  fit_a_damage <- hazard_headline_spec(
    "log_commits", ~ ge_est + log_population + log_damage_climate_lag1,
    did_panel_damage_ok)
  fit_c_ma3 <- hazard_headline_spec(
    "log_commits", ~ ge_est + log_population + log_affected_climate_ma3,
    did_panel_hazard_ok %>% filter(!is.na(log_affected_climate_ma3)))
  fit_variant_all <- hazard_headline_spec(
    "log_commits", ~ ge_est + log_population + log_affected_all_lag1,
    did_panel_hazard_ok %>% filter(!is.na(log_affected_all_lag1)))

  # --- (b) residualised-outcome spec (primary hazard var: affected, climate) -
  never_treated_sub <- did_panel_hazard_ok %>% filter(cohort_year == 0)
  resid_fit <- fixest::feols(
    log_commits ~ log_affected_climate_lag1 | country_id + year, data = never_treated_sub
  )
  beta_hazard <- coef(resid_fit)[["log_affected_climate_lag1"]]
  beta_hazard_se <- fixest::se(resid_fit)[["log_affected_climate_lag1"]]
  message(sprintf(
    "Never-treated hazard slope (country + year FE): beta = %.5f (SE = %.5f, n = %d units).",
    beta_hazard, beta_hazard_se, n_distinct(never_treated_sub$country_id)
  ))
  did_panel_resid <- did_panel_hazard_ok %>%
    mutate(log_commits_hazard_resid = log_commits - beta_hazard * log_affected_climate_lag1)
  fit_b_resid <- hazard_headline_spec("log_commits_hazard_resid", ~ ge_est + log_population,
                                      did_panel_resid)

  # --- Aligned-IF difference tests vs. baseline (§6 helper) ------------------
  diff_a_affected <- compute_att_difference(
    fit_baseline$agg_s$overall.att, fit_a_affected$agg_s$overall.att,
    fit_baseline$agg_s$inf.function$simple.att, fit_a_affected$agg_s$inf.function$simple.att,
    fit_baseline$agg_s$overall.se, fit_a_affected$agg_s$overall.se,
    fit_baseline$ids, fit_a_affected$ids)
  diff_b_resid <- compute_att_difference(
    fit_baseline$agg_s$overall.att, fit_b_resid$agg_s$overall.att,
    fit_baseline$agg_s$inf.function$simple.att, fit_b_resid$agg_s$inf.function$simple.att,
    fit_baseline$agg_s$overall.se, fit_b_resid$agg_s$overall.se,
    fit_baseline$ids, fit_b_resid$ids)
  # Row (c) is estimated on a DIFFERENT time window (2011+, three lags
  # needed) than baseline (first_year+). Its unit-id SET typically still
  # matches baseline's (dropping early years drops observations, not whole
  # units), so compute_att_difference()'s own same_units check -- which only
  # compares unit ids, not time coverage -- would wrongly take the aligned-IF
  # branch. Rather than modify the verbatim-copied helper (05:428-475) to add
  # time-window awareness it was never designed for, the independence
  # approximation is computed directly here (replicating that helper's own
  # fallback formula) and forced for this one comparison, with the reason
  # stated explicitly.
  indep_diff <- function(fit_a, fit_b, reason) {
    d_att <- fit_a$agg_s$overall.att - fit_b$agg_s$overall.att
    d_se  <- sqrt(fit_a$agg_s$overall.se^2 + fit_b$agg_s$overall.se^2)
    d_z   <- if (d_se > 1e-12) d_att / d_se else NA_real_
    list(diff = d_att, se = d_se, z = d_z,
         pval = if (is.na(d_z)) NA_real_ else 2 * pnorm(-abs(d_z)),
         method = paste0("independence approximation (forced: ", reason, ")"),
         n_units = NA_integer_)
  }
  # Damage row: its estimation sample differs from the
  # baseline by every recipient-year without a recorded damage value, so the
  # aligned-IF test (which requires the SAME estimation sample) does not apply.
  diff_a_damage <- indep_diff(fit_baseline, fit_a_damage, paste0(
    "the damage row drops ", n_dropped_damage_na, " observations with events but no ",
    "recorded damage value; aligned-IF requires the same estimation sample"))
  # Same-sample baseline vs. Table 2 headline: different samples.
  diff_ref_baseline <- indep_diff(fit_headline_ref, fit_baseline, paste0(
    "the hazard sample drops ", n_dropped_hazard_na, " observations and ",
    n_units_lost_hazard, " recipient(s) relative to the headline sample"))
  diff_c_ma3 <- {
    d_att   <- fit_baseline$agg_s$overall.att - fit_c_ma3$agg_s$overall.att
    d_se    <- sqrt(fit_baseline$agg_s$overall.se^2 + fit_c_ma3$agg_s$overall.se^2)
    d_z     <- if (d_se > 1e-12) d_att / d_se else NA_real_
    list(diff = d_att, se = d_se, z = d_z,
        pval = if (is.na(d_z)) NA_real_ else 2 * pnorm(-abs(d_z)),
        method = paste0("independence approximation (forced: row (c) is estimated on a ",
                        "2011--", last_year, " window vs. baseline's ", first_year, "--",
                        last_year, "; aligned-IF requires the same estimation sample, not ",
                        "just overlapping unit ids)"),
        n_units = NA_integer_)
  }
  diff_variant_all <- if (fit_variant_all$n_obs != fit_baseline$n_obs) {
    indep_diff(fit_baseline, fit_variant_all, paste0(
      "the all-natural variant's sample differs from the baseline by ",
      abs(fit_variant_all$n_obs - fit_baseline$n_obs), " observation(s)"))
  } else compute_att_difference(
    fit_baseline$agg_s$overall.att, fit_variant_all$agg_s$overall.att,
    fit_baseline$agg_s$inf.function$simple.att, fit_variant_all$agg_s$inf.function$simple.att,
    fit_baseline$agg_s$overall.se, fit_variant_all$agg_s$overall.se,
    fit_baseline$ids, fit_variant_all$ids)

  # --- Propensity-score overlap diagnostic for row (a), affected/climate ----
  # Fitting the logit on a SINGLE cross-sectional year would drop any unit
  # missing a covariate in exactly that year BEFORE the overlap trim, so
  # "trimmed N of N" would reflect data availability, not overlap. The logit
  # is therefore fit on RECIPIENT-LEVEL PRE-TREATMENT MEANS (mean
  # of each covariate over a unit's own pre-treatment years -- all years for
  # never-treated units, years < cohort_year for treated units), so every
  # main-sample recipient with at least one valid pre-treatment observation
  # gets a pscore. Units whose pre-treatment window is ENTIRELY missing one
  # covariate (e.g. no EM-DAT entity) cannot get a pscore either way
  # and are reported separately, not silently folded into "trimmed".
  pretrend_means <- did_panel_hazard_ok %>%
    mutate(pretrend_cutoff = if_else(cohort_year > 0, cohort_year - 1L, last_year)) %>%
    filter(year <= pretrend_cutoff) %>%
    group_by(country_id, recipient_name, cohort_year) %>%
    summarise(
      ge_est_mean                     = mean(ge_est, na.rm = TRUE),
      log_population_mean             = mean(log_population, na.rm = TRUE),
      log_affected_climate_lag1_mean  = mean(log_affected_climate_lag1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(is_treated = as.integer(cohort_year > 0))
  n_main_sample_units <- n_units_main   # 127 main-sample recipients (Table 2), not the hazard sample
  n_no_pscore <- sum(!complete.cases(
    pretrend_means[, c("ge_est_mean", "log_population_mean", "log_affected_climate_lag1_mean")]))
  cross_pscore <- pretrend_means %>%
    filter(is.finite(ge_est_mean), is.finite(log_population_mean),
          is.finite(log_affected_climate_lag1_mean))
  pscore_model <- glm(is_treated ~ ge_est_mean + log_population_mean +
                      log_affected_climate_lag1_mean,
                      data = cross_pscore, family = binomial())
  cross_pscore$pscore <- predict(pscore_model, type = "response")
  share_treated_out <- mean(cross_pscore$pscore[cross_pscore$is_treated == 1] < 0.05 |
                            cross_pscore$pscore[cross_pscore$is_treated == 1] > 0.95)
  share_control_out <- mean(cross_pscore$pscore[cross_pscore$is_treated == 0] < 0.05 |
                            cross_pscore$pscore[cross_pscore$is_treated == 0] > 0.95)
  trimmed_ids <- cross_pscore$country_id[cross_pscore$pscore >= 0.05 & cross_pscore$pscore <= 0.95]
  n_trimmed_out <- nrow(cross_pscore) - length(trimmed_ids)
  message(sprintf(
    paste0("Pscore overlap (row a, affected/climate, logit on RECIPIENT-LEVEL PRE-TREATMENT ",
          "MEANS, n = %d of %d main-sample units; %d lack a pscore for missing pre-treatment ",
          "covariate data): %.1f%% of treated and %.1f%% of control units have pscore outside ",
          "[0.05, 0.95]; %d unit(s) trimmed on overlap, %d retained."),
    nrow(cross_pscore), n_main_sample_units, n_no_pscore,
    100 * share_treated_out, 100 * share_control_out, n_trimmed_out, length(trimmed_ids)
  ))
  no_unit_trimmed <- (n_trimmed_out == 0)
  fit_a_trimmed <- if (no_unit_trimmed) {
    message("No unit trimmed on overlap; row (a) trimmed re-estimate is identical to row (a).")
    fit_a_affected
  } else {
    hazard_headline_spec(
      "log_commits", ~ ge_est + log_population + log_affected_climate_lag1,
      did_panel_hazard_ok %>% filter(country_id %in% trimmed_ids))
  }

  # --- Warning roll-up ---------------------------------------------------
  all_fits_for_warn <- list(fit_headline_ref, fit_baseline, fit_a_affected, fit_a_damage,
                            fit_c_ma3, fit_variant_all, fit_b_resid)
  if (!no_unit_trimmed) all_fits_for_warn <- c(all_fits_for_warn, list(fit_a_trimmed))
  n_warnings_total <- sum(vapply(all_fits_for_warn, function(f) length(f$warnings), integer(1)))

  hazard_rows <- list(
    list(label = "Headline (full main panel)",    set = "--", fit = fit_headline_ref, diff = NULL),
    list(label = "Baseline, hazard sample (no hazard control)", set = "--",
        fit = fit_baseline, diff = diff_ref_baseline),
    list(label = "(a) Affected, $g-1$",          set = "Climate-related",
        fit = fit_a_affected, diff = diff_a_affected),
    list(label = "(a) Damage, $g-1$",            set = "Climate-related",
        fit = fit_a_damage, diff = diff_a_damage),
    list(label = "(b) Residualised outcome",     set = "Climate-related",
        fit = fit_b_resid, diff = diff_b_resid),
    list(label = "(c) 3-yr moving avg., affected", set = "Climate-related",
        fit = fit_c_ma3, diff = diff_c_ma3),
    list(label = "Variant: affected, $g-1$",     set = "All natural",
        fit = fit_variant_all, diff = diff_variant_all)
  )

  hazard_tab <- bind_rows(lapply(hazard_rows, function(r) {
    if (is.null(r$fit) || is.null(r$fit$agg_s)) {
      return(data.frame(spec = r$label, set = r$set, ATT = NA_real_, SE = NA_real_,
                        p = NA_real_, n_obs = NA_integer_, n_treated = NA_integer_,
                        diff = NA_real_, diff_se = NA_real_, diff_p = NA_real_))
    }
    att <- r$fit$agg_s$overall.att; se <- r$fit$agg_s$overall.se
    z <- att / se; p <- 2 * pnorm(-abs(z))
    data.frame(spec = r$label, set = r$set, ATT = att, SE = se, p = p,
              n_obs = r$fit$n_obs, n_treated = r$fit$n_treated,
              diff = if (is.null(r$diff)) NA_real_ else r$diff$diff,
              diff_se = if (is.null(r$diff)) NA_real_ else r$diff$se,
              diff_p  = if (is.null(r$diff)) NA_real_ else r$diff$pval)
  })) %>% mutate(stars = stars_of(p), diff_stars = stars_of(diff_p))

  message("\n--- Hazard-controls results ---")
  print(hazard_tab, row.names = FALSE)

  hazard_lines <- c(
    "\\begin{tabular}{llccccccc}",
    "\\toprule",
    "Specification & Hazard set & ATT & SE & $N$ & $N_{treated}$ & $\\Delta$ vs.\\ reference & SE($\\Delta$) & $p(\\Delta)$ \\\\",
    "\\midrule",
    paste0(hazard_tab$spec, " & ", hazard_tab$set, " & ", fmt4(hazard_tab$ATT), hazard_tab$stars,
           " & (", fmt4(hazard_tab$SE), ") & ",
           ifelse(is.na(hazard_tab$n_obs), "--", format(hazard_tab$n_obs, big.mark = ",")), " & ",
           ifelse(is.na(hazard_tab$n_treated), "--", hazard_tab$n_treated), " & ",
           fmt4(hazard_tab$diff), hazard_tab$diff_stars, " & (", fmt4(hazard_tab$diff_se),
           ") & ", fmt3(hazard_tab$diff_p), " \\\\"),
    "\\bottomrule",
    "\\end{tabular}"
  )

  all_warn_texts <- unlist(lapply(all_fits_for_warn, function(f) f$warnings))
  warn_tab <- if (length(all_warn_texts) > 0) sort(table(gsub("\\s+", " ", trimws(all_warn_texts))), decreasing = TRUE) else NULL
  # Singular-design warnings (one per failed (g,t) cell) are collapsed into a
  # single entry listing the cells; every other distinct warning is kept.
  warn_listing <- if (!is.null(warn_tab)) {
    nm <- names(warn_tab); is_sing <- grepl("singular", nm, ignore.case = TRUE)
    cells <- regmatches(nm[is_sing], regexpr("\\(g, t\\) = \\([0-9]+, [0-9]+\\)", nm[is_sing]))
    cells <- sub("\\(g, t\\) = ", "", cells)
    sing_g_all <- rep(sub("\\((\\d+), \\d+\\)", "\\1", cells), times = as.integer(warn_tab[is_sing]))
    sing_by_g  <- table(sing_g_all)
    sing_txt <- if (any(is_sing)) paste0("``singular design, $ATT(g,t)$ set NA'' for ",
      sum(warn_tab[is_sing]), " cells (",
      esc_tex(paste(paste0(names(sing_by_g), ":", as.integer(sing_by_g)), collapse = ", ")),
      "; cell list in run log)") else NULL
    # Non-singular warnings are grouped by their first 55 characters so near-
    # duplicate messages (same complaint, different group list) collapse into
    # one note entry with a summed count, rather than one entry each.
    other_txt <- if (any(!is_sing)) {
      prefixes <- substr(nm[!is_sing], 1, 32)
      by_prefix <- tapply(as.integer(warn_tab[!is_sing]), prefixes, sum)
      vapply(seq_along(by_prefix), function(i) paste0("``", esc_tex(names(by_prefix)[i]),
             "\\dots'' ($\\times$", by_prefix[[i]], ")"), character(1))
    } else NULL
    paste(c(sing_txt, other_txt), collapse = "; ")
  } else "none"

  write_tex_float(
    out_path      = here("output", "tables", "hazard", "att_hazard_controls.tex"),
    caption_title = "NAP effect on log adaptation commitments with EM-DAT hazard-realisation controls",
    label         = "tab:hazard_controls",
    tabular_lines = hazard_lines,
    notes_text    = paste0(
      "CS\\,(2021) DR, never-treated controls, WGI GE + log population, cohorts $\\geq 5$, ",
      "multiplier-bootstrap SE (999 reps, seed 1242), \\texttt{did} ", DID_VERSION,
      ". Row 1: headline fit ($\\Delta$ reference). Row 2: hazard sample (missing EM-DAT ",
      "lag dropped). Hazard set CLIMATE-RELATED unless noted (last row ALL-",
      "NATURAL). (a) adds the lagged hazard measure to \\texttt{xformla}; (b) residualises the ",
      "outcome on it via never-treated FE, re-estimates on the residual; (c) uses a 3-yr ",
      "moving average. $N_{treated}$: recipients in cohorts with an estimable post-treatment ",
      "cell. $p(\\Delta)$: independence-approximation upper bound"
    ),
    source_text = paste0("OECD CRS (", CRS_VINTAGE, "); UNFCCC NAP Central; World Bank ",
                         "WGI/WDI; ", EMDAT_VINTAGE, ", EM-DAT, CRED / UCLouvain, Brussels, ",
                         "Belgium -- www.emdat.be")
  )
}

##############################################################################
# SECTION 4. NAP-timing orthogonality checks (two separate tables)
#
# Discrete-time hazard/duration setup common to both: each adopter-country-
# year is "at risk" of NAP submission until it either adopts (event = 1,
# then exits the risk set -- later years are not included) or is right-
# censored (never-treated countries, observed every year, event always 0).
# For each hazard channel this runs: a two-lag LPM (feols, country+year FE,
# clustered SE), single-lag LPMs, a two-lag LPM dropping 2020-2022 (COVID
# confound in both aid flows and disaster reporting), a two-lag cloglog
# discrete-time hazard model (feglm, same FE/clustering -- the correct
# functional form for a rare-event hazard, unlike the LPM used for
# comparison), and a joint Wald test of both lags for the two-lag LPM and
# the cloglog model.
#
# Table A: "NAP timing vs. humanitarian aid" -- the CRS emergency-response-
#   AID proxy (01 §23). This is a check that NAP adoption timing does not
#   simply follow recent humanitarian-aid surges -- NOT a disaster/hazard
#   test (the regressor measures aid, not disaster realisations).
# Table B: "NAP timing vs. EM-DAT hazard" -- CLIMATE-RELATED log(1+affected)
#   at lag 1/2 (headline); the n_events version is reported compactly in the
#   note (joint Wald only) rather than as a full extra table to keep this
#   exhibit a manageable size.
##############################################################################

message("\n=== SECTION 4: NAP-timing orthogonality checks ===\n")

at_risk_panel <- did_panel_full %>%
  filter(cohort_year == 0 | (cohort_year > 0 & year <= cohort_year)) %>%
  mutate(nap_submit = as.integer(cohort_year > 0 & year == cohort_year))

message(sprintf(
  "At-risk sample (shared by both orthogonality tables): %d country-year rows (%d countries), ",
  nrow(at_risk_panel), n_distinct(at_risk_panel$recipient_name)))
message(sprintf("  %d NAP-submission events.", sum(at_risk_panel$nap_submit)))

#' Runs the full orthogonality battery for one hazard channel and writes the
#' table. lag1_var/lag2_var must be columns of at_risk_panel.
run_orthogonality_battery <- function(lag1_var, lag2_var, out_path, caption_title,
                                      label, hazard_desc, source_text, extra_note = "") {
  two_lag <- at_risk_panel %>% filter(!is.na(.data[[lag1_var]]), !is.na(.data[[lag2_var]]))
  n_dropped <- nrow(at_risk_panel) - nrow(two_lag)
  two_lag_nocovid <- two_lag %>% filter(!year %in% 2020:2022)

  fm2_str <- sprintf("nap_submit ~ %s + %s | country_id + year", lag1_var, lag2_var)
  fm1_str <- sprintf("nap_submit ~ %s | country_id + year", lag1_var)
  fm2b_str <- sprintf("nap_submit ~ %s | country_id + year", lag2_var)

  m_lpm      <- fixest::feols(as.formula(fm2_str), data = two_lag, cluster = ~country_id)
  m_lag1     <- fixest::feols(as.formula(fm1_str), data = two_lag, cluster = ~country_id)
  m_lag2     <- fixest::feols(as.formula(fm2b_str), data = two_lag, cluster = ~country_id)
  m_nocovid  <- fixest::feols(as.formula(fm2_str), data = two_lag_nocovid, cluster = ~country_id)
  cloglog_error_msg <- NULL
  m_cloglog  <- tryCatch(
    fixest::feglm(as.formula(fm2_str), data = two_lag,
                  family = binomial(link = "cloglog"), cluster = ~country_id),
    error = function(e) {
      # Sanctioned exception to the `<<-` rule (condition capture; see §2 note).
      cloglog_error_msg <<- conditionMessage(e)
      message("  cloglog failed: ", cloglog_error_msg)
      NULL
    }
  )

  wald_lpm     <- tryCatch(fixest::wald(m_lpm, keep = c(lag1_var, lag2_var)),
                           error = function(e) list(stat = NA_real_, p = NA_real_))
  wald_nocovid <- tryCatch(fixest::wald(m_nocovid, keep = c(lag1_var, lag2_var)),
                           error = function(e) list(stat = NA_real_, p = NA_real_))
  wald_cloglog <- if (!is.null(m_cloglog)) tryCatch(
    fixest::wald(m_cloglog, keep = c(lag1_var, lag2_var)),
    error = function(e) list(stat = NA_real_, p = NA_real_)) else list(stat = NA_real_, p = NA_real_)

  extract <- function(model, term) {
    if (is.null(model)) return(list(est = NA_real_, se = NA_real_, p = NA_real_))
    tt <- tryCatch(broom::tidy(model), error = function(e) NULL)
    if (is.null(tt) || !(term %in% tt$term)) return(list(est = NA_real_, se = NA_real_, p = NA_real_))
    r <- tt[tt$term == term, ]
    list(est = r$estimate[1], se = r$std.error[1], p = r$p.value[1])
  }

  models <- list(
    list(lbl = "Two-lag (LPM)",              m = m_lpm,     has1 = TRUE,  has2 = TRUE),
    list(lbl = "Lag-1 only (LPM)",            m = m_lag1,    has1 = TRUE,  has2 = FALSE),
    list(lbl = "Lag-2 only (LPM)",            m = m_lag2,    has1 = FALSE, has2 = TRUE),
    list(lbl = "Two-lag, drop 2020--22 (LPM)", m = m_nocovid, has1 = TRUE,  has2 = TRUE),
    list(lbl = "Two-lag (cloglog)",           m = m_cloglog, has1 = TRUE,  has2 = TRUE)
  )
  wald_p <- c(wald_lpm$p, NA_real_, NA_real_, wald_nocovid$p, wald_cloglog$p)
  n_rows <- c(nrow(two_lag), nrow(two_lag), nrow(two_lag), nrow(two_lag_nocovid),
              if (is.null(m_cloglog)) NA_integer_ else nrow(two_lag))

  col1 <- vapply(seq_along(models), function(i) {
    if (!models[[i]]$has1) return("--")
    e <- extract(models[[i]]$m, lag1_var)
    paste0(fmt4(e$est), stars_of(e$p))
  }, character(1))
  col1_se <- vapply(seq_along(models), function(i) {
    if (!models[[i]]$has1) return("")
    e <- extract(models[[i]]$m, lag1_var)
    paste0("(", fmt4(e$se), ")")
  }, character(1))
  col2 <- vapply(seq_along(models), function(i) {
    if (!models[[i]]$has2) return("--")
    e <- extract(models[[i]]$m, lag2_var)
    paste0(fmt4(e$est), stars_of(e$p))
  }, character(1))
  col2_se <- vapply(seq_along(models), function(i) {
    if (!models[[i]]$has2) return("")
    e <- extract(models[[i]]$m, lag2_var)
    paste0("(", fmt4(e$se), ")")
  }, character(1))

  hdr <- paste0(" & ", paste(vapply(models, `[[`, character(1), "lbl"), collapse = " & "), " \\\\")
  row_lag1    <- paste0("Lag-1 & ", paste(col1, collapse = " & "), " \\\\")
  row_lag1_se <- paste0(" & ", paste(col1_se, collapse = " & "), " \\\\")
  row_lag2    <- paste0("Lag-2 & ", paste(col2, collapse = " & "), " \\\\")
  row_lag2_se <- paste0(" & ", paste(col2_se, collapse = " & "), " \\\\")
  row_wald    <- paste0("Joint Wald $p$ (both lags) & ", paste(fmt3(wald_p), collapse = " & "), " \\\\")
  row_n       <- paste0("$N$ & ", paste(ifelse(is.na(n_rows), "--", format(n_rows, big.mark = ",")), collapse = " & "), " \\\\")
  row_events  <- paste0("NAP-submission events & ",
                        paste(rep(sum(two_lag$nap_submit), length(models)), collapse = " & "),
                        " \\\\")
  row_fe      <- "Country + year FE & \\multicolumn{5}{c}{Yes} \\\\"

  ortho_lines <- c(
    paste0("\\begin{tabular}{l", strrep("c", length(models)), "}"),
    "\\toprule",
    hdr, "\\midrule",
    row_lag1, row_lag1_se, row_lag2, row_lag2_se, "\\midrule",
    row_wald, row_n, row_events, row_fe,
    "\\bottomrule", "\\end{tabular}"
  )

  # Auto-generated joint-rejection / COVID-sensitivity sentence: states, with the actual computed numbers, whether
  # the two-lag joint test rejects, and whether that rejection survives
  # dropping 2020-2022 -- rather than a hand-written claim that can drift out
  # of sync with the numbers.
  reject_full    <- !is.na(wald_lpm$p) && wald_lpm$p < 0.05
  reject_nocovid <- !is.na(wald_nocovid$p) && wald_nocovid$p < 0.05
  lag1_sig <- !is.na(extract(m_lpm, lag1_var)$p) && extract(m_lpm, lag1_var)$p < 0.05
  lag2_sig <- !is.na(extract(m_lpm, lag2_var)$p) && extract(m_lpm, lag2_var)$p < 0.05
  driver <- if (lag2_sig && !lag1_sig) "lag-2" else if (lag1_sig && !lag2_sig) "lag-1" else
    if (lag1_sig && lag2_sig) "both lags" else "neither lag individually"
  rejection_sentence <- if (is.na(wald_lpm$p)) "" else paste0(
    "Joint two-lag test ", if (reject_full) "rejects" else "does not reject",
    " the null at 5\\% ($p = ", fmt3(wald_lpm$p), "$), driven by ", driver,
    "; dropping 2020--2022 ",
    if (is.na(wald_nocovid$p)) "could not be tested" else paste0(
      "gives $p = ", fmt3(wald_nocovid$p), "$ (",
      if (reject_full && !reject_nocovid) "no longer significant -- COVID-sensitive"
      else if (reject_full && reject_nocovid) "still significant -- not just COVID"
      else "still not significant", ")"
    ), ". "
  )

  n_countries_two_lag <- n_distinct(two_lag$recipient_name)
  n_years_two_lag     <- n_distinct(two_lag$year)

  write_tex_float(
    out_path = out_path, caption_title = caption_title, label = label,
    tabular_lines = ortho_lines,
    notes_text = paste0(
      "LPM (cols 1--4), cloglog hazard (col 5); DV $= 1$ in the NAP-submission year, else 0; ",
      "post-adoption years excluded, never-treated countries censored. Country + year FE, SE ",
      "clustered by recipient. Regressor: ", hazard_desc, ". Col.\\ 4 drops 2020--2022. Joint ",
      "Wald tests both lags jointly nonzero. Base rate ",
      sprintf("%.1f", 100 * mean(two_lag$nap_submit)), "\\% per country-year (",
      sum(two_lag$nap_submit), " events, ", n_countries_two_lag, " countries). ", n_dropped,
      " rows lacking a valid two-year lag dropped. Full cohort set (not $\\geq 5$). ",
      "\\texttt{did} ", DID_VERSION,
      if (is.null(m_cloglog))
        "; col.\\ 5 not estimable (fixed-effects singleton), shown as ``--''" else "",
      if (nzchar(extra_note)) paste0(". ", extra_note) else ""
    ),
    source_text = source_text
  )

  list(m_lpm = m_lpm, m_cloglog = m_cloglog, wald_lpm = wald_lpm, wald_cloglog = wald_cloglog,
      lag1 = extract(m_lpm, lag1_var), lag2 = extract(m_lpm, lag2_var), n = nrow(two_lag))
}

# --- Table A: NAP timing vs. humanitarian aid (CRS proxy) -------------------
res_aid <- run_orthogonality_battery(
  lag1_var = "log_aid_lag1", lag2_var = "log_aid_lag2",
  out_path = here("output", "tables", "hazard", "nap_timing_vs_humanitarian_aid.tex"),
  caption_title = "NAP submission timing and lagged emergency-response aid",
  label = "tab:nap_timing_humanitarian_aid",
  hazard_desc = paste0("one- and two-year-lagged log(1+CRS emergency-response aid), purpose ",
                       "codes 72010/72040/72050/73010/74020; a ",
                       "humanitarian-aid-flow proxy, not a hazard measure"),
  source_text = paste0("OECD CRS (", CRS_VINTAGE, "); UNFCCC NAP Central")
)
message(sprintf("NAP timing vs. humanitarian aid: lag-1 = %.4f (p=%.3f), lag-2 = %.4f (p=%.3f), joint Wald p = %s.",
                res_aid$lag1$est, res_aid$lag1$p, res_aid$lag2$est, res_aid$lag2$p,
                fmt3(res_aid$wald_lpm$p)))

# --- Table B: NAP timing vs. EM-DAT hazard (climate-related) ---------------
if (emdat_available) {
  res_emdat <- run_orthogonality_battery(
    lag1_var = "log_affected_climate_lag1", lag2_var = "log_affected_climate_lag2",
    out_path = here("output", "tables", "hazard", "nap_timing_vs_emdat_hazard.tex"),
    caption_title = "NAP submission timing and lagged EM-DAT natural-hazard realisations",
    label = "tab:nap_timing_emdat_hazard",
    hazard_desc = paste0("one- and two-year-lagged log(1+Total Affected), climate-related ",
                         "EM-DAT hazard set (Hydrological/Meteorological/Climatological)"),
    source_text = paste0("UNFCCC NAP Central; ", EMDAT_VINTAGE,
                         ", EM-DAT, CRED / UCLouvain, Brussels, Belgium -- www.emdat.be"),
    extra_note = ""
  )
  message(sprintf("NAP timing vs. EM-DAT hazard (affected): lag-1 = %.4f (p=%.3f), lag-2 = %.4f (p=%.3f), joint Wald p = %s.",
                  res_emdat$lag1$est, res_emdat$lag1$p, res_emdat$lag2$est, res_emdat$lag2$p,
                  fmt3(res_emdat$wald_lpm$p)))

  # Compact n_events alternative (joint Wald only, reported inline -- not a
  # separate table, to keep this exhibit a manageable size).
  events_2lag <- at_risk_panel %>%
    filter(!is.na(n_events_climate_lag1), !is.na(n_events_climate_lag2))
  m_events <- fixest::feols(
    nap_submit ~ n_events_climate_lag1 + n_events_climate_lag2 | country_id + year,
    data = events_2lag, cluster = ~country_id)
  wald_events <- tryCatch(
    fixest::wald(m_events, keep = c("n_events_climate_lag1", "n_events_climate_lag2")),
    error = function(e) list(p = NA_real_))
  message(sprintf(
    "NAP timing vs. EM-DAT hazard (n_events, compact alternative, N=%d): joint Wald p = %s.",
    nrow(events_2lag), fmt3(wald_events$p)))
} else {
  message("SKIPPED (EM-DAT not found): output/tables/hazard/nap_timing_vs_emdat_hazard.tex ",
          "NOT regenerated -- see data/raw/emdat/README.md.")
}

##############################################################################
# SECTION 5. Prior-NAPA list construction
#
# Preferred source: data/raw/napa/napa_list_unfccc.csv, the UNFCCC "NAPAs
# received" list (51 NAPAs, latest South Sudan Feb. 2017), transcribed by
# the authors on 2026-09-15 from the live UNFCCC page (which blocks scripted
# access -- see data/raw/napa/README.md) and cross-checked 51/51 against the
# archived page PDF (data/raw/napa/Submitted_NAPAs_UNFCCC.pdf). This script
# does NOT attempt to scrape the page itself: a live curl/system() call in
# an analysis script is fragile (environment-dependent, network-dependent,
# silently stale) and unnecessary now that a verified, provenance-documented
# source file exists. If the file is absent the script stops: there is no
# fallback list. See data/raw/napa/README.md for the primary-source
# retrieval method.
##############################################################################

message("\n=== SECTION 5: NAPA list construction ===\n")

napa_src  <- here("data", "raw", "napa", "napa_list_unfccc.csv")
napa_list <- NULL

if (file.exists(napa_src)) {
  napa_list <- read.csv(napa_src, stringsAsFactors = FALSE) %>%
    select(iso3, country, napa_year, napa_month, source)
  message("Verified NAPA list found: data/raw/napa/napa_list_unfccc.csv (",
          nrow(napa_list), " countries; source = ", unique(napa_list$source), ").")
} else {
  stop("data/raw/napa/napa_list_unfccc.csv not found. This file ships with the ",
       "replication package; see data/raw/napa/README.md for its source and how it ",
       "was transcribed from the UNFCCC 'NAPAs received' page.")
}

stopifnot(!any(is.na(napa_list$iso3)))
n_napa <- nrow(napa_list)  # computed from the file, never hardcoded

# Provenance string reused in every NAPA table note.
napa_source_note <- paste0(
  "UNFCCC `NAPAs received' page (", n_napa, " NAPAs, latest South Sudan, Feb.\\ 2017), ",
  "transcribed by the authors 2026-09-15, cross-checked ", n_napa, "/", n_napa,
  " against the archived page PDF; data/raw/napa/napa\\_list\\_unfccc.csv")
# Short label for Source: fields (the full provenance sentence belongs in
# Notes only; repeating it in Source would be redundant).
napa_source_short <- "UNFCCC `NAPAs received' list (transcribed by the authors, verified -- see Notes)"

write.csv(napa_list, here("data", "processed", "napa_list.csv"), row.names = FALSE)
message("Wrote: data/processed/napa_list.csv  (", n_napa, " countries; source = ",
        paste(unique(napa_list$source), collapse = ", "), ")")

##############################################################################
# SECTION 6. LDC classification (for the prior-NAPA x LDC cross-tab in §7)
# Copied VERBATIM from 05_heterogeneity.R's ldc_iso3_2013 vintage (49
# economies, UN list as of 1 January 2013 -- the pre-treatment vintage 05
# itself settled on for its own LDC split; 05
# itself is NOT modified or sourced by this script).
##############################################################################

ldc_iso3_2013 <- c(
  "AFG", "AGO", "BGD", "BEN", "BTN", "BFA", "BDI", "KHM", "CAF", "TCD",
  "COM", "COD", "DJI", "GNQ", "ERI", "ETH", "GMB", "GIN", "GNB", "HTI",
  "KIR", "LAO", "LSO", "LBR", "MDG", "MWI", "MLI", "MRT", "MOZ", "MMR",
  "NPL", "NER", "RWA", "WSM", "STP", "SEN", "SLE", "SLB", "SOM", "SSD",
  "SDN", "TLS", "TGO", "TUV", "UGA", "TZA", "VUT", "YEM", "ZMB"
)
stopifnot(length(ldc_iso3_2013) == 49L, !anyDuplicated(ldc_iso3_2013))

##############################################################################
# SECTION 7. Prior-NAPA split
# Among NAP adopters in the main sample (cohorts >= 5 treated units), split
# by whether the country had submitted a NAPA before its NAP, and estimate
# the NAP effect in each cell (same spec as the LDC split in 05: CS (2021),
# doubly-robust, never-treated control, analytical SE -- see 05 §22 for the
# split machinery this mirrors; 05 itself was not edited).
#
# Because every
# verified NAPA year is <= 2017 and every main-sample cohort is 2021-2024,
# "prior NAPA" is essentially "was an LDC in 2013" (NAPAs were an LDC-only
# UNFCCC instrument) on the current panel -- see the cross-tab below for the
# exact per-cell counts (nearly all prior-NAPA adopters are LDCs, and every
# no-prior-NAPA adopter is a non-LDC). This split is close to the LDC split
# under another name, not an independent test; read it accordingly.
##############################################################################

message("\n=== SECTION 7: Prior-NAPA split ===\n")

did_panel_napa <- did_panel_main %>%
  left_join(napa_list %>% select(recipient_iso = iso3, napa_year), by = "recipient_iso") %>%
  mutate(
    prior_napa = as.integer(cohort_year > 0 & !is.na(napa_year) & napa_year < cohort_year),
    is_ldc     = as.integer(recipient_iso %in% ldc_iso3_2013)
  )

adopters_prior_napa <- did_panel_napa %>%
  filter(cohort_year > 0) %>%
  distinct(recipient_name, cohort_year, prior_napa, is_ldc)
message("Adopters in the main sample by prior-NAPA status:")
print(count(adopters_prior_napa, prior_napa), row.names = FALSE)

# --- Cross-tab: prior_napa x LDC status, with per-cell treated-unit counts -
cross_tab <- adopters_prior_napa %>%
  count(prior_napa, is_ldc, name = "n_treated") %>%
  mutate(
    prior_napa_lbl = ifelse(prior_napa == 1L, "Prior NAPA", "No prior NAPA"),
    is_ldc_lbl     = ifelse(is_ldc == 1L, "LDC (2013 list)", "Non-LDC")
  )
message("\nPrior-NAPA x LDC-status cross-tab (treated-unit counts, main sample):")
print(cross_tab %>% select(prior_napa_lbl, is_ldc_lbl, n_treated), row.names = FALSE)

cross_tab_full <- expand.grid(prior_napa = c(0L, 1L), is_ldc = c(0L, 1L)) %>%
  left_join(cross_tab, by = c("prior_napa", "is_ldc")) %>%
  mutate(
    n_treated      = coalesce(n_treated, 0L),
    prior_napa_lbl = ifelse(prior_napa == 1L, "Prior NAPA", "No prior NAPA"),
    is_ldc_lbl     = ifelse(is_ldc == 1L, "LDC (2013 list)", "Non-LDC")
  ) %>%
  arrange(prior_napa, is_ldc)

cross_tab_lines <- c(
  "\\begin{tabular}{lcc}",
  "\\toprule",
  " & LDC (2013 list) & Non-LDC \\\\",
  "\\midrule",
  paste0("Prior NAPA & ",
         cross_tab_full$n_treated[cross_tab_full$prior_napa == 1 & cross_tab_full$is_ldc == 1],
         " & ",
         cross_tab_full$n_treated[cross_tab_full$prior_napa == 1 & cross_tab_full$is_ldc == 0],
         " \\\\"),
  paste0("No prior NAPA & ",
         cross_tab_full$n_treated[cross_tab_full$prior_napa == 0 & cross_tab_full$is_ldc == 1],
         " & ",
         cross_tab_full$n_treated[cross_tab_full$prior_napa == 0 & cross_tab_full$is_ldc == 0],
         " \\\\"),
  "\\bottomrule",
  "\\end{tabular}"
)
write_tex_float(
  out_path      = here("output", "tables", "napa", "prior_napa_ldc_crosstab.tex"),
  caption_title = "Prior-NAPA status by LDC classification, main-sample NAP adopters",
  label         = "tab:prior_napa_ldc_crosstab",
  tabular_lines = cross_tab_lines,
  notes_text    = paste0(
    "Treated-unit counts (main-sample NAP adopters, cohorts $\\geq 5$) by prior-NAPA status ",
    "and LDC classification (UN 2013 list, matching 05\\_heterogeneity.R's split). See main ",
    "text: the two partitions nearly coincide, so this split is the LDC split under another ",
    "name. NAPA list: ", napa_source_note
  ),
  source_text = paste0("UNFCCC NAP Central; ", napa_source_short)
)

napa_split_spec <- function(panel_in, keep_treated_flag) {
  panel_split <- panel_in %>%
    filter(cohort_year == 0 | prior_napa == keep_treated_flag)
  set.seed(1242)
  cap <- run_att_gt_captured(
    yname = "log_commits", tname = "year", idname = "country_id", gname = "cohort_year",
    xformla = ~ ge_est + log_population, data = panel_split, est_method = "dr",
    bstrap = FALSE, cband = FALSE, control_group = "nevertreated", anticipation = 0,
    base_period = "universal", panel = TRUE, allow_unbalanced_panel = TRUE
  )
  if (is.null(cap$gt_obj)) return(list(gt_obj = NULL, agg_s = NULL, n_units = NA_integer_,
                                       ids = integer(0)))
  agg_s <- tryCatch(aggte(cap$gt_obj, type = "simple", na.rm = TRUE), error = function(e) NULL)
  list(gt_obj = cap$gt_obj, agg_s = agg_s, n_units = n_distinct(panel_split$country_id),
      ids = sort(unique(panel_split$country_id)))
}

fit_prior_yes <- napa_split_spec(did_panel_napa, 1L)
fit_prior_no  <- napa_split_spec(did_panel_napa, 0L)

# Per-cell treated-unit counts (reported in the split table itself, not just
# the cross-tab).
n_treated_yes <- sum(adopters_prior_napa$prior_napa == 1L)
n_treated_no  <- sum(adopters_prior_napa$prior_napa == 0L)

napa_split_row <- function(label, fit, n_treated_cell) {
  if (is.null(fit) || is.null(fit$agg_s)) {
    return(data.frame(group = label, ATT = NA_real_, SE = NA_real_, p = NA_real_,
                      n_units = NA_integer_, n_treated = n_treated_cell))
  }
  att <- fit$agg_s$overall.att; se <- fit$agg_s$overall.se
  z <- att / se; p <- 2 * pnorm(-abs(z))
  data.frame(group = label, ATT = att, SE = se, p = p, n_units = fit$n_units,
            n_treated = n_treated_cell)
}
napa_split_tab <- bind_rows(
  napa_split_row("Prior NAPA", fit_prior_yes, n_treated_yes),
  napa_split_row("No prior NAPA", fit_prior_no, n_treated_no)
) %>% mutate(stars = stars_of(p))

# The two subgroup fits' unit-id
# INTERSECTION is exactly the 87 shared never-treated controls (fit_prior_yes
# and fit_prior_no both include every never-treated unit; their TREATED
# units are mutually exclusive by construction). Aligning influence functions
# on a set that contains zero treated units from EITHER fit is invalid --
# it would estimate Var(ATT_yes - ATT_no) from the control-side IF alone,
# discarding the very treated-unit variation the ATT difference depends on.
# did has no single-call mechanism for a joint fit of both subgroups against
# one shared control group with a group indicator, so this uses the independence approximation directly
# (replicating compute_att_difference()'s own fallback formula) rather than
# routing through that helper's id-intersection branch, with the reason
# stated explicitly in the table note.
diff_napa <- if (!is.null(fit_prior_yes$agg_s) && !is.null(fit_prior_no$agg_s)) {
  d_att <- fit_prior_yes$agg_s$overall.att - fit_prior_no$agg_s$overall.att
  d_se  <- sqrt(fit_prior_yes$agg_s$overall.se^2 + fit_prior_no$agg_s$overall.se^2)
  d_z   <- if (d_se > 1e-12) d_att / d_se else NA_real_
  list(diff = d_att, se = d_se, z = d_z,
      pval = if (is.na(d_z)) NA_real_ else 2 * pnorm(-abs(d_z)),
      method = paste0("independence approximation (the subgroups share all never-treated ",
                      "controls, so aligned-IF is not meaningful; conservative, $p$ is an upper bound)"),
      n_units = NA_integer_)
} else {
  list(diff = NA_real_, se = NA_real_, z = NA_real_, pval = NA_real_, method = "not estimable")
}
message(sprintf("Difference (prior NAPA - no prior NAPA) = %.4f, SE = %.4f, p = %.3f [%s]",
                diff_napa$diff, diff_napa$se, diff_napa$pval, diff_napa$method))

napa_split_lines <- c(
  "\\begin{tabular}{lccccc}",
  "\\toprule",
  "Group & ATT & SE & $p$-value & $N$ units & $N_{treated}$ \\\\",
  "\\midrule",
  paste0(napa_split_tab$group, " & ", fmt4(napa_split_tab$ATT), napa_split_tab$stars, " & (",
         fmt4(napa_split_tab$SE), ") & ", fmt3(napa_split_tab$p), " & ",
         ifelse(is.na(napa_split_tab$n_units), "--", napa_split_tab$n_units), " & ",
         napa_split_tab$n_treated, " \\\\"),
  "\\midrule",
  paste0("Difference (Prior $-$ No prior) & ", fmt4(diff_napa$diff), " & (",
         fmt4(diff_napa$se), ") & ", fmt3(diff_napa$pval), " & & \\\\"),
  "\\bottomrule",
  "\\end{tabular}"
)


write_tex_float(
  out_path      = here("output", "tables", "napa", "att_prior_napa_split.tex"),
  caption_title = "NAP effect by prior NAPA submission status",
  label         = "tab:prior_napa_split",
  tabular_lines = napa_split_lines,
  notes_text    = paste0(
    "CS\\,(2021) DR, never-treated controls, WGI GE + log population, cohorts $\\geq 5$, ",
    "analytical (IF) SE, \\texttt{did} ", DID_VERSION,
    ". ``Prior NAPA'': submitted a NAPA ",
    "to the UNFCCC before the NAP; NAPA list: ", napa_source_note,
    ". Difference row: ", diff_napa$method, ".",
    " Sample: ", first_year, "--", last_year
  ),
  source_text = paste0("OECD CRS (", CRS_VINTAGE, "); UNFCCC NAP Central; ", napa_source_short)
)

##############################################################################
# SECTION 8. NAPA-adoption falsification check
# Relabelled from "placebo test" (NAPA is the
# NAP's institutional predecessor, not an unrelated planning event, so this
# is a LOW-POWER FALSIFICATION CHECK, not a clean placebo -- a significant
# result here is informative but does not by itself indict the main
# estimate). Treats napa_year as G_i on the same outcome (log_commits).
#
# Design: on the full 2009-2024 panel, real NAP adopters (main-sample cohorts
# 2021-2024) would sit in the "control" pool while already treated by the
# REAL NAP, and any NAPA country that later adopted a real NAP inside the
# window would have its post-NAPA years contaminated by the real treatment.
# This is avoided by restricting the ENTIRE placebo panel to years <=
# (earliest real NAP adoption year - 1), computed from the data (not
# hardcoded): no real NAP adoption -- treated or control-side -- occurs
# within the restricted window, verified by the stopifnot() below, so no
# separate "remove real adopters from the control pool" step is needed on
# top of the truncation itself. Countries whose NAPA predates first_year are
# coded napa_cohort_year = 0 (pseudo-never-treated), IDENTICAL to how the
# headline spec's own cohort_year construction (Section 1) handles NAP
# adoptions before first_year -- not a new convention invented here.
##############################################################################

message("\n=== SECTION 8: NAPA-adoption falsification check ===\n")

earliest_real_nap <- min(did_panel_full$nap_year[did_panel_full$nap_year >= first_year],
                         na.rm = TRUE)
placebo_cutoff_year <- earliest_real_nap - 1L
message(sprintf(
  paste0("Earliest real NAP adoption in the panel: %d. Falsification panel restricted to ",
        "%d-%d so no real NAP adoption (treated or control side) occurs within the ",
        "estimation window."),
  earliest_real_nap, first_year, placebo_cutoff_year))

napa_cohort_lookup <- did_panel_full %>%
  distinct(recipient_name, recipient_iso) %>%
  left_join(napa_list %>% select(recipient_iso = iso3, napa_year), by = "recipient_iso") %>%
  mutate(
    napa_cohort_year = case_when(
      is.na(napa_year)                       ~ 0,
      napa_year < first_year                  ~ 0,
      napa_year > placebo_cutoff_year         ~ 0,
      TRUE                                    ~ napa_year
    )
  )

n_late_napa <- sum(napa_cohort_lookup$napa_cohort_year > 0)
napa_treated_names <- napa_cohort_lookup %>% filter(napa_cohort_year > 0) %>%
  arrange(napa_cohort_year, recipient_name)
message(sprintf(
  paste0("NAPA cohorts falling within the falsification window (%d-%d): %d countries (out of ",
        "%d NAPA-listed countries; the rest predate first_year or fall after the cutoff and ",
        "are pooled with never-treated)."),
  first_year, placebo_cutoff_year, n_late_napa, sum(!is.na(napa_cohort_lookup$napa_year))))
print(napa_treated_names %>% select(recipient_name, napa_cohort_year), row.names = FALSE)

did_panel_placebo <- did_panel_full %>%
  select(-any_of("napa_cohort_year")) %>%
  left_join(napa_cohort_lookup %>% select(recipient_name, napa_cohort_year),
            by = "recipient_name") %>%
  filter(year <= placebo_cutoff_year) %>%
  mutate(cohort_year = napa_cohort_year)  # overwrite: NAPA is now the falsification "treatment"

# Verify the design fix: no real NAP adoption falls inside the restricted window.
stopifnot(all(is.na(did_panel_placebo$nap_year) |
             did_panel_placebo$nap_year > placebo_cutoff_year |
             did_panel_placebo$nap_year < first_year))

if (n_late_napa < 5L) {

  message(sprintf(
    paste0("Falsification check INFEASIBLE: only %d NAPA cohorts fall within the restricted ",
          "window (need >= 5 treated units). Writing a placeholder table documenting the ",
          "infeasibility instead of a spurious estimate."),
    n_late_napa))
  placebo_lines <- c(
    "\\begin{tabular}{l}",
    "\\toprule",
    "Result \\\\",
    "\\midrule",
    paste0("Falsification check not estimated: only ", n_late_napa,
           " NAPA-adoption cohorts fall within the ", first_year, "--", placebo_cutoff_year,
           " restricted window (below the project's $\\geq 5$ treated-unit threshold). \\\\"),
    "\\bottomrule",
    "\\end{tabular}"
  )
  write_tex_float(
    out_path      = here("output", "tables", "napa", "att_napa_falsification.tex"),
    caption_title = "NAPA submission as a low-power falsification check (infeasible)",
    label         = "tab:napa_falsification",
    tabular_lines = placebo_lines,
    notes_text    = paste0(
      "NAPA list: ", napa_source_note, ". Falsification panel restricted to ", first_year,
      "--", placebo_cutoff_year, " (before the earliest real NAP adoption, ",
      earliest_real_nap, ") so no real NAP treatment contaminates either side. Only ",
      n_late_napa, " NAPA cohorts fall inside that window, below the $\\geq 5$ treated-unit ",
      "threshold used throughout this project, so no estimate is reported"
    ),
    source_text = paste0(napa_source_short, "; OECD CRS")
  )

} else {

  napa_cohort_sizes <- did_panel_placebo %>%
    filter(cohort_year > 0) %>% distinct(recipient_name, cohort_year) %>%
    count(cohort_year, name = "n_treated")
  napa_thin <- napa_cohort_sizes$cohort_year[napa_cohort_sizes$n_treated < thin_threshold]
  message("NAPA-cohort sizes (falsification treatment, restricted window):")
  print(napa_cohort_sizes, row.names = FALSE)

  # With treated units spread thinly across several NAPA-adoption years, no
  # individual cohort may meet the project's per-cohort thin_threshold = 5.
  # Applying the headline's own "drop cohorts < 5" rule here could drop every
  # treated cohort and leave att_gt() with zero treated groups. Instead this
  # uses the SAME fallback the headline spec itself uses for its own small
  # 2015-2020 cohorts (03_main_results.R's "cohorts_retained" spec): retain
  # every NAPA cohort, outcome-regression estimator, analytical SE.
  use_retained_fallback <- length(napa_thin) > 0
  if (use_retained_fallback) {
    message(sprintf(
      paste0("All or some NAPA cohorts (%s) fall below the per-cohort thin_threshold = %d ",
            "(total treated units = %d, which is >= 5, so the check is still run -- using ",
            "the 'cohorts_retained' style fallback: est_method = \"reg\", bstrap = FALSE, ",
            "analytical SE, matching 03_main_results.R's own small-cohort convention)."),
      paste(napa_thin, collapse = ", "), thin_threshold, n_late_napa))
  }
  use_dr_placebo     <- if (use_retained_fallback) "reg"  else "dr"
  use_bstrap_placebo <- if (use_retained_fallback) FALSE else TRUE

  set.seed(1242)
  cap_placebo <- run_att_gt_captured(
    yname = "log_commits", tname = "year", idname = "country_id", gname = "cohort_year",
    xformla = ~ ge_est + log_population, data = did_panel_placebo, est_method = use_dr_placebo,
    bstrap = use_bstrap_placebo, biters = 999L, cband = FALSE, control_group = "nevertreated",
    anticipation = 0, base_period = "universal", panel = TRUE, allow_unbalanced_panel = TRUE
  )
  gt_placebo <- cap_placebo$gt_obj

  agg_placebo <- if (!is.null(gt_placebo)) {
    tryCatch(aggte(gt_placebo, type = "simple", na.rm = TRUE), error = function(e) NULL)
  } else NULL

  # --- Everything below is DERIVED FROM THE FITTED OBJECT (every count and
  #     country name in the note), not pre-computed/hardcoded. unique(gt_placebo$
  #     group) is did's own record of which treated cohorts survived
  #     estimation (a cohort exactly at first_year is dropped internally,
  #     "already treated in first period" -- confirmed via did's own warning
  #     text and gt_obj$group in a controlled test); gt_placebo$n is did's
  #     own count of units used.
  if (is.null(gt_placebo) || is.null(agg_placebo)) {
    surviving_cohorts   <- integer(0)
    surviving_countries <- character(0)
    n_units_used        <- NA_integer_
    placebo_att <- placebo_se <- placebo_p <- NA_real_
    ci_lo <- ci_hi <- mde <- NA_real_
  } else {
    surviving_cohorts <- sort(unique(gt_placebo$group))
    surviving_countries <- napa_treated_names %>%
      filter(napa_cohort_year %in% surviving_cohorts) %>%
      pull(recipient_name)
    n_units_used <- gt_placebo$n
    placebo_att  <- agg_placebo$overall.att
    placebo_se   <- agg_placebo$overall.se
    placebo_p    <- 2 * pnorm(-abs(placebo_att / placebo_se))
    ci_lo <- placebo_att - 1.96 * placebo_se
    ci_hi <- placebo_att + 1.96 * placebo_se
    mde   <- 2.8 * placebo_se  # standard approx. MDE, 80% power at 5% significance
  }
  n_dropped_cohorts <- setdiff(unique(did_panel_placebo$cohort_year[did_panel_placebo$cohort_year > 0]),
                               surviving_cohorts)
  message(sprintf("Falsification-check ATT = %s, SE = %s, p = %s, 95%% CI = [%s, %s], MDE = %s.",
                  fmt4(placebo_att), fmt4(placebo_se), fmt3(placebo_p), fmt4(ci_lo), fmt4(ci_hi),
                  fmt4(mde)))
  message(sprintf("Surviving cohorts (from gt_obj$group): %s. Contributing countries: %s.",
                  if (length(surviving_cohorts) == 0) "(none)" else paste(surviving_cohorts, collapse = ", "),
                  if (length(surviving_countries) == 0) "(none)" else paste(surviving_countries, collapse = ", ")))

  # --- Power-diagnostic stats for the note: treated count, per-cohort sizes, post-periods per surviving
  #     cohort, and the composition of the control pool.
  n_treated_used <- length(surviving_countries)
  surviving_cohort_sizes <- napa_cohort_sizes %>%
    filter(cohort_year %in% surviving_cohorts) %>%
    mutate(post_periods = placebo_cutoff_year - cohort_year + 1L) %>%
    arrange(cohort_year)
  n_pre_first_year_napa <- sum(napa_cohort_lookup$napa_year < first_year, na.rm = TRUE)
  n_future_nap_in_control <- did_panel_placebo %>%
    filter(cohort_year == 0, !is.na(nap_year), nap_year > placebo_cutoff_year) %>%
    distinct(recipient_name) %>%
    nrow()
  message(sprintf(
    paste0("Power diagnostics: %d treated units in cohorts %s (sizes %s; post-periods %s); ",
          "%d NAPA countries predate %d and are pooled with never-treated; %d control-pool ",
          "countries go on to adopt a REAL NAP after %d (fine within this window, since none ",
          "of them are treated by it before %d)."),
    n_treated_used, paste(surviving_cohort_sizes$cohort_year, collapse = "/"),
    paste(surviving_cohort_sizes$n_treated, collapse = "/"),
    paste(surviving_cohort_sizes$post_periods, collapse = "/"),
    n_pre_first_year_napa, first_year, n_future_nap_in_control, placebo_cutoff_year,
    placebo_cutoff_year + 1L
  ))

  placebo_lines <- c(
    "\\begin{tabular}{lccccc}",
    "\\toprule",
    "Falsification treatment & ATT & SE & 95\\% CI & $p$-value & $N$ cohorts \\\\",
    "\\midrule",
    paste0("NAPA submission ($G_i$ = napa\\_year) & ", fmt4(placebo_att), stars_of(placebo_p),
           " & (", fmt4(placebo_se), ") & [", fmt4(ci_lo), ", ", fmt4(ci_hi), "] & ",
           fmt3(placebo_p), " & ", length(surviving_cohorts), " \\\\"),
    "\\bottomrule",
    "\\end{tabular}"
  )
  write_tex_float(
    out_path      = here("output", "tables", "napa", "att_napa_falsification.tex"),
    caption_title = "NAPA submission as a low-power falsification check",
    label         = "tab:napa_falsification",
    tabular_lines = placebo_lines,
    notes_text    = paste0(
      "CS\\,(2021), ", if (use_retained_fallback) "OR (analytical SE)"
      else "DR (multiplier-bootstrap SE, 999 reps, seed 1242)",
      " fallback, never-treated controls; low-power falsification, not a placebo (see main ",
      "text for the point estimate, cohort composition, and MDE). $G_i$ = NAPA-submission ",
      "year; panel restricted to ", first_year, "--", placebo_cutoff_year,
      ", verified free of real-NAP treatment; out-of-window NAPAs are pseudo-never-treated. ",
      "NAPA list: ", napa_source_note, ".", " \\texttt{did} ", DID_VERSION
    ),
    source_text = paste0(napa_source_short, "; OECD CRS (", CRS_VINTAGE, ")")
  )
}

message("\n=== 14_hazard_napa.R COMPLETE ===")
message("Outputs:")
message("  data/processed/napa_list.csv")
message("  data/processed/emdat_panel.csv (if EM-DAT present)")
message("  output/tables/hazard/att_hazard_controls.tex (if EM-DAT present)")
message("  output/tables/hazard/nap_timing_vs_humanitarian_aid.tex")
message("  output/tables/hazard/nap_timing_vs_emdat_hazard.tex (if EM-DAT present)")
message("  output/tables/napa/prior_napa_ldc_crosstab.tex")
message("  output/tables/napa/att_prior_napa_split.tex")
message("  output/tables/napa/att_napa_falsification.tex")
