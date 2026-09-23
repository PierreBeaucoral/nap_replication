##############################################################################
# 01_prepare_data.R
# Author : Pierre Beaucoral
# Date   : 2026-06-10
# Purpose: All data construction (plus three scope tables):
#            CRS → adaptation panel → NAP merge;
#            WDI/WGI merge, IMF adjustments, C/V orientation, derived vars
#              → simple_panel_wgi.csv;
#            descriptives panel → adaptationNAP.csv;
#            mitigation panel → mitigation_panel.csv
#
# Inputs : data/raw/CRS/CRS <year> Data.txt   (2007–2024)
#          data/raw/shared_nap_data/nap_information.csv
#          data/raw/PVCCI.csv
#          WDI API  (pinned via data/raw/wdi_cache/*.rds)
#
# Outputs (data/processed/):
#   adaptationNAP_donortype_wgi.csv  — donor-type panel with WGI/WDI controls
#   simple_panel_wgi.csv             — collapsed panel (one row per country-year)
#   adaptationNAP.csv                — descriptives panel
#   mitigation_panel.csv             — recipient-year mitigation commitments
#   donor_recipient_year_adaptation.csv — donor x recipient x year adaptation
#                                     commitments panel (§3c;
#                                     not used by the paper's exhibits)
#   adaptation_panel_oda_only.csv    — ODA-only variant of simple_panel_wgi.csv
#                                     (§22; ODA-only robustness input)
#   emergency_response_panel.csv     — recipient x year emergency-response AID
#                                     commitments (§23).
#                                     Used ONLY by 14_hazard_napa.R's NAP-
#                                     timing-vs-humanitarian-aid orthogonality
#                                     check (§4) -- NOT a hazard-realisation
#                                     measure and not used as a covariate in
#                                     any hazard-control specification (same
#                                     CRS, same donors, same recipient
#                                     envelopes as the outcome, so it would
#                                     be a bad control). The hazard-
#                                     control specifications use EM-DAT
#                                     (data/raw/emdat/emdat.csv), built by
#                                     14_hazard_napa.R itself, not by this file.
#   output/tables/scope/flow_type_shares.tex, regional_exclusion.tex,
#   sample_funnel.tex                — scope statements (§22)
##############################################################################

##############################################################################
# §0. PACKAGES (all library() calls at top)
##############################################################################

library(data.table)
library(dplyr)
library(tidyr)
library(countrycode)
library(lubridate)
library(ggplot2)
library(xtable)
library(WDI)
library(janitor)
library(broom)
library(stringr)
library(here)

# Single global seed at the top of the script
set.seed(20240601)

##############################################################################
# §0b. OUTPUT DIRECTORIES
##############################################################################

dir.create(here("data", "processed"),      recursive = TRUE, showWarnings = FALSE)
dir.create(here("data", "raw", "wdi_cache"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "figures"),      recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "tables"),       recursive = TRUE, showWarnings = FALSE)

##############################################################################
# §0c. HELPER: LaTeX-escape strings for column headers / table text
##############################################################################

esc_header <- function(x) {
  x <- gsub("%",  "\\\\%",  x)
  x <- gsub("_",  "\\\\_",  x)
  x <- gsub("#",  "\\\\#",  x)
  x
}

esc_tex <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([&%$#_{}])", "\\\\\\1", x)
  x
}

##############################################################################
# §0c-2. HELPER: RecipientName -> ISO3, with the two custom matches this
# project needs. Consolidates what used to be six separate
# `countrycode(...) + if_else(RecipientName == "Kosovo", ...)` call sites
# into one function.
#   - Kosovo: not in ISO 3166-1, CRS reports it as "Kosovo" -- XKX is the
#     provisional user-assigned code this project has always used for it.
#   - "Micronesia" (CRS RecipientCode 860): countrycode() cannot resolve this
#     unambiguously because "Micronesia" is also a subregion name. CRS code
#     860 is the Federated States of Micronesia specifically (verified against
#     the CRS RecipientName/RecipientCode crosswalk: 860 = "Micronesia", a
#     single country-level entry, distinct from 1034 = "Micronesia, regional"
#     in `regionalflows`) -- FSM, not the subregion.
##############################################################################

recipient_iso_from_name <- function(recipient_name) {
  # suppressWarnings(): "Kosovo" and "Micronesia" are EXPECTED to be
  # ambiguous to countrycode() -- that is exactly why the two explicit
  # overrides below exist; the warning would otherwise fire on every call.
  iso <- suppressWarnings(
    countrycode(recipient_name, origin = "country.name", destination = "iso3c")
  )
  iso[recipient_name %in% "Kosovo"]     <- "XKX"   # %in% is NA-safe
  iso[recipient_name %in% "Micronesia"] <- "FSM"
  iso
}

##############################################################################
# §0c-3. RECIPIENT CODES EXCLUDED AS REGIONAL/UNSPECIFIED (single definition;
# defined once here and referenced everywhere).
# 2026-09-15: removed code 860 ("Micronesia"). Verified against the
# CRS RecipientCode/RecipientName crosswalk (§0c-2 above): 860 is the
# Federated States of Micronesia, a genuine country-level recipient, distinct
# from 1034 ("Micronesia, regional"), which correctly stays in this list. 860
# was previously misclassified as regional and its adaptation-marked
# commitments ($92.0M, 0.02% of the 2009-2024 total per the prior
# regional_exclusion.tex) were dropped from the panel entirely; FSM now
# enters the panel as recipient 145 (see §22's effect-on-panel message below).
##############################################################################

regionalflows <- c(88, 89, 189, 237, 289, 298, 389, 489, 498, 589, 619,
                   679, 689, 789, 798, 889, 1027:1035, 9998)

##############################################################################
# §0d. WDI CACHE HELPER
# If a cached .rds exists, read it; otherwise pull from the API and cache.
# This pins the API data so reruns use the same World Bank vintage.
##############################################################################

wdi_cached <- function(cache_name, indicator, start, end, extra = FALSE,
                       wdi_cache = NULL) {
  path <- here("data", "raw", "wdi_cache", paste0(cache_name, ".rds"))
  if (file.exists(path)) {
    message("WDI cache hit: ", cache_name)
    return(readRDS(path))
  }
  message("WDI cache miss — pulling: ", cache_name)
  if (is.null(wdi_cache)) {
    dat <- WDI(country = "all", indicator = indicator,
               start = start, end = end, extra = extra)
  } else {
    dat <- WDI(country = "all", indicator = indicator,
               start = start, end = end, extra = extra, cache = wdi_cache)
  }
  saveRDS(dat, path)
  dat
}

##############################################################################
# §1. ANALYSIS PERIOD
##############################################################################

years        <- 2007:2024
period_label <- paste0(min(years), "-", max(years))
annual_years <- years

# Resolve a CRS file path tolerating the two capitalizations OECD has used
# ("CRS <year> Data.txt" vs "CRS <year> data.txt"). macOS file systems are
# case-insensitive so both work locally; on Linux (coauthor reruns, CI) the
# exact name matters.
crs_path <- function(y) {
  for (nm in c(paste0("CRS ", y, " Data.txt"), paste0("CRS ", y, " data.txt"))) {
    p <- here("data", "raw", "CRS", nm)
    if (file.exists(p)) return(p)
  }
  stop("CRS file for year ", y, " not found in data/raw/CRS/ ",
       "(looked for 'CRS ", y, " Data.txt' and 'CRS ", y, " data.txt').")
}

##############################################################################
# §2. LOAD RAW CRS DATA
##############################################################################

for (year in annual_years) {
  file_name <- paste0("CRS ", year, " Data")
  assign(file_name, fread(crs_path(year), encoding = "UTF-8"))
}

CRS <- lapply(annual_years, function(year) {
  get(paste0("CRS ", year, " Data"))
})
CRS <- lapply(CRS, as.data.frame)
gc()

##############################################################################
# §3. PROCESS CLIMATE ADAPTATION DATA
# Exhibit: feeds into all descriptive figures and regression panels
##############################################################################

process_rio_adaptation <- function(data_list) {
  processed_list <- lapply(data_list, function(df) {
    required_cols <- c("ClimateAdaptation", "Year", "DonorName", "DonorCode",
                       "RecipientName", "RecipientCode",
                       "USD_Commitment_Defl", "USD_Disbursement_Defl")
    if (!all(required_cols %in% names(df))) return(NULL)

    df %>%
      filter(ClimateAdaptation %in% c(1, 2)) %>%
      group_by(Year, DonorName, RecipientName, RecipientCode, DonorCode) %>%
      summarise(
        USD_Commitment_Defl  = sum(USD_Commitment_Defl,  na.rm = TRUE),
        USD_Disbursement_Defl = sum(USD_Disbursement_Defl, na.rm = TRUE),
        .groups = "drop"
      )
  })
  processed_list <- processed_list[!sapply(processed_list, is.null)]
  return(processed_list)
}

rio_data_adaptation_list <- process_rio_adaptation(CRS)

# Donor code lists
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

# Add donor type classification
rio_data_adaptation_list <- lapply(rio_data_adaptation_list, function(df) {
  mutate(df, DonorType = case_when(
    DonorCode %in% donor_list$DAC_members        ~ "Bilateral_members",
    DonorCode %in% donor_list$Multilateral_donors ~ "Multilateral_donors",
    TRUE                                          ~ "Other"
  ))
})

rio_data_adaptation <- bind_rows(rio_data_adaptation_list)

# Sort for reproducibility
df <- rio_data_adaptation[order(rio_data_adaptation$DonorType,
                                rio_data_adaptation$DonorName), ]
df_unique <- unique(df[, c("DonorType", "DonorName")])
df_unique <- df_unique[order(df_unique$DonorType, df_unique$DonorName), ]

# Persist the donor list so 02_descriptive_stats.R (list.tex) never has to
# touch raw CRS.
write.csv(df_unique, here("data", "processed", "donor_list.csv"),
          row.names = FALSE)
message("Wrote: data/processed/donor_list.csv  (", nrow(df_unique), " donors)")

# Persist total adaptation commitments by DonorName (+DonorType) so 02 §8
# (top_donors figure) never has to touch raw CRS on the happy path.
# Aggregation matches the computation in 02 §8's raw-CRS fallback branch.
donor_totals_01 <- rio_data_adaptation %>%
  group_by(DonorName, DonorType) %>%
  summarise(Total_Commitments = sum(USD_Commitment_Defl, na.rm = TRUE) / 1e3,
            .groups = "drop")
write.csv(donor_totals_01,
          here("data", "processed", "donor_totals.csv"),
          row.names = FALSE)
message("Wrote: data/processed/donor_totals.csv  (", nrow(donor_totals_01), " donor rows)")

##############################################################################
# §3c. DONOR x RECIPIENT x YEAR ADAPTATION PANEL
# Persists a donor x recipient x year adaptation-commitments panel derived
# from the SAME rio_data_adaptation frame already built in §3, so
# 06_reallocation_spillover.R (leave-one-out reallocation spillover test)
# never has to reload the raw CRS. Purely additive: does not touch or
# reorder any other section of this script.
#
# Recipient identity is standardized to ISO3 and regional/multi-country
# flows are excluded, using the single `regionalflows` constant (§0c-3 above)
# and the single `recipient_iso_from_name()` helper (§0c-2 above).
##############################################################################

donor_recipient_year_adaptation <- rio_data_adaptation %>%
  filter(!RecipientCode %in% regionalflows) %>%
  mutate(
    RecipientISO = recipient_iso_from_name(RecipientName)
  ) %>%
  filter(!is.na(RecipientISO)) %>%
  group_by(DonorCode, DonorName, DonorType,
           RecipientCode, RecipientName, RecipientISO, Year) %>%
  summarise(
    adaptation_commitments  = sum(USD_Commitment_Defl,  na.rm = TRUE),
    adaptation_disbursements = sum(USD_Disbursement_Defl, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(year = Year) %>%
  arrange(DonorCode, RecipientISO, year)

write.csv(donor_recipient_year_adaptation,
          here("data", "processed", "donor_recipient_year_adaptation.csv"),
          row.names = FALSE)
message("Wrote: data/processed/donor_recipient_year_adaptation.csv  (",
        nrow(donor_recipient_year_adaptation), " rows; ",
        n_distinct(donor_recipient_year_adaptation$DonorCode), " donors x ",
        n_distinct(donor_recipient_year_adaptation$RecipientISO), " recipients)")

##############################################################################
# §4. CREATE ADAPTATION AID PANEL — DONOR-TYPE LEVEL
# Exhibit: adaptationNAP_donortype_wgi.csv
# regionalflows and the ISO3 helper are both defined once in §0c-2/§0c-3.
##############################################################################

adaptation_aid <- rio_data_adaptation %>%
  filter(!RecipientCode %in% regionalflows) %>%
  group_by(RecipientName, Year, DonorType) %>%
  summarise(
    Commitments  = sum(USD_Commitment_Defl,  na.rm = TRUE),
    Disbursements = sum(USD_Disbursement_Defl, na.rm = TRUE),
    .groups = "drop"
  )

adaptation_aid$RecipientISO <- recipient_iso_from_name(adaptation_aid$RecipientName)

adaptation_aid_panel <- expand.grid(
  RecipientName = unique(adaptation_aid$RecipientName),
  Year          = min(adaptation_aid$Year):max(adaptation_aid$Year),
  DonorType     = unique(adaptation_aid$DonorType)
) %>%
  left_join(adaptation_aid, by = c("RecipientName", "Year", "DonorType")) %>%
  mutate(
    Commitments   = coalesce(Commitments,   0),
    Disbursements = coalesce(Disbursements, 0)
  ) %>%
  arrange(RecipientName, Year, DonorType)

# Panel balance check
panel_check <- adaptation_aid_panel %>%
  group_by(RecipientName, DonorType) %>%
  summarise(
    n_years     = n(),
    min_year    = min(Year),
    max_year    = max(Year),
    is_balanced = n_years == (max(Year) - min(Year) + 1),
    .groups = "drop"
  )
message("Number of countries: ", length(unique(adaptation_aid_panel$RecipientName)))
message("Number of years: ",      length(unique(adaptation_aid_panel$Year)))
message("Total observations: ",   nrow(adaptation_aid_panel))
message("Is panel balanced? ",    all(panel_check$is_balanced))

##############################################################################
# §5. MERGE NAP DATA
##############################################################################

nap_data <- fread(here("data", "raw", "shared_nap_data", "nap_information.csv"))
nap_data$RecipientISO <- countrycode(
  nap_data$Country, origin = "country.name", destination = "iso3c"
)
adaptation_aid_panel <- left_join(adaptation_aid_panel, nap_data, by = "RecipientISO")

##############################################################################
# §6. DATE VARIABLES & TREATMENT TIMING
##############################################################################

adaptation_aid_panel <- adaptation_aid_panel %>%
  mutate(
    date_posted = parse_date_time(`Date Posted`, orders = c("dmy", "mdy", "ymd")),
    NAP_Year    = year(date_posted)
  )

adaptation_aid_panel <- adaptation_aid_panel %>%
  mutate(
    time_to_nap = if_else(
      !is.na(NAP_Year) & NAP_Year > 0, Year - NAP_Year, NA_real_
    )
  )

##############################################################################
# §7. MERGE WDI POPULATION & GDP
# WDI pulls are cached in data/raw/wdi_cache/
##############################################################################

pop_data <- wdi_cached(
  "wdi_population",
  indicator = "SP.POP.TOTL",
  start = min(adaptation_aid_panel$Year),
  end   = max(adaptation_aid_panel$Year)
)
pop_clean <- pop_data %>%
  select(iso3c, year, SP.POP.TOTL) %>%
  rename(RecipientISO = iso3c, Year = year, Population = SP.POP.TOTL) %>%
  filter(!is.na(Population))

gdp_data <- wdi_cached(
  "wdi_gdp",
  indicator = "NY.GDP.MKTP.KD",
  start = min(adaptation_aid_panel$Year),
  end   = max(adaptation_aid_panel$Year)
)
gdp_clean <- gdp_data %>%
  select(iso3c, year, NY.GDP.MKTP.KD) %>%
  rename(RecipientISO = iso3c, Year = year, GDP = NY.GDP.MKTP.KD) %>%
  filter(!is.na(GDP))

adaptation_aid_panel$RecipientISO <- recipient_iso_from_name(adaptation_aid_panel$RecipientName)

adaptation_aid_panel <- adaptation_aid_panel %>%
  left_join(pop_clean, by = c("RecipientISO", "Year")) %>%
  left_join(gdp_clean, by = c("RecipientISO", "Year")) %>%
  mutate(
    Commitments_pc  = if_else(Population > 0, Commitments  / Population, NA_real_),
    Disbursements_pc = if_else(Population > 0, Disbursements / Population, NA_real_)
  )

##############################################################################
# §8. MERGE WGI GOVERNMENT EFFECTIVENESS
##############################################################################

wgi_ge_wdi <- wdi_cached(
  "wdi_wgi_ge",
  indicator = "GOV_WGI_GE.EST",
  start = min(adaptation_aid_panel$Year, na.rm = TRUE),
  end   = max(adaptation_aid_panel$Year, na.rm = TRUE)
)

wgi_ge_clean <- wgi_ge_wdi %>%
  select(iso3c, year, GOV_WGI_GE.EST) %>%
  rename(RecipientISO = iso3c, Year = year) %>%
  mutate(ge_est = GOV_WGI_GE.EST) %>%
  select(RecipientISO, Year, ge_est) %>%
  filter(!is.na(ge_est))

adaptation_aid_panel <- adaptation_aid_panel %>%
  left_join(wgi_ge_clean, by = c("RecipientISO", "Year"))

##############################################################################
# §9. IMF GDP ADJUSTMENTS
##############################################################################

adaptation_aid_panel <- adaptation_aid_panel %>%
  mutate(GDP = if_else(RecipientName == "Bhutan" & Year == 2023, 2.87e+09, GDP))

imf_adjustments_djibouti <- data.frame(
  RecipientName = rep("Djibouti", 4),
  Year  = c(2009, 2010, 2011, 2012),
  value = c("1,43", "1,54", "1,74", "1,74"),
  stringsAsFactors = FALSE
)
imf_adjustments_eritrea <- data.frame(
  RecipientName = rep("Eritrea", 12),
  Year  = 2012:2023,
  value = c("2,255", "1,958", "2,604", "2,016", "2,213", "1,904",
            "2,006", "1,982", NA, NA, NA, NA),
  stringsAsFactors = FALSE
)
imf_adjustments_lebanon <- data.frame(
  RecipientName = "Lebanon",
  Year  = 2023,
  value = "24,023",
  stringsAsFactors = FALSE
)
imf_adjustments_southsudan <- data.frame(
  RecipientName = rep("South Sudan", 8),
  Year  = c(2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023),
  value = c("12,42673", "13,21288", "14,32292", "14,03651",
            "12,88751", "14,62216", "14,47636", "16"),
  stringsAsFactors = FALSE
)
imf_adjustments_sar <- data.frame(
  RecipientName = "Syrian Arab Republic",
  Year  = 2023,
  value = "16",
  stringsAsFactors = FALSE
)
imf_adjustments_tonga <- data.frame(
  RecipientName = "Tonga",
  Year  = 2023,
  value = "0,52",
  stringsAsFactors = FALSE
)
imf_adjustments_venezuela <- data.frame(
  RecipientName = rep("Venezuela", 15),
  Year  = 2009:2023,
  value = c(
    "268,624", "318,281", "316,482", "372,592", "258,931", "214,69",
    "125,449", "112,915", "115,883", "102,021", "73,003", "43,788",
    "57,666", "92,104", "99,203"
  ),
  stringsAsFactors = FALSE
)
imf_adjustments_yemen <- data.frame(
  RecipientName = rep("Yemen", 4),
  Year  = 2019:2022,
  value = c("21,888", "20,22", "19,394", "23,534"),
  stringsAsFactors = FALSE
)

imf_adjustments <- rbind(
  imf_adjustments_djibouti, imf_adjustments_eritrea, imf_adjustments_lebanon,
  imf_adjustments_southsudan, imf_adjustments_sar, imf_adjustments_tonga,
  imf_adjustments_venezuela, imf_adjustments_yemen
)
imf_adjustments <- imf_adjustments %>%
  mutate(
    clean_value = as.numeric(gsub(",", ".", value)),
    GDP         = clean_value * 1e9
  )

adaptation_aid_panel <- adaptation_aid_panel %>%
  left_join(
    imf_adjustments %>% select(RecipientName, Year, GDP),
    by = c("RecipientName", "Year"), suffix = c("", "_new")
  ) %>%
  mutate(GDP = if_else(!is.na(GDP_new), GDP_new, GDP)) %>%
  select(-GDP_new)

##############################################################################
# §10. C- vs V-ORIENTED DONORS (donor sensitivities)
##############################################################################

yvar    <- if ("Commitments_pc" %in% names(adaptation_aid_panel)) "Commitments_pc" else "Commitments"
cap_var <- "ge_est"

stopifnot("ge_est" %in% names(adaptation_aid_panel))

##############################################################################
# §10a. LOAD & MERGE PVCCI
##############################################################################

df_pvcci <- suppressWarnings(
  fread(here("data", "raw", "PVCCI.csv"), encoding = "UTF-8")
)
if (!all(c("recipient_name", "PVCCI") %in% names(df_pvcci))) {
  stop("PVCCI.csv must contain columns 'recipient_name' and 'PVCCI'.")
}

# recipient_iso_from_name() (this site previously duplicated a Kosovo-only countrycode() call instead of using the
# shared helper -- §0c-2 -- which also carries the FSM custom match; PVCCI.csv
# does not contain a Micronesia/FSM row, so FSM correctly gets an honest NA
# PVCCI value either way, but the merge now goes through one canonical path).
pvcci_iso <- df_pvcci %>%
  mutate(recipientiso = recipient_iso_from_name(recipient_name)) %>%
  filter(!is.na(recipientiso) & !is.na(PVCCI)) %>%
  group_by(recipientiso) %>%
  summarise(PVCCI = mean(PVCCI, na.rm = TRUE), .groups = "drop")

adaptation_aid_panel <- adaptation_aid_panel %>%
  mutate(
    RecipientISO = if_else(
      is.na(RecipientISO),
      recipient_iso_from_name(RecipientName),
      RecipientISO
    )
  ) %>%
  left_join(pvcci_iso, by = c("RecipientISO" = "recipientiso")) %>%
  mutate(PVCCI_std = if (all(is.na(PVCCI))) NA_real_ else as.numeric(scale(PVCCI)))

vul_var <- "PVCCI_std"

message("PVCCI merged for ", sum(!is.na(adaptation_aid_panel$PVCCI)),
        " / ", nrow(adaptation_aid_panel), " rows; using ", vul_var, " as vulnerability.")

##############################################################################
# §10b. ESTIMATE DONOR SENSITIVITIES
##############################################################################

adapt_pre <- adaptation_aid_panel %>%
  filter(!is.na(NAP_Year), Year < NAP_Year) %>%
  filter(
    !is.na(.data[[yvar]]),
    !is.na(.data[[cap_var]]),
    !is.na(.data[[vul_var]]),
    !is.na(RecipientISO),
    !is.na(DonorType)
  )

donor_sensitivities <- adapt_pre %>%
  group_by(DonorType) %>%
  do({
    fm  <- reformulate(
      c(cap_var, vul_var, "factor(RecipientISO)", "factor(Year)"),
      response = yvar
    )
    mod <- lm(fm, data = .)
    tidy(mod)
  }) %>%
  ungroup() %>%
  filter(term %in% c(cap_var, vul_var)) %>%
  select(DonorType, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate, names_prefix = "phi_")

donor_sensitivities <- donor_sensitivities %>%
  mutate(
    C_score    = rank(phi_ge_est,    ties.method = "average"),
    V_score    = rank(phi_PVCCI_std, ties.method = "average"),
    C_oriented = if_else(C_score > V_score, 1L, 0L)
  ) %>%
  rename(
    phi_capacity      = phi_ge_est,
    phi_vulnerability = phi_PVCCI_std
  )

adaptation_aid_panel <- adaptation_aid_panel %>%
  left_join(
    donor_sensitivities %>%
      select(DonorType, phi_capacity, phi_vulnerability, C_score, V_score, C_oriented),
    by = "DonorType"
  )

message("C_oriented distribution:")
print(table(adaptation_aid_panel$C_oriented, useNA = "ifany"))

##############################################################################
# §11. EXPORT DONOR-TYPE PANEL
# (only adaptationNAP_donortype_wgi.csv is written; an Africa-only panel is
#  not built because nothing downstream uses it)
##############################################################################

write.csv(
  adaptation_aid_panel,
  here("data", "processed", "adaptationNAP_donortype_wgi.csv"),
  row.names = FALSE
)
message("Wrote: data/processed/adaptationNAP_donortype_wgi.csv  (",
        nrow(adaptation_aid_panel), " rows)")

##############################################################################
# §12. FURTHER DATA PREPARATION — SIMPLE BALANCED PANEL
# Output: simple_panel_wgi.csv (the estimation panel used by stages 03-14)
##############################################################################

data <- fread(here("data", "processed", "adaptationNAP_donortype_wgi.csv"))
data <- clean_names(data)

# Donor-type-specific adaptation commitments (heterogeneity analysis)
donor_commits <- data %>%
  group_by(recipient_name, year, donor_type) %>%
  summarise(commits_dt = sum(commitments, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = donor_type, values_from = commits_dt, values_fill = 0) %>%
  mutate(
    commits_dac   = if ("Bilateral_members"   %in% names(.)) Bilateral_members   else 0,
    commits_multi = if ("Multilateral_donors" %in% names(.)) Multilateral_donors else 0,
    commits_other = if ("Other"               %in% names(.)) Other               else 0,
    log_commits_dac   = log1p(commits_dac),
    log_commits_multi = log1p(commits_multi),
    log_commits_other = log1p(commits_other)
  ) %>%
  select(recipient_name, year,
         commits_dac, commits_multi, commits_other,
         log_commits_dac, log_commits_multi, log_commits_other)

vars_to_include <- names(data)[!(names(data) %in% c("donor_type", "lcommitments_pc"))]
df_keep <- data %>% select(any_of(vars_to_include))

aggregated <- df_keep %>%
  group_by(recipient_name, year) %>%
  summarise(
    commitments    = sum(commitments,  na.rm = TRUE),
    disbursements  = sum(disbursements, na.rm = TRUE),
    lcommitments   = log1p(commitments),
    ldisbursements = log1p(disbursements),
    across(-c(commitments, disbursements, lcommitments, ldisbursements), first),
    .groups = "drop"
  )

aggregated <- aggregated %>%
  left_join(donor_commits, by = c("recipient_name", "year"))
message("Donor-type columns added: log_commits_dac, log_commits_multi, log_commits_other")

# Generate dynamic NAP treatment dummies
for (i in 1:5) {
  varname   <- paste0("NAP_minus", i)
  aggregated <- aggregated %>%
    mutate(
      !!varname := if_else(
        !is.na(nap_year) & (year >= (nap_year - i)) & (year < nap_year), 1, 0
      )
    )
}
aggregated <- aggregated %>%
  mutate(
    NAP_0  = if_else(!is.na(nap_year) & year >= nap_year,       1, 0),
    NAP_01 = if_else(!is.na(nap_year) & year >= (nap_year + 1), 1, 0)
  )

aggregated <- aggregated %>%
  group_by(recipient_name) %>%
  mutate(
    treated = if_else(!is.na(nap_year) & nap_year >= 2013 & nap_year <= 2024, 1, 0)
  ) %>%
  ungroup()

setDT(aggregated)
aggregated[, log_commits := log1p(commitments)]
aggregated[, year_factor      := as.factor(year)]
aggregated[, recipient_factor := as.factor(recipient_iso)]

y_col        <- "log_commits"
covariates   <- c("population", "gdp")
cluster_cols <- c("recipient_factor", "year_factor")

all_columns  <- c(y_col, covariates, cluster_cols)
missing_cols <- setdiff(all_columns, names(aggregated))
if (length(missing_cols) > 0) {
  stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
}

date_cols <- names(aggregated)[sapply(aggregated, function(x) inherits(x, "IDate"))]
aggregated[, (date_cols) := lapply(.SD, as.Date),      .SDcols = date_cols]
aggregated[, (date_cols) := lapply(.SD, as.character), .SDcols = date_cols]

##############################################################################
# §13. ADDITIONAL WDI INDICATORS (WGI governance + macro/social)
##############################################################################

wdi_cache_obj <- WDIcache()

wgi_codes <- c(
  "GOV_WGI_GE.EST",
  "GOV_WGI_CC.EST",
  "GOV_WGI_RL.EST",
  "GOV_WGI_RQ.EST",
  "GOV_WGI_PV.EST",
  "GOV_WGI_VA.EST"
)

other_codes <- c(
  "NY.GDP.MKTP.KD.ZG",
  "SL.UEM.TOTL.ZS",
  "SP.DYN.LE00.IN",
  "SH.XPD.CHEX.PC.CD",
  "SE.PRM.ENRR",
  "SE.SEC.ENRR"
)

all_codes <- unique(c(wgi_codes, other_codes))

wdi_raw <- wdi_cached(
  "wdi_additional",
  indicator = all_codes,
  start     = 2009,
  end       = 2024,
  extra     = FALSE,
  wdi_cache = wdi_cache_obj
)

wdi_tidy <- wdi_raw %>%
  rename(
    recipient_name = country,
    iso2c          = iso2c,
    year           = year
  )

aggregated <- aggregated %>%
  left_join(wdi_tidy, by = c("recipient_name", "year"))

# --- PVCCI merge ---------------------------------------------------------------
# The panel already carries a correctly ISO3-merged `pvcci` column: §10a above
# builds `pvcci_iso` via countrycode(..., "country.name", "iso3c") and joins
# it onto adaptation_aid_panel by RecipientISO (57/57 PVCCI.csv countries
# matched — full coverage of the source index, which covers 57 of the panel's
# 144 recipients by construction). That column survives §11's write/read
# round-trip and §12's clean_names() as lowercase `pvcci`.
#
# PVCCI is NOT re-merged by RecipientName: a name match covers only 51/57
# countries (string mismatches on China, "Cote d'Ivoire", "Congo, Dem. Rep." /
# "Congo, Rep", "Cape Verde" and "Swaziland"). §17's PVCCI interaction terms
# (PVCCI_GE, PVCCI_sq, PVCCI_th) and 02's balance table read the ISO3-merged
# `pvcci` column. 02 prints the PVCCI merge-coverage audit.
message("PVCCI: using the ISO3-merged `pvcci` column from Section 10a (",
        sum(!is.na(aggregated$pvcci)), " / ", nrow(aggregated),
        " row-obs non-missing; deleted the redundant name-based re-merge).")

##############################################################################
# §14. MISSING VALUE AUDIT
# (audit messages only; missing values are handled by listwise deletion
#  below, no imputation)
##############################################################################

setDT(aggregated)

all_num  <- names(aggregated)[vapply(aggregated, is.numeric, logical(1))]
to_drop  <- c(
  "recipient_name", "year", "commitments", "disbursements",
  "lcommitments", "ldisbursements", "recipient_iso", "no",
  "country", "region", "ldc_sids", "date_posted", "language_1",
  "language_2", "date_posted_2", "nap_year", "time_to_nap",
  "commitments_pc", "disbursements_pc",
  "NAP_minus1", "NAP_minus2", "NAP_minus3", "NAP_minus4",
  "NAP_minus5", "NAP_0", "NAP_01", "treated", "treatment",
  "log_commits", "iso2c", "iso3c"
)
num_cols <- setdiff(all_num, to_drop)
z_data   <- aggregated[, ..num_cols]

message("\n=== MISSING VALUE AUDIT (pre-imputation) ===\n")

n_obs     <- nrow(z_data)
n_cols    <- ncol(z_data)
total_cells   <- n_obs * n_cols
total_missing <- sum(is.na(z_data))

message(sprintf("Dataset: %d rows x %d columns (%d cells total)", n_obs, n_cols, total_cells))
message(sprintf("Total missing cells: %d (%.1f%%)",
                total_missing, 100 * total_missing / total_cells))

miss_by_col <- data.frame(
  variable    = names(z_data),
  n_missing   = sapply(z_data, function(x) sum(is.na(x))),
  pct_missing = round(100 * sapply(z_data, function(x) mean(is.na(x))), 1),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(n_missing))

message("\n--- Missingness by variable (sorted) ---")
print(miss_by_col, row.names = FALSE)

miss_by_country <- aggregated %>%
  select(recipient_name, all_of(num_cols)) %>%
  group_by(recipient_name) %>%
  summarise(
    n_rows      = n(),
    n_missing   = sum(is.na(across(all_of(num_cols)))),
    pct_missing = round(100 * n_missing / (n_rows * length(num_cols)), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_missing))

message("\n--- Missingness by country (sorted, top 30) ---")
print(head(miss_by_country, 30), row.names = FALSE)

miss_by_year <- aggregated %>%
  select(year, all_of(num_cols)) %>%
  group_by(year) %>%
  summarise(
    n_rows      = n(),
    n_missing   = sum(is.na(across(all_of(num_cols)))),
    pct_missing = round(100 * n_missing / (n_rows * length(num_cols)), 1),
    .groups = "drop"
  ) %>%
  arrange(year)

message("\n--- Missingness by year ---")
print(miss_by_year, row.names = FALSE)

complete_vars <- miss_by_col$variable[miss_by_col$n_missing == 0]
message(sprintf("\n%d variable(s) are complete (no imputation needed):",
                length(complete_vars)))
if (length(complete_vars) > 0) message(paste(complete_vars, collapse = ", "))

high_miss_vars <- miss_by_col$variable[miss_by_col$pct_missing > 50]
message(sprintf("\n%d variable(s) exceed 50%% missing (imputation may be unreliable):",
                length(high_miss_vars)))
if (length(high_miss_vars) > 0) message(paste(high_miss_vars, collapse = ", "))

message("\n--- Decision guidance ---")
message("If total missingness is low (<5%) and no variable exceeds ~20%,")
message("consider complete-case analysis instead of missForest.")
message("missForest is most valuable when missingness is scattered and MCAR/MAR.")
message("=== END OF MISSING VALUE AUDIT ===\n")

##############################################################################
# §15. LISTWISE DELETION ON CONTROLS
##############################################################################

n_before   <- nrow(aggregated)
aggregated <- aggregated[!is.na(ge_est) & !is.na(population)]
n_after    <- nrow(aggregated)
message(sprintf(
  "Listwise deletion on ge_est + population: dropped %d rows (%.1f%%), kept %d.",
  n_before - n_after,
  100 * (n_before - n_after) / n_before,
  n_after
))

aggregated[, (cluster_cols) := lapply(.SD, as.factor), .SDcols = cluster_cols]

##############################################################################
# §16. COUNTRY-LEVEL VARIABLES: CONTINENT, WB REGION
##############################################################################

aggregated$continent <- countrycode(
  sourcevar = aggregated$iso3c, origin = "iso3c", destination = "continent"
)
aggregated$WB_region <- countrycode(
  sourcevar = aggregated$iso3c, origin = "iso3c", destination = "region"
)
aggregated$continent[aggregated$iso3c == "XKX"] <- "Europe"
aggregated$WB_region[aggregated$iso3c == "XKX"] <- "Europe & Central Asia"

##############################################################################
# §17. DERIVED VARIABLES
# 2026-09-15: PVCCI_GE/PVCCI_sq/PVCCI_th now read the ISO3-merged
# `pvcci` column (57/144 coverage) instead of the deleted name-based `PVCCI`
# column (51/144 coverage) -- see the §13 comment above.
##############################################################################

aggregated[, `:=`(
  PVCCI_GE       = pvcci * ge_est,
  PVCCI_sq       = pvcci^2,
  GE_sq          = ge_est^2,
  PVCCI_th       = pvcci^3,
  GE_th          = ge_est^3,
  population_gdp = population * gdp,
  population_sq  = population^2,
  gdp_sq         = gdp^2,
  population_th  = population^3,
  gdp_th         = gdp^3
)]

z_col <- c("pvcci", "ge_est", "PVCCI_GE", "PVCCI_sq", "GE_sq", "PVCCI_th", "GE_th")
x_col <- c("population", "gdp", "population_gdp", "population_sq", "gdp_sq",
           "population_th", "gdp_th", z_col)

message("Balanced panel dimensions: ", nrow(aggregated), " observations.")

##############################################################################
# §18. ALL-COMMITMENTS & PRINCIPAL-ONLY & FLOW-TYPE SUPPLEMENTS
# Builds: commitments_all, lcommitments_nonadapt, lcommitments_principal,
#         commitments_oda/oof/private, share_adapt
##############################################################################

all_commitments_ry <- lapply(CRS, function(df) {
  df %>%
    filter(!RecipientCode %in% regionalflows) %>%
    group_by(RecipientName, Year) %>%
    summarise(commitments_all = sum(USD_Commitment_Defl, na.rm = TRUE),
              .groups = "drop")
}) %>%
  bind_rows() %>%
  rename(recipient_name = RecipientName, year = Year) %>%
  group_by(recipient_name, year) %>%
  summarise(commitments_all = sum(commitments_all, na.rm = TRUE), .groups = "drop")

aggregated <- aggregated %>%
  left_join(all_commitments_ry, by = c("recipient_name", "year")) %>%
  mutate(
    commitments_all       = if_else(is.na(commitments_all), 0, commitments_all),
    commitments_nonadapt  = pmax(commitments_all - commitments, 0),
    lcommitments_all      = log1p(commitments_all),
    lcommitments_nonadapt = log1p(commitments_nonadapt)
  )

stopifnot(all(aggregated$commitments_all >= aggregated$commitments))

# Principal-only adaptation finance (Rio marker = 2)
principal_ry <- lapply(CRS, function(df) {
  df %>%
    filter(ClimateAdaptation == 2, !RecipientCode %in% regionalflows) %>%
    group_by(RecipientName, Year) %>%
    summarise(commitments_principal = sum(USD_Commitment_Defl, na.rm = TRUE),
              .groups = "drop")
}) %>%
  bind_rows() %>%
  rename(recipient_name = RecipientName, year = Year) %>%
  group_by(recipient_name, year) %>%
  summarise(commitments_principal = sum(commitments_principal, na.rm = TRUE),
            .groups = "drop")

aggregated <- aggregated %>%
  left_join(principal_ry, by = c("recipient_name", "year")) %>%
  mutate(
    commitments_principal  = if_else(is.na(commitments_principal), 0, commitments_principal),
    lcommitments_principal = log1p(commitments_principal)
  )

# Flow type: ODA / OOF / private. Aggregates BOTH commitments and
# disbursements by flow type in the same pass (the ODA-only panel below needs disbursements_oda too, so it can
# rebuild disbursements/ldisbursements consistently with the ODA-only
# commitments figure instead of dropping or leaving them stale).
flow_ry <- lapply(CRS, function(df) {
  df %>%
    filter(ClimateAdaptation %in% c(1, 2), !RecipientCode %in% regionalflows) %>%
    mutate(flow_grp = case_when(
      grepl("^ODA", FlowName)               ~ "oda",
      grepl("Other Official", FlowName)     ~ "oof",
      grepl("Private|Equity|PSI", FlowName) ~ "private",
      TRUE                                  ~ "other"
    )) %>%
    group_by(RecipientName, Year, flow_grp) %>%
    summarise(v = sum(USD_Commitment_Defl,  na.rm = TRUE),
             d = sum(USD_Disbursement_Defl, na.rm = TRUE),
             .groups = "drop")
}) %>%
  bind_rows() %>%
  rename(recipient_name = RecipientName, year = Year) %>%
  group_by(recipient_name, year, flow_grp) %>%
  summarise(v = sum(v, na.rm = TRUE), d = sum(d, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = flow_grp, values_from = c(v, d),
                     values_fill = 0, names_glue = "{.value}_{flow_grp}") %>%
  rename_with(~ sub("^v_", "commitments_", .x), starts_with("v_")) %>%
  rename_with(~ sub("^d_", "disbursements_", .x), starts_with("d_"))

aggregated <- aggregated %>%
  left_join(flow_ry, by = c("recipient_name", "year")) %>%
  mutate(
    commitments_oda      = if_else(is.na(commitments_oda),     0, commitments_oda),
    commitments_oof      = if_else(is.na(commitments_oof),     0, commitments_oof),
    commitments_private  = if_else(is.na(commitments_private), 0, commitments_private),
    lcommitments_oda     = log1p(commitments_oda),
    lcommitments_oof     = log1p(commitments_oof),
    lcommitments_private = log1p(commitments_private),
    disbursements_oda     = if_else(is.na(disbursements_oda),     0, disbursements_oda),
    disbursements_oof     = if_else(is.na(disbursements_oof),     0, disbursements_oof),
    disbursements_private = if_else(is.na(disbursements_private), 0, disbursements_private)
  )

# Adaptation share of global total
setDT(aggregated)
aggregated[, global_adapt_t := sum(commitments, na.rm = TRUE), by = year]
aggregated[, share_adapt    := fifelse(
  global_adapt_t > 0, commitments / global_adapt_t * 100, NA_real_
)]
aggregated[, global_adapt_t := NULL]

##############################################################################
# §19. EXPORT simple_panel_wgi.csv
##############################################################################

write.csv(aggregated, here("data", "processed", "simple_panel_wgi.csv"),
          row.names = FALSE)
message("Wrote: data/processed/simple_panel_wgi.csv  (", nrow(aggregated), " rows)")

##############################################################################
# §20. DESCRIPTIVES PANEL (adaptationNAP.csv)
# This is the panel used by 02_descriptive_stats.R: no WGI/WDI, no donor-type
# split, just recipient x year aggregation for descriptive tables and figures.
# Exhibit: 02_descriptive_stats.R (all 6 figures + stats_des + nap_regional
#          + finance_change + list.tex)
##############################################################################

# Build CRS extract directly from per-year files
crs_years <- 2007:2024
crs_cols  <- c("ClimateAdaptation", "Year", "DonorName", "DonorCode",
               "RecipientName", "RecipientCode",
               "USD_Commitment_Defl", "USD_Disbursement_Defl")

rio_data_04 <- rbindlist(lapply(crs_years, function(y) {
  dt <- fread(crs_path(y), encoding = "UTF-8")
  dt[, intersect(crs_cols, names(dt)), with = FALSE]
}), fill = TRUE)

# Filter and aggregate
rio_data_adaptation_04 <- rio_data_04 %>%
  filter(ClimateAdaptation %in% c(1, 2)) %>%
  group_by(Year, DonorName, RecipientName, RecipientCode) %>%
  summarise(
    USD_Commitment_Defl   = sum(USD_Commitment_Defl,   na.rm = TRUE),
    USD_Disbursement_Defl = sum(USD_Disbursement_Defl, na.rm = TRUE),
    .groups = "drop"
  )

# Exclude regional flows; single `regionalflows`
# constant (§0c-3) and `recipient_iso_from_name()` helper (§0c-2).
adaptation_aid_04 <- rio_data_adaptation_04 %>%
  filter(!RecipientCode %in% regionalflows) %>%
  group_by(RecipientName, Year) %>%
  summarise(
    Commitments   = sum(USD_Commitment_Defl,   na.rm = TRUE),
    Disbursements = sum(USD_Disbursement_Defl, na.rm = TRUE),
    .groups = "drop"
  )

adaptation_aid_04$RecipientISO <- recipient_iso_from_name(adaptation_aid_04$RecipientName)

# Merge NAP data
nap_data_04 <- fread(here("data", "raw", "shared_nap_data", "nap_information.csv"))
nap_data_04$RecipientISO <- countrycode(
  nap_data_04$Country, origin = "country.name", destination = "iso3c"
)
adaptation_aid_04 <- left_join(adaptation_aid_04, nap_data_04, by = "RecipientISO")

# Add date and NAP year
adaptation_aid_04 <- adaptation_aid_04 %>%
  mutate(
    date_posted = parse_date_time(`Date Posted`, orders = c("dmy", "mdy", "ymd")),
    NAP_Year    = year(date_posted)
  )

# NAP treatment variable for finance_change table
adaptation_aid_04 <- adaptation_aid_04 %>%
  mutate(
    NAP = case_when(
      is.na(NAP_Year) ~ 0,
      Year >= NAP_Year ~ 1,
      TRUE             ~ 0
    ),
    NAP_Year = case_when(
      is.na(NAP_Year) ~ 0,
      TRUE            ~ NAP_Year
    )
  )

# Drop columns that don't exist across all year-files before write
cols_to_drop <- intersect(c("No.", "Country", "Region", "LDC/SIDS"), names(adaptation_aid_04))
if (length(cols_to_drop) > 0) {
  adaptation_aid_04 <- adaptation_aid_04 %>% select(-all_of(cols_to_drop))
}

write.csv(adaptation_aid_04, here("data", "processed", "adaptationNAP.csv"),
          row.names = FALSE)
message("Wrote: data/processed/adaptationNAP.csv  (", nrow(adaptation_aid_04), " rows)")

##############################################################################
# §21. MITIGATION PANEL (mitigation_panel.csv)
# Data build only; the mitigation falsification test itself is estimated in
# 04_robustness.R.
##############################################################################

message("\n=== SECTION 21: BUILD MITIGATION PANEL ===")

mitigation_ry <- lapply(CRS, function(df) {
  required_cols <- c("ClimateMitigation", "Year",
                     "RecipientName", "RecipientCode",
                     "USD_Commitment_Defl")
  if (!all(required_cols %in% names(df))) return(NULL)
  df %>%
    filter(ClimateMitigation %in% c(1, 2)) %>%
    group_by(Year, RecipientName) %>%
    summarise(commits_mitigation = sum(USD_Commitment_Defl, na.rm = TRUE),
              .groups = "drop")
}) %>%
  bind_rows() %>%
  rename(recipient_name = RecipientName, year = Year) %>%
  group_by(recipient_name, year) %>%
  summarise(commits_mitigation = sum(commits_mitigation, na.rm = TRUE),
            .groups = "drop")

message(sprintf("Mitigation panel: %d recipient-year rows (%d unique recipients)",
                nrow(mitigation_ry),
                n_distinct(mitigation_ry$recipient_name)))

write.csv(mitigation_ry, here("data", "processed", "mitigation_panel.csv"),
          row.names = FALSE)
message("Wrote: data/processed/mitigation_panel.csv  (", nrow(mitigation_ry), " rows)")

##############################################################################
# §22. SCOPE STATEMENTS
# Documents (a) which CRS flow types the
# adaptation panel includes (ODA-only vs. ODA + OOF + private/equity) and
# (b) which recipient codes the regional-flow exclusion drops and what share
# of adaptation-marked commitments they represent. Does NOT change the main panel (03/04/05 depend on `adaptation_aid`
# / `aggregated` exactly as built above) — exports scope tables plus an
# ODA-only panel variant for a future robustness run. Window: 2009-2024,
# matching the main estimation sample (2009 is the first full year of Rio
# adaptation-marker reporting; see §2's stats_des note in 02).
##############################################################################

message("\n=== SECTION 22: SCOPE STATEMENTS ===")

dir.create(here("output", "tables", "scope"), recursive = TRUE, showWarnings = FALSE)

fmt2 <- function(x) formatC(x, format = "f", digits = 2, big.mark = ",")

## --- Flow-type composition of the adaptation panel ------------------------
# Same filter as adaptation_aid (§4): ClimateAdaptation in {1,2}, regional
# flows excluded, 2009-2024. FlowName grouping matches §18's flow_ry
# classification (ODA / OOF / private) verbatim.
flow_type_ry <- lapply(CRS, function(df) {
  req <- c("ClimateAdaptation", "Year", "RecipientCode", "FlowName", "USD_Commitment_Defl")
  if (!all(req %in% names(df))) return(NULL)
  df %>%
    filter(ClimateAdaptation %in% c(1, 2), !RecipientCode %in% regionalflows,
           Year >= 2009, Year <= 2024) %>%
    group_by(FlowName) %>%
    summarise(commitments = sum(USD_Commitment_Defl, na.rm = TRUE), .groups = "drop")
}) %>%
  bind_rows() %>%
  group_by(FlowName) %>%
  summarise(commitments = sum(commitments, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    flow_group = case_when(
      grepl("^ODA", FlowName)              ~ "ODA (grants + loans)",
      grepl("Other Official", FlowName)    ~ "Other Official Flows (OOF)",
      grepl("Private|Equity|PSI", FlowName) ~ "Private / equity",
      TRUE                                  ~ "Other"
    )
  ) %>%
  arrange(desc(commitments))

total_adapt_commit_scope <- sum(flow_type_ry$commitments, na.rm = TRUE)
flow_type_ry <- flow_type_ry %>%
  mutate(share_pct = 100 * commitments / total_adapt_commit_scope)

flow_group_shares <- flow_type_ry %>%
  group_by(flow_group) %>%
  summarise(commitments = sum(commitments, na.rm = TRUE), .groups = "drop") %>%
  mutate(share_pct = 100 * commitments / total_adapt_commit_scope) %>%
  arrange(desc(commitments))

message("Flow-type composition of the adaptation panel (2009-2024):")
print(flow_group_shares, row.names = FALSE)
oda_share_scope <- flow_group_shares$share_pct[flow_group_shares$flow_group == "ODA (grants + loans)"]
message(sprintf(
  "  -> Main panel is %.2f%% ODA; %.2f%% non-ODA (OOF + private/equity + other). Main panel = ODA + OOF + private/equity (all FlowName categories are retained; the panel is NOT restricted to ODA).",
  oda_share_scope, 100 - oda_share_scope
))

flow_lines <- c(
  "\\begin{tabular}{llrr}",
  "\\toprule",
  "Flow group & CRS FlowName & Commitments (USD M) & Share (\\%) \\\\",
  "\\midrule",
  paste0(esc_tex(flow_type_ry$flow_group), " & ", esc_tex(flow_type_ry$FlowName), " & ",
         fmt2(flow_type_ry$commitments), " & ", fmt2(flow_type_ry$share_pct), " \\\\"),
  "\\midrule",
  paste0("\\multicolumn{2}{l}{Total} & ", fmt2(total_adapt_commit_scope), " & 100.00 \\\\"),
  "\\bottomrule",
  "\\end{tabular}"
)
writeLines(flow_lines, here("output", "tables", "scope", "flow_type_shares.tex"))
message("Wrote: output/tables/scope/flow_type_shares.tex")

## ODA-only panel variant (does NOT replace the main panel; for a future
## robustness run). Same row structure as `aggregated` (simple_panel_wgi.csv),
## with commitments/lcommitments/log_commits AND disbursements/ldisbursements
## recomputed from commitments_oda/disbursements_oda (built in §18 above,
## flow_ry now aggregates both); the original all-flow-type totals are kept
## as commitments_all_flow_types / disbursements_all_flow_types for reference.
## Columns defined relative to ALL-flow-type commitments would be
## inconsistent with the ODA-only figure this file represents, so
## disbursements/ldisbursements are rebuilt from disbursements_oda, and
## commitments_all, lcommitments_all,
## commitments_nonadapt, lcommitments_nonadapt, commitments_pc,
## disbursements_pc, and share_adapt still cannot be rebuilt without a much
## larger CRS re-aggregation ("all commitments across every purpose code",
## not just the adaptation marker) this file does not otherwise need, so
## those remain DROPPED (not silently stale).
stopifnot("commitments_oda" %in% names(aggregated), "disbursements_oda" %in% names(aggregated))
oda_inconsistent_cols <- c("commitments_all", "lcommitments_all", "commitments_nonadapt",
                           "lcommitments_nonadapt", "commitments_pc", "disbursements_pc",
                           "share_adapt")
adaptation_panel_oda_only <- aggregated %>%
  mutate(
    commitments_all_flow_types   = commitments,
    disbursements_all_flow_types = disbursements,
    commitments    = commitments_oda,
    lcommitments   = log1p(commitments_oda),
    log_commits    = log1p(commitments_oda),
    disbursements  = disbursements_oda,
    ldisbursements = log1p(disbursements_oda)
  ) %>%
  select(-any_of(oda_inconsistent_cols))
write.csv(adaptation_panel_oda_only,
          here("data", "processed", "adaptation_panel_oda_only.csv"),
          row.names = FALSE)
message("Wrote: data/processed/adaptation_panel_oda_only.csv  (",
        nrow(adaptation_panel_oda_only), " rows; ODA-only commitments AND disbursements ",
        "(rebuilt from disbursements_oda); dropped ",
        length(intersect(oda_inconsistent_cols, names(aggregated))),
        " flow-type-dependent columns that cannot be rebuilt without a larger CRS pass: ",
        paste(oda_inconsistent_cols, collapse = ", "))

## --- Regional-flow exclusion documentation --------------------------------
# regionalflows (§4) drops recipient codes that identify regional/unspecified
# aggregates rather than single countries; these commitments are DROPPED, not
# apportioned to member countries. Uses rio_data_adaptation (§3: ClimateAdaptation in {1,2},
# all recipient codes, before the §4 exclusion), 2009-2024.
rio_scope_window <- rio_data_adaptation %>% filter(Year >= 2009, Year <= 2024)
total_pre_exclusion <- sum(rio_scope_window$USD_Commitment_Defl, na.rm = TRUE)

regional_excl_tab <- rio_scope_window %>%
  filter(RecipientCode %in% regionalflows) %>%
  group_by(RecipientCode, RecipientName) %>%
  summarise(commitments = sum(USD_Commitment_Defl, na.rm = TRUE), .groups = "drop") %>%
  mutate(share_pct = 100 * commitments / total_pre_exclusion) %>%
  arrange(desc(commitments))

regional_total_share <- 100 * sum(regional_excl_tab$commitments) / total_pre_exclusion
message(sprintf(
  paste0("Regional-flow exclusion: %d recipient codes with adaptation-marked ",
        "commitments dropped, %.2f%% of 2009-2024 adaptation-marked commitments ",
        "(dropped from the panel, not apportioned to member countries)."),
  nrow(regional_excl_tab), regional_total_share
))

regional_lines <- c(
  "\\begin{tabular}{llrr}",
  "\\toprule",
  "Recipient code & Recipient name & Commitments (USD M) & Share (\\%) \\\\",
  "\\midrule",
  paste0(regional_excl_tab$RecipientCode, " & ", esc_tex(regional_excl_tab$RecipientName), " & ",
         fmt2(regional_excl_tab$commitments), " & ", fmt2(regional_excl_tab$share_pct), " \\\\"),
  "\\midrule",
  paste0("\\multicolumn{2}{l}{Total excluded} & ", fmt2(sum(regional_excl_tab$commitments)),
         " & ", fmt2(regional_total_share), " \\\\"),
  "\\bottomrule",
  "\\end{tabular}"
)
writeLines(regional_lines, here("output", "tables", "scope", "regional_exclusion.tex"))
message("Wrote: output/tables/scope/regional_exclusion.tex")

## --- Effect of the code-860 (FSM) fix on the panel -------------------------
fsm_commit <- sum(rio_scope_window$USD_Commitment_Defl[rio_scope_window$RecipientCode == 860],
                  na.rm = TRUE)
fsm_in_panel <- "FSM" %in% unique(adaptation_aid$RecipientISO)
n_adaptation_aid_countries <- n_distinct(adaptation_aid$RecipientISO)
message(sprintf(
  paste0("Code-860 (FSM) fix: recipient code 860 ('Micronesia' = Federated States of ",
        "Micronesia) is a genuine country, no longer excluded as regional. Its 2009-2024 ",
        "adaptation-marked commitments ($%.2fM, %.4f%% of the pre-exclusion total) now ",
        "enter the panel. FSM present in adaptation_aid post-fix: %s. adaptation_aid ",
        "country count (broader than the final estimation panel -- see the sample funnel ",
        "below): %d (would be %d without this fix)."),
  fsm_commit, 100 * fsm_commit / total_pre_exclusion, fsm_in_panel,
  n_adaptation_aid_countries, n_adaptation_aid_countries - 1L
))

## --- Sample funnel (documents the full construction
## chain from raw adaptation-marked CRS rows to the final estimation panel).
n_regional_codes_used <- n_distinct(regional_excl_tab$RecipientCode)
funnel <- data.frame(
  stage = c(
    "Raw adaptation-marked commitments (Rio marker 1 or 2), 2009-2024, all recipient codes",
    paste0("After dropping regional/unspecified recipient codes (regionalflows, n = ",
           n_regional_codes_used, " codes with nonzero commitments, of ",
           length(regionalflows), " defined)"),
    "After ISO3 mapping (countrycode + Kosovo/FSM custom matches; drops any residual NA)",
    "Recipient-year panel after listwise deletion on controls (ge\\_est, population) -- \\S15"
  ),
  commitments_usd_m = c(
    total_pre_exclusion,
    total_pre_exclusion - sum(regional_excl_tab$commitments),
    NA_real_,  # commitments unaffected by ISO3 mapping once regional codes are already dropped
    sum(aggregated$commitments, na.rm = TRUE)
  ),
  n_countries = c(
    NA_integer_,
    n_distinct(rio_scope_window$RecipientName[!rio_scope_window$RecipientCode %in% regionalflows]),
    n_adaptation_aid_countries,
    n_distinct(aggregated$recipient_iso)
  )
)
message("\n--- Sample funnel (regional exclusion -> ISO3 mapping -> listwise deletion) ---")
print(funnel, row.names = FALSE)

funnel_lines <- c(
  "\\begin{tabular}{lrr}",
  "\\toprule",
  "Stage & Commitments (USD M) & Countries \\\\",
  "\\midrule",
  paste0(funnel$stage, " & ",
         ifelse(is.na(funnel$commitments_usd_m), "---", fmt2(funnel$commitments_usd_m)), " & ",
         ifelse(is.na(funnel$n_countries), "---", funnel$n_countries), " \\\\"),
  "\\bottomrule",
  "\\end{tabular}"
)
writeLines(funnel_lines, here("output", "tables", "scope", "sample_funnel.tex"))
message("Wrote: output/tables/scope/sample_funnel.tex")

##############################################################################
# §23. EMERGENCY-RESPONSE AID COMMITMENTS PANEL
# This is a CRS proxy for humanitarian/emergency-response AID -- NOT a
# hazard-realisation measure. Design decision: lagged emergency-response aid is a bad control for the hazard-control
# specifications (same CRS, same donors, same recipient envelopes as the
# adaptation-commitments outcome -- a baseline-outcome proxy, not an
# independent hazard measure) and must not appear as a covariate in any
# hazard-control exhibit. It is retained ONLY as an input to
# code/14_hazard_napa.R's separate "NAP timing vs. lagged emergency-response
# aid" orthogonality check (§4 there) -- a check that NAP adoption timing does
# not simply follow recent humanitarian-aid surges, distinct from the
# EM-DAT-based hazard-realisation check (data/raw/emdat/emdat.csv), which
# code/14_hazard_napa.R builds itself, not this file.
#
# Purpose codes (OECD CRS purpose-code list, confirmed present 2007-2024):
#   72010  Material relief assistance and services
#   72040  Emergency food assistance
#   72050  Relief co-ordination and support services
#   73010  Immediate post-emergency reconstruction and rehabilitation
#   74020  Multi-hazard response preparedness
# All donors, all purpose-code-matching commitments regardless of the Rio
# adaptation marker (disaster response is not itself adaptation-marked).
# Same recipient standardisation (ISO3/Kosovo/FSM, `recipient_iso_from_name()`)
# and regional-flow exclusion (`regionalflows`, §0c-3) as the adaptation
# panel. Balanced to a full recipient x year grid (2008-2024) over the
# recipients retained in the final estimation panel (`aggregated`, post §15
# listwise deletion), zero-filled -- absence of emergency-response aid is a
# meaningful "no reported humanitarian response," not missingness.
##############################################################################

message("\n=== SECTION 23: EMERGENCY-RESPONSE COMMITMENTS PANEL ===")

emergency_purpose_codes <- c(72010, 72040, 72050, 73010, 74020)

emergency_response_ry <- lapply(CRS, function(df) {
  req <- c("PurposeCode", "Year", "RecipientName", "RecipientCode", "USD_Commitment_Defl")
  if (!all(req %in% names(df))) return(NULL)
  df %>%
    filter(PurposeCode %in% emergency_purpose_codes,
           !RecipientCode %in% regionalflows,
           Year >= 2008, Year <= 2024) %>%
    group_by(RecipientName, Year) %>%
    summarise(emergency_commitments = sum(USD_Commitment_Defl, na.rm = TRUE),
              .groups = "drop")
}) %>%
  bind_rows() %>%
  mutate(
    RecipientISO = recipient_iso_from_name(RecipientName)
  ) %>%
  filter(!is.na(RecipientISO)) %>%
  group_by(RecipientISO, Year) %>%
  summarise(emergency_commitments = sum(emergency_commitments, na.rm = TRUE),
            .groups = "drop") %>%
  rename(recipient_iso = RecipientISO, year = Year)

recipients_scope <- unique(aggregated$recipient_iso)
recipients_scope <- recipients_scope[!is.na(recipients_scope)]

emergency_response_panel <- expand.grid(
  recipient_iso    = recipients_scope,
  year             = 2008:2024,
  stringsAsFactors = FALSE
) %>%
  left_join(emergency_response_ry, by = c("recipient_iso", "year")) %>%
  mutate(emergency_commitments = coalesce(emergency_commitments, 0)) %>%
  arrange(recipient_iso, year)

message(sprintf(
  paste0("Emergency-response panel: %d recipient-year rows (%d recipients x %d years, ",
        "2008-2024); %d rows with commitments > 0."),
  nrow(emergency_response_panel), length(recipients_scope), length(2008:2024),
  sum(emergency_response_panel$emergency_commitments > 0)
))

write.csv(emergency_response_panel,
          here("data", "processed", "emergency_response_panel.csv"),
          row.names = FALSE)
message("Wrote: data/processed/emergency_response_panel.csv  (",
        nrow(emergency_response_panel), " rows)")

message("\n=== 01_prepare_data.R COMPLETE ===")
message("Outputs:")
message("  data/processed/adaptationNAP_donortype_wgi.csv")
message("  data/processed/simple_panel_wgi.csv")
message("  data/processed/adaptationNAP.csv")
message("  data/processed/mitigation_panel.csv")
message("  data/processed/donor_list.csv")
message("  data/processed/donor_totals.csv")
message("  data/processed/donor_recipient_year_adaptation.csv")
message("  data/processed/adaptation_panel_oda_only.csv")
message("  data/processed/emergency_response_panel.csv")
message("  output/tables/scope/flow_type_shares.tex")
message("  output/tables/scope/regional_exclusion.tex")
