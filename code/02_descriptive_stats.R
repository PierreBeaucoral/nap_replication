##############################################################################
# 02_descriptive_stats.R
# Author : Pierre Beaucoral
# Date   : 2026-06-10
# Purpose: All descriptive tables and figures.
#          Reads processed intermediates produced by 01_prepare_data.R.
#          Falls back to a minimal raw-CRS read ONLY when
#          data/processed/donor_list.csv or data/processed/donor_totals.csv
#          are absent (i.e. 01_prepare_data.R was not run first).
#          Tables: stats_des, list, nap_regional, finance_change,
#          balance_adopters; 6 descriptive figures.
#
# Inputs (primary — written by 01_prepare_data.R):
#          data/processed/adaptationNAP.csv
#          data/processed/adaptationNAP_donortype_wgi.csv
#          data/processed/donor_list.csv      (for list.tex)
#          data/processed/donor_totals.csv    (for top_donors figure)
#          data/processed/simple_panel_wgi.csv (for balance_adopters.tex, §8)
# Inputs (fallback only, if processed intermediates missing):
#          data/raw/CRS/CRS <year> Data.txt  (minimal 2–3 column read)
#
# Outputs:
#   output/tables/stats_des.tex
#   output/tables/nap_regional.tex
#   output/tables/finance_change.tex
#   output/tables/list.tex
#   output/tables/balance_adopters.tex
#   output/figures/climate_finance_evolution.png
#   output/figures/top_donors.png
#   output/figures/recipient_map.png
#   output/figures/nap_status_map.png
#   output/figures/nap_adoption_timeline.png
#   output/figures/cumulative_nap_adoption.png
##############################################################################

##############################################################################
# §0. PACKAGES
##############################################################################

library(data.table)
library(dplyr)
library(countrycode)
library(lubridate)
library(ggplot2)
library(xtable)
library(scales)
library(here)

# Single global seed
set.seed(20240601)

##############################################################################
# §0b. Helper: write_tex_float()
# Wraps a bare tabular block in the project's standard complete float.
# Caption on top, \adjustbox, bottom Notes + Source minipage. No bold/italic
# in caption, notes, or source lines.
##############################################################################

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

##############################################################################
# §0c. OUTPUT DIRECTORIES
##############################################################################

dir.create(here("output", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "tables"),  recursive = TRUE, showWarnings = FALSE)

##############################################################################
# §1. LOAD PROCESSED DATA
# Two processed frames are needed:
#
#  adaptation_aid       — recipient × year (adaptationNAP.csv, 2133 rows)
#                         Used for: totals, N_Recipients, nap_regional, figures
#
#  adaptation_dtype_nz  — recipient × year × DonorType, NON-ZERO rows only
#                         (adaptationNAP_donortype_wgi.csv filtered to nonzero)
#                         Used for: Mean/Median/SD in stats_des (computed
#                         over the sparse recipient × year × DonorType frame)
#
# Note on stats_des universe:
#   The sparse frame is built from the CRS raw files via
#   group_by(RecipientName, Year, DonorType) with sum(na.rm=TRUE).  For 2009,
#   Bangladesh (BGD) and Jordan (JOR) had a single NA-commitment row each;
#   sum(..., na.rm=TRUE) = 0 created a zero-row in the sparse frame, so they
#   appear in n_distinct(RecipientISO) = 10. adaptationNAP_donortype_wgi.csv
#   (the expanded balanced panel) also carries these rows as zeros, but the
#   nonzero filter excludes them. Hence N_Recipients comes from adaptation_aid
#   (recipient × year) which correctly gives 10 for 2009.
##############################################################################

adaptation_aid <- fread(here("data", "processed", "adaptationNAP.csv"))

# Restore date column (was coerced to character on write)
if (!inherits(adaptation_aid$date_posted, "POSIXct")) {
  adaptation_aid$date_posted <- parse_date_time(
    adaptation_aid$date_posted, orders = c("ymd HMS", "ymd", "dmy", "mdy")
  )
}

# Ensure NAP_Year is numeric
adaptation_aid$NAP_Year <- suppressWarnings(as.numeric(adaptation_aid$NAP_Year))

# Donor-type panel: keep only rows with actual (non-zero) flows.
# This is the sparse frame used for Mean/Median/SD.
adaptation_dtype_nz <- fread(
  here("data", "processed", "adaptationNAP_donortype_wgi.csv")
) %>%
  filter(Commitments > 0 | Disbursements > 0)

##############################################################################
# §2. STATS_DES TABLE
# Exhibit: paper/Tables/stats_des  (\input{Tables/stats_des})
# Headers are in USD Bn (values are /1e3); adjustbox wrapper for the wide table.
##############################################################################

# Build summary stats by year
# adaptation_aid$Commitments is in USD millions (raw CRS unit); /1e3 = USD Bn
#
# The statistics are defined on the SPARSE
# recipient × year × DonorType frame (3-way grouping, observed rows only).
# N_Recipients = n_distinct(RecipientISO) over that frame gives the correct
# counts (e.g. 10 for 2009, including BGD/JOR whose raw records had NA
# commitments that summed to 0).  Mean/Median/SD computed across the sparse
# donor-type rows reproduce the published cell values.
#
# We approximate the sparse frame by:
#   - N_Recipients and Totals from adaptation_aid (recipient × year), which has
#     the same country-year universe as the sparse frame.
#   - Mean/Median/SD from adaptation_dtype_nz (sparse nonzero donor-type rows),
#     which gives the correct within-year distributions.

totals_n <- adaptation_aid %>%
  group_by(Year) %>%
  summarise(
    Total_Commitments   = sum(Commitments,  na.rm = TRUE) / 1e3,
    Total_Disbursements = sum(Disbursements, na.rm = TRUE) / 1e3,
    N_Recipients        = n_distinct(RecipientISO),
    .groups = "drop"
  )

distrib <- adaptation_dtype_nz %>%
  group_by(Year) %>%
  summarise(
    Mean_Commitments   = mean(Commitments,   na.rm = TRUE) / 1e3,
    Median_Commitments = median(Commitments, na.rm = TRUE) / 1e3,
    SD_Commitments     = sd(Commitments,     na.rm = TRUE) / 1e3,
    .groups = "drop"
  )

summary_stats <- totals_n %>%
  left_join(distrib, by = "Year") %>%
  select(Year, Total_Commitments, Total_Disbursements, N_Recipients,
         Mean_Commitments, Median_Commitments, SD_Commitments)

latex_table <- xtable(
  summary_stats,
  label  = "tab:summary_stats",
  digits = c(0, 0, 2, 2, 0, 2, 2, 2)
)

# Fix 1 + Fix 2: all value columns are /1e3 → "(USD Bn)" not "(USD Mn)"
colnames(latex_table) <- c(
  "Year",
  "Total Commitments (USD Bn)",
  "Total Disbursements (USD Bn)",
  "Number of Recipients",
  "Mean Commitments (USD Bn)",
  "Median Commitments (USD Bn)",
  "SD Commitments (USD Bn)"
)

stats_des_path <- here("output", "tables", "stats_des.tex")

# Capture bare tabular (floating = FALSE) then wrap with write_tex_float
stats_raw <- capture.output(
  print(latex_table,
        include.rownames       = FALSE,
        floating               = FALSE,
        booktabs               = TRUE,
        size                   = "small",
        sanitize.text.function = identity)
)

write_tex_float(
  out_path      = stats_des_path,
  caption_title = "Summary Statistics of Climate Adaptation Finance by Year",
  label         = "tab:summary_stats",
  tabular_lines = stats_raw,
  notes_text    = paste0(
    "Adaptation-related finance identified by the Rio adaptation markers ",
    "(principal or significant), OECD CRS, constant USD. Totals sum all ",
    "flows per year; recipients counts distinct countries with ",
    "adaptation-marked flows; mean, median, SD computed across ",
    "recipient--donor-type cells per year. The 2009 row reflects the ",
    "partial first year of marker reporting"
  ),
  source_text   = "OECD CRS"
)

##############################################################################
# §3. LIST.TEX — Donor list by donor type
# Exhibit: paper/Tables/list  (\input{Tables/list})
# Longtable column spec {cc} (no vertical rule; booktabs rules only).
##############################################################################

# Preferred path: 01_prepare_data.R persists the classified donor list as
# data/processed/donor_list.csv. Fall back to a minimal raw-CRS read (two
# columns) only if that intermediate is missing.
donor_list_path <- here("data", "processed", "donor_list.csv")
if (file.exists(donor_list_path)) {
  message("Building list.tex from data/processed/donor_list.csv ...")
  df_unique_list <- read.csv(donor_list_path, stringsAsFactors = FALSE) %>%
    distinct(DonorType, DonorName) %>%
    arrange(DonorType, DonorName)
} else {
  message("donor_list.csv not found — minimal raw CRS read for list.tex ...")

  crs_years_list <- 2007:2024
  donor_codes_raw <- rbindlist(lapply(crs_years_list, function(y) {
    fp <- here("data", "raw", "CRS", paste0("CRS ", y, " Data.txt"))
    if (!file.exists(fp)) fp <- here("data", "raw", "CRS", paste0("CRS ", y, " data.txt"))
    if (!file.exists(fp)) return(NULL)
    dt <- fread(fp, encoding = "UTF-8",
                select = c("ClimateAdaptation", "DonorCode", "DonorName"),
                showProgress = FALSE)
    dt[ClimateAdaptation %in% c(1, 2), .(DonorCode, DonorName)]
  }), fill = TRUE)

  # Same classification as 01_prepare_data.R §3
  dac_members <- c(
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 18, 20, 21, 22, 40, 50, 61, 68,
    69, 75, 76, 301, 302, 701, 742, 801, 820, 918, 30, 45, 55, 62, 70, 72,
    77, 82, 83, 84, 87, 130, 133, 358, 543, 546, 552, 561, 566, 576, 611,
    613, 732, 764, 765
  )
  multilateral_donors <- c(
    104, 807, 811, 812, 901, 902, 903, 905, 906, 907, 909, 913, 914, 915,
    921, 923, 926, 928, 932, 940, 944, 948, 951, 952, 953, 954, 956, 958,
    959, 960, 963, 964, 966, 967, 971, 974, 976, 978, 979, 980, 981, 982,
    983, 988, 990, 992, 997, 1011, 1012, 1013, 1014, 1015, 1016, 1017,
    1018, 1019, 1020, 1023, 1024, 1025, 1037, 1038, 1039, 1058,
    1311, 1312, 1403
  )

  df_unique_list <- donor_codes_raw %>%
    mutate(DonorType = case_when(
      DonorCode %in% dac_members         ~ "Bilateral_members",
      DonorCode %in% multilateral_donors ~ "Multilateral_donors",
      TRUE                               ~ "Other"
    )) %>%
    distinct(DonorType, DonorName) %>%
    arrange(DonorType, DonorName)
}

esc_tex_fn <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([&%$#_{}])", "\\\\\\1", x)
  x
}

list_body_lines <- vapply(seq_len(nrow(df_unique_list)), function(i) {
  paste0(
    esc_tex_fn(trimws(df_unique_list$DonorType[i])), " & ",
    esc_tex_fn(trimws(df_unique_list$DonorName[i])), " \\\\"
  )
}, character(1))

# Fix 3: {cc} — no vertical rule; booktabs rules only
# Fix 5: caption on top, plain text; notes + source as plain rows at end;
#   no \textit{} in caption or notes.
list_lines <- c(
  "\\begin{longtable}{cc}",
  paste0(
    "\\caption{List of donors by donor type} \\label{tab:list}\\\\"
  ),
  "",
  "\\toprule",
  "Donor Type & Donor Name \\\\",
  "\\midrule",
  "\\endfirsthead",
  "",
  "\\multicolumn{2}{c}%",
  "{{\\tablename\\ \\thetable{} -- continued from previous page}} \\\\",
  "\\toprule",
  "Donor Type & Donor Name \\\\",
  "\\midrule",
  "\\endhead",
  "",
  "\\midrule \\multicolumn{2}{r}{{Continued on next page}} \\\\",
  "\\endfoot",
  "",
  "\\bottomrule",
  "\\endlastfoot",
  list_body_lines,
  "\\multicolumn{2}{p{\\dimexpr\\linewidth-4\\tabcolsep}}{\\footnotesize Notes: All providers of adaptation-related finance observed in the estimation sample, grouped into the donor-type categories used in the heterogeneity analysis (DAC and non-DAC bilateral providers, multilateral institutions, and private/other providers).} \\\\",
  "\\multicolumn{2}{p{\\dimexpr\\linewidth-4\\tabcolsep}}{\\footnotesize Source: OECD CRS.} \\\\",
  "\\end{longtable}"
)

list_path <- here("output", "tables", "list.tex")
writeLines(list_lines, list_path)
message("Wrote: output/tables/list.tex")

##############################################################################
# §4. NAP_REGIONAL TABLE
# Exhibit: paper/Tables/nap_regional  (\input{Tables/nap_regional})
# Human-readable column headers.
##############################################################################

nap_regional <- adaptation_aid %>%
  filter(!is.na(date_posted)) %>%
  distinct(RecipientName, date_posted) %>%
  mutate(
    Region = countrycode(RecipientName, origin = "country.name",
                         destination = "region"),
    Year   = year(date_posted)
  ) %>%
  group_by(Region) %>%
  summarise(
    N_Countries    = n(),
    First_Adoption = min(Year),
    Last_Adoption  = max(Year),
    Median_Year    = median(Year),
    .groups = "drop"
  ) %>%
  arrange(First_Adoption) %>%
  # Escape literal '&' in region names (e.g. "Latin America & Caribbean") —
  # the table is printed with sanitization off, so unescaped '&' breaks the
  # tabular alignment.
  mutate(Region = gsub("&", "\\\\&", Region, fixed = FALSE))

xt_regional <- xtable(
  nap_regional,
  label  = "tab:nap_regional",
  align  = c("l", "l", "c", "c", "c", "c"),
  # digits=0 for all columns: year columns must print as integers (e.g. 2022, not 2022.0).
  # median(Year) returns a double, so digits=0 is needed to print "2022", not "2022.0".
  digits = c(0, 0, 0, 0, 0, 0)
)

# Human-readable column headers replace
# the raw summarise() names (N_Countries, First_Adoption, ...).
colnames(xt_regional) <- c("Region", "Countries", "First adoption",
                           "Last adoption", "Median year")

nap_regional_path <- here("output", "tables", "nap_regional.tex")

nap_raw <- capture.output(
  print(xt_regional,
        include.rownames       = FALSE,
        booktabs               = TRUE,
        floating               = FALSE,
        sanitize.text.function = identity)
)

write_tex_float(
  out_path      = nap_regional_path,
  caption_title = "Regional Distribution of NAP Adoption",
  label         = "tab:nap_regional",
  tabular_lines = nap_raw,
  notes_text    = paste0(
    "Number of countries having submitted a National Adaptation Plan to the UNFCCC ",
    "by World Bank region, with the first, last, and median year of submission in each region"
  ),
  source_text   = "UNFCCC NAP Central tracking tool"
)

##############################################################################
# §5. FINANCE_CHANGE TABLE
# Exhibit: paper/Tables/finance_change  (\input{Tables/finance_change})
# Human-readable column headers; notes go in the caption, not in a body row.
# Fully descriptive, computed from the current panel (recipient × year means,
# matching the caption wording); none of its values are quoted in the text.
##############################################################################

finance_analysis <- adaptation_aid %>%
  filter(!is.na(NAP)) %>%
  group_by(RecipientName) %>%
  summarise(
    Pre_NAP_Mean  = mean(Commitments[NAP == 0], na.rm = TRUE) / 1000,
    Post_NAP_Mean = mean(Commitments[NAP == 1], na.rm = TRUE) / 1000,
    .groups = "drop"
  ) %>%
  filter(!is.na(Post_NAP_Mean)) %>%
  # Guard against division by zero: only compute Change_Pct when Pre_NAP_Mean > 0
  filter(Pre_NAP_Mean > 0) %>%
  mutate(Change_Pct = ((Post_NAP_Mean - Pre_NAP_Mean) / Pre_NAP_Mean) * 100) %>%
  arrange(desc(Change_Pct))

xt_finance <- xtable(
  head(finance_analysis, 50),
  label  = "tab:finance_changes",
  digits = c(0, 0, 2, 2, 1)
)

# Fix 2: human-readable column headers; Fix 4: no embedded note row (table-only)
colnames(xt_finance) <- c(
  "Countries", "Pre-NAP mean (USD Bn)", "Post-NAP mean (USD Bn)", "Change (\\%)"
)

finance_path <- here("output", "tables", "finance_change.tex")

fin_raw <- capture.output(
  print(xt_finance,
        include.rownames       = FALSE,
        booktabs               = TRUE,
        floating               = FALSE,
        sanitize.text.function = identity)
)

write_tex_float(
  out_path      = finance_path,
  caption_title = "Changes in Adaptation Finance After NAP Adoption",
  label         = "tab:finance_changes",
  tabular_lines = fin_raw,
  notes_text    = paste0(
    "Mean annual adaptation-related commitments ",
    "(USD billion, constant prices) received by each NAP-adopting country ",
    "before and after its adoption year; the final column is the percentage ",
    "change of the post-adoption mean relative to the pre-adoption mean. ",
    "Descriptive only --- these raw changes do not adjust for time trends or composition"
  ),
  source_text   = "OECD CRS (Rio adaptation markers); UNFCCC NAP Central"
)
message("Wrote: output/tables/finance_change.tex")

##############################################################################
# §6. FIGURE HELPERS
# labs(title = NULL, subtitle = NULL) on all 6 figures.
# Figure titles and notes go in LaTeX \caption{} only.
##############################################################################

theme_paper <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      text              = element_text(family = "serif", size = base_size),
      axis.title        = element_text(size = base_size, face = "bold"),
      axis.text         = element_text(size = base_size - 2),
      legend.position   = "bottom",
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(color = "gray90"),
      plot.margin       = margin(1, 1, 1, 1, "cm"),
      plot.caption      = element_text(size = 10, hjust = 0,
                                       margin = margin(t = 20))
    )
}

##############################################################################
# §7. FIGURE — climate_finance_evolution.png
# Exhibit: paper/Figures/climate_finance_evolution.png
# No in-figure title/subtitle/caption (they go in the LaTeX caption).
##############################################################################

yearly_totals <- adaptation_aid %>%
  group_by(Year) %>%
  summarise(
    Total_Commitments   = sum(Commitments,  na.rm = TRUE) / 1000,
    Total_Disbursements = sum(Disbursements, na.rm = TRUE) / 1000,
    .groups = "drop"
  )

p_finance <- ggplot(yearly_totals, aes(x = Year)) +
  geom_line(aes(y = Total_Commitments,   color = "Commitments"),   linewidth = 1) +
  geom_line(aes(y = Total_Disbursements, color = "Disbursements"),  linewidth = 1) +
  geom_point(aes(y = Total_Commitments,  color = "Commitments"),   size = 3) +
  geom_point(aes(y = Total_Disbursements, color = "Disbursements"), size = 3) +
  scale_color_manual(
    values = c("Commitments" = "#2E86C1", "Disbursements" = "#E67E22")
  ) +
  scale_x_continuous(
    breaks = seq(min(yearly_totals$Year), max(yearly_totals$Year), by = 1)
  ) +
  scale_y_continuous(labels = comma) +
  # No title, subtitle, or caption — those go in LaTeX \caption{}
  labs(title = NULL, subtitle = NULL, caption = NULL,
       x = "Year", y = "USD (Billions)", color = "Type") +
  theme_paper()

ggsave(here("output", "figures", "climate_finance_evolution.png"),
       p_finance, width = 10, height = 6, dpi = 300)
message("Wrote: output/figures/climate_finance_evolution.png")
# Appendix A (\ref{fig:adaptation_trends_app}) inputs the same figure under the
# name climate_adaptation_flows.png; emit it so the manuscript compiles from the pipeline.
ggsave(here("output", "figures", "climate_adaptation_flows.png"),
       p_finance, width = 10, height = 6, dpi = 300)
message("Wrote: output/figures/climate_adaptation_flows.png")

##############################################################################
# §8. FIGURE — top_donors.png
# Exhibit: paper/Figures/top_donors.png
# No in-figure title/subtitle/caption (they go in the LaTeX caption).
#
# Read donor_totals.csv when present (written by
# 01_prepare_data.R §3); fall back to a minimal raw-CRS read ONLY if absent.
# All raw-CRS objects (crs_years_list) are defined locally in the else branch
# so that no unconditional reference to objects from §3 exists.
##############################################################################

donor_totals_path <- here("data", "processed", "donor_totals.csv")

if (file.exists(donor_totals_path)) {
  message("Building top_donors figure from data/processed/donor_totals.csv ...")
  donor_totals <- read.csv(donor_totals_path, stringsAsFactors = FALSE) %>%
    group_by(DonorName) %>%
    summarise(Total_Commitments = sum(Total_Commitments, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(desc(Total_Commitments)) %>%
    slice_head(n = 10)
} else {
  message("donor_totals.csv not found — minimal raw CRS read for top_donors figure ...")
  crs_years_list <- 2007:2024
  donor_totals_usd <- rbindlist(lapply(crs_years_list, function(y) {
    fp <- here("data", "raw", "CRS", paste0("CRS ", y, " Data.txt"))
    if (!file.exists(fp)) fp <- here("data", "raw", "CRS", paste0("CRS ", y, " data.txt"))
    if (!file.exists(fp)) return(NULL)
    dt <- fread(fp, encoding = "UTF-8",
                select = c("ClimateAdaptation", "DonorName", "USD_Commitment_Defl"),
                showProgress = FALSE)
    dt[ClimateAdaptation %in% c(1, 2), .(DonorName, USD_Commitment_Defl)]
  }), fill = TRUE)

  donor_totals <- donor_totals_usd %>%
    group_by(DonorName) %>%
    summarise(Total_Commitments = sum(USD_Commitment_Defl, na.rm = TRUE) / 1e3,
              .groups = "drop") %>%
    arrange(desc(Total_Commitments)) %>%
    slice_head(n = 10)
}

p_donors <- ggplot(
  donor_totals,
  aes(x = reorder(DonorName, Total_Commitments), y = Total_Commitments)
) +
  geom_bar(stat = "identity", fill = "#2E86C1", alpha = 0.8) +
  coord_flip() +
  # No title, subtitle, or caption — those go in LaTeX \caption{}
  labs(title = NULL, subtitle = NULL, caption = NULL,
       x = "", y = "USD (Billions)") +
  theme_paper()

ggsave(here("output", "figures", "top_donors.png"),
       p_donors, width = 10, height = 6, dpi = 300)
message("Wrote: output/figures/top_donors.png")

##############################################################################
# §9. FIGURE — recipient_map.png
# Exhibit: paper/Figures/recipient_map.png
# No in-figure title/subtitle/caption (they go in the LaTeX caption).
##############################################################################

world <- map_data("world")

recipient_totals <- adaptation_aid %>%
  group_by(RecipientName) %>%
  summarise(Total_Received = sum(Commitments, na.rm = TRUE) / 1000,
            .groups = "drop")

recipient_totals$region <- countrycode(
  recipient_totals$RecipientName,
  origin = "country.name", destination = "country.name",
  custom_match = c("Kosovo" = "Kosovo")
)

world_data_recv <- left_join(world, recipient_totals, by = "region")

p_map <- ggplot(world_data_recv,
                aes(x = long, y = lat, group = group, fill = Total_Received)) +
  geom_polygon(color = "white", linewidth = 0.1) +
  scale_fill_viridis_c(option = "magma", direction = -1,
                       na.value = "grey80", name = "USD (Billions)") +
  coord_fixed(1.3) +
  # No title, subtitle, or caption — those go in LaTeX \caption{}
  labs(title = NULL, subtitle = NULL, caption = NULL) +
  theme_minimal() +
  theme(
    text          = element_text(family = "serif", size = 12),
    axis.text     = element_blank(),
    axis.title    = element_blank(),
    panel.grid    = element_blank()
  )

ggsave(here("output", "figures", "recipient_map.png"),
       p_map, width = 12, height = 8, dpi = 300)
message("Wrote: output/figures/recipient_map.png")

##############################################################################
# §10. FIGURE — nap_status_map.png
# Exhibit: paper/Figures/nap_status_map.png
# No in-figure title/subtitle/caption (they go in the LaTeX caption).
##############################################################################

nap_status <- adaptation_aid %>%
  group_by(RecipientName, RecipientISO) %>%
  summarise(has_nap = !all(is.na(date_posted)), .groups = "drop") %>%
  distinct()

nap_status$region <- countrycode(
  nap_status$RecipientName,
  origin = "country.name", destination = "country.name",
  custom_match = c("Kosovo" = "Kosovo")
)

world_data_nap <- left_join(world, nap_status, by = "region")

p_nap_map <- ggplot(world_data_nap,
                    aes(x = long, y = lat, group = group, fill = has_nap)) +
  geom_polygon(color = "white", linewidth = 0.1) +
  scale_fill_manual(
    values = c("TRUE" = "#2E86C1", "FALSE" = "#E67E22"),
    labels = c("TRUE" = "Has NAP", "FALSE" = "No NAP"),
    name   = "National Adaptation Plan Status",
    na.value = "grey80"
  ) +
  coord_fixed(1.3) +
  # No title, subtitle, or caption — those go in LaTeX \caption{}
  labs(title = NULL, subtitle = NULL, caption = NULL) +
  theme_minimal() +
  theme(
    text           = element_text(family = "serif", size = 12),
    axis.text      = element_blank(),
    axis.title     = element_blank(),
    panel.grid     = element_blank(),
    legend.position = "bottom",
    legend.title   = element_text(size = 12),
    plot.margin    = margin(1, 1, 1, 1, "cm")
  )

ggsave(here("output", "figures", "nap_status_map.png"),
       p_nap_map, width = 12, height = 8, dpi = 300)
message("Wrote: output/figures/nap_status_map.png")

##############################################################################
# §11. FIGURE — nap_adoption_timeline.png
# Exhibit: paper/Figures/nap_adoption_timeline.png
# No in-figure title/subtitle/caption (they go in the LaTeX caption).
##############################################################################

nap_timeline <- adaptation_aid %>%
  filter(!is.na(date_posted)) %>%
  distinct(RecipientName, date_posted) %>%
  arrange(date_posted) %>%
  mutate(
    adoption_order = row_number(),
    year           = year(date_posted),
    month          = month(date_posted),
    date           = as.Date(date_posted)
  )

p_timeline <- ggplot(
  nap_timeline,
  aes(x = date, y = reorder(RecipientName, adoption_order))
) +
  geom_point(color = "#2E86C1", size = 3) +
  geom_segment(
    aes(x = min(date), xend = date,
        yend = reorder(RecipientName, adoption_order)),
    color = "gray80", linewidth = 0.5
  ) +
  # No title, subtitle, or caption — those go in LaTeX \caption{}
  labs(title = NULL, subtitle = NULL, caption = NULL,
       x = "Adoption Date", y = NULL) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme_minimal() +
  theme(
    text               = element_text(family = "serif", size = 12),
    axis.title.x       = element_text(size = 12, face = "bold"),
    axis.text.y        = element_text(size = 12),
    axis.text.x        = element_text(size = 12, angle = 45, hjust = 1),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin        = margin(1, 1, 1, 1, "cm")
  )

ggsave(here("output", "figures", "nap_adoption_timeline.png"),
       p_timeline,
       width  = 12,
       height = 0.25 * nrow(nap_timeline),
       dpi    = 300)
message("Wrote: output/figures/nap_adoption_timeline.png")

##############################################################################
# §12. FIGURE — cumulative_nap_adoption.png
# Exhibit: paper/Figures/cumulative_nap_adoption.png
# No in-figure title/subtitle/caption (they go in the LaTeX caption).
##############################################################################

cumulative_adoption <- nap_timeline %>%
  arrange(date) %>%
  mutate(
    cumulative_count = row_number(),
    year_month       = floor_date(date, "month")
  )

p_cumulative <- ggplot(cumulative_adoption, aes(x = date, y = cumulative_count)) +
  geom_step(color = "#2E86C1", linewidth = 1) +
  geom_point(size = 2, color = "#2E86C1") +
  # No title, subtitle, or caption — those go in LaTeX \caption{}
  labs(title = NULL, subtitle = NULL, caption = NULL,
       x = "Year", y = "Number of Countries with NAP") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(
    breaks = seq(0, max(cumulative_adoption$cumulative_count), by = 5)
  ) +
  theme_minimal() +
  theme(
    text            = element_text(family = "serif", size = 12),
    axis.title      = element_text(size = 12, face = "bold"),
    axis.text       = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.margin     = margin(1, 1, 1, 1, "cm")
  )

ggsave(here("output", "figures", "cumulative_nap_adoption.png"),
       p_cumulative, width = 10, height = 6, dpi = 300)
message("Wrote: output/figures/cumulative_nap_adoption.png")

##############################################################################
# §8. BALANCE TABLE: adopters vs never-adopters, pre-adoption period
# Documents pre-adoption comparability of the treatment and comparison groups
# (app:descriptives).
# The pre-period window and its
# "before the first NAP cohort" justification used to hardcode 2009-2012 and
# "first cohort 2013" -- the actual first cohort on the current panel is
# 2015 (14_hazard_napa.R computes this the same way and gets the same
# answer), so 2009-2012 was neither the true pre-treatment window nor
# internally consistent with 14's own number. The window is now derived from
# the panel itself: 2009 (panel start) through (first real cohort - 1), i.e.
# every year before ANY unit is treated -- no hardcoded cutoff year anywhere
# below. Country-level means over that window, Welch two-sample t-tests on
# the group difference. Reads the estimation panel (single source of truth
# for the DiD sample).
##############################################################################

message("\n=== Section 8: Balance table (adopters vs never-adopters) ===\n")

panel_bal <- fread(here("data", "processed", "simple_panel_wgi.csv"))
panel_bal <- as.data.frame(panel_bal)

## --- PVCCI merge coverage audit ------------------------------------------
## The `pvcci` column read below is merged in 01_prepare_data.R §10a by ISO3
## (countrycode(recipient_name, "country.name", "iso3c")), NOT by raw
## recipient-name string match. 01 previously ALSO carried a second,
## name-based re-merge that fed the derived PVCCI_GE/PVCCI_sq/PVCCI_th
## interaction terms (fixed in the same 01 edit; those terms are not consumed
## here). This block re-derives the name-based match count directly from the
## raw inputs, cheaply (no CRS read), purely to document the coverage gain
## the ISO3 merge buys over a naive name join.
pvcci_raw_check <- suppressWarnings(
  fread(here("data", "raw", "PVCCI.csv"), encoding = "UTF-8")
)
recipients_in_panel <- unique(panel_bal$recipient_name)
n_recipients_panel  <- length(recipients_in_panel)
name_matched   <- intersect(recipients_in_panel, unique(pvcci_raw_check$recipient_name))
n_name_match   <- length(name_matched)
name_unmatched <- setdiff(unique(pvcci_raw_check$recipient_name), recipients_in_panel)
n_iso3_match   <- sum(!is.na(panel_bal$pvcci[!duplicated(panel_bal$recipient_name)]))

# Classify each name-unmatched PVCCI country programmatically (hardcoded prose
# would silently break if PVCCI.csv, the panel's recipient set, or
# countrycode's matching ever changed). A country
# is "recoverable by ISO3" if its countrycode()-derived ISO3 appears anywhere
# in the panel's recipient_iso column; otherwise it is genuinely absent from
# the estimation panel (not a re-merge failure).
panel_iso3_set <- unique(panel_bal$recipient_iso)
unmatched_iso3 <- suppressWarnings(
  countrycode(name_unmatched, origin = "country.name", destination = "iso3c")
)
recoverable   <- name_unmatched[!is.na(unmatched_iso3) & unmatched_iso3 %in% panel_iso3_set]
truly_absent  <- setdiff(name_unmatched, recoverable)

message(sprintf(
  paste0("PVCCI merge coverage: naive name-based match = %d / %d panel recipients. ",
        "%d PVCCI.csv countries do not exact-string-match a panel recipient name (%s). ",
        "Of these, %d are genuinely absent from the panel (not an adaptation-finance ",
        "recipient): %s. The remaining %d are name-spelling variants recoverable by ",
        "ISO3: %s. ISO3-based match (as merged in 01_prepare_data.R §10a, feeding ",
        "the `pvcci` column read below) = %d / %d panel recipients -- the full coverage ",
        "of the %d-country PVCCI source index."),
  n_name_match, n_recipients_panel,
  length(name_unmatched), paste(name_unmatched, collapse = ", "),
  length(truly_absent), if (length(truly_absent) == 0) "(none)" else paste(truly_absent, collapse = ", "),
  length(recoverable), if (length(recoverable) == 0) "(none)" else paste(recoverable, collapse = ", "),
  n_iso3_match, n_recipients_panel, n_distinct(pvcci_raw_check$recipient_name)
))

panel_start_year <- min(panel_bal$year, na.rm = TRUE)
first_cohort_year <- min(panel_bal$nap_year[panel_bal$nap_year >= panel_start_year],
                         na.rm = TRUE)
pre_period_end <- first_cohort_year - 1L
message(sprintf(
  "Balance-table pre-period: %d-%d (panel start through first real cohort [%d] - 1).",
  panel_start_year, pre_period_end, first_cohort_year))

pre_bal <- panel_bal %>%
  filter(year >= panel_start_year, year <= pre_period_end) %>%
  group_by(recipient_name) %>%
  summarise(
    adopter        = as.integer(any(!is.na(nap_year) &
                                      nap_year >= panel_start_year &
                                      nap_year <= max(panel_bal$year, na.rm = TRUE))),
    commitments    = mean(commitments,     na.rm = TRUE),
    disbursements  = mean(disbursements,   na.rm = TRUE),
    lcommitments   = mean(lcommitments,    na.rm = TRUE),
    ldisbursements = mean(ldisbursements,  na.rm = TRUE),
    share_adapt    = mean(share_adapt,     na.rm = TRUE),
    ge_est         = mean(ge_est,          na.rm = TRUE),
    log_population = mean(log(population), na.rm = TRUE),
    pvcci          = mean(pvcci,           na.rm = TRUE),
    .groups = "drop"
  )

n_adopt <- sum(pre_bal$adopter == 1L)
n_never <- sum(pre_bal$adopter == 0L)
message(sprintf("  Pre-period country means: %d adopters | %d never-adopters",
                n_adopt, n_never))

bal_vars <- c(
  commitments    = "Adaptation commitments (USD M)",
  lcommitments   = "Log(1+adaptation commitments)",
  disbursements  = "Adaptation disbursements (USD M)",
  ldisbursements = "Log(1+adaptation disbursements)",
  share_adapt    = "Share of global adaptation finance (\\%)",
  ge_est         = "Government Effectiveness (WGI)",
  log_population = "Log population",
  pvcci          = "PVCCI"
)

bal_rows <- vector("list", length(bal_vars))
for (i in seq_along(bal_vars)) {
  v  <- names(bal_vars)[i]
  x1 <- pre_bal[[v]][pre_bal$adopter == 1L]
  x0 <- pre_bal[[v]][pre_bal$adopter == 0L]
  x1 <- x1[is.finite(x1)]
  x0 <- x0[is.finite(x0)]
  stopifnot(length(x1) >= 2L, length(x0) >= 2L)
  tt <- t.test(x1, x0)   # Welch: unequal variances by default
  # Normalized difference (Imbens-Rubin): |mean diff| / sqrt((s1^2 + s0^2)/2);
  # values above 0.25 conventionally flag meaningful imbalance.
  nd <- abs(mean(x1) - mean(x0)) / sqrt((var(x1) + var(x0)) / 2)
  bal_rows[[i]] <- data.frame(
    variable   = bal_vars[[i]],
    n1         = length(x1),
    n0         = length(x0),
    adopters   = mean(x1),
    never      = mean(x0),
    difference = mean(x1) - mean(x0),
    se         = tt$stderr,
    p          = tt$p.value,
    norm_diff  = nd
  )
}
bal_tab <- do.call(rbind, bal_rows)

fmt3 <- function(x) formatC(x, format = "f", digits = 3)

bal_lines <- c(
  "\\begin{tabular}{lccccccc}",
  "  \\toprule",
  paste0("  Variable & $N$ (A/N) & Adopters & Never-adopters & Difference & ",
         "SE & $p$-value & Norm.\\ diff. \\\\"),
  "  \\midrule",
  paste0("  ", bal_tab$variable, " & ", bal_tab$n1, "/", bal_tab$n0, " & ",
         fmt3(bal_tab$adopters), " & ", fmt3(bal_tab$never), " & ",
         fmt3(bal_tab$difference), " & (", fmt3(bal_tab$se), ") & ",
         fmt3(bal_tab$p), " & ", fmt3(bal_tab$norm_diff), " \\\\"),
  "  \\bottomrule",
  "\\end{tabular}"
)

write_tex_float(
  out_path      = here("output", "tables", "balance_adopters.tex"),
  caption_title = "Pre-adoption comparability of adopters and never-adopters",
  label         = "tab:balance_adopters",
  tabular_lines = bal_lines,
  notes_text    = paste0(
    "Country-level means, pre-adoption period ", panel_start_year, "--",
    pre_period_end, " (through the year before cohort ", first_cohort_year,
    "). Adopters: ", n_adopt, "; never-adopters: ", n_never,
    ". $N$(A/N): non-missing per variable. Commitments/disbursements: annual ",
    "USD-million averages (Rio markers); log rows use $\\log(1+\\cdot)$. PVCCI = ",
    "Physical Vulnerability to Climate Change Index. Differences: Welch ",
    "$t$-tests, SE of the difference; Norm.\\ diff.\\ = Imbens--Rubin normalized ",
    "difference (threshold in text)"),
  source_text   = paste0("OECD CRS, UNFCCC NAP Central, World Bank WGI, ",
                         "FERDI PVCCI")
)

message("\n=== 02_descriptive_stats.R COMPLETE ===")
message("Tables: stats_des.tex, nap_regional.tex, finance_change.tex, list.tex")
message("Figures: climate_finance_evolution.png, top_donors.png, recipient_map.png,")
message("         nap_status_map.png, nap_adoption_timeline.png, cumulative_nap_adoption.png")
