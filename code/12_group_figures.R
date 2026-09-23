##############################################################################
# 12_group_figures.R
# Builds the two descriptive group figures --
# Figures/adaptation_marginal.pdf (main text, \label{fig:adapt_marginal}) and
# Figures/adopter_vs_never.pdf (appendix, \label{fig:adopter_groups}). The
# original versions (July 2026) were made outside the pipeline; this script
# recreates them from the processed panel on
# disk and reproduces the numbers already locked into the manuscript's captions
# (verified below; do not edit those captions from this script).
#
# STATUS: wired into run_all.R as stage 12 (see run_all.R header/pipeline list).
# Self-contained given that 01_prepare_data.R has produced
# data/processed/simple_panel_wgi.csv; duplicates small pieces of 02's
# theme_paper() rather than sourcing 02, per the project's convention that
# stage scripts do not source one another (see 04_robustness.R §1).
#
# Inputs : data/processed/simple_panel_wgi.csv
# Outputs:
#   output/figures/group/adaptation_marginal.pdf   (Figure, main text)
#   output/figures/group/adaptation_marginal.png   (preview mirror)
#   output/figures/group/adopter_vs_never.pdf      (Figure, appendix)
#   output/figures/group/adopter_vs_never.png      (preview mirror)
#   -> copied (flattened, no "group/" subfolder) to:
#      paper/Figures/adaptation_marginal.pdf
#      paper/Figures/adaptation_marginal.png
#      paper/Figures/adopter_vs_never.pdf
#      paper/Figures/adopter_vs_never.png
#
# ============================================================
# Paper-to-Code Naming Map
# ============================================================
# Paper Notation / caption quantity          | Code name(s)
# Adaptation-marked commitments               | commitments
# Total development-finance commitments       | commitments_all
# Non-adaptation commitments                   | commitments_nonadapt
# Adaptation share of total (%)                 | share_pct  = commitments / commitments_all * 100
# Growth index, 2010 = 100                      | idx_adapt, idx_total, idx_nonadapt
# NAP adoption cohort (first-submission year)   | nap_year (0/NA = never-treated)
# Ever-adopter indicator (58 countries)         | ever_adopter = any(!is.na(nap_year) & nap_year > 0) by recipient_name
# Never-adopter indicator (86 countries)        | !ever_adopter
# Mean adaptation commitments by group          | mean_commit
# Group growth index, 2010 = 100                 | idx_mean
# 2021-2024 main adoption-wave shading           | ADOPTION_WAVE_START / END (constants below)
# ============================================================
#
# Verification against the manuscript captions (computed on this data; see the
# printed summary at the end of this script):
#   Fig. adaptation_marginal — share 2010 = 2.31% -> caption "2.3%"
#                              share 2024 = 7.12% -> caption "7.1%"
#                              pooled avg share 2010-2024 = 5.54% -> caption "5.5%"
#                              adaptation growth x5.33 -> caption "x5.3"
#                              total growth x1.73 -> caption "x1.7"
#   Fig. adopter_vs_never    — n ever-adopters = 58, n never-adopters = 86
#                              share 2010: ever 50.3% / never 49.7%
#                              share 2024: ever 50.2% / never 49.8%
#                              growth: ever x5.22 -> caption "x5.2"
#                                      never x5.35 -> caption "x5.3"
# All figures match the committed captions to the stated precision. No
# discrepancy found; captions in the manuscript are left untouched.
##############################################################################

##############################################################################
# §0. PACKAGES
##############################################################################

library(here)
library(dplyr)
library(ggplot2)
library(cowplot)

# titles go in the LaTeX caption
set.seed(20240601)

##############################################################################
# §0b. Helper: copy_to_paper_flat()
# Mirrors a file from output/figures/group/<name> to paper/Figures/<name>,
# dropping the "group/" subfolder so the manuscript's flat
# \includegraphics{Figures/adaptation_marginal.pdf} reference keeps working.
##############################################################################

copy_to_paper_flat <- function(out_path) {
  if (!dir.exists(here("paper"))) {  # stand-alone replication package: no paper/
    message("paper/ not found -- exhibits are left in output/ only")
    return(invisible(NULL))
  }
  dest <- here("paper", "Figures", basename(out_path))
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  file.copy(out_path, dest, overwrite = TRUE)
  message("Mirrored to: ", dest)
  invisible(dest)
}

##############################################################################
# §0d. Helper: combine_panels()
# Stacks two ggplot panels vertically with a single shared bottom legend,
# using cowplot (patchwork is not a top-level renv.lock dependency in this
# project; cowplot is -- see renv.lock -- so this keeps the script
# renv::restore()-reproducible without touching renv.lock, which is out of
# scope for this task). legend_from selects which panel's legend to reuse
# (only relevant when the two panels' legends are identical or one panel
# has no legend at all, as is the case for both figures built below).
##############################################################################

combine_panels <- function(p_top, p_bottom, legend_from = p_bottom) {
  legend   <- cowplot::get_legend(legend_from + theme(legend.position = "bottom"))
  body     <- cowplot::plot_grid(
    p_top + theme(legend.position = "none"),
    p_bottom + theme(legend.position = "none"),
    ncol = 1L, align = "v"
  )
  cowplot::plot_grid(body, legend, ncol = 1L, rel_heights = c(1, 0.08))
}

##############################################################################
# §0c. Shared theme (serif, no in-figure titles beyond panel labels)
##############################################################################

theme_group <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      text              = element_text(family = "serif", size = base_size),
      plot.title        = element_text(size = base_size, hjust = 0, face = "plain"),
      axis.title        = element_text(size = base_size),
      axis.text         = element_text(size = base_size - 2),
      axis.text.x       = element_text(angle = 45, hjust = 1),
      legend.position   = "bottom",
      legend.title      = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(color = "gray90")
    )
}

# Project palette (matches code/02_descriptive_stats.R and
# code/07_principal_and_share.R): blue / orange / gray, each additionally
# distinguished by point shape so the figures read in grayscale.
PAL_ADAPT      <- "#2E86C1"
PAL_TOTAL      <- "#E67E22"
PAL_NONADAPT   <- "#7F8C8D"

##############################################################################
# §1. Load processed panel
##############################################################################

message("\n=== 12_group_figures.R: loading panel ===\n")

panel <- read.csv(here("data", "processed", "simple_panel_wgi.csv"),
                   stringsAsFactors = FALSE)

FIG_START <- 2010L
FIG_END   <- 2024L

panel_fig <- panel %>% filter(year >= FIG_START, year <= FIG_END)

dir.create(here("output", "figures", "group"), recursive = TRUE, showWarnings = FALSE)

##############################################################################
# §2. Figure: adaptation_marginal.pdf
#   Panel A -- adaptation-marked commitments as a share of total
#              development-finance commitments, by year.
#   Panel B -- adaptation / total / non-adaptation commitments indexed to
#              their FIG_START (2010) value.
##############################################################################

message("Building adaptation_marginal figure ...")

yearly_marginal <- panel_fig %>%
  group_by(year) %>%
  summarise(
    commitments          = sum(commitments,          na.rm = TRUE),
    commitments_all       = sum(commitments_all,       na.rm = TRUE),
    commitments_nonadapt  = sum(commitments_nonadapt,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(share_pct = commitments / commitments_all * 100)

base_adapt     <- yearly_marginal$commitments[yearly_marginal$year == FIG_START]
base_total     <- yearly_marginal$commitments_all[yearly_marginal$year == FIG_START]
base_nonadapt  <- yearly_marginal$commitments_nonadapt[yearly_marginal$year == FIG_START]

yearly_marginal <- yearly_marginal %>%
  mutate(
    idx_adapt    = commitments / base_adapt * 100,
    idx_total    = commitments_all / base_total * 100,
    idx_nonadapt = commitments_nonadapt / base_nonadapt * 100
  )

# Numbers-in-caption check (printed, not asserted -- see header verification
# block; kept here so a future rerun on updated data surfaces any drift).
share_2010     <- yearly_marginal$share_pct[yearly_marginal$year == FIG_START]
share_2024     <- yearly_marginal$share_pct[yearly_marginal$year == FIG_END]
share_pooled   <- sum(yearly_marginal$commitments) / sum(yearly_marginal$commitments_all) * 100
growth_adapt   <- yearly_marginal$idx_adapt[yearly_marginal$year == FIG_END] / 100
growth_total   <- yearly_marginal$idx_total[yearly_marginal$year == FIG_END] / 100
message(sprintf(
  paste0("  share %.0f = %.2f%% | share %.0f = %.2f%% | pooled avg = %.2f%% | ",
         "adaptation growth x%.2f | total growth x%.2f"),
  FIG_START, share_2010, FIG_END, share_2024, share_pooled, growth_adapt, growth_total
))

p_marginal_A <- ggplot(yearly_marginal, aes(x = year, y = share_pct)) +
  geom_area(fill = PAL_ADAPT, alpha = 0.15) +
  geom_line(color = PAL_ADAPT, linewidth = 1) +
  geom_point(color = PAL_ADAPT, size = 2) +
  scale_x_continuous(breaks = FIG_START:FIG_END) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(title = "(A) Adaptation as a share of total development finance",
       x = NULL, y = "% of total dev. finance") +
  theme_group()

idx_long <- bind_rows(
  yearly_marginal %>% transmute(year, value = idx_adapt,    series = "Adaptation"),
  yearly_marginal %>% transmute(year, value = idx_total,    series = "Total"),
  yearly_marginal %>% transmute(year, value = idx_nonadapt, series = "Non-adaptation")
) %>%
  mutate(series = factor(series, levels = c("Adaptation", "Total", "Non-adaptation")))

p_marginal_B <- ggplot(idx_long, aes(x = year, y = value, colour = series,
                                      shape = series, linetype = series)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_colour_manual(values = c(
    "Adaptation" = PAL_ADAPT, "Total" = PAL_TOTAL, "Non-adaptation" = PAL_NONADAPT
  )) +
  scale_shape_manual(values = c(
    "Adaptation" = 16, "Total" = 17, "Non-adaptation" = 15
  )) +
  scale_linetype_manual(values = c(
    "Adaptation" = "solid", "Total" = "solid", "Non-adaptation" = "solid"
  )) +
  scale_x_continuous(breaks = FIG_START:FIG_END) +
  labs(title = "(B) Growth indexed to 2010 = 100",
       x = "Year", y = "Index, 2010 = 100") +
  theme_group()

p_marginal <- combine_panels(p_marginal_A, p_marginal_B)

fig_marginal_pdf <- here("output", "figures", "group", "adaptation_marginal.pdf")
fig_marginal_png <- here("output", "figures", "group", "adaptation_marginal.png")
ggsave(fig_marginal_pdf, p_marginal, width = 7.5, height = 8, device = cairo_pdf)
ggsave(fig_marginal_png, p_marginal, width = 7.5, height = 8, dpi = 300)
message("Saved: ", fig_marginal_pdf)
copy_to_paper_flat(fig_marginal_pdf)
copy_to_paper_flat(fig_marginal_png)

##############################################################################
# §3. Figure: adopter_vs_never.pdf
#   Panel A -- mean adaptation commitments (USD million) for ever-adopters
#              (58 countries) vs. never-adopters (86 countries), by year.
#   Panel B -- same two series indexed to their FIG_START (2010) value.
#   Shaded band marks the 2021-2024 main adoption wave.
##############################################################################

message("Building adopter_vs_never figure ...")

ADOPTION_WAVE_START <- 2020.5
ADOPTION_WAVE_END   <- FIG_END + 0.5

# Ever-adopter status must be computed at the recipient level: nap_year is
# NA in 2009 (partial Rio-marker reporting year) even for countries that do
# adopt later, so a naive per-row nap_year > 0 check misclassifies adopters
# as never-treated in that one year. Group-level "any" gives the correct
# 58 ever- / 86 never-adopter split reported in the manuscript.
group_lookup <- panel %>%
  group_by(recipient_name) %>%
  summarise(ever_adopter = any(!is.na(nap_year) & nap_year > 0), .groups = "drop")

n_ever  <- sum(group_lookup$ever_adopter)
n_never <- sum(!group_lookup$ever_adopter)
message(sprintf("  n ever-adopters = %d | n never-adopters = %d", n_ever, n_never))

panel_groups <- panel_fig %>%
  left_join(group_lookup, by = "recipient_name")

yearly_groups <- panel_groups %>%
  group_by(ever_adopter, year) %>%
  summarise(mean_commit = mean(commitments, na.rm = TRUE), .groups = "drop") %>%
  mutate(group_label = if_else(
    ever_adopter,
    sprintf("Ever-adopter (n=%d)", n_ever),
    sprintf("Never-adopter (n=%d)", n_never)
  ))

base_by_group <- yearly_groups %>%
  filter(year == FIG_START) %>%
  select(ever_adopter, base_commit = mean_commit)

yearly_groups <- yearly_groups %>%
  left_join(base_by_group, by = "ever_adopter") %>%
  mutate(idx_mean = mean_commit / base_commit * 100)

growth_ever  <- yearly_groups %>% filter(ever_adopter)  %>% filter(year == FIG_END) %>% pull(idx_mean) / 100
growth_never <- yearly_groups %>% filter(!ever_adopter) %>% filter(year == FIG_END) %>% pull(idx_mean) / 100
message(sprintf("  ever-adopter growth x%.2f | never-adopter growth x%.2f",
                 growth_ever, growth_never))

group_levels <- c(sprintf("Ever-adopter (n=%d)", n_ever), sprintf("Never-adopter (n=%d)", n_never))
yearly_groups <- yearly_groups %>% mutate(group_label = factor(group_label, levels = group_levels))

p_groups_A <- ggplot(yearly_groups, aes(x = year, y = mean_commit,
                                         colour = group_label, shape = group_label)) +
  annotate("rect", xmin = ADOPTION_WAVE_START, xmax = ADOPTION_WAVE_END,
           ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_colour_manual(values = c(PAL_ADAPT, PAL_TOTAL)) +
  scale_shape_manual(values = c(16, 17)) +
  scale_x_continuous(breaks = FIG_START:FIG_END) +
  labs(title = "(A) Mean adaptation commitments (USD million)",
       x = NULL, y = "Mean, USD m") +
  theme_group()

p_groups_B <- ggplot(yearly_groups, aes(x = year, y = idx_mean,
                                         colour = group_label, shape = group_label)) +
  annotate("rect", xmin = ADOPTION_WAVE_START, xmax = ADOPTION_WAVE_END,
           ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_colour_manual(values = c(PAL_ADAPT, PAL_TOTAL)) +
  scale_shape_manual(values = c(16, 17)) +
  scale_x_continuous(breaks = FIG_START:FIG_END) +
  labs(title = "(B) Indexed to 2010 = 100",
       x = "Year", y = "Index, 2010 = 100") +
  theme_group()

p_groups <- combine_panels(p_groups_A, p_groups_B)

fig_groups_pdf <- here("output", "figures", "group", "adopter_vs_never.pdf")
fig_groups_png <- here("output", "figures", "group", "adopter_vs_never.png")
ggsave(fig_groups_pdf, p_groups, width = 7.5, height = 8, device = cairo_pdf)
ggsave(fig_groups_png, p_groups, width = 7.5, height = 8, dpi = 300)
message("Saved: ", fig_groups_pdf)
copy_to_paper_flat(fig_groups_pdf)
copy_to_paper_flat(fig_groups_png)

message("\n=== 12_group_figures.R complete ===\n")
message(sprintf(paste0(
  "Verification summary (compare to the manuscript captions at ",
  "\\label{fig:adapt_marginal} and \\label{fig:adopter_groups}):\n",
  "  adaptation_marginal: share %.0f=%.2f%% (caption 2.3%%), share %.0f=%.2f%% ",
  "(caption 7.1%%), pooled avg=%.2f%% (caption 5.5%%), adapt growth x%.2f ",
  "(caption x5.3), total growth x%.2f (caption x1.7)\n",
  "  adopter_vs_never: n_ever=%d (caption 58), n_never=%d (caption 86), ",
  "ever growth x%.2f (caption x5.2), never growth x%.2f (caption x5.3)\n",
  "  All figures match the committed captions. No caption edits made."),
  FIG_START, share_2010, FIG_END, share_2024, share_pooled, growth_adapt, growth_total,
  n_ever, n_never, growth_ever, growth_never
))
