# Replication package for "A Larger Slice, Not a Larger Pie: National Adaptation Plans and the Composition of Adaptation Finance"

Pierre Beaucoral, Michaël Goujon and Sébastien Marchand (Université Clermont Auvergne, CNRS, IRD, CERDI).

## Overview

The paper estimates the effect of submitting a National Adaptation Plan (NAP) on the adaptation finance a developing country receives. It uses a panel of 145 recipient countries over 2009–2024, built from OECD Creditor Reporting System (CRS) activity data, and the Callaway and Sant'Anna (2021) staggered difference-in-differences estimator.

The code is written in R. One master script, `run_all.R`, runs 13 stage scripts in `code/` and writes every table (`.tex`) and figure (`.png`/`.pdf`) of the paper to `output/`.

There are two ways to run the package:

- **Default (about 10 minutes).** Start from the CRS-derived panels shipped in `data/processed/`. No external download is needed. Every exhibit is regenerated except the three scope tables built from the raw CRS files (they are shipped) and the two tables that need EM-DAT data, which cannot be redistributed (see below).
- **Full rebuild from raw data (`NAP_FROM_RAW=1`).** Download the raw OECD CRS files (about 4.8 GB). `run_all.R` checks their SHA-256 checksums against the vintage used in the paper, rebuilds `data/processed/` from scratch, compares every rebuilt file with the shipped one, and then runs the whole analysis on the rebuilt data.

The exhibits shipped in `output/` are the ones in the paper. A default run reproduces every regenerated `.tex` table exactly; the only difference is the date-stamp comment line (e.g. `% Wed Sep 23 10:52:40 2026`) that the `xtable` package writes at the top of some tables.

## Data availability and provenance

### Statement about rights

- [x] I certify that the author(s) of the manuscript have legitimate access to and permission to use the data used in this manuscript.
- [x] I certify that the author(s) of the manuscript have documented permission to redistribute/publish the data contained within this replication package. Every source except EM-DAT may be redistributed with attribution. EM-DAT is not included (see below).

### Summary of availability

- [ ] All data **are** publicly available.
- [x] **Some** data **cannot be made** publicly available: EM-DAT (free for registered non-commercial users, but not redistributable). It enters only two appendix tables.
- [ ] **No data can be made** publicly available.

### Data sources

| Source | Files in this package | Shipped | Licence / terms | Access |
|---|---|:-:|---|---|
| OECD Creditor Reporting System (CRS), bulk files 2007–2024, downloaded 1 April 2026 | Raw: none (`data/raw/CRS/README.md`, `crs_checksums.csv`). Derived: `data/processed/*.csv`, `data/processed/crs_adaptation_activities/*.csv.gz` | Derived panels: yes. Raw files: no (4.8 GB) | OECD terms and conditions: free reuse with attribution | [OECD Data Explorer](https://data-explorer.oecd.org), see `data/raw/CRS/README.md` |
| UNFCCC NAP Central, list of submitted NAPs (snapshot 4 June 2026) | `data/raw/shared_nap_data/nap_information.csv` | Yes | UNFCCC public information; compiled by the authors, CC BY 4.0 | <https://napcentral.org/submitted-naps> |
| UNFCCC "NAPAs received" list, transcribed 15 September 2026 | `data/raw/napa/napa_list_unfccc.csv`, `Submitted_NAPAs_UNFCCC.pdf` (archived page) | Yes | UNFCCC public information; transcription by the authors, CC BY 4.0 | <https://unfccc.int/topics/resilience/workstreams/national-adaptation-programmes-of-action/napas-received>, see `data/raw/napa/README.md` |
| World Bank World Development Indicators (population `SP.POP.TOTL`, GDP `NY.GDP.MKTP.KD` and six other series), pulled 10 June 2026 | `data/raw/wdi_cache/*.rds` | Yes | CC BY 4.0 | WDI API via the `WDI` R package |
| Worldwide Governance Indicators (government effectiveness `GOV_WGI_GE.EST` and five other WGI estimates), pulled 10 June 2026 | `data/raw/wdi_cache/wdi_wgi_ge.rds`, `wdi_additional.rds` | Yes | CC BY 4.0 | WDI API via the `WDI` R package |
| FERDI Physical Vulnerability to Climate Change Index (PVCCI), copy obtained 4 June 2026 (FERDI does not version the file; SHA-256 below) | `data/raw/PVCCI.csv` | Yes | FERDI, free for research use with attribution; not covered by this package's CC BY 4.0 grant | <https://ferdi.fr/en/indicators/the-physical-vulnerability-to-climate-change-index-pvcci> |
| World Bank historical income classifications (OGHIST), FY2013 column, downloaded 15 September 2026 | `data/raw/oghist/OGHIST.xlsx` | Yes | CC BY 4.0 | <https://datacatalogfiles.worldbank.org/ddh-published/0037712/DR0090754/OGHIST.xlsx> |
| UN list of Least Developed Countries, 2013 and 2024 | Written into `code/05_heterogeneity.R` (§ LDC classification) and copied in `code/14_hazard_napa.R` §6 | Yes (in code) | UN public information | UN DESA / Committee for Development Policy, <https://www.un.org/development/desa/dpad/least-developed-country-category.html> |
| EM-DAT international disaster database, version 2026-09-11 | None | **No** | EM-DAT terms of use: no redistribution, no derivative databases | Free registration at <https://public.emdat.be>; exact query in `data/raw/emdat/README.md` |

Details for each source:

**OECD CRS.** The adaptation outcome counts every CRS activity whose Rio adaptation marker is principal (2) or significant (1) (`ClimateAdaptation %in% c(1, 2)`). Amounts are commitments in constant US dollars (`USD_Commitment_Defl`). Regional and unspecified recipient codes are excluded. Stage 01 aggregates the raw files to the panels in `data/processed/`. Stage 09 needs activity-level records and reads them from `data/processed/crs_adaptation_activities/` (one gzip-compressed CSV per year, 2009–2024, all activities for the 145 panel recipients). The OECD revises past CRS years, so a new download will usually not match the April 2026 vintage. `data/raw/CRS/crs_checksums.csv` gives the size and SHA-256 of each file used.

**NAP and NAPA lists.** NAP submission dates come from NAP Central. The 51 NAPAs were transcribed by hand from the UNFCCC page, which blocks scripted access. The transcription was checked 51/51 against the archived PDF of the page. SHA-256: `nap_information.csv` `95699a5675a58aaf9cef290c11dee234e3af8d49d1ae5e6e7847e137f294c525`; `napa_list_unfccc.csv` `bcf34444fc37e00cc2be8ea47cc325488d01bad0a629526ccac817a45359952d`.

**WDI and WGI.** Stage 01 pulls these series from the World Bank API and caches them in `data/raw/wdi_cache/`. The shipped caches pin the 10 June 2026 vintage. Stage 01 reads the caches and does not pull again as long as they are present.

**PVCCI.** SHA-256 `74200d420fe4af010cd99ed59d23f7d3be15aeae5ca9af8dd7ba2ed6d6f425ad`. Two columns: `recipient_name`, `PVCCI`.

**OGHIST.** `code/05_heterogeneity.R` uses the FY2013 column of the sheet "Country Analytical History", the classification announced in July 2012 and therefore fixed before any NAP in the sample. It checks the file's SHA-256 (`17eb9e67b2eaf7d489ceb303f39f53697ce5b6238ed13c0b11c0846c33024d7b`) and warns if it differs. It downloads the file if it is missing.

**UN LDC lists.** The heterogeneity split uses the UN list in force on 1 January 2013 (49 countries), a pre-treatment vintage. The 2024 list (45 countries) is used only to report which countries graduated. Both lists are written in the code, with the graduation dates they rely on.

**EM-DAT.** Two appendix tables control for, or test NAP timing against, natural-hazard realisations from EM-DAT. EM-DAT's terms of use forbid redistributing the data or derived databases, so the package includes neither the extract nor the recipient-year panel built from it. `data/raw/emdat/README.md` gives the exact query (Disaster Group = Natural, all countries, Start Year 2000–2026, version 2026-09-11, 10,896 records), the conversion to `data/raw/emdat/emdat.csv`, and the checksums of the authors' files. Without it, stage 14 prints which two tables it cannot regenerate and leaves the shipped copies in place. The rest of the pipeline runs normally.

### Citations for the data

- OECD (2026). *Creditor Reporting System (CRS)*. OECD Data Explorer. Bulk files downloaded 1 April 2026.
- UNFCCC (2026). *NAP Central: Submitted National Adaptation Plans*. <https://napcentral.org/submitted-naps>. Accessed 4 June 2026.
- UNFCCC (2026). *NAPAs received*. <https://unfccc.int/topics/resilience/workstreams/national-adaptation-programmes-of-action/napas-received>. Accessed 15 September 2026.
- World Bank (2026). *World Development Indicators*. Washington, DC: The World Bank. Accessed 10 June 2026.
- Kaufmann, D. and A. Kraay (2026). *Worldwide Governance Indicators*. Washington, DC: The World Bank. Accessed via WDI, 10 June 2026.
- World Bank (2026). *World Bank Country and Lending Groups: Historical Classification by Income* (OGHIST). Accessed 15 September 2026.
- Feindouno, S., P. Guillaumont and C. Simonet (2020). "The Physical Vulnerability to Climate Change Index: An Index to Be Used for International Policy." *Ecological Economics* 176: 106752. Data: FERDI.
- United Nations, Committee for Development Policy and UN DESA (2013, 2024). *List of Least Developed Countries*.
- EM-DAT, CRED / UCLouvain, Brussels, Belgium — www.emdat.be (version 2026-09-11). Delforge, D. et al. (2025). "EM-DAT: the Emergency Events Database." *International Journal of Disaster Risk Reduction* 124: 105509.

## Computational requirements

### Software

- R 4.6.0. Every package version is pinned in `renv.lock` (did 2.5.0, HonestDiD 0.2.8, data.table 1.18.4, dplyr 1.2.1, ggplot2 4.0.3, fixest 0.14.1, DIDmultiplegtDYN 2.3.3, polars 1.11.0 from r-multiverse, bacondecomp 0.1.1, cowplot 1.2.0, WDI 2.7.10, readxl 1.4.5, digest 0.6.39, here 1.0.2). Install them once, from this folder:

  ```r
  install.packages("renv")
  renv::restore()
  ```

  When `renv::restore()` asks, choose "Activate the project". This installs the pinned versions into a project library and creates `.Rprofile` and `renv/`, so that every later R session started in this folder (including `Rscript run_all.R`) uses them. Then run `Rscript run_all.R` from this folder.

  Several stages stop with an explanatory error if `did` is older than 2.5.0 (see "Package-version note" below).
- The randomization-inference stage (08) uses `parallel::mclapply()`, which forks and does not parallelise on Windows. On Windows it runs on one core; only a full recompute of the draws (not needed by default) is affected.

### Hardware and memory

The results were produced on an Apple M2 Pro (10 cores), 32 GB RAM, macOS 26 (Darwin 25.6.0). No GPU is used.

- Default mode peaks at about 6 GB of RAM (measured peak resident memory 5.9 GB in stage 09, which reads the 3.65 million-row activity extract; 6.55 GB for the whole run). **8 GB of free RAM is recommended.**
- Full rebuild: stage 01 reads each raw CRS year (up to 520 MB of text) with `data.table::fread()`. Keep **at least 8 GB of free RAM** (16 GB recommended), and run the stages one after another, as `run_all.R` does. Do not run two CRS-reading R sessions at the same time.
- Disk: the package takes about 115 MB. The raw CRS files add 4.8 GB.

### Runtime

Measured on the machine above, stages run one at a time:

| Stage | Default mode | From-raw mode | Notes |
|---|--:|--:|---|
| Checksums of raw CRS files | — | < 0.5 min | SHA-256 of 4.8 GB |
| `01_prepare_data.R` | skipped | 1.0 min | Reads 18 CRS years; can take several minutes when the files are not already in the operating system's disk cache |
| `02_descriptive_stats.R` | 0.1 min | 0.1 min | |
| `03_main_results.R` | 0.1 min | 0.1 min | |
| `04_robustness.R` | 6.5 min | 6.6 min | HonestDiD sensitivity, dCDH, Goodman-Bacon |
| `05_heterogeneity.R` | 0.1 min | 0.1 min | |
| `13_cohort_anticipation.R` | 0.1 min | 0.1 min | |
| `14_hazard_napa.R` | < 0.1 min | < 0.1 min | |
| `07_principal_and_share.R` | 0.1 min | 0.1 min | |
| `08_randomization_inference.R` | < 0.1 min | < 0.1 min | With the shipped draws. Full recompute (draws deleted): about 29 min on 8 cores |
| `09_remarking_decomposition.R` | 2.2 min | 2.3 min | From-raw: rebuilds the activity extract from 16 CRS years (0.4–10 min depending on disk cache) |
| `10_model_tests.R` | < 0.1 min | < 0.1 min | |
| `11_base_year_sensitivity.R` | 0.1 min | 0.1 min | |
| `12_group_figures.R` | < 0.1 min | < 0.1 min | |
| **Total** | **about 9.5 min** | **about 11 min** | Add about 29 min in either mode to recompute the permutation draws |

(Measured on 23 September 2026. Stage 08 timings are for a separate run of the stage with the shipped draws; the full runs above used `NAP_SKIP_RI=1`.)

## Instructions for replicators

1. Install R 4.6.0 and restore the packages: open R in this folder and run `renv::restore()`. When it asks, choose "Activate the project"; then run `Rscript run_all.R` from this folder.
2. From this folder (it contains a `.here` file, so all paths are resolved relative to it), run:

   ```bash
   Rscript run_all.R
   ```

   This is the default mode: it starts at stage 02 from the shipped panels in `data/processed/`. Stage 08 reuses the shipped permutation draws (`output/tables/randomization/ri_draws.csv`) and only rebuilds its table and figure.
3. Compare the regenerated `output/tables/**/*.tex` with the shipped versions, for example with `git diff` if you put the folder under version control before running, or against `MANIFEST.csv` (SHA-256 of every shipped file). Tables match the shipped ones except for the `xtable` date-stamp comment line. On the authors' machine the PNG figures are also byte-identical and the two PDF figures differ only in embedded metadata (they render identically). On other systems, fonts and graphics libraries can change figure files at the byte level without changing their content.

Options (environment variables, combinable):

| Variable | Effect |
|---|---|
| `NAP_SKIP_RI=1` | Skip stage 08 entirely (its table and figure are then the shipped ones). |
| `NAP_START_AT=<stage file>` | Resume at a stage, e.g. `NAP_START_AT=09_remarking_decomposition.R Rscript run_all.R`. Earlier outputs are used as they are. |
| `NAP_FROM_RAW=1` | Full rebuild from raw CRS files (below). |

To recompute the randomization-inference draws from scratch (about 29 minutes on 8 cores), delete `output/tables/randomization/ri_draws.csv` before running. The draws are deterministic (see "Seeds and determinism"), so the recomputed file and table match the shipped ones.

### Full rebuild from the raw CRS files

1. Download the 18 CRS bulk files for 2007–2024 as described in `data/raw/CRS/README.md` and put them in `data/raw/CRS/`.
2. Make sure the machine is online: stage 01 calls the World Bank API catalogue (`WDI::WDIcache()`) even though the data series themselves are read from the shipped caches.
3. Run:

   ```bash
   NAP_FROM_RAW=1 Rscript run_all.R
   ```

   `run_all.R` then:
   - verifies the SHA-256 of each raw CRS file against `data/raw/CRS/crs_checksums.csv` and prints `MATCH` or `DIFFERENT VINTAGE` for each year. A different vintage triggers a warning, not an error, because the OECD revises past years;
   - moves the shipped panels to `data/processed_shipped/` and rebuilds `data/processed/` with stage 01. Stage 09 rebuilds the activity-level extract from the raw files (about 10 minutes);
   - after each stage, compares every newly rebuilt data file with its shipped counterpart and prints `PASS (byte-identical)`, `PASS (equal within 1e-10 after sorting)` or `FAIL`, followed by a summary;
   - runs the rest of the analysis on the rebuilt data. With the April 2026 CRS vintage every comparison passes and the tables match the shipped ones.

   A second from-raw run keeps `data/processed_shipped/` as the reference and rebuilds `data/processed/` again. To go back to the default mode, delete `data/processed/` and rename `data/processed_shipped/` to `data/processed/`.

### Seeds and determinism

Every stage script sets a global seed at the top (`set.seed(20240601)`). The scripts also call `set.seed(1242)` immediately before each bootstrap estimator, so several seeds appear in each script. This is deliberate and needed to reproduce the published standard errors. In the `did` package, `aggte()` draws its own multiplier bootstrap, so an aggregated standard error depends on the random-number state when `aggte()` is called, not only when `att_gt()` is called. Reseeding before each estimator makes each reported standard error independent of the code that runs before it. The headline fits are computed once in stage 03, saved in `output/fits/`, and read by later stages instead of being re-estimated, so each specification has exactly one standard error in the paper.

The randomization-inference draws are generated sequentially in the parent R process (`set.seed(1242)` once per design) before any forking. The forked workers draw no random numbers. The draws therefore do not depend on the number of cores.

`ri_draws.csv` stores each double both as readable decimal text and as an exact hexadecimal column (`*_hex`), which stage 08 reads, so the cached draws reload bit for bit.

### Package-version note

The published numbers use `did` 2.5.0. Earlier versions (for example 2.3.0) compute the multiplier bootstrap differently on unbalanced panels, and the panel is unbalanced (South Sudan is observed in 14 of 16 years). Across versions:

- point estimates of the doubly robust specifications are identical;
- standard errors and pre-trend p-values differ on unbalanced panels (balanced subsamples, such as the non-LDC split, match exactly);
- the placebo test, which uses the regression-adjustment estimator (`est_method = "reg"`) on a truncated panel, gives a **different point estimate** as well as a different standard error.

So match `did` 2.5.0 through `renv::restore()` to reproduce the tables exactly.

## List of programs

Each stage runs in a fresh R process, in the order below (stage 06 is an internal diagnostic that is not part of the paper and is not included).

| Order | Script | Reads | Writes |
|---|---|---|---|
| 1 | `code/01_prepare_data.R` (from-raw mode only) | `data/raw/CRS/`, `data/raw/shared_nap_data/`, `data/raw/PVCCI.csv`, `data/raw/wdi_cache/` | `data/processed/*.csv`; `output/tables/scope/` |
| 2 | `code/02_descriptive_stats.R` | `data/processed/adaptationNAP.csv`, `adaptationNAP_donortype_wgi.csv`, `simple_panel_wgi.csv`, `donor_list.csv`, `donor_totals.csv`, `data/raw/PVCCI.csv` | descriptive tables and figures at the top level of `output/tables/`, `output/figures/` |
| 3 | `code/03_main_results.R` | `simple_panel_wgi.csv` | `output/tables/cohorts_dropped/att_combined_wide.tex`, `extensive_margin/`, `nap_cohorts.tex`, main figures; headline fits in `output/fits/` |
| 4 | `code/04_robustness.R` | `simple_panel_wgi.csv`, `mitigation_panel.csv`, `output/fits/` | `cohorts_retained/`, `notyettreated/`, `outlier_india/`, `units_zeros/`, `placebo/`, `mitigation/`, `bacon/`, `dcdh/`, HonestDiD tables in `cohorts_dropped/` |
| 5 | `code/05_heterogeneity.R` | `simple_panel_wgi.csv`, `data/raw/oghist/OGHIST.xlsx`, `output/fits/` | `heterogeneity/` |
| 6 | `code/13_cohort_anticipation.R` | `simple_panel_wgi.csv`, `output/fits/` | `cohort_battery/`, `anticipation/` |
| 7 | `code/14_hazard_napa.R` | `simple_panel_wgi.csv`, `emergency_response_panel.csv`, `data/raw/napa/`, `data/raw/emdat/emdat.csv` (if present) | `hazard/`, `napa/`; `data/processed/napa_list.csv` |
| 8 | `code/07_principal_and_share.R` | `simple_panel_wgi.csv`, `output/fits/` | `principal_share/` |
| 9 | `code/08_randomization_inference.R` | `simple_panel_wgi.csv`, `output/fits/`, `output/tables/randomization/ri_draws.csv` (cache) | `randomization/` |
| 10 | `code/09_remarking_decomposition.R` | `simple_panel_wgi.csv`, `data/processed/crs_adaptation_activities/` (or raw CRS if absent), `output/fits/` | `remarking/`; run log in `output/logs/` |
| 11 | `code/10_model_tests.R` | `simple_panel_wgi.csv`, `output/fits/` | `model_tests/` |
| 12 | `code/11_base_year_sensitivity.R` | `simple_panel_wgi.csv`, `output/fits/` | `base_year/` |
| 13 | `code/12_group_figures.R` | `simple_panel_wgi.csv` | `output/figures/group/` |

After the last stage, `run_all.R` wraps each table's `tabular` in `\adjustbox{max width=\textwidth}` and sets floats to `[H]` (idempotent). Folder names under `output/tables/` and `output/figures/` are the same.

Processed data files (`data/processed/`): `simple_panel_wgi.csv` (estimation panel, one row per recipient-year), `adaptationNAP.csv` (descriptive panel), `adaptationNAP_donortype_wgi.csv` (recipient-year-donor type panel), `mitigation_panel.csv`, `emergency_response_panel.csv` (CRS emergency-response aid, used only for a NAP-timing check), `adaptation_panel_oda_only.csv` (ODA-only variant, not used in the paper), `donor_list.csv`, `donor_totals.csv`, `donor_recipient_year_adaptation.csv` (not used in the paper), `napa_list.csv` (copy of the NAPA list written by stage 14), and `crs_adaptation_activities/` (activity-level extract for stage 09).

## List of tables and figures

Numbers refer to the paper as compiled on 23 September 2026 (Appendix A: model and descriptive appendices; Appendix B: supplementary robustness exhibits). All files are in `output/`.

| Exhibit | Location | File in `output/` | Produced by |
|---|---|---|---|
| Figure 1 | Main text | `figures/group/adaptation_marginal.pdf` | `code/12_group_figures.R` |
| Table 2 | Main text | `tables/cohorts_dropped/att_combined_wide.tex` | `code/03_main_results.R` |
| Table 3 | Main text | `tables/principal_share/within_share_wide.tex` | `code/07_principal_and_share.R` |
| Table 4 | Main text | `tables/principal_share/share_reconciliation.tex` | `code/07_principal_and_share.R` |
| Figure 2 | Main text | `figures/cohorts_dropped/did_combined_es_wgi.png` | `code/03_main_results.R` |
| Table 5 | Main text | `tables/placebo/att_placebo.tex` | `code/04_robustness.R` |
| Table 6 | Main text | `tables/mitigation/att_mitigation.tex` | `code/04_robustness.R` |
| Table 7 | Main text | `tables/randomization/ri_pvalues.tex` | `code/08_randomization_inference.R` |
| Table 8 | Main text | `tables/cohorts_dropped/honestdid_rm.tex` | `code/04_robustness.R` |
| Table 9 | Main text | `tables/base_year/pretrend_tests_full.tex` | `code/11_base_year_sensitivity.R` |
| Table 10 | Main text | `tables/base_year/base_year_sensitivity.tex` | `code/11_base_year_sensitivity.R` |
| Table 11 | Main text | `tables/notyettreated/att_notyettreated_wide.tex` | `code/04_robustness.R` |
| Table 12 | Main text | `tables/cohorts_retained/att_combined_wide.tex` | `code/04_robustness.R` |
| Table 13 | Main text | `tables/units_zeros/diagnostic_units_zeros.tex` | `code/04_robustness.R` |
| Table 14 | Main text | `tables/outlier_india/att_combined_wide.tex` | `code/04_robustness.R` |
| Table 15 | Main text | `tables/principal_share/principal_wide.tex` | `code/07_principal_and_share.R` |
| Table 16 | Main text | `tables/remarking/att_remarking_margins.tex` | `code/09_remarking_decomposition.R` |
| Table 17 | Main text | `tables/remarking/att_remarking_margins_unlinked.tex` | `code/09_remarking_decomposition.R` |
| Table 18 | Main text | `tables/cohort_battery/att_loco.tex` | `code/13_cohort_anticipation.R` |
| Table 19 | Main text | `tables/anticipation/att_anticipation.tex` | `code/13_cohort_anticipation.R` |
| Table 20 | Main text | `tables/anticipation/att_eventdate.tex` | `code/13_cohort_anticipation.R` |
| Table 21 | Main text | `tables/extensive_margin/att_extensive.tex` | `code/03_main_results.R` |
| Table 22 | Main text | `tables/heterogeneity/het_mde.tex` | `code/05_heterogeneity.R` |
| Table 23 | Main text | `tables/heterogeneity/donor_type/att_combined_wide.tex` | `code/05_heterogeneity.R` |
| Table 24 | Main text | `tables/heterogeneity/governance/het_gov_wide.tex` | `code/05_heterogeneity.R` |
| Table 25 | Main text | `tables/heterogeneity/ldc/het_ldc_wide.tex` | `code/05_heterogeneity.R` |
| Table 26 | Main text | `tables/heterogeneity/income_group/het_income_wide.tex` | `code/05_heterogeneity.R` |
| Table 27 | Main text | `tables/heterogeneity/het_difference_tests.tex` | `code/05_heterogeneity.R` |
| Table A.1 | Appendix | `tables/model_tests/lemma2_size.tex` | `code/10_model_tests.R` |
| Figure A.1 | Appendix | `figures/model_tests/fig_lemma2_size.png` | `code/10_model_tests.R` |
| Table A.2 | Appendix | `tables/model_tests/alpha_half.tex` | `code/10_model_tests.R` |
| Figure A.2 | Appendix | `figures/group/adopter_vs_never.pdf` | `code/12_group_figures.R` |
| Table A.3 | Appendix | `tables/stats_des.tex` | `code/02_descriptive_stats.R` |
| Figure A.3 | Appendix | `figures/top_donors.png` | `code/02_descriptive_stats.R` |
| Figure A.4 | Appendix | `figures/recipient_map.png` | `code/02_descriptive_stats.R` |
| Figure A.5 | Appendix | `figures/nap_status_map.png` | `code/02_descriptive_stats.R` |
| Table A.4 | Appendix | `tables/nap_cohorts.tex` | `code/03_main_results.R` |
| Table A.5 | Appendix | `tables/balance_adopters.tex` | `code/02_descriptive_stats.R` |
| Table A.6 | Appendix | `tables/cohort_battery/cohort_contributions.tex` | `code/13_cohort_anticipation.R` |
| Table A.7 | Appendix | `tables/anticipation/cohort_redating_crosstab.tex` | `code/13_cohort_anticipation.R` |
| Table A.8 | Appendix | `tables/cohort_battery/att_cohort_battery.tex` | `code/13_cohort_anticipation.R` |
| Table A.9 | Appendix | `tables/cohort_battery/listwise_losses_13.tex` | `code/13_cohort_anticipation.R` |
| Table A.10 | Appendix | `tables/heterogeneity/listwise_losses.tex` | `code/05_heterogeneity.R` |
| Table A.11 | Appendix | `tables/napa/prior_napa_ldc_crosstab.tex` | `code/14_hazard_napa.R` |
| Table A.12 | Appendix | `tables/scope/flow_type_shares.tex` | `code/01_prepare_data.R` |
| Table A.13 | Appendix | `tables/scope/sample_funnel.tex` | `code/01_prepare_data.R` |
| Table A.14 | Appendix | `tables/scope/regional_exclusion.tex` | `code/01_prepare_data.R` |
| Table A.15 | Appendix | `tables/dcdh/att_dcdh.tex` | `code/04_robustness.R` |
| Figure A.6 | Appendix | `figures/dcdh/did_dcdh_es.png` | `code/04_robustness.R` |
| Table A.16 | Appendix | `tables/bacon/bacon_decomp.tex` | `code/04_robustness.R` |
| Figure A.7 | Appendix | `figures/bacon/bacon_scatter.png` | `code/04_robustness.R` |
| Table A.17 | Appendix | `tables/remarking/att_remarking_counts.tex` | `code/09_remarking_decomposition.R` |
| Table A.18 | Appendix | `tables/remarking/att_remarking_exclusions.tex` | `code/09_remarking_decomposition.R` |
| Table A.19 | Appendix | `tables/remarking/tab_remarking_sectors.tex` | `code/09_remarking_decomposition.R` |
| Figure A.8 | Appendix | `figures/remarking/fig_remarking_sectors.png` | `code/09_remarking_decomposition.R` |
| Table A.20 | Appendix | `tables/remarking/remarking_flag_shares_by_year.tex` | `code/09_remarking_decomposition.R` |
| Table A.21 | Appendix | `tables/remarking/remarking_flag_shares_by_group.tex` | `code/09_remarking_decomposition.R` |
| Table A.22 | Appendix | `tables/remarking/att_remarking_exclusions_regex_comparison.tex` | `code/09_remarking_decomposition.R` |
| Figure B.1 | Appendix | `figures/nap_adoption_timeline.png` | `code/02_descriptive_stats.R` |
| Figure B.2 | Appendix | `figures/principal_share/fig_within_share_es.png` | `code/07_principal_and_share.R` |
| Figure B.3 | Appendix | `figures/principal_share/fig_share_reconciliation.png` | `code/07_principal_and_share.R` |
| Figure B.4 | Appendix | `figures/cohorts_dropped/did_combined_cohort_wgi.png` | `code/03_main_results.R` |
| Figure B.5 | Appendix | `figures/placebo/did_placebo_es.png` | `code/04_robustness.R` |
| Figure B.6 | Appendix | `figures/mitigation/did_mitigation_es.png` | `code/04_robustness.R` |
| Table B.1 | Appendix | `tables/hazard/nap_timing_vs_emdat_hazard.tex` | `code/14_hazard_napa.R` (needs EM-DAT) |
| Table B.2 | Appendix | `tables/hazard/nap_timing_vs_humanitarian_aid.tex` | `code/14_hazard_napa.R` |
| Table B.3 | Appendix | `tables/hazard/att_hazard_controls.tex` | `code/14_hazard_napa.R` (needs EM-DAT) |
| Table B.4 | Appendix | `tables/napa/att_napa_falsification.tex` | `code/14_hazard_napa.R` |
| Table B.5 | Appendix | `tables/napa/att_prior_napa_split.tex` | `code/14_hazard_napa.R` |
| Figure B.7 | Appendix | `figures/randomization/fig_ri_distributions.png` | `code/08_randomization_inference.R` |
| Table B.6 | Appendix | `tables/cohorts_dropped/honestdid_prepriods.tex` | `code/04_robustness.R` |
| Table B.7 | Appendix | `tables/cohorts_dropped/honestdid_sd.tex` | `code/04_robustness.R` |
| Table B.8 | Appendix | `tables/base_year/pretrend_cells_2021.tex` | `code/11_base_year_sensitivity.R` |
| Figure B.8 | Appendix | `figures/base_year/fig_pretrend_cells_by_cohort.png` | `code/11_base_year_sensitivity.R` |
| Figure B.9 | Appendix | `figures/base_year/fig_es_full_window.png` | `code/11_base_year_sensitivity.R` |
| Figure B.10 | Appendix | `figures/notyettreated/did_notyettreated_es.png` | `code/04_robustness.R` |
| Figure B.11 | Appendix | `figures/cohorts_retained/did_combined_cohort_wgi.png` | `code/04_robustness.R` |
| Table B.9 | Appendix | `tables/principal_share/principal_retained.tex` | `code/07_principal_and_share.R` |
| Table B.10 | Appendix | `tables/cohort_battery/att_drop2024.tex` | `code/13_cohort_anticipation.R` |
| Table B.11 | Appendix | `tables/cohort_battery/att_balance.tex` | `code/13_cohort_anticipation.R` |
| Table B.12 | Appendix | `tables/cohort_battery/att_2x2.tex` | `code/13_cohort_anticipation.R` |
| Table B.13 | Appendix | `tables/cohort_battery/att_conditioning.tex` | `code/13_cohort_anticipation.R` |
| Table B.14 | Appendix | `tables/anticipation/att_placebo_ladder.tex` | `code/13_cohort_anticipation.R` |
| Table B.15 | Appendix | `tables/heterogeneity/donor_type/zero_shares.tex` | `code/05_heterogeneity.R` |
| Figure B.12 | Appendix | `figures/heterogeneity/donor_type/did_donor_type_es.png` | `code/05_heterogeneity.R` |
| Figure B.13 | Appendix | `figures/heterogeneity/governance/did_governance_es.png` | `code/05_heterogeneity.R` |
| Figure B.14 | Appendix | `figures/heterogeneity/ldc/did_ldc_es.png` | `code/05_heterogeneity.R` |
| Figure B.15 | Appendix | `figures/heterogeneity/income_group/did_income_es.png` | `code/05_heterogeneity.R` |

Table 1 (mapping of the IPCC AR6 risk components onto aid-allocation roles) is typed directly in the manuscript and is not produced by code. Every other table and figure is listed above.

Exhibits in `output/` that the paper does not use are regenerated as well (for example `output/figures/climate_finance_evolution.png`, `output/tables/finance_change.tex`, `output/tables/nap_regional.tex`, `output/tables/list.tex`).

## Licence

Code (`run_all.R`, `code/`): MIT License. Data compiled by the authors, the processed panels and the exhibits: Creative Commons Attribution 4.0 International (CC BY 4.0). Third-party data keep their own terms (table above). In particular, the PVCCI values in `data/raw/PVCCI.csv` and in the processed panels (`pvcci` and the columns derived from it) remain under FERDI's terms and are excluded from the CC BY 4.0 grant. EM-DAT data are not included and may not be redistributed. See `LICENSE`.

## How to cite

Beaucoral, P., M. Goujon and S. Marchand (2026). "A Larger Slice, Not a Larger Pie: National Adaptation Plans and the Composition of Adaptation Finance." Working paper, CERDI, Université Clermont Auvergne.

Replication package: Beaucoral, P., M. Goujon and S. Marchand (2026). Replication package for "A Larger Slice, Not a Larger Pie: National Adaptation Plans and the Composition of Adaptation Finance." **[REPOSITORY DOI: TO BE ADDED ON DEPOSIT]**

## Acknowledgements

This work was supported by the Agence Nationale de la Recherche of the French government through the program "Investissements d'avenir" (ANR-10-LABX-14-01).
