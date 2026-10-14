# Replication package for "Do National Adaptation Plans Attract Adaptation Finance? Evidence from the Composition of Aid"

Pierre Beaucoral, Michaël Goujon and Sébastien Marchand (Université Clermont Auvergne, CNRS, IRD, CERDI).

## Overview

The paper estimates the effect of submitting a National Adaptation Plan (NAP) on the adaptation finance a developing country receives. It uses a panel of 145 recipient countries over 2009–2024, built from OECD Creditor Reporting System (CRS) activity data, and the Callaway and Sant'Anna (2021) staggered difference-in-differences estimator.

The code is written in R. One master script, `run_all.R`, runs 13 stage scripts in `code/` and writes every table (`.tex`) and figure (`.png`/`.pdf`) of the paper to `output/`.

There are two ways to run the package:

- **Default (about 12 minutes with the shipped permutation draws).** Start from the CRS-derived panels shipped in `data/processed/`. No external download is needed. Every exhibit is regenerated except the three scope tables built from the raw CRS files (they are shipped) and the two tables that need EM-DAT data, which cannot be redistributed (see below).
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

**OECD CRS.** The adaptation outcome counts every CRS activity whose Rio adaptation marker is principal (2) or significant (1) (`ClimateAdaptation %in% c(1, 2)`). Amounts are commitments in constant US dollars (`USD_Commitment_Defl`). Regional and unspecified recipient codes are excluded. A recipient-year in which the country has no CRS record of any kind (after it leaves the DAC List of ODA recipients) is missing, not a year of zero adaptation finance. Stage 01 aggregates the raw files to the panels in `data/processed/`. Stage 09 needs activity-level records and reads them from `data/processed/crs_adaptation_activities/` (one gzip-compressed CSV per year, 2009–2024, all activities for the 145 panel recipients). The OECD revises past CRS years, so a new download will usually not match the April 2026 vintage. `data/raw/CRS/crs_checksums.csv` gives the size and SHA-256 of each file used.

**NAP and NAPA lists.** NAP submission dates come from NAP Central. Treatment is the first submission; the one entry that lists two postings in a single cell (Paraguay: May 2020 and July 2022) is dated to the first. The 51 NAPAs were transcribed by hand from the UNFCCC page, which blocks scripted access. The transcription was checked 51/51 against the archived PDF of the page. SHA-256: `nap_information.csv` `95699a5675a58aaf9cef290c11dee234e3af8d49d1ae5e6e7847e137f294c525`; `napa_list_unfccc.csv` `bcf34444fc37e00cc2be8ea47cc325488d01bad0a629526ccac817a45359952d`.

**WDI and WGI.** Stage 01 pulls these series from the World Bank API and caches them in `data/raw/wdi_cache/`. The shipped caches pin the 10 June 2026 vintage. Stage 01 reads the caches and does not pull again as long as they are present.

**PVCCI.** SHA-256 `74200d420fe4af010cd99ed59d23f7d3be15aeae5ca9af8dd7ba2ed6d6f425ad`. Two columns: `recipient_name`, `PVCCI`.

**OGHIST.** `code/05_heterogeneity.R` uses the FY2013 column of the sheet "Country Analytical History", the classification announced in July 2012 and therefore fixed before any NAP in the sample. It checks the file's SHA-256 (`17eb9e67b2eaf7d489ceb303f39f53697ce5b6238ed13c0b11c0846c33024d7b`) and warns if it differs. It downloads the file if it is missing.

**UN LDC lists.** The heterogeneity split uses the UN list in force on 1 January 2013 (49 countries), a pre-treatment vintage. The 2024 list (45 countries) is used only to report which countries graduated. Both lists are written in the code, with the graduation dates they rely on.

**EM-DAT.** Two appendix tables control for, or test NAP timing against, natural-hazard realisations from EM-DAT. EM-DAT's terms of use forbid redistributing the data or derived databases, so the package includes neither the extract nor the recipient-year panel built from it. `data/raw/emdat/README.md` gives the exact query (Disaster Group = Natural, all countries, Start Year 2000–2026, version 2026-09-11, 10,896 records), the conversion to `data/raw/emdat/emdat.csv`, and the checksums of the authors' files. Without it, stage 14 prints which two tables it cannot regenerate and leaves the shipped copies in place (these two tables are the only exhibits a stage keeps without regenerating them). The rest of the pipeline runs normally. To regenerate the two tables after adding `emdat.csv`, rerun the pipeline (`Rscript run_all.R`), or, if stage 03 has already run in this copy, only stage 14 (`NAP_START_AT=14_hazard_napa.R Rscript run_all.R`): stage 14 reads the headline fit that stage 03 stores in `output/fits/`, which is not shipped, so on a fresh copy a stage-14-only run stops.

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

- R 4.6.1. Every package version is pinned in `renv.lock` (did 2.5.0, HonestDiD 0.2.8, data.table 1.18.4, dplyr 1.2.1, ggplot2 4.0.3, fixest 0.14.1, DIDmultiplegtDYN 2.3.3, polars 1.11.0 from r-multiverse, bacondecomp 0.1.1, cowplot 1.2.0, WDI 2.7.10, readxl 1.4.5, digest 0.6.39, here 1.0.2, stringi 1.8.7). Install them once, from this folder:

  ```r
  install.packages("renv")
  renv::restore()
  ```

  When `renv::restore()` asks, choose "Activate the project". This installs the pinned versions into a project library and creates `.Rprofile` and `renv/`, so that every later R session started in this folder (including `Rscript run_all.R`) uses them. Then run `Rscript run_all.R` from this folder.

  Several stages stop with an explanatory error if `did` is older than 2.5.0 (see "Package-version note" below).
- The randomization-inference stage (08) uses `parallel::mclapply()`, which forks and does not parallelise on Windows. On Windows stage 08 sets itself to one core; only a full recompute of the draws (not needed by default) is affected, and the draws do not depend on the number of cores.
- On Windows, where R cannot fork, the doubly-robust cells that stages 08 and 13 isolate in a forked process (because `att_gt(est_method = "dr")` can segfault on thin panels) run in the main process instead, so such a segfault stops the stage rather than being recorded as a failed cell (stage 08: rerun on macOS or Linux, or keep the shipped draws).
- `run_all.R` runs each stage with the `Rscript` of the R installation that runs it (`R.home("bin")`), not whichever `Rscript` is first on the `PATH`.
- Locale. The results do not depend on the session's collation. The unit identifiers that `did` sorts on come from `code/functions/make_country_id.R`, which ranks recipient names with the ICU collation rules for `en_US` (package `stringi`) instead of the session's collation, under which names such as "Côte d'Ivoire" and "Türkiye" sort differently (C collation puts them after "Cuba" and "Tuvalu", language locales before "Croatia" and "Turkmenistan"); a machine-dependent order would change the identifiers. The ICU order is the one of language locales such as `en_US.UTF-8` and `fr_FR.UTF-8`, under which the published results were produced, and a session with C collation (for example `LC_COLLATE=C`) gives the same identifiers and results. The session must use a UTF-8 character encoding (any UTF-8 locale, including `C.UTF-8`): under the ASCII-only `C`/`POSIX` locale (`LC_ALL=C`), R cannot match the non-ASCII recipient names across files, and stage 01 stops on its consistency checks.

### System requirements

- GLPK (the GNU Linear Programming Kit), needed to build `Rglpk`, on which `HonestDiD` depends: `brew install glpk` (macOS), `apt-get install libglpk-dev` (Debian/Ubuntu); the CRAN Windows binary includes it.
- A Rust toolchain (<https://rustup.rs>, `rustc` >= 1.70 and `cargo`) when `polars` or `clarabel` is built from source: `polars` if the r-multiverse repository recorded in `renv.lock` has no binary for your platform, `clarabel` whenever CRAN has no binary (Linux).
- On Linux, where `renv::restore()` builds packages from source unless it is pointed at a binary repository (for example Posit Public Package Manager), also: C, C++ and Fortran compilers, `make` and `cmake` (Debian/Ubuntu: `build-essential gfortran cmake`), and the development headers of GMP, FreeType, Fontconfig, HarfBuzz, FriBidi, libpng, libjpeg, libtiff, libwebp, zlib and libuv (`libgmp-dev libfreetype6-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libpng-dev libjpeg-dev libtiff-dev libwebp-dev zlib1g-dev libuv1-dev`), plus OpenGL/GLU and X11 headers for `rgl` (`libgl1-mesa-dev libglu1-mesa-dev libx11-dev`; the pipeline itself never opens an rgl window and sets `RGL_USE_NULL=TRUE`).
- `bash`: on macOS and Linux, `run_all.R` runs each stage through `bash` (to copy its output to the log while keeping Rscript's exit status).

### Hardware and memory

The results were produced on an Apple M2 Pro (10 cores), 32 GB RAM, macOS 27 (Darwin 27.0.0). No GPU is used.

- Default mode peaks at about 6 GB of RAM (measured peak resident memory 5.9 GB in stage 09, which reads the 3.65 million-row activity extract; 6.55 GB for the whole run). **8 GB of free RAM is recommended.**
- Full rebuild: stage 01 reads each raw CRS year (up to 520 MB of text) with `data.table::fread()`. Keep **at least 8 GB of free RAM** (16 GB recommended), and run the stages one after another, as `run_all.R` does. Do not run two CRS-reading R sessions at the same time.
- Disk: the package takes about 115 MB. The raw CRS files add 4.8 GB.

### Runtime

Measured on the machine above, stages run one at a time:

| Stage | Default mode | From-raw mode | Notes |
|---|--:|--:|---|
| Checksums of raw CRS files | — | < 0.5 min | SHA-256 of 4.8 GB |
| `01_prepare_data.R` | skipped | 1.5 min | Reads 18 CRS years; can take several minutes when the files are not already in the operating system's disk cache |
| `02_descriptive_stats.R` | 0.1 min | 0.1 min | |
| `03_main_results.R` | 0.1 min | 0.1 min | |
| `04_robustness.R` | 9.0 min | 9.1 min | HonestDiD sensitivity, dCDH, Goodman-Bacon |
| `05_heterogeneity.R` | 0.1 min | 0.1 min | |
| `13_cohort_anticipation.R` | 0.2 min | 0.2 min | |
| `14_hazard_napa.R` | < 0.1 min | < 0.1 min | |
| `07_principal_and_share.R` | 0.1 min | 0.1 min | |
| `08_randomization_inference.R` | < 0.1 min | < 0.1 min | With the shipped draws. Full recompute (both draw files deleted): about 15 min on 8 cores (10.9 min for the level and global share, 3.6 min for the within-country share), longer under load |
| `09_remarking_decomposition.R` | 2.2 min | 3.4 min | From-raw: rebuilds the activity extract from 16 CRS years (0.4–10 min depending on disk cache) |
| `10_model_tests.R` | < 0.1 min | < 0.1 min | |
| `11_base_year_sensitivity.R` | 0.1 min | 0.1 min | |
| `12_group_figures.R` | < 0.1 min | < 0.1 min | |
| **Total** | **about 12 min** | **about 15 min** | Add about 15 min in either mode to recompute the permutation draws. A clean-room default run on 13 October 2026 that recomputed both draw files took 20.3 min (stage 04: 6.9 min; stage 08: 10.6 min; stage 09: 2.1 min; peak resident memory 6.53 GB) |

(Per-stage times measured on 29 September 2026: a default-mode run from stage 03 with the shipped draws, and a from-raw run in a fresh copy of this folder with the draws deleted. The 13 October 2026 clean-room total is in the last row.)

## Instructions for replicators

1. Install R 4.6.1 and restore the packages: open R in this folder and run `renv::restore()`. When it asks, choose "Activate the project"; then run `Rscript run_all.R` from this folder.
2. From this folder (it contains a `.here` file, so all paths are resolved relative to it), run:

   ```bash
   Rscript run_all.R
   ```

   This is the default mode: it starts at stage 02 from the shipped panels in `data/processed/`. Stage 08 reuses the shipped permutation draws (`output/tables/randomization/ri_draws.csv` for the level and the global share, `ri_draws_within.csv` for the within-country share) and only rebuilds its tables and figures. It reuses a file only if it was built on the current estimation sample: each file stores a SHA-256 hash of every recipient's (unit, adoption cohort, World Bank region) triple, and stage 08 also checks the sample size and the observed ATTs. If any of these differ, it recomputes those draws (see below).
3. Compare the regenerated `output/tables/**/*.tex` with the shipped versions, for example with `git diff` if you put the folder under version control before running, or against `MANIFEST.csv` (SHA-256 of every shipped file). Tables match the shipped ones except for the `xtable` date-stamp comment line. On the authors' machine the PNG figures are also byte-identical and the two PDF figures differ only in embedded metadata (they render identically). On other systems, fonts and graphics libraries can change figure files at the byte level without changing their content.

Options (environment variables, combinable):

| Variable | Effect |
|---|---|
| `NAP_SKIP_RI=1` | Skip stage 08, only if both cached draw files, `output/tables/randomization/ri_draws.csv` and `ri_draws_within.csv`, exist (its tables and figures are then the shipped ones). Without either file stage 08 runs anyway. |
| `NAP_START_AT=<stage file>` | Resume at a stage, e.g. `NAP_START_AT=09_remarking_decomposition.R Rscript run_all.R`. Earlier outputs are used as they are. Every stage after 03 reads the fits that stages 03, 04 and 07 write to `output/fits/`. That folder is not shipped (about 26 MB, regenerated with the stages), so on a fresh copy run the pipeline once from stage 03 or earlier (the default mode starts at 02) before resuming at a later stage. Cannot be combined with `NAP_FROM_RAW=1` (the rebuild starts at stage 01); `run_all.R` stops if both are set. |
| `NAP_FROM_RAW=1` | Full rebuild from raw CRS files (below). |

To recompute the randomization-inference draws from scratch (about 15 minutes on 8 cores), delete `output/tables/randomization/ri_draws.csv` and `ri_draws_within.csv` before running. The draws are deterministic (see "Seeds and determinism"), so the recomputed file and table match the shipped ones.

### Full rebuild from the raw CRS files

1. Download the 18 CRS bulk files for 2007–2024 as described in `data/raw/CRS/README.md` and put them in `data/raw/CRS/`.
2. No internet connection is needed: stage 01 reads the World Bank series from the shipped caches in `data/raw/wdi_cache/` and contacts the World Bank API only if one of those files is missing.
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

`ri_draws.csv` and `ri_draws_within.csv` store each double both as readable decimal text and as an exact hexadecimal column (`*_hex`), which stage 08 reads, so the cached draws reload bit for bit. The within-country share uses the same permuted assignments as the other two outcomes (each design re-seeds with 1242 before drawing).

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
| 2 | `code/02_descriptive_stats.R` | `data/processed/adaptationNAP.csv`, `adaptationNAP_donortype_wgi.csv`, `simple_panel_wgi.csv`, `donor_totals.csv`, `data/raw/PVCCI.csv` | descriptive tables and figures at the top level of `output/tables/`, `output/figures/` |
| 3 | `code/03_main_results.R` | `simple_panel_wgi.csv` | `output/tables/cohorts_dropped/att_combined_wide.tex`, `extensive_margin/`, `nap_cohorts.tex`, main figures; headline fits in `output/fits/` |
| 4 | `code/07_principal_and_share.R` | `simple_panel_wgi.csv`, `output/fits/` | `principal_share/`; within-country share fits in `output/fits/` |
| 5 | `code/04_robustness.R` | `simple_panel_wgi.csv`, `mitigation_panel.csv`, `output/fits/` | `cohorts_retained/`, `balanced_panel/`, `panel_2010/`, `notyettreated/`, `outlier_india/`, `units_zeros/`, `placebo/`, `mitigation/`, `bacon/`, `dcdh/`, HonestDiD tables in `cohorts_dropped/`; placebo and mitigation fits in `output/fits/` |
| 6 | `code/05_heterogeneity.R` | `simple_panel_wgi.csv`, `data/raw/oghist/OGHIST.xlsx`, `output/fits/` | `heterogeneity/` |
| 7 | `code/13_cohort_anticipation.R` | `simple_panel_wgi.csv`, `output/fits/` | `cohort_battery/`, `anticipation/`, `timing/` |
| 8 | `code/14_hazard_napa.R` | `simple_panel_wgi.csv`, `emergency_response_panel.csv`, `data/raw/napa/`, `data/raw/emdat/emdat.csv` (if present), `output/fits/` | `hazard/`, `napa/`; `data/processed/napa_list.csv` |
| 9 | `code/08_randomization_inference.R` | `simple_panel_wgi.csv`, `output/fits/`, `output/tables/randomization/ri_draws.csv` and `ri_draws_within.csv` (caches) | `randomization/` |
| 10 | `code/09_remarking_decomposition.R` | `simple_panel_wgi.csv`, `data/processed/crs_adaptation_activities/` (or raw CRS if absent), `output/fits/` | `remarking/`; run log in `output/logs/` |
| 11 | `code/10_model_tests.R` | `simple_panel_wgi.csv`, `output/fits/` | `model_tests/` |
| 12 | `code/11_base_year_sensitivity.R` | `simple_panel_wgi.csv`, `output/fits/` | `base_year/` |
| 13 | `code/12_group_figures.R` | `simple_panel_wgi.csv` | `output/figures/group/` |

`code/diagnostics/crs_commitment_date_check.py` (Python 3, standard library; from-raw only, about 15 s) is not part of `run_all.R`. It scans the raw 2021–2024 CRS files and writes `output/logs/diag_crs_commitment_dates.log`, the source of the commitment-date figures in the timing appendix (Appendix C, "Timing: adopters-only comparisons and the 2020–21 window").

`code/functions/pretrend_test.R` holds the joint Wald pre-trend test (`compute_pretrend_test()`). Stages 03, 04, 05, 07, 09, 11 and 13 source it, so every pre-trend statistic in the paper comes from one function. It inverts the pre-treatment covariance block with `solve()` at full rank, and with an SVD pseudo-inverse otherwise, in which case the test's degrees of freedom are the block's numerical rank. `code/functions/read_headline_fit.R` reads the headline fits that stage 03 stores in `output/fits/`; stages 04, 05, 07, 08, 09, 10, 11, 13 and 14 report the headline specification (ATT, SE, pre-trend test) from them instead of re-estimating it. Stage 04 stores the placebo and mitigation fits (read by 05 and 13) and stage 07 the within-country share fits (read by 04, 08, 11 and 13), which is why 07 runs right after 03. `code/functions/simple_att_event_weights.R` gives the event-time weights that reproduce a simple ATT, the target of the HonestDiD analyses (stages 03 and 07). `code/functions/mde.R` holds the one minimum-detectable-effect formula, $(z_{0.975} + z_{0.80}) \times$ SE $= 2.8016 \times$ SE, used by stages 03, 05, 13 and 14. `code/functions/make_wide_table.R` builds the wide tables of simple ATTs across outcomes (Table 2 in stage 03, the robustness variants in stage 04, the donor-type table in stage 05). `code/functions/make_country_id.R` builds the unit identifiers (see "Locale" above) and `code/functions/crs_positive.R` the zero/positive split of summed CRS amounts (stages 01, 02, 03, 04, 05, 07, 08, 11 and 13), with a tolerance for the residue that offsetting negative CRS records can leave.

**Failed runs.** Stages 04, 05, 07 to 11, 13 and 14 delete the exhibits they own before re-estimating them, and stages 03, 04 and 07 delete the fits they store in `output/fits/`; stages 01, 02, 03 and 12 overwrite their exhibits in place. Every stage stops (so `run_all.R` stops) if an estimation step fails, including the analytical fit and aggregation behind a pre-trend test, so a failed run cannot report success with an old exhibit or fit, or with an empty cell in place of an estimate. The exceptions, each disclosed in the exhibit concerned:

- the two EM-DAT tables are kept when EM-DAT is absent (see "EM-DAT" above);
- stage 13 runs the doubly-robust cell of the cohort-rule-by-estimator grid on the retained-cohort panel in a forked process, because `att_gt(est_method = "dr")` can crash on 2-unit cohorts; a crash is reported as unavailable ("---") in that table and in the cohort-battery summary, and both notes say so;
- stage 13, event-date row (a2) (all re-dated cohorts retained): if the doubly-robust fit cannot be produced in its forked process, the row uses outcome regression, and its label and the notes of the event-date and cohort-battery tables say so;
- stage 08: a permutation draw whose doubly-robust fit fails is re-estimated with outcome regression, and a draw that still fails is dropped; the table reports the number of valid draws and of fallbacks for each design (none with the shipped data);
- stage 14, the two NAP-timing tables: the complementary log-log hazard model is dropped from a table, and the note says so, when it cannot be estimated (with the shipped data it cannot, in both tables: every observation is a fixed-effects singleton or perfectly explained by the fixed effects);
- `did`'s own pre-test ("`did` pre-test $p$") is shown as "---" when `did` does not compute it, with the reason in the table note.

With the shipped data, the only exception that applies is the cloglog column of the NAP-timing tables.

After the last stage, `run_all.R` wraps each table's `tabular` in `\adjustbox{max width=\textwidth}` and sets floats to `[H]` (idempotent). Folder names under `output/tables/` and `output/figures/` are the same.

**Table format.** Stages 02 to 14 write each table as a complete LaTeX float (`table` environment with caption, label, notes and source), which the paper includes with `\input` as it stands. The three scope tables written by stage 01 (`output/tables/scope/`) are the exception: they contain only the `tabular` (wrapped in `\adjustbox` by `run_all.R`, like every table), and the manuscript supplies their float, caption and notes.

**Logs.** `run_all.R` copies each stage's console output, messages included, to `output/logs/<stage>.log` (e.g. `output/logs/03_main_results.log`), overwriting it on every run. Numbers that the paper quotes but that appear in no table or figure are printed there by the stage that computes them: for example the recipients that leave the DAC List or enter late and the NAP-adopter reconciliation (stage 01), the number of reporting donors and the recipients of the 2009 summary-statistics row (stage 02), the cohort-level ATTs (stages 03 and 04), the exact zeros in 2009 and the 2009 reporters (stage 04), the income-classification changes (stage 05), the listwise-deletion losses of each specification (stages 05 and 13), the propensity-score range (stage 13), the prior-NAPA by LDC-status counts of the NAP adopters (stage 14), the recipients without an EM-DAT entity, the NAPA submission years, the NAPA check's cohorts, control pool, estimate to eight significant digits (ATT 0.49002138, SE 0.30103212) and minimum detectable effect (stage 14), the extensive-margin minimum detectable effect (stage 03), the event-study coefficients of the governance split (stage 05) and of the within-country share (stage 07), and every cell of the within-share panels of the base-year table (stage 11). Stage 09 also writes its own detailed log, `output/logs/09_remarking_decomposition_log.txt` (rewritten on every run), which includes the activity-linkage rates. `output/logs/` is not shipped.

Processed data files (`data/processed/`): `simple_panel_wgi.csv` (estimation panel, one row per recipient-year), `adaptationNAP.csv` (descriptive panel), `adaptationNAP_donortype_wgi.csv` (recipient-year-donor type panel), `mitigation_panel.csv`, `emergency_response_panel.csv` (CRS emergency-response aid, used only for a NAP-timing check), `donor_totals.csv`, `donor_recipient_year_adaptation.csv` (used only for the donor count printed by stage 02), `napa_list.csv` (copy of the NAPA list written by stage 14), and `crs_adaptation_activities/` (activity-level extract for stage 09).

## List of tables and figures

Numbers refer to the paper as compiled on 13 October 2026 (Appendix A: model; Appendix B: data; Appendix C: supplementary estimation results). All files are in `output/`.

| Exhibit | Location | File in `output/` | Produced by |
|---|---|---|---|
| Figure 1 | Main text | `figures/group/adaptation_marginal.pdf` | `code/12_group_figures.R` |
| Table 2 | Main text | `tables/cohorts_dropped/att_combined_wide.tex` | `code/03_main_results.R` |
| Table 3 | Main text | `tables/principal_share/within_share_wide.tex` | `code/07_principal_and_share.R` |
| Figure 2 | Main text | `figures/cohorts_dropped/did_combined_es_wgi.png` | `code/03_main_results.R` |
| Table 4 | Main text | `tables/placebo/att_placebo.tex` | `code/04_robustness.R` |
| Table 5 | Main text | `tables/mitigation/att_mitigation.tex` | `code/04_robustness.R` |
| Table 6 | Main text | `tables/randomization/ri_pvalues.tex` | `code/08_randomization_inference.R` |
| Table 7 | Main text | `tables/cohorts_dropped/honestdid_rm.tex` | `code/04_robustness.R` |
| Table 8 | Main text | `tables/base_year/base_year_sensitivity.tex` | `code/11_base_year_sensitivity.R` |
| Table 9 | Main text | `tables/cohorts_retained/att_combined_wide.tex` | `code/04_robustness.R` |
| Table 10 | Main text | `tables/principal_share/principal_wide.tex` | `code/07_principal_and_share.R` |
| Table 11 | Main text | `tables/remarking/att_remarking_margins.tex` | `code/09_remarking_decomposition.R` |
| Table 12 | Main text | `tables/heterogeneity/donor_type/att_combined_wide.tex` | `code/05_heterogeneity.R` |
| Table 13 | Main text | `tables/heterogeneity/het_capacity.tex` | `code/05_heterogeneity.R` |
| Table 14 | Main text | `tables/heterogeneity/het_difference_tests.tex` | `code/05_heterogeneity.R` |
| Table A.1 | Appendix | `tables/model_tests/lemma2_size.tex` | `code/10_model_tests.R` |
| Table A.2 | Appendix | `tables/model_tests/alpha_half.tex` | `code/10_model_tests.R` |
| Table B.1 | Appendix | `tables/scope/flow_type_shares.tex` | `code/01_prepare_data.R` |
| Table B.2 | Appendix | `tables/scope/sample_funnel.tex` | `code/01_prepare_data.R` |
| Table B.3 | Appendix | `tables/scope/regional_exclusion.tex` | `code/01_prepare_data.R` |
| Table B.4 | Appendix | `tables/stats_des.tex` | `code/02_descriptive_stats.R` |
| Figure B.1 | Appendix | `figures/nap_status_map.png` | `code/02_descriptive_stats.R` |
| Table B.5 | Appendix | `tables/nap_cohorts.tex` | `code/03_main_results.R` |
| Figure B.2 | Appendix | `figures/cumulative_nap_adoption.png` | `code/02_descriptive_stats.R` |
| Table B.6 | Appendix | `tables/balance_adopters.tex` | `code/02_descriptive_stats.R` |
| Figure B.3 | Appendix | `figures/group/adopter_vs_never.pdf` | `code/12_group_figures.R` |
| Table B.7 | Appendix | `tables/principal_share/share_reconciliation.tex` | `code/07_principal_and_share.R` |
| Figure B.4 | Appendix | `figures/top_donors.png` | `code/02_descriptive_stats.R` |
| Figure B.5 | Appendix | `figures/recipient_map.png` | `code/02_descriptive_stats.R` |
| Table C.1 | Appendix | `tables/balanced_panel/att_combined_wide.tex` | `code/04_robustness.R` |
| Table C.2 | Appendix | `tables/panel_2010/att_combined_wide.tex` | `code/04_robustness.R` |
| Table C.3 | Appendix | `tables/outlier_india/att_combined_wide.tex` | `code/04_robustness.R` |
| Figure C.1 | Appendix | `figures/principal_share/fig_within_share_es.png` | `code/07_principal_and_share.R` |
| Figure C.2 | Appendix | `figures/principal_share/fig_share_reconciliation.png` | `code/07_principal_and_share.R` |
| Figure C.3 | Appendix | `figures/cohorts_dropped/did_combined_cohort_wgi.png` | `code/03_main_results.R` |
| Table C.4 | Appendix | `tables/extensive_margin/att_extensive.tex` | `code/03_main_results.R` |
| Table C.5 | Appendix | `tables/anticipation/att_anticipation.tex` | `code/13_cohort_anticipation.R` |
| Table C.6 | Appendix | `tables/cohort_battery/cohort_contributions.tex` | `code/13_cohort_anticipation.R` |
| Table C.7 | Appendix | `tables/anticipation/att_eventdate.tex` | `code/13_cohort_anticipation.R` |
| Table C.8 | Appendix | `tables/anticipation/cohort_redating_crosstab.tex` | `code/13_cohort_anticipation.R` |
| Table C.9 | Appendix | `tables/cohort_battery/att_cohort_battery.tex` | `code/13_cohort_anticipation.R` |
| Table C.10 | Appendix | `tables/hazard/nap_timing_vs_emdat_hazard.tex` | `code/14_hazard_napa.R` (needs EM-DAT) |
| Table C.11 | Appendix | `tables/hazard/nap_timing_vs_humanitarian_aid.tex` | `code/14_hazard_napa.R` |
| Table C.12 | Appendix | `tables/hazard/att_hazard_controls.tex` | `code/14_hazard_napa.R` (needs EM-DAT) |
| Table C.13 | Appendix | `tables/napa/att_napa_falsification.tex` | `code/14_hazard_napa.R` |
| Table C.14 | Appendix | `tables/napa/att_prior_napa_split.tex` | `code/14_hazard_napa.R` |
| Figure C.4 | Appendix | `figures/placebo/did_placebo_es.png` | `code/04_robustness.R` |
| Figure C.5 | Appendix | `figures/mitigation/did_mitigation_es.png` | `code/04_robustness.R` |
| Figure C.6 | Appendix | `figures/randomization/fig_ri_distributions.png` | `code/08_randomization_inference.R` |
| Table C.15 | Appendix | `tables/cohorts_dropped/honestdid_prepriods.tex` | `code/04_robustness.R` |
| Table C.16 | Appendix | `tables/cohorts_dropped/honestdid_sd.tex` | `code/04_robustness.R` |
| Table C.17 | Appendix | `tables/randomization/ri_pvalues_within.tex` | `code/08_randomization_inference.R` |
| Figure C.7 | Appendix | `figures/randomization/fig_ri_distributions_within.png` | `code/08_randomization_inference.R` |
| Table C.18 | Appendix | `tables/cohorts_dropped/honestdid_within.tex` | `code/04_robustness.R` |
| Table C.19 | Appendix | `tables/base_year/pretrend_tests_full.tex` | `code/11_base_year_sensitivity.R` |
| Table C.20 | Appendix | `tables/base_year/pretrend_cells_2021.tex` | `code/11_base_year_sensitivity.R` |
| Figure C.8 | Appendix | `figures/base_year/fig_pretrend_cells_by_cohort.png` | `code/11_base_year_sensitivity.R` |
| Figure C.9 | Appendix | `figures/base_year/fig_es_full_window.png` | `code/11_base_year_sensitivity.R` |
| Table C.21 | Appendix | `tables/base_year/base_year_sensitivity_logit.tex` | `code/11_base_year_sensitivity.R` |
| Table C.22 | Appendix | `tables/base_year/att_premean_baseline.tex` | `code/11_base_year_sensitivity.R` |
| Table C.23 | Appendix | `tables/timing/att_adopters_only.tex` | `code/13_cohort_anticipation.R` |
| Table C.24 | Appendix | `tables/timing/cohort_submission_months.tex` | `code/13_cohort_anticipation.R` |
| Table C.25 | Appendix | `tables/notyettreated/att_notyettreated_wide.tex` | `code/04_robustness.R` |
| Figure C.10 | Appendix | `figures/notyettreated/did_notyettreated_es.png` | `code/04_robustness.R` |
| Figure C.11 | Appendix | `figures/cohorts_retained/did_combined_cohort_wgi.png` | `code/04_robustness.R` |
| Table C.26 | Appendix | `tables/cohorts_retained/att_group_retained.tex` | `code/04_robustness.R` |
| Table C.27 | Appendix | `tables/dcdh/att_dcdh.tex` | `code/04_robustness.R` |
| Figure C.12 | Appendix | `figures/dcdh/did_dcdh_es.png` | `code/04_robustness.R` |
| Table C.28 | Appendix | `tables/bacon/bacon_decomp.tex` | `code/04_robustness.R` |
| Table C.29 | Appendix | `tables/units_zeros/diagnostic_units_zeros.tex` | `code/04_robustness.R` |
| Table C.30 | Appendix | `tables/cohort_battery/att_2x2.tex` | `code/13_cohort_anticipation.R` |
| Table C.31 | Appendix | `tables/anticipation/att_placebo_ladder.tex` | `code/13_cohort_anticipation.R` |
| Table C.32 | Appendix | `tables/principal_share/principal_retained.tex` | `code/07_principal_and_share.R` |
| Table C.33 | Appendix | `tables/remarking/att_remarking_counts.tex` | `code/09_remarking_decomposition.R` |
| Table C.34 | Appendix | `tables/remarking/att_remarking_exclusions.tex` | `code/09_remarking_decomposition.R` |
| Table C.35 | Appendix | `tables/remarking/tab_remarking_sectors.tex` | `code/09_remarking_decomposition.R` |
| Table C.36 | Appendix | `tables/remarking/remarking_flag_shares_by_group.tex` | `code/09_remarking_decomposition.R` |
| Table C.37 | Appendix | `tables/remarking/att_remarking_exclusions_regex_comparison.tex` | `code/09_remarking_decomposition.R` |
| Table C.38 | Appendix | `tables/heterogeneity/donor_type/zero_shares.tex` | `code/05_heterogeneity.R` |
| Table C.39 | Appendix | `tables/heterogeneity/governance/het_gov_wide.tex` | `code/05_heterogeneity.R` |
| Table C.40 | Appendix | `tables/heterogeneity/ldc/het_ldc_wide.tex` | `code/05_heterogeneity.R` |
| Table C.41 | Appendix | `tables/heterogeneity/income_group/het_income_wide.tex` | `code/05_heterogeneity.R` |
| Figure C.13 | Appendix | `figures/heterogeneity/fig_het_es_panel.png` | `code/05_heterogeneity.R` |

Table 1 (mapping of the IPCC AR6 risk components onto aid-allocation roles) is typed directly in the manuscript and is not produced by code. Every other table and figure is listed above.

## Licence

Code (`run_all.R`, `code/`): MIT License. Data compiled by the authors, the processed panels and the exhibits: Creative Commons Attribution 4.0 International (CC BY 4.0). Third-party data keep their own terms (table above). In particular, the PVCCI values in `data/raw/PVCCI.csv` and in the processed panels (`pvcci` and the columns derived from it) remain under FERDI's terms and are excluded from the CC BY 4.0 grant. EM-DAT data are not included and may not be redistributed. See `LICENSE`.

## How to cite

Beaucoral, P., M. Goujon and S. Marchand (2026). "Do National Adaptation Plans Attract Adaptation Finance? Evidence from the Composition of Aid." Working paper, CERDI, Université Clermont Auvergne.

Replication package: Beaucoral, P., M. Goujon and S. Marchand (2026). Replication package for "Do National Adaptation Plans Attract Adaptation Finance? Evidence from the Composition of Aid." **[REPOSITORY DOI: TO BE ADDED ON DEPOSIT]**

## Acknowledgements

This work was supported by the Agence Nationale de la Recherche of the French government through the program "Investissements d'avenir" (ANR-10-LABX-14-01).
