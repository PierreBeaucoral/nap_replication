# EM-DAT (CRED / UCLouvain) — not included, obtain it yourself

EM-DAT's terms of use (<https://doc.emdat.be/docs/legal/terms-of-use/>) do not
allow redistributing the database or any substantial part of it, nor
distributing derivative databases. The replication package therefore contains
**no EM-DAT data**: neither the extract used in the paper nor the
recipient-year panel built from it (`data/processed/emdat_panel.csv`).

## What needs EM-DAT

Only two appendix exhibits, both produced by `code/14_hazard_napa.R`:

- `output/tables/hazard/att_hazard_controls.tex` (hazard-control specifications)
- `output/tables/hazard/nap_timing_vs_emdat_hazard.tex` (NAP timing vs. lagged EM-DAT hazard)

Without `emdat.csv` the pipeline still runs end to end. Stage 14 prints a
notice naming these two tables and leaves the authors' copies in
`output/tables/hazard/` untouched (they are then not regenerated).

## How to obtain the extract used in the paper

1. Register (free, non-commercial academic use) at <https://public.emdat.be>.
2. Open **Data** > custom request and select:
   - Disaster classification: **Natural** (Disaster Group = Natural; all subgroups),
   - Location: **all countries / all regions**,
   - Time period: **Start Year 2000 to 2026**,
   - include historical events: default setting.
3. Download the Excel file. The paper uses the file created on
   15 September 2026 (EM-DAT **version 2026-09-11**, table type
   `public_emdat_custom_request`, **10,896 records**; the "EM-DAT Info" sheet
   of the download shows version and record count), saved as
   `public_emdat_custom_request_2026-09-15.xlsx`.
4. Convert the "EM-DAT Data" sheet to `data/raw/emdat/emdat.csv`, keeping these
   20 columns in this order (R, from the project root):

   ```r
   library(readxl)
   cols <- c("DisNo.", "ISO", "Country", "Region", "Disaster Group",
             "Disaster Subgroup", "Disaster Type", "Disaster Subtype",
             "Start Year", "Start Month", "End Year", "Total Deaths",
             "No. Injured", "No. Affected", "No. Homeless", "Total Affected",
             "Total Damage ('000 US$)", "Total Damage, Adjusted ('000 US$)",
             "Magnitude", "Magnitude Scale")
   x <- read_excel("data/raw/emdat/public_emdat_custom_request_2026-09-15.xlsx",
                   sheet = "EM-DAT Data")
   write.csv(x[, cols], "data/raw/emdat/emdat.csv", row.names = FALSE, na = "")
   ```

5. Rerun stage 14: `NAP_START_AT=14_hazard_napa.R Rscript run_all.R`
   (or `Rscript code/14_hazard_napa.R`).

EM-DAT is updated continuously and past events are revised, so a download made
later than the paper's will not be identical, and the two tables can differ
slightly. For reference, the SHA-256 checksums of the authors' files are:

| File | SHA-256 |
|---|---|
| `public_emdat_custom_request_2026-09-15.xlsx` | `6f0e3bd7bad683a4622f107dd7bf6c2aab80d5b2da089a47f495ab5ffb3357a6` |
| `emdat.csv` | `e1091041540e082ff87ffa44a5e209d230ad2d4c8a91f88b0c302ced75e11b17` |

A CSV written with the snippet above contains the same data but is not
byte-identical to the authors' `emdat.csv` (number formatting and quoting
differ, e.g. `4` vs `4.0`), so compare row counts and values, not checksums.

## How stage 14 uses it

`code/14_hazard_napa.R` §2 aggregates the event-level extract to a
recipient-year panel (2008-2024) for two hazard sets: all natural events, and
climate-related events (Hydrological, Meteorological and Climatological
subgroups). Events are assigned to their Start Year.

- `n_events` is a genuine 0 when EM-DAT records no event for that
  country-year.
- `Total Affected` and `Total Damage, Adjusted` are 0 only when there was no
  event, and NA (not 0) when events occurred but none reported the value.
- Kosovo (XKX) and Nauru (NRU) have no EM-DAT entity at all (Kosovo's events
  are booked under Serbia), so their whole hazard series is NA rather than 0.

## Citation

EM-DAT, CRED / UCLouvain, Brussels, Belgium — www.emdat.be (version
2026-09-11, accessed 15 September 2026). Delforge, D., Wathelet, V., Below, R.,
Lanfredi Sofia, C., Tonnelier, M., van Loenhout, J. A. F., and Speybroeck, N.
(2025). EM-DAT: the Emergency Events Database. *International Journal of
Disaster Risk Reduction*, 124, 105509.
