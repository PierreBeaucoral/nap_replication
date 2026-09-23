# Raw OECD CRS files (needed only for a full rebuild)

The default run (`Rscript run_all.R`) does **not** need these files: it starts
from the CRS-derived panels shipped in `data/processed/`. They are needed only
for a full rebuild from raw data (`NAP_FROM_RAW=1 Rscript run_all.R`).

## What to download

The per-year OECD Creditor Reporting System (CRS) bulk text files for
**2007 to 2024** (18 files, about 4.8 GB in total). Stage 01 reads all 18
years; stage 09 reads 2009-2024 when it has to rebuild its activity-level
extract. Earlier years (1973-2006) are not used.

Source: OECD Data Explorer, "Creditor Reporting System (flows)", full
bulk download of the CRS micro-data by year:
<https://data-explorer.oecd.org/vis?lc=en&df[ds]=DcdDisseminateFinalDMZ&df[id]=DSD_CRS%40DF_CRS&df[ag]=OECD.DCD.FSD&dq=DAC..1000.100._T._T.D.Q._T..&lom=LASTNPERIODS&lo=5&to[TIME_PERIOD]=false>

Place the files in this folder, named exactly:

    CRS 2007 data.txt
    CRS 2008 data.txt
    ...
    CRS 2024 data.txt

(`CRS <year> Data.txt`, with a capital D, is also accepted.)

## Vintage and checksums

The paper uses the files downloaded on 1 April 2026. `crs_checksums.csv`
lists the size and SHA-256 of each of the 18 files. In from-raw mode
`run_all.R` verifies every file against this table before rebuilding and
prints `MATCH` or `DIFFERENT VINTAGE` per year. The OECD revises past years of
the CRS, so a later download will usually differ: the rebuild then still
runs, but `run_all.R` warns and its panel-by-panel comparison with the shipped
`data/processed/` files shows how much the results move.

To check by hand: `shasum -a 256 "CRS 2007 data.txt"` (macOS/Linux) or
`certutil -hashfile "CRS 2007 data.txt" SHA256` (Windows).

## Terms of use and citation

OECD data are free to reuse with attribution (OECD Terms and Conditions,
<https://www.oecd.org/en/about/terms-conditions.html>). The raw files are not
included in the replication package because of their size; the CRS-derived
panels in `data/processed/` are.

Citation: OECD (2026), *Creditor Reporting System (CRS)*, OECD Data Explorer,
<https://data-explorer.oecd.org> (bulk files downloaded 1 April 2026).
