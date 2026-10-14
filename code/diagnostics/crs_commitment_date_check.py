"""Check whether raw CRS CommitmentDate supports a within-year timing test.

For adaptation-marked records (ClimateAdaptation 1 or 2) in the 2021-2024 raw
CRS files, report the share with a filled CommitmentDate, the share of filled
dates on 31 December or 1 January (year-boundary placeholders), and the share
falling in the reporting year. Cited in the timing appendix (app:timing).

Usage (from the project root): python3 code/diagnostics/crs_commitment_date_check.py
Writes output/logs/diag_crs_commitment_dates.log. Reads the raw CRS files only.
"""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
YEARS = range(2021, 2025)
LOG = ROOT / "output" / "logs" / "diag_crs_commitment_dates.log"

csv.field_size_limit(sys.maxsize)


def crs_file(year):
    """Return the raw CRS file for a year (the case of 'Data' varies)."""
    for name in (f"CRS {year} data.txt", f"CRS {year} Data.txt"):
        p = ROOT / "data" / "raw" / "CRS" / name
        if p.exists():
            return p
    raise FileNotFoundError(f"raw CRS file for {year} not found")


def scan(year):
    """Count adaptation-marked records and their CommitmentDate patterns."""
    n = filled = boundary = same_year = 0
    with open(crs_file(year), encoding="utf-8", errors="replace", newline="") as fh:
        rd = csv.reader(fh, delimiter="|")
        head = next(rd)
        i_ad, i_cd, i_yr = head.index("ClimateAdaptation"), head.index("CommitmentDate"), head.index("Year")
        for row in rd:
            if len(row) <= max(i_ad, i_cd) or row[i_ad].strip() not in ("1", "2"):
                continue
            n += 1
            d = row[i_cd].strip()
            if not d:
                continue
            filled += 1
            md = d[5:10]  # ISO date YYYY-MM-DD...
            boundary += md in ("12-31", "01-01")
            same_year += d[:4] == row[i_yr].strip()
    if filled == 0:
        raise ValueError(f"{year}: no filled CommitmentDate; check the date format")
    return n, filled, boundary, same_year


def main():
    """Scan each year and write the log."""
    lines = ["year | adaptation-marked records | filled (%) | 31 Dec or 1 Jan, of filled (%) | in reporting year, of filled (%)"]
    for y in YEARS:
        n, f, b, s = scan(y)
        lines.append(f"{y} | {n} | {100 * f / n:.1f} | {100 * b / f:.1f} | {100 * s / f:.1f}")
    LOG.parent.mkdir(parents=True, exist_ok=True)
    LOG.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
