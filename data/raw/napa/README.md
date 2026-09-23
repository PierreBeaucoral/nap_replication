# NAPA submissions (UNFCCC)

`napa_list_unfccc.csv`: the 51 National Adaptation Programmes of Action (NAPAs) posted on the UNFCCC "NAPAs received" page (<https://unfccc.int/topics/resilience/workstreams/national-adaptation-programmes-of-action/napas-received>; latest: South Sudan, February 2017). The authors transcribed the list by hand on 2026-09-15 because the page blocks scripted access. The page as retrieved is archived here as `Submitted_NAPAs_UNFCCC.pdf`, and the CSV was cross-checked against it (51/51 country–month–year matches). Bangladesh: first submission November 2005 (updated June 2009).

Columns: `iso3`, `country`, `napa_year`, `napa_month`, `source`.

`code/14_hazard_napa.R` §5 reads this file and stops with an error if it is missing; there is no fallback list.

SHA-256 (`napa_list_unfccc.csv`): `bcf34444fc37e00cc2be8ea47cc325488d01bad0a629526ccac817a45359952d`
