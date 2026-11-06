# Planning vertiport infrastructure from existing helipad sites

**Published research article (2026)** · *Case Studies on Transport Policy*, Article 102000 · https://doi.org/10.1016/j.cstp.2026.102000

Rodrigo Mollo Furlan, Lucas Orbolato Carvalho, Marcelo Xavier Guterres, and Anderson Ribeiro Correia

**Publication page:** https://orbolato.github.io/publications/helipad-based-vertiport-planning/

**Scholarly discovery:** This repository links to the [canonical publication page](https://orbolato.github.io/publications/helipad-based-vertiport-planning/), which is listed in the [author website sitemap](https://orbolato.github.io/sitemap.xml) and covered by the site's [robots.txt](https://orbolato.github.io/robots.txt). Search engines index the public repository separately under GitHub's own crawling rules. These links enable discovery but do not guarantee search-engine indexing.

## Study

The study examines whether registered helipads can serve as candidates for strategic vertiport planning across three São Paulo planning contexts. The capacitated facility-location formulation incorporates population-weighted planning demand, candidate capacities, geometric access, candidate-reachable coverage, and a soft spacing preference. Exact optimization benchmarks the two smaller instances, while variable neighborhood search (VNS) and an SMA-inspired binary search/VNS relay support statewide screening.

| Planning context | Candidates | Selected sites |
| --- | ---: | --- |
| São Paulo municipality | 185 | Exact optimum: 4 |
| Pooled São Paulo, Campinas, and Guarulhos | 222 | Exact optimum: 13 |
| São Paulo State | 609 | Feasible relay solutions: 166–169 |

The statewide scenario contains 250,827 one-kilometer grid cells. Feasible relay configurations reach approximately 80% of **candidate-reachable** territory, equivalent to about **14.34% of total state grid area**. Of 45 demand-intensity configurations, 44 met both coverage floors. These are conditional screening results, not implementation-ready site recommendations.

## Lightweight research package

This is a **compact, CSV-based research archive**, rather than the full original computational snapshot.

| Directory or file | Contents |
| --- | --- |
| `code/` | Original research Python modules for optimization, sensitivity analysis, evaluation, and report generation |
| `data/` | Derived candidate/cell CSV inputs and scenario parameters |
| `results/` | Tabular experiment and spatial-sensitivity outputs |
| `requirements.txt`, `software_versions.json` | Original environment specifications |
| `PROVENANCE.md` | Scope of the public package and reproduction limitations |
| `CITATION.cff` | Citation metadata |

**Upload status:** `code/` and `data/` are already on GitHub. The separate `results/`, `requirements.txt`, and `software_versions.json` files may not be present until the lightweight ZIP is uploaded.

### What can be inspected

The optional results CSV package documents 147 analysis records, 45 demand-intensity runs, two exact local configurations, and spatial-sensitivity summaries. As of this update, the `results/` directory has not yet been uploaded to the public repository. The original scripts are retained for transparency, but **several scripts require research assets omitted from the lightweight release** (GeoPackages, `.npz` selection/certificate files, and some saved JSON records). Thus, commands such as `run_intensity.py`, `run_spatial.py`, `generate_figures.py`, `generate_tables.py`, and `validate_results.py` are **not fully executable from this package alone**.

Once the additional results files have been uploaded, they can be inspected directly with pandas:

```python
import pandas as pd

summary = pd.read_csv("results/intensity_summary.csv")
records = pd.read_csv("results/all_analysis_records.csv")
print(summary)
print(f"{len(records)} documented analysis records")
```

To reproduce the original end-to-end experiments, additional archived computational inputs and, for exact optimization, a licensed Gurobi installation are necessary. See [PROVENANCE.md](PROVENANCE.md).

## Data provenance and licensing

Candidate records derive from Brazil's ANAC civil-aerodrome inventory; spatial population and boundaries derive from IBGE data. The underlying source providers retain their rights. The public package contains no publisher-produced proof PDF and no general software/data license is asserted.

## Citation

Cite the published journal article for its scientific conclusions: **Furlan, R. M., Carvalho, L. O., Guterres, M. X., and Correia, A. R. (2026).** *Planning vertiport infrastructure from existing helipad sites: A multi-scale case study in São Paulo State.* *Case Studies on Transport Policy*, Article 102000. https://doi.org/10.1016/j.cstp.2026.102000. The volume and complete publication date will be added once independently confirmed. See [CITATION.cff](CITATION.cff).

**Related:** [Author publications](https://orbolato.github.io/publications/) · [São Paulo UAM network planning](https://github.com/Orbolato/UAM-Network-Planning-Sao-Paulo)
