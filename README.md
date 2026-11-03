# Planning vertiport infrastructure from existing helipad sites

**Research article (accepted, 2026)** · *Case Studies on Transport Policy*

**Rodrigo Mollo Furlan, Lucas Orbolato Carvalho, Marcelo Xavier Guterres, and Anderson Ribeiro Correia**

Dedicated paper page: https://orbolato.github.io/publications/helipad-based-vertiport-planning/

## Overview

This study tests whether São Paulo's registered helipad inventory can support **strategic screening** of potential urban air mobility (UAM) vertiport locations at three geographic scales. A capacitated location model combines population-weighted planning demand, type-dependent site capacity, geometric access, candidate-reachable territorial coverage, and a soft site-spacing preference.

Exact optimization benchmarks municipal and pooled multi-city cases; stochastic search, including variable neighborhood search (VNS) and an SMA-inspired binary search/VNS relay, supports statewide screening. The selected sites are **planning candidates, not approved or implementation-ready vertiports**.

## Key results

| Planning scale | Candidate inventory | Result |
| --- | ---: | --- |
| São Paulo municipality | 185 | Exact optimum: 4 selected sites |
| São Paulo, Campinas, and Guarulhos (pooled) | 222 | Exact optimum: 13 sites |
| São Paulo State | 609 | Feasible relay configurations: 166–169 sites across five seeds |

The statewide case contains **250,827** 1-km grid cells. Feasible relay configurations covered approximately **80% of the territory reachable from the candidate inventory**, but only **14.34% of total state grid area**. In demand-intensity sensitivity tests, **44 of 45** configurations met both coverage floors. These are conditional, model-based screening results, not forecasts of actual UAM demand.

## Computational supplement

The **author-provided computational supplement** is organized under `research/` so its README and file manifest remain separate from this repository's publication metadata. The clean upload package contains **497 files**, with no compiled Python caches or local-only markers.

Once the research package has been imported, its directory includes:

- `research/code/`: capacitated location model, SMA-inspired search, VNS/relay procedures, sensitivity experiments, validation, and figure/table generation.
- `research/data/`: derived ANAC/IBGE grid and candidate inputs, parameters, original selections, and statewide traces.
- `research/results/`: 45 demand-intensity runs and saved composition, boundary, resolution, and geometric-envelope experiments, including paired JSON/NPZ allocation certificates.
- `research/figures/` and `research/tables/`: generated figures and manuscript table fragments.
- `research/Supplementary_Information.pdf` and `.tex`: detailed computational results supplied by the authors, not an official publisher-hosted supplement.
- `research/README.md`, `research/requirements.txt`, `research/software_versions.json`, `research/validation_report.json`, and `research/MANIFEST_SHA256.txt`: usage, environment, reported validation, and file-integrity information.

**Availability:** If the `research/` folder is absent on GitHub, the computational ZIP has not yet been imported; the documentation and publication page alone do not constitute an uploaded code/data release.

### Quick start

With Python 3.12, from the repository root:

```bash
python -m pip install -r research/requirements.txt
python research/code/generate_tables.py
python research/code/generate_figures.py
python research/code/validate_results.py
```

To confirm file integrity **before** running scripts (which can rewrite generated files), run `sha256sum -c MANIFEST_SHA256.txt` from the `research/` directory on platforms that provide `sha256sum`.

The supplied `research/validation_report.json` records 147 checked allocation records and 44 feasible outcomes among 45 demand-intensity runs. It is a saved report; it is not a guarantee of cross-platform bitwise reproduction. Some exact reoptimization commands require a **Gurobi license**. See `research/README.md` for the optional `run_intensity.py` and `run_spatial.py` commands, including their overwrite behavior.

## Data and software provenance

The candidate coordinates/types and cell population/geometry are derived from official **ANAC** and **IBGE** sources. Original providers retain their rights. Raw registry workbooks are not needed to rerun the supplied analysis from the derived inputs, but original source-data extraction is not reconstructed here. No blanket open-source license is asserted for all materials. See `PROVENANCE.md`.

## Citation

Please cite the journal article when using its results. Publication DOI, volume, and article number will be added when assigned. See [CITATION.cff](CITATION.cff) for provisional bibliographic metadata.

**Related work:** [Author publications](https://orbolato.github.io/publications/) · [UAM network planning in São Paulo](https://github.com/Orbolato/UAM-Network-Planning-Sao-Paulo)
