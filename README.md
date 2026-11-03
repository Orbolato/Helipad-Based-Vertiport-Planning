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

The author-provided computational package includes the Python optimization implementation, derived scenario inputs, original solution records, 45 demand-intensity runs, spatial sensitivity and feasibility checks, tables, figures, and supplementary information. It was originally prepared as a **local research/verification record**, not as an officially published Elsevier supplementary file.

**Publication status of files:** The website and citation metadata are available here. The 39-MB compressed research package requires a separate Git upload for its binary inputs; until that upload is completed, not all directories described below are available in this repository.

When present, the package contains:

- `code/`: optimization, stochastic search, experiment drivers, result validation, and figure/table generation.
- `data/`: derived ANAC/IBGE candidate/grid inputs and saved original configurations.
- `results/`: seed-level selection, allocation, and sensitivity records.
- `figures/` and `tables/`: generated visualizations and manuscript table fragments.
- `Supplementary_Information.pdf`: descriptive computational results, not the journal's published article.
- `MANIFEST_SHA256.txt`: integrity hashes for research files.

With Python 3.12 and the dependencies in `requirements.txt`:

```bash
python -m pip install -r requirements.txt
python code/generate_tables.py
python code/generate_figures.py
python code/validate_results.py
```

Fresh stochastic experiments: `python code/run_intensity.py --rerun`. See the computational package's original instructions and `PROVENANCE.md` for coverage, data sourcing, environment, and Gurobi licensing requirements.

## Data and software provenance

The candidate coordinates/types and cell population/geometry are derived from official **ANAC** and **IBGE** sources. Original providers retain their rights. Raw registry workbooks are not needed to rerun the supplied analysis from the derived inputs, but original source-data extraction is not reconstructed here. No blanket open-source license is asserted for all materials. See `PROVENANCE.md`.

## Citation

Please cite the journal article when using its results. Publication DOI, volume, and article number will be added when assigned. See [CITATION.cff](CITATION.cff) for provisional bibliographic metadata.

**Related work:** [Author publications](https://orbolato.github.io/publications/) · [UAM network planning in São Paulo](https://github.com/Orbolato/UAM-Network-Planning-Sao-Paulo)
