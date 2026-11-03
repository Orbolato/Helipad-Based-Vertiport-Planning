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

The author-provided computational package is organized at the **repository root**. This avoids duplicating the source already published under `code/` and `data/`. The research files comprise eight Python source modules, derived spatial input tables, saved run records and allocation certificates, spatial sensitivity results, reproducible figures and tables, and a SHA-256 manifest.

| Location | Contents |
| --- | --- |
| `code/` | Capacitated location model, SMA-inspired/VNS/relay algorithms, sensitivity experiments, validation, and figure/table generation |
| `data/` | Derived ANAC/IBGE candidate and grid inputs, parameters, original site selections, and statewide traces |
| `results/` | Saved seed-level intensity, spatial-demand, grid-resolution, boundary, and envelope experiments (JSON, NPZ, CSV) |
| `figures/`, `tables/` | Generated graphics and manuscript table fragments |
| `Supplementary_Information.tex`, `.pdf` | Detailed author-provided results, not the publisher's formally hosted supplement |
| `requirements.txt`, `software_versions.json` | Pinned computational environment |
| `validation_report.json`, `MANIFEST_SHA256.txt` | Supplied validation report and file-integrity checks |

If some of these directories are absent on GitHub, the full computational archive has not yet been uploaded. Do not interpret a partially populated repository as an independently reproducible package.

### Quick start

Use Python 3.12. From the repository root:

```bash
python -m pip install -r requirements.txt
sha256sum -c MANIFEST_SHA256.txt   # On systems with sha256sum; before regeneration
python code/generate_tables.py
python code/generate_figures.py
python code/validate_results.py
```

The supplied validation report records **147** checked continuous-allocation results and **44 of 45** feasible demand-intensity runs. This is a saved report, not a claim of independent cross-platform rerunning.

To check the saved intensity selections, run `python code/run_intensity.py`. Fresh stochastic runs require `--rerun` and can replace checkpoint files, so perform them in a separate copy. Controlled spatial experiments are available through `python code/run_spatial.py`. Some exact reoptimization procedures require **Gurobi** and an active solver license.

## Data and software provenance

The candidate coordinates/types and cell population/geometry are derived from official **ANAC** and **IBGE** sources. Original providers retain their rights. Raw registry workbooks are not needed to rerun the supplied analysis from the derived inputs, but original source-data extraction is not reconstructed here. No blanket open-source license is asserted for all materials. See `PROVENANCE.md`.

## Citation

Please cite the journal article when using its results. Publication DOI, volume, and article number will be added when assigned. See [CITATION.cff](CITATION.cff) for provisional bibliographic metadata.

**Related work:** [Author publications](https://orbolato.github.io/publications/) · [UAM network planning in São Paulo](https://github.com/Orbolato/UAM-Network-Planning-Sao-Paulo)
