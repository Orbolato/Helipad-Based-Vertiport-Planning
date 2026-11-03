# Provenance and scope of the lightweight research package

This repository accompanies the accepted 2026 *Case Studies on Transport Policy* article **Planning vertiport infrastructure from existing helipad sites: A multi-scale case study in São Paulo State**, by Furlan, Carvalho, Guterres, and Correia.

## Included

- Original Python research modules under `code/`. These are preserved for methodological inspection.
- Derived ANAC candidate-site and IBGE population/grid CSV tables for municipal (`light`), pooled three-city (`medium`), and statewide (`heavy`) scenarios under `data/`.
- Aggregated numerical results under `results/`, including `all_analysis_records.csv` (147 records), `intensity_seed_results.csv` (45 runs), spatial-demand redistribution tables, and exact municipal/regional allocation summaries.
- Library version specifications under `requirements.txt` and `software_versions.json`.

## Intentionally omitted

- Geospatial intermediates such as `*.gpkg`, including the large boundary/grid GeoPackages.
- NumPy binary archives `*.npz` containing original solutions, fixed preference/distance matrices, and sparse allocation certificates.
- Pickled convergence traces and original per-run JSON records.
- Generated graphics and LaTeX output, temporary caches, and locally generated working directories.

These omissions are intentional to keep the repository lightweight. **Do not describe the archive as sufficient to rerun all original simulations or regenerate every paper figure/table.** Some original scripts still refer to the omitted assets. The existing CSV summaries document results but are not equivalent to independent replication from raw data.

## Data use

ANAC registration is a candidate-inventory source, not evidence that a site can be converted to a vertiport. IBGE-derived cells and populations are planning inputs. Demand and capacity values are synthetic assumptions, not observed commercial demand or measured aircraft throughput.

Original ANAC/IBGE sources retain their own usage conditions. This repository does not redistribute the publisher's final proof PDF or assert that the publisher formally hosts this computational material. No blanket third-party data or software license is asserted.

## Reproduction

A full reproduction requires omitted intermediate geospatial and solution files, source-data preparation, the specified Python environment, and a licensed Gurobi optimizer for exact methods. Retaining the original scripts here documents the method, but the repository provides **results inspection and source transparency**, not fully automated end-to-end reproducibility.
