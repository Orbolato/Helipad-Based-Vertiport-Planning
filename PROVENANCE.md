# Computational artifact provenance and limitations

## Article and scope

This repository accompanies the **accepted** journal article *Planning vertiport infrastructure from existing helipad sites: A multi-scale case study in São Paulo State* (Furlan, Carvalho, Guterres, and Correia; *Case Studies on Transport Policy*, 2026). It is not a claim that the publisher released these files as formal supplemental material.

The original computational supplement (`revision_supplement.zip`) contains three derived scenarios: `light` = São Paulo municipality, `medium` = pooled São Paulo/Campinas/Guarulhos, and `heavy` = São Paulo State. The analyses use an ANAC-derived inventory of 609 registered helipad candidates statewide and an IBGE-derived 1-km grid of 250,827 statewide cells.

## Included evidence and reproduction levels

1. **Artifact integrity.** The research inputs, selections, and results have a SHA-256 manifest. Two stale entries in the supplied archive were identified for `code/generate_tables.py` and `tables/intensity_results.tex`; the public upload bundle corrects those manifest values. The research file bytes are otherwise preserved.
2. **Saved-analysis validation.** The provided `validation_report.json` records 147 continuous-allocation validations and status `passed`. This is a *supplied report*, not a claim that this repository's maintainers have rerun every calculation on every platform.
3. **Figures and tables.** The scripts regenerate figures and table fragments from included derived inputs and saved results.
4. **Fresh stochastic experiments.** `code/run_intensity.py --rerun` runs the sensitivity experiments and overwrites the corresponding checkpoint files in that working copy. Stochastic details and runtime may vary with hardware and dependencies.
5. **Exact reoptimization.** Some procedures require `gurobipy` and an active Gurobi license, including exact mixed-integer reoptimization of the refined scenarios.

## Data boundaries

- **ANAC:** derived coordinates and candidate types from registered Brazilian civil aerodromes. Registration does not establish site suitability or current operational status.
- **IBGE:** derived grid populations, geometry and boundary data from Brazilian census/statistical products. The original source provider retains its rights.
- **No individually identifiable passenger trajectory records** are needed for the supplied scenario analyses.
- The package starts from prepared derived scenarios; it does not reproduce every step from original registry workbooks or source downloads.
- Municipal and regional boundary grid-support cells are not identical to census municipal totals due to whole-cell boundary inclusion.
- Synthetic capacities and population-weighted planning demand must not be presented as observed service demand, aircraft throughput, or adoption forecasts.

## Integrity and publication notes

The research files were originally marked *local-only* and were not submitted as publisher-hosted supplementary files. Publishing them here makes them an author-provided open research artifact, but does not grant rights over third-party data or remove proprietary dependencies.

The final Elsevier proof PDF is **not** included; do not redistribute the production proof as an authorized published version. Link the published article or DOI after formal publication.

The ZIP that was provided for this repository includes a 55-MiB GeoPackage (`data/heavy_grid.gpkg`). It is below GitHub's 100-MiB single-file hard limit, but may trigger the 50-MiB file-size warning on direct Git pushes.

No blanket license is applied to code or data until the respective rightsholders and redistribution terms have been confirmed.
