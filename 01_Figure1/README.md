# Figure 1 — Antibody waning

Code to reproduce Figure 1 of *Myeloma Precursors Erode Durable Immunity*.

| Panel | Script | Description |
|-------|--------|-------------|
| 1B | `Figure1B.R` | Spike IgG titers 2–8 weeks post-2nd dose (HD / MGUS / SMM / MM) |
| 1C | `Figure1C.R` | Spike IgG titers 2–4 months post-2nd dose (waning timepoint) |
| 1D | `Figure1D.R` | Per-individual serial titer trajectories |
| 1E | `Figure1E.R` | Linear mixed-effects model of waning slope (MGUS vs SMM-untreated) |
| 1F | `Figure1F.R` | Spike IgG titers after the 3rd dose (SMM split into Untreated / Treated) |

## Statistical framework

Cross-sectional disease-vs-HD contrasts use age- and sex-adjusted rank-based ANCOVA with Benjamini–Hochberg correction. Ordered-trend statistics use the Jonckheere–Terpstra test on age+sex residuals (helper `jt_test_residuals_age_sex()` in `../config.R`). The significance threshold throughout is q (or p) < 0.1.

## Inputs

- `data/elisa/elisa_spike_post2nd.csv` — de-identified post-2nd-dose ELISA titers (panels 1B/1C).
- `data/elisa/elisa_spike_post3rd.csv` — de-identified post-3rd-dose ELISA titers (panel 1F).
- `data/elisa/elisa_serial_titers_all.csv` and `elisa_serial_titers_filtered.csv` — serial titer trajectories (panels 1D/1E).

All inputs are hosted on Zenodo (DOI [10.5281/zenodo.18989223](https://doi.org/10.5281/zenodo.18989223)). Edit `SCRNA_DIR` in `../config.py` / `../config.R` to point at your local copy. If you mirror the deposit to an internal Google Cloud Storage bucket (e.g. `gs://your-bucket/impact_data/`), `gsutil cp -r gs://your-bucket/impact_data/elisa data/` reproduces the layout expected by these scripts.

## Run

```bash
cd Figure1
Rscript Figure1B.R   # outputs ../figures/Figure1B.png
Rscript Figure1C.R
Rscript Figure1D.R
Rscript Figure1E.R
Rscript Figure1F.R
```
