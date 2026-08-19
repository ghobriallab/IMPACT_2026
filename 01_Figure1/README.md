# Figure 1: Antibody waning

| Panel | Script | Description |
|-------|--------|-------------|
| 1B | `01_Figure1B.R` | Spike IgG titers 2-8 weeks post-2nd dose (HD / MGUS / SMM / MM) |
| 1C | `02_Figure1C.R` | Spike IgG titers 2-4 months post-2nd dose |
| 1D | `03_Figure1D.R` | Per-individual serial titer trajectories |
| 1E | `04_Figure1E.R` | Linear mixed-effects model of waning slope (MGUS vs untreated SMM) |
| 1F | `05_Figure1F.R` | Spike IgG titers after the 3rd dose, SMM split into untreated and treated |

## Statistics

Cross-sectional disease-vs-HD contrasts use age- and sex-adjusted rank-based ANCOVA with
Benjamini-Hochberg correction. Ordered trends use the Jonckheere-Terpstra test on age+sex
residuals (`jt_test_residuals_age_sex()` in `../config.R`). Significance threshold q < 0.1.

## Inputs

| File | Panels |
|------|--------|
| `data/elisa/elisa_spike_post2nd.csv` | 1B, 1C |
| `data/elisa/elisa_spike_post3rd.csv` | 1F |
| `data/elisa/elisa_serial_titers_all.csv`, `elisa_serial_titers_filtered.csv` | 1D, 1E |

All inputs come from the Zenodo deposit (DOI [10.5281/zenodo.18989222](https://doi.org/10.5281/zenodo.18989222));
see the root README for the expected `data/` layout.

## Run

```bash
cd 01_Figure1
Rscript 01_Figure1B.R    # writes ../figures/Figure1B.{png,pdf,svg}
Rscript 02_Figure1C.R
Rscript 03_Figure1D.R
Rscript 04_Figure1E.R
Rscript 05_Figure1F.R
```
