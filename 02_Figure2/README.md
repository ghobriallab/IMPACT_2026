# Figure 2 — MMR titers and APRIL depletion

Code to reproduce Figure 2.

| Panel | Script | Description |
|-------|--------|-------------|
| 2B | `Figure2B.R` | MMR (Measles / Mumps / Rubella) antibody titers, HD vs SMM |
| 2D | `Figure2D.R` | Plasma APRIL (TNFSF13) Olink levels, paired pre/post-vaccination |
| 2E | `Figure2E.R` | Post-vaccination APRIL vs serum M-spike Spearman correlation, treatment-naive SMM |

## Statistical framework

Cross-sectional disease-vs-HD contrasts use age- and sex-adjusted rank-based ANCOVA (Age modeled via Age_Range midpoint when ages are binned). Paired pre/post comparisons use Wilcoxon signed-rank. Panel 2E reports Spearman ρ with a nonparametric bootstrap CI (10,000 resamples) and a log10 Pearson sensitivity.

Treatment-naive Olink cohort: the single MGUS patient with prior systemic therapy is excluded from panels 2D, 5A/B and 5D-F (R2 revision).

## Inputs

- `data/elisa/elisa_mmr.csv` — de-identified MMR ELISA (panel 2B).
- `data/olink/olink_paired_prepost.csv` — paired pre/post Olink panel (panel 2D).
- `data/olink/olink_summary_tumor_burden.csv` — Olink with linked tumor-burden metrics (panel 2E).

All inputs are hosted on Zenodo (DOI [10.5281/zenodo.18989223](https://doi.org/10.5281/zenodo.18989223)). To mirror the deposit to an internal Google Cloud Storage bucket, copy the Zenodo file tree into your bucket and use `gsutil cp -r gs://your-bucket/impact_data/olink data/` to populate the layout the scripts expect.

## Run

```bash
cd Figure2
Rscript Figure2B.R
Rscript Figure2D.R
Rscript Figure2E.R
```
