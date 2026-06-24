# Figure 2 — MMR titers and APRIL depletion

Code to reproduce Figure 2.

| Panel | Script | Description |
|-------|--------|-------------|
| 2B | `01_Figure2B.R` | MMR (Measles / Mumps / Rubella) antibody titers, HD vs SMM |
| 2D | `02_Figure2D.R` | Plasma APRIL (TNFSF13) Olink levels, paired pre/post-vaccination |
| 2E | `03_Figure2E.R` | Post-vaccination APRIL vs serum M-spike Spearman correlation, treatment-naive SMM |
| 2F | `04_Figure2F.py` | Paired pre/post-teclistamab plasma Olink for APRIL / BAFF / sBCMA (n=10 HRSMM on the teclistamab arm of Immuno-PRISM, NCT05469893); composite 3-protein decoy-sink readout |

## Statistical framework

Cross-sectional disease-vs-HD contrasts use age- and sex-adjusted rank-based ANCOVA (Age modeled via Age_Range midpoint when ages are binned). Paired pre/post comparisons use Wilcoxon signed-rank. Panel 2E reports Spearman ρ with a nonparametric bootstrap CI (10,000 resamples) and a log10 Pearson sensitivity. Panel 2F reports paired Wilcoxon signed-rank q-values with Benjamini-Hochberg correction across the three proteins, plus median delta-NPX as the magnitude annotation (NPX is the Olink log2-scale unit, so delta-NPX is the directly-measured paired difference).

Treatment-naive Olink cohort: the single MGUS patient with prior systemic therapy is excluded from panels 2D, 5A/B and 5D-F (R2 revision).

## Inputs

- `data/elisa/elisa_mmr.csv` — de-identified MMR ELISA (panel 2B).
- `data/olink/olink_paired_prepost.csv` — paired pre/post Olink panel (panel 2D).
- `data/olink/olink_summary_tumor_burden.csv` — Olink with linked tumor-burden metrics (panel 2E).
- `data/external/PrePostTEC_olink_deid.csv` — de-identified paired plasma Olink (APRIL / BAFF / sBCMA) pre- and post-teclistamab in 10 HRSMM participants from Immuno-PRISM (panel 2F).

All inputs are hosted on Zenodo (DOI [10.5281/zenodo.18989222](https://doi.org/10.5281/zenodo.18989222)). To mirror the deposit to an internal Google Cloud Storage bucket, copy the Zenodo file tree into your bucket and use `gsutil cp -r gs://your-bucket/impact_data/olink data/` to populate the layout the scripts expect.

## Run

```bash
cd 02_Figure2
Rscript 01_Figure2B.R
Rscript 02_Figure2D.R
Rscript 03_Figure2E.R
python3 04_Figure2F.py
```
