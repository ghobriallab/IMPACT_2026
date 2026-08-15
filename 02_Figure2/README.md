# Figure 2: MMR titers and APRIL depletion

Code to reproduce Figure 2.

| Panel | Script | Description |
|-------|--------|-------------|
| 2B | `01_Figure2B.R` | MMR (Measles / Mumps / Rubella) antibody titers, HD vs SMM |
| 2D | `02_Figure2D.R` | Panel-wide Olink screen: pre- vs post-vaccination q-values for every protein, showing APRIL as the top hit rather than a pre-selected candidate |
| 2E | `03_Figure2E.R` | Plasma APRIL (TNFSF13) Olink levels, paired pre/post-vaccination |
| 2F | `04_Figure2F.R` | Post-vaccination APRIL vs serum M-spike Spearman correlation, treatment-naive SMM |
| 2G | `05_Figure2G.py` | Paired pre/post-teclistamab plasma Olink for APRIL / BAFF / sBCMA (n=10 HRSMM on the teclistamab arm of Immuno-PRISM, NCT05469893); four sub-panels: serum M-spike, sBCMA, APRIL and BAFF |

## Statistical framework

Panel 2D screens all 52 panel proteins against HD at both timepoints using age- and sex-adjusted rank-based ANCOVA, Benjamini-Hochberg corrected across all 156 tests (52 proteins x 3 pairwise contrasts) within a timepoint. Age and sex come from `data/elisa/elisa_cohort_demographics.csv`, not from the Olink file, whose `Age` column carries a single fill value for every healthy donor. The full ranked results for both timepoints are written to `tables/` as figure source data, carrying the adjusted and unadjusted q side by side.

Cross-sectional disease-vs-HD contrasts use age- and sex-adjusted rank-based ANCOVA (Age modeled via Age_Range midpoint when ages are binned). Paired pre/post comparisons use Wilcoxon signed-rank. Panel 2E draws the APRIL contrasts from the same 156-test family as panel 2D. Panel 2F reports Spearman ρ with a nonparametric bootstrap CI (10,000 resamples, percentile method) and a leave-one-out sensitivity analysis. Panel 2G reports paired Wilcoxon signed-rank q-values with Benjamini-Hochberg correction across its four sub-panels (M-spike, sBCMA, APRIL, BAFF), annotated with the median linear fold-change.

Treatment-naive Olink cohort: the single MGUS patient with prior systemic therapy is excluded from panels 2D, 2E, 5A/B and 5D-F (R2 revision).

## Inputs

- `data/elisa/elisa_mmr.csv`: de-identified MMR ELISA (panel 2B).
- `data/olink/olink_paired_prepost.csv`: paired pre/post Olink panel (panels 2D and 2E).
- `data/elisa/elisa_cohort_demographics.csv`: de-identified per-participant age and sex (panel 2D covariates).
- `data/olink/olink_summary_tumor_burden.csv`: Olink with linked tumor-burden metrics (panel 2F).
- `data/external/PrePostTEC_olink_deid.csv`: de-identified paired plasma Olink (APRIL / BAFF / sBCMA) pre- and post-teclistamab in 10 HRSMM participants from Immuno-PRISM (panel 2G).

All inputs are hosted on Zenodo (DOI [10.5281/zenodo.18989222](https://doi.org/10.5281/zenodo.18989222)). To mirror the deposit to an internal Google Cloud Storage bucket, copy the Zenodo file tree into your bucket and use `gsutil cp -r gs://your-bucket/impact_data/olink data/` to populate the layout the scripts expect.

## Run

```bash
cd 02_Figure2
Rscript 01_Figure2B.R
Rscript 02_Figure2D.R
Rscript 03_Figure2E.R
Rscript 04_Figure2F.R
python3 05_Figure2G.py
```
