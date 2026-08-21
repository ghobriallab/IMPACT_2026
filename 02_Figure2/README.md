# Figure 2: MMR titers and APRIL depletion

| Panel | Script | Description |
|-------|--------|-------------|
| 2B | `01_Figure2B.R` | MMR (measles / mumps / rubella) antibody titers, HD vs SMM |
| 2D | `02_Figure2D.R` | Panel-wide Olink screen: q-values for every protein at both timepoints, with APRIL as the top hit |
| 2E | `03_Figure2E.R` | Plasma APRIL (TNFSF13), paired pre vs post vaccination |
| 2F | `04_Figure2F.R` | Post-vaccination APRIL vs serum M-spike, treatment-naive SMM |
| 2G | `05_Figure2G.py` | Paired pre vs post-teclistamab measurements in 10 participants with HR-SMM on the teclistamab arm of Immuno-PRISM (NCT05469893): serum M-spike, sBCMA, APRIL and BAFF |

## Statistics

Panel 2D screens all 52 panel proteins against HD at both timepoints with age- and sex-adjusted
rank-based ANCOVA, Benjamini-Hochberg corrected across all 156 tests (52 proteins x 3 contrasts)
within a timepoint. Age and sex come from `data/elisa/elisa_cohort_demographics.csv`, not from the
Olink file, whose `Age` column carries one fill value for every healthy donor. Panel 2E draws the
APRIL contrasts from that same 156-test family. The full ranked screen for both timepoints is
written to `tables/` as figure source data, carrying adjusted and unadjusted q side by side.

Panel 2F reports Spearman rho with a nonparametric bootstrap CI (10,000 resamples, percentile
method) and a leave-one-out sensitivity analysis. Panel 2G reports paired Wilcoxon signed-rank
q-values, BH corrected across its four sub-panels, annotated with the median linear fold-change.
Significance threshold q (or p) < 0.1.

## Inputs

| File | Panels |
|------|--------|
| `data/elisa/elisa_mmr.csv` | 2B |
| `data/olink/olink_paired_prepost.csv` | 2D, 2E |
| `data/elisa/elisa_cohort_demographics.csv` | 2D, 2E covariates |
| `data/olink/olink_summary_tumor_burden.csv` | 2F |
| `data/external/PrePostTEC_olink_deid.csv` | 2G |

All inputs come from the Zenodo deposit (DOI [10.5281/zenodo.18989222](https://doi.org/10.5281/zenodo.18989222));
see the root README for the expected `data/` layout.

## Run

```bash
cd 02_Figure2
Rscript 01_Figure2B.R    # writes ../figures/Figure2B.{png,pdf,svg}
Rscript 02_Figure2D.R
Rscript 03_Figure2E.R
Rscript 04_Figure2F.R
python3 05_Figure2G.py
```
