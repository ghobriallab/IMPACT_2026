# Figure 4 — T cell responses

Code to reproduce Figure 4.

| Panel | Script | Description |
|-------|--------|-------------|
| 4B | `Figure4B.R` | SARS-CoV-2 spike-specific TCR clonotype proportions, paired pre vs post (HD/MGUS/SMM) |
| 4C | `Figure4C.R` | CEF (CMV/EBV/Flu) recall TCR clonotype proportions, paired pre vs post (control for recall vs de novo) |
| 4E | `Figure4E.R` | IFN-gamma ELISPOT for SARS spike and CEF pools, HD vs SMM, age+sex adjusted |

## Statistical framework

Paired pre vs post within group: Wilcoxon signed-rank. Cross-group differences: age- and sex-adjusted rank-based ANCOVA (panel 4E) and Mann–Whitney with BH correction (panels 4B/4C). ELISPOT values are DMSO-normalized and technical-duplicate-averaged before group comparison.

The TCR clonotype data are computed against curated spike-specific and CEF-specific TCR reference panels (panel of curated reference TCRs deposited in the Supplementary Tables).

## Inputs

- `data/tcr/tcr_clonotype_proportions.rds` — per-patient ClusTCR proportions, pre and post vaccination (panels 4B/4C).
- `data/elisa/elispot_spike_cef.csv` — well-level ELISPOT counts (panel 4E).

All inputs are hosted on Zenodo (DOI [10.5281/zenodo.18989223](https://doi.org/10.5281/zenodo.18989223)). If you mirror to an internal Google Cloud Storage bucket, `gsutil cp -r gs://your-bucket/impact_data/tcr data/` and `gsutil cp -r gs://your-bucket/impact_data/elisa data/` reproduce the expected layout.

## Run

```bash
cd Figure4
Rscript Figure4B.R
Rscript Figure4C.R
Rscript Figure4E.R
```
