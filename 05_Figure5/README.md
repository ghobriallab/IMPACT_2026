# Figure 5 — Blunted innate / cytokine responses

Code to reproduce Figure 5.

| Panel | Script | Description |
|-------|--------|-------------|
| 5A / 5B | `Figure5AB.R` | IL-1B and IL-18 plasma protein levels (Olink), paired pre vs post vaccination |
| 5C | `Figure5C.py` | IL-1B response gene signature score (scRNA-seq), paired pre vs post |
| 5D / 5E / 5F | `Figure5DEF.R` | DDX58 (RIG-I), NUB1 and MMP7 Olink levels, paired pre vs post |

`Figure5DEF.R` `source()`s `Figure5AB.R` to reuse its data prep and the `plot_cytokine()` helper. Run order: 5AB before 5DEF.

## Statistical framework

Paired pre vs post within group: Wilcoxon signed-rank. Significance threshold p < 0.1; paired tests are inherently controlled for age and sex (each individual is their own baseline). The IL-1B signature score is re-normalized from the counts layer on a 2,678-HVG control pool before `score_genes`.

Treatment-naive Olink cohort: the single MGUS patient with prior systemic therapy is excluded from panels 5A/B and 5D-F (R2 revision).

## Inputs

- `data/olink/olink_cytokines.csv` — paired Olink cytokine matrix (panels 5A/B/D/E/F).
- `scRNAseq_IMPACT_Zenodo.h5ad` + `data/hvg_2678_genes.txt` (control pool) + `data/il1b_response_genes_human.csv` (gene set) — for panel 5C.

All inputs are hosted on Zenodo (DOI [10.5281/zenodo.18989222](https://doi.org/10.5281/zenodo.18989222)). To mirror the deposit to an internal Google Cloud Storage bucket: `gsutil cp gs://your-bucket/impact_data/scRNAseq_IMPACT_Zenodo.h5ad .` and `gsutil cp -r gs://your-bucket/impact_data/olink data/`.

## Run

```bash
cd Figure5
Rscript Figure5AB.R    # outputs Figure5A.png, Figure5B.png
python  Figure5C.py    # outputs Figure5C.png
Rscript Figure5DEF.R   # outputs Figure5D/E/F.png (sources Figure5AB.R)
```
