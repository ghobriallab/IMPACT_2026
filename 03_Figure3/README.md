# Figure 3 — Impaired APRIL signaling in B cells (scRNA-seq)

Code to reproduce Figure 3 from the integrated peripheral-blood scRNA-seq (1,109,633 clean cells; QC_removed, doublets, Platelets and CLL are excluded by these scripts; the deposit retains them with explicit labels for traceability).

| Panel | Script | Description |
|-------|--------|-------------|
| 3B (UMAP) | `Figure3C.py` | Full UMAP with Annotation_Level_2 cell-type labels |
| 3C (myeloid APRIL) | `Figure3D.py` | Myeloid TNFSF13 expression across HD/MGUS/SMM, pre vs post vaccination, age+sex adjusted |
| 3D (B-cell APRIL module) | `Figure3E.py` | APRIL-responsive module score in B cells, treatment-naive HD/MGUS/SMM |
| 3E (BM external validation) | `Figure3F.py` | Boiarsky et al. GSE193531: APRIL module in normal vs malignant PCs across NBM/MGUS/SMM/NDMM |

*Note on panel-letter drift:* the script names are historical (when Fig 3 had a different layout); the table above maps each script to the final manuscript panel.

## Statistical framework

Per-timepoint, age- and sex-adjusted rank-based ANCOVA + BH for the myeloid panel. Jonckheere–Terpstra ordered-trend test across the disease continuum for the B-cell module score and the external normal-PC validation. Pairwise Mann–Whitney + BH correction for selected contrasts. Significance threshold q (or p) < 0.1.

## Inputs

- `scRNAseq_IMPACT_Zenodo.h5ad` — the de-identified comprehensive scRNA-seq object (1.4M cells, 42k genes, normalized log1p X + int64 `counts` layer + `obsm['X_umap']` + obs columns `Annotation_Level_1`, `Annotation_Level_2`, `Diagnosis`, `Timepoint`, `TreatmentStatus`). Download from Zenodo (DOI [10.5281/zenodo.18989223](https://doi.org/10.5281/zenodo.18989223)) and point `SCRNA_DIR` in `../config.py` at the containing directory. If you mirror the deposit to an internal Google Cloud Storage bucket, `gsutil cp gs://your-bucket/impact_data/scRNAseq_IMPACT_Zenodo.h5ad .` puts it where the scripts expect.
- External validation (panel 3E): `data/external/GSE193531_cell-level-metadata.csv` (companion cell-level malignant/normal calls from inferCNV) shipped in the Zenodo deposit; the UMI matrix is fetched from GSE193531 directly.

## Run

```bash
cd Figure3
python Figure3C.py
python Figure3D.py
python Figure3E.py
python Figure3F.py
```
