# Figure 5: Blunted innate and cytokine responses

| Panel | Script | Description |
|-------|--------|-------------|
| 5A, 5B | `01_Figure5AB.R` | Plasma IL-1B and IL-18 (Olink), paired pre vs post vaccination |
| 5C | `02_Figure5C.py` | IL-1B response gene signature score (scRNA-seq), paired pre vs post |
| 5D, 5E, 5F | `03_Figure5DEF.R` | Plasma DDX58 (RIG-I), NUB1 and MMP7 (Olink), paired pre vs post |

`03_Figure5DEF.R` sources `01_Figure5AB.R` for its data prep and the `plot_cytokine()` helper,
so run 5AB first.

## Statistics

Paired pre vs post uses a two-sided Wilcoxon signed-rank test. Paired designs are inherently
controlled for age and sex, since each individual is their own baseline. Olink q-values come from
the panel-wide Benjamini-Hochberg family of 156 tests (52 analytes x 3 disease groups) defined in
`01_Figure5AB.R` and inherited by `03_Figure5DEF.R`, matching the Methods. The IL-1B signature is
scored on cells re-normalized from the counts layer against a 2,678-HVG control pool, after
excluding QC-failed cells, doublets, platelets and CLL. Significance threshold q (or p) < 0.1.

The single MGUS participant with prior systemic therapy is excluded from the Olink panels, which
are restricted to treatment-naive individuals.

## Inputs

| File | Panels |
|------|--------|
| `data/olink/olink_cytokines.csv` | 5A, 5B, 5D, 5E, 5F |
| `scRNAseq_IMPACT_Zenodo.h5ad`, `data/hvg_2678_genes.txt`, `data/il1b_response_genes_human.csv` | 5C |

All inputs come from the Zenodo deposit (DOI [10.5281/zenodo.18989222](https://doi.org/10.5281/zenodo.18989222));
see the root README for the expected `data/` layout and for where to point `SCRNA_DIR`.

## Run

```bash
cd 05_Figure5
Rscript 01_Figure5AB.R    # writes ../figures/Figure5A.* and Figure5B.*
python3 02_Figure5C.py
Rscript 03_Figure5DEF.R   # sources 01_Figure5AB.R
```
