# Figure 4: T cell responses

| Panel | Script | Description |
|-------|--------|-------------|
| 4B | `01_Figure4B.R` | Spike-specific TCR clonotype proportions, paired pre vs post |
| 4C | `02_Figure4C.R` | CEF (CMV / EBV / influenza) recall clonotype proportions, the recall control for 4B |
| 4E | `03_Figure4E.R` | IFN-gamma ELISPOT for spike and CERI pools, HD vs SMM |

Clonotypes are called against the curated spike- and CEF-specific reference panels listed in
the Supplementary Tables.

## Statistics

Panels 4B and 4C compare pre with post within each group using a two-sided paired Wilcoxon
signed-rank test, Benjamini-Hochberg corrected across the groups of the panel. Panel 4E compares
HD with SMM using age- and sex-adjusted rank-based ANCOVA, run separately for each peptide pool;
counts are DMSO-normalized and averaged over technical duplicates first. Significance threshold
q (or p) < 0.1.

## Inputs

| File | Panels |
|------|--------|
| `data/tcr/tcr_clonotype_proportions.rds` | 4B, 4C |
| `data/elisa/elispot_spike_cef.csv` | 4E |

All inputs come from the Zenodo deposit (DOI [10.5281/zenodo.18989222](https://doi.org/10.5281/zenodo.18989222));
see the root README for the expected `data/` layout.

## Run

```bash
cd 04_Figure4
Rscript 01_Figure4B.R    # writes ../figures/Figure4B.{png,pdf,svg}
Rscript 02_Figure4C.R
Rscript 03_Figure4E.R
```
