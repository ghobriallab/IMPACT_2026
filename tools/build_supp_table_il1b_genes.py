#!/usr/bin/env python3
# ============================================================================
# Purpose:      Build the IL-1beta-response gene table underlying Figure 5C.
#               Reproduces exactly the gene selection performed by
#               05_Figure5/02_Figure5C.py, which scores the signature on the
#               2,678-gene high-variance analysis universe:
#                   available_genes = [g for g in il1b_genes if g in adata.var_names]
#               after `adata` has been restricted to that HVG list. Genes are
#               therefore reported with two separate flags: whether they are
#               present in the sequenced object at all, and whether they survive
#               the restriction to the analysis universe and are actually scored.
# Inputs:       data/il1b_response_genes_human.csv   (1,475 human-converted genes)
#               data/hvg_2678_genes.txt              (2,678-gene analysis universe)
#               H5AD_IL1B from config.py (or $IMPACT_H5AD), optional: only used to
#               flag which genes are present in the object. Blank column if absent.
# Outputs:      tables/Supplementary_Table_7_IL1B-response-genes.csv
# Dependencies: Python + pandas (h5py only when the h5ad is available); reads config.py
# ============================================================================
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import DATA_DIR, H5AD_IL1B

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
OUT_DIR = REPO / "tables"
OUT = OUT_DIR / "Supplementary_Table_7_IL1B-response-genes.csv"

il1b = pd.read_csv(DATA_DIR / "il1b_response_genes_human.csv")["gene"].tolist()
hvg = set(l.strip() for l in open(DATA_DIR / "hvg_2678_genes.txt") if l.strip())

# Presence in the sequenced object. Read from the var index only, so the 9.6 GB
# matrix is never loaded. Left blank if the deposit is not available locally.
# Allow pointing at a downloaded Zenodo copy without editing config.py.
H5 = Path(os.environ.get("IMPACT_H5AD", str(H5AD_IL1B)))
detected = None
if H5.exists():
    import h5py
    with h5py.File(H5, "r") as f:
        vi = f["var"]["_index"][:]
    detected = set(x.decode() if isinstance(x, bytes) else str(x) for x in vi)
else:
    print(f"NOTE: {H5} not found; 'Detected_in_scRNAseq_dataset' left blank.")

df = pd.DataFrame([{
    "Gene": g,
    "Detected_in_scRNAseq_dataset": "" if detected is None else ("Yes" if g in detected else "No"),
    "Retained_in_analysis_universe_and_scored": "Yes" if g in hvg else "No",
} for g in il1b]).sort_values("Gene", kind="mergesort").reset_index(drop=True)

OUT_DIR.mkdir(exist_ok=True)
df.to_csv(OUT, index=False)

n_scored = int((df["Retained_in_analysis_universe_and_scored"] == "Yes").sum())
print(f"signature genes (total)              : {len(df)}")
if detected is not None:
    print(f"detected in the dataset              : {int((df['Detected_in_scRNAseq_dataset']=='Yes').sum())}")
print(f"retained in analysis universe, scored: {n_scored}")
assert len(df) == 1475, len(df)
assert n_scored == 300, n_scored
print(f"wrote {OUT}")
