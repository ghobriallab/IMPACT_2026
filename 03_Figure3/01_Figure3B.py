#!/usr/bin/env python3
"""Figure 3B (manuscript Figure 3B): full UMAP with Annotation_Level_2 labels.

Purpose:      Figure 3B: full UMAP of the integrated scRNA-seq object (~1.1M clean cells) colored by Annotation_Level_2 cell-type labels.

Inputs:       H5AD_ANNOTATED (the de-identified comprehensive deposit, scRNAseq_IMPACT_Zenodo.h5ad).

Outputs:      figures/Figure3B.png (and PDF + SVG).

Dependencies: Python + scanpy, matplotlib, numpy, pandas; reads config.py.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import *

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib.patheffects import withStroke

# Every figure in this paper is set in Arial. Register the metrically identical Liberation Sans
# as a fallback and rewrite the SVG font-family afterwards; svg.fonttype="none" keeps text editable.
import glob as _glob
from matplotlib import font_manager as _fm
FONT = "Arial"
_av = {f.name for f in _fm.fontManager.ttflist}
if FONT not in _av:
    for _p in _glob.glob("/usr/share/fonts/**/LiberationSans-*.ttf", recursive=True):
        _fm.fontManager.addfont(_p)
    _av = {f.name for f in _fm.fontManager.ttflist}
PLOT_FONT = FONT if FONT in _av else ("Liberation Sans" if "Liberation Sans" in _av else "sans-serif")
plt.rcParams["font.family"] = PLOT_FONT
plt.rcParams["svg.fonttype"] = "none"


def save_figure(basename, dpi=300):
    """Write PNG + PDF + SVG; rewrite the SVG font-family to Arial (see note above)."""
    for ext in ("png", "pdf", "svg"):
        out = FIGURES_DIR / f"{basename}.{ext}"
        plt.savefig(out, dpi=dpi, bbox_inches="tight", facecolor="white")
        if ext == "svg" and PLOT_FONT != FONT:
            t = Path(out).read_text(encoding="utf-8")
            for q in ('"', "'"):
                t = t.replace(f"font-family: {q}{PLOT_FONT}{q}", f"font-family: {FONT}")
            t = t.replace(f"font-family: {PLOT_FONT}", f"font-family: {FONT}")
            Path(out).write_text(t, encoding="utf-8")
        print(f"Saved: {out.name}")
plt.rcParams['font.size'] = 14

# Load data
adata = sc.read_h5ad(H5AD_ANNOTATED)
print(f"Loaded {adata.n_obs:,} cells")

# Remove QC_removed and doublets
mask_qc = adata.obs['Annotation_Level_2'] != 'QC_removed'
mask_doublets = ~adata.obs['Annotation_Level_2'].str.startswith('db:')
mask_cll = adata.obs['Annotation_Level_2'] != 'CLL'
mask_clean = mask_qc & mask_doublets & mask_cll

n_cells_clean = mask_clean.sum()
print(f"Cells to plot: {n_cells_clean:,}")

# Extract coordinates and annotations
umap_coords = adata.obsm['X_umap'][mask_clean.values]
annotations = adata.obs.loc[mask_clean, 'Annotation_Level_2'].values

unique_cats = pd.Series(annotations).unique()
n_cats = len(unique_cats)

# Assign colors
if n_cats <= 20:
    cmap = colormaps['tab20']
    color_dict = {cat: cmap(i % 20) for i, cat in enumerate(unique_cats)}
else:
    colors1 = [colormaps['tab20'](i) for i in range(20)]
    colors2 = [colormaps['tab20b'](i) for i in range(20)]
    colors3 = [colormaps['tab20c'](i) for i in range(20)]
    all_colors = colors1 + colors2 + colors3
    color_dict = {cat: all_colors[i % 60] for i, cat in enumerate(unique_cats)}

# Cluster centroids for label placement
centroids = {}
for label in np.unique(annotations):
    mask = annotations == label
    centroids[label] = np.median(umap_coords[mask], axis=0)

# Manual label offsets for crowded regions
manual_offsets = {
    'cDC1': (2, -0.8), 'cDC2': (1.5, 0.8), 'acDC2': (2, 1.5),
    'MoDC': (2, 2.2), 'Cytokine+ MoDC': (2.5, 4), 'Cytokine+ cDC2': (2.5, -1.2),
    'ccDC2': (2, -2), 'IFN+ cDC2': (2.5, -0.3),
    'NBC': (1.5, -1.5), 'TBC': (0, -2), 'aNBC': (-1.5, -1.5),
    'IFN+ BC': (-2, -0.8), 'aMBC': (-1.8, 0.8),
    'aNCSMBC': (1.8, 0.8), 'MBC': (2, 0), 'NCSMBC': (2, -0.8),
    'ABC': (1.5, 1.2),
    'HLA-DR-high CD14+ Mono': (1.5, 1), 'IFN+ CD14+ Mono': (2, -0.5),
    'cNK': (-1.5, 0.5), 'IFN+ NK': (-1.5, 1),
    'AS-DC': (1.5, 0.5),
}

# Plot. Match the originally-published Figure 3B proportions (optimized for small panel
# in a large figure): rectangular ~10x11.2, auto aspect, font sizes unchanged.
fig, ax = plt.subplots(figsize=(10, 11.2))

for cat in unique_cats:
    mask = annotations == cat
    ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
               c=[color_dict[cat]], s=2, alpha=0.6, label=cat, rasterized=True)

for cluster, (cx, cy) in centroids.items():
    if cluster in manual_offsets:
        ox, oy = manual_offsets[cluster]
        tx, ty = cx + ox, cy + oy
        ax.plot([cx, tx], [cy, ty], 'k-', lw=1, alpha=0.6)
        ax.text(tx, ty, cluster, fontsize=14, fontweight='bold',
                ha='center', va='center', color='black',
                path_effects=[withStroke(linewidth=5, foreground='white')])
    else:
        ax.text(cx, cy, cluster, fontsize=14, fontweight='bold',
                ha='center', va='center', color='black',
                path_effects=[withStroke(linewidth=5, foreground='white')])

ax.set_title(f'PBMCs (n={n_cells_clean:,})', fontsize=26, fontweight='bold', pad=20)
ax.set_xlabel('UMAP1', fontsize=18)
ax.set_ylabel('UMAP2', fontsize=18)
ax.tick_params(labelsize=16)
ax.set_aspect('auto')

plt.tight_layout()
save_figure("Figure3B")
plt.close()

