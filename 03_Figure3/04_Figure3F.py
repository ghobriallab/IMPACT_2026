#!/usr/bin/env python3
"""Figure 3F (= manuscript Figure 3E): APRIL-responsive gene module score in plasma cells from
external bone-marrow scRNA-seq (Boiarsky et al. 2022, GSE193531), with malignant ('neoplastic')
and non-malignant ('normal') plasma cells classified by inferCNV by the original authors.

Side-by-side normal vs neoplastic PCs at each disease stage along the tumor-mass axis
(NBM normal | MGUS normal+neoplastic | SMM normal+neoplastic | NDMM normal+neoplastic).
Statistic of record (shown on the figure): Jonckheere-Terpstra ordered trend across normal PCs
(NBM -> NDMM) and across neoplastic PCs (MGUS -> NDMM). Both are reported as on the figure;
JT(normal) reaches significance, JT(neoplastic) is a directionally consistent non-significant
trend (small samples in MGUS-neoplastic and NDMM-normal -- residual normal PCs are rare in
MM bone marrows). Pairwise within-stage normal-vs-neoplastic comparisons are computed for the
figure caption but not annotated on the plot (all BH-q > 0.4 at the available sample sizes).

Purpose:      Figure 3F (= manuscript Figure 3E): external-validation APRIL-responsive gene module score in plasma cells from Boiarsky et al. (Nat Commun 2022, GSE193531), side-by-side normal vs malignant PCs across NBM/MGUS/SMM/NDMM, with JT ordered-trend test on normal PCs.

Inputs:       External GSE193531 data + companion cell-level metadata (CNV-based malignant calls).

Outputs:      figures/Figure3F.png (and SVG).

Dependencies: Python + scanpy, numpy, pandas, scipy, matplotlib; reads config.py.
"""
import sys
import warnings
warnings.filterwarnings('ignore')
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import *

import pandas as pd
import numpy as np
import scanpy as sc
from scipy import stats
from scipy.sparse import csr_matrix
from statsmodels.stats.multitest import multipletests
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
# Make mathtext bold ($\mathbf{...}$) render in DejaVu Sans so it matches the surrounding text
mpl.rcParams['mathtext.fontset'] = 'dejavusans'

np.random.seed(42)

# Same 15-gene module as Figure 3D / SupFig 3
APRIL_TARGETS = {
    'Survival': ['MCL1', 'BCL2'],
    'NF-kB_pathway': ['NFKB2', 'RELB', 'NFKB1', 'NFKBIA'],
    'Adhesion': ['CD44', 'ICAM1'],
    'Immunomodulation': ['CCL3', 'CCL4', 'VEGFA', 'IL10', 'CD274', 'TGFB1', 'CXCL8'],
}
ALL_APRIL_GENES = [g for genes in APRIL_TARGETS.values() for g in genes]

EXT_DATA_DIR = DATA_DIR / "external"
MIN_CELLS_PER_SAMPLE = 10   # per-sample mean only included when >=10 PCs of that population in the sample

print("Loading GSE193531 (Boiarsky 2022)...")
metadata = pd.read_csv(EXT_DATA_DIR / "GSE193531_cell-level-metadata.csv", index_col=0)
print(f"  metadata: {len(metadata):,} cells")
counts = pd.read_csv(EXT_DATA_DIR / "GSE193531_umi-count-matrix.csv.gz", index_col=0)
print(f"  count matrix: {counts.shape[0]:,} genes x {counts.shape[1]:,} cells")

counts_T = counts.T
common = counts_T.index.intersection(metadata.index)
counts_T = counts_T.loc[common]; metadata = metadata.loc[common]
adata = sc.AnnData(X=csr_matrix(counts_T.values), obs=metadata,
                   var=pd.DataFrame(index=counts_T.columns))
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
available = [g for g in ALL_APRIL_GENES if g in adata.var_names]
print(f"  APRIL genes present: {len(available)}/{len(ALL_APRIL_GENES)}")
sc.tl.score_genes(adata, gene_list=available, score_name='APRIL_score', random_state=42)

obs = adata.obs[['sample_ID', 'disease_stage', 'normal_or_neoplastic']].copy()
obs['score'] = adata.obs['APRIL_score'].values
obs = obs.dropna(subset=['normal_or_neoplastic'])

# Per-(sample, neoplastic-status) aggregate; require >= MIN_CELLS_PER_SAMPLE
agg = (obs.groupby(['sample_ID', 'disease_stage', 'normal_or_neoplastic'])
          .agg(n_cells=('score', 'size'), mean=('score', 'mean'))
          .reset_index())
agg = agg[agg['n_cells'] >= MIN_CELLS_PER_SAMPLE].copy()

# Box order, left -> right along increasing plasma-cell tumor mass:
# NBM normal | MGUS normal | MGUS neo | SMM normal | SMM neo
# (NDMM dropped to match the cohort scope of Fig 3D, which has no treatment-naive MM. With only
# two neoplastic groups remaining (MGUS, SMM), the Jonckheere-Terpstra ordered-trend test cannot
# be computed for the neoplastic axis -- the test requires >=3 groups -- so only the JT across
# normal PCs (NBM -> MGUS -> SMM, 3 groups) is reported.)
# Box colors encode the population (blue shades = non-malignant normal PCs going light-to-dark
# along NBM -> SMM; red shades = malignant neoplastic PCs going light-to-dark along MGUS -> SMM).
NORMAL_SHADES     = ['#B5CEE1', '#6CA1CC', '#3577B1']   # NBM, MGUS, SMM (normal)
NEOPLASTIC_SHADES = ['#F4B6A5', '#D86A52']              # MGUS, SMM (neoplastic)
SPECS = [
    ('NBM\nnorm',  'NBM',  'normal',     NORMAL_SHADES[0], False),
    ('MGUS\nnorm', 'MGUS', 'normal',     NORMAL_SHADES[1], False),
    ('MGUS\nmalignant',  'MGUS', 'neoplastic', NEOPLASTIC_SHADES[0], True),
    ('SMM\nnorm',  'SMM',  'normal',     NORMAL_SHADES[2], False),
    ('SMM\nmalignant',   'SMM',  'neoplastic', NEOPLASTIC_SHADES[1], True),
]
NORMAL_BRACKET_COLOR     = '#0B2A4A'   # dark blue for cross-stage normal-PC comparisons
NEOPLASTIC_BRACKET_COLOR = '#5C1812'   # dark red  for cross-stage neoplastic-PC comparisons
data, colors, hatches, labels = [], [], [], []
for lab, stage, status, color, hatch in SPECS:
    vals = agg.loc[(agg.disease_stage == stage) & (agg.normal_or_neoplastic == status), 'mean'].values
    data.append(vals); colors.append(color); hatches.append('///' if hatch else None); labels.append(lab)

# Statistics ----
def jonckheere_terpstra(groups_in_order):
    """Two-sided JT z, p; normal approximation, mid-rank ties."""
    g = [np.asarray(x) for x in groups_in_order]; k = len(g); J = 0.0
    for i in range(k):
        for j in range(i + 1, k):
            xi, xj = g[i], g[j]
            for v in xi:
                J += float(np.sum(v < xj)) + 0.5 * float(np.sum(v == xj))
    n = [len(x) for x in g]; N = sum(n)
    mu = (N * N - sum(ni * ni for ni in n)) / 4.0
    var = (N * N * (2 * N + 3) - sum(ni * ni * (2 * ni + 3) for ni in n)) / 72.0
    z = (J - mu) / np.sqrt(var); p = 2 * (1 - stats.norm.cdf(abs(z)))
    return z, p

norm_groups = [data[0], data[1], data[3]]   # NBM_n, MGUS_n, SMM_n  (3 groups for JT)
z_n, p_n = jonckheere_terpstra(norm_groups)
# Neoplastic JT not reported: only 2 groups remain (MGUS_neo, SMM_neo) after dropping NDMM,
# and the Jonckheere-Terpstra test requires at least 3 ordered groups.

# Pairwise within-stage normal vs neoplastic (computed for figure as subtle gray brackets)
pair_idx = [('MGUS', 1, 2), ('SMM', 3, 4)]
pair_p = []
for stage, i_norm, i_neo in pair_idx:
    a, b = data[i_norm], data[i_neo]
    p = stats.mannwhitneyu(a, b, alternative='two-sided').pvalue if (len(a) >= 2 and len(b) >= 2) else np.nan
    pair_p.append(p)
pair_q = multipletests([p for p in pair_p if not np.isnan(p)], method='fdr_bh')[1]
qi = 0; qmap = {}
for (stage, _, _), p in zip(pair_idx, pair_p):
    qmap[stage] = (pair_q[qi] if not np.isnan(p) else np.nan)
    if not np.isnan(p): qi += 1

# Cross-stage pairwise comparisons within each population:
# Normal PCs: NBM-vs-MGUS, NBM-vs-SMM (2 comparisons, BH/2).
# Neoplastic PCs: MGUS-vs-SMM only (single comparison; reported as raw p, no BH needed).
# Indices into `data` (NBM_n=0, MGUS_n=1, MGUS_neo=2, SMM_n=3, SMM_neo=4).
norm_cross = [('NBM-MGUS', 0, 1), ('NBM-SMM', 0, 3)]
neo_cross  = [('MGUS-SMM', 2, 4)]

def cross_qs(pairs):
    ps = []
    for _, i, j in pairs:
        a, b = data[i], data[j]
        ps.append(stats.mannwhitneyu(a, b, alternative='two-sided').pvalue
                  if (len(a) >= 2 and len(b) >= 2) else np.nan)
    valid_p = [p for p in ps if not np.isnan(p)]
    qs_valid = multipletests(valid_p, method='fdr_bh')[1] if valid_p else []
    qmap, k = {}, 0
    for (lbl, _, _), p in zip(pairs, ps):
        qmap[lbl] = (qs_valid[k] if not np.isnan(p) else np.nan)
        if not np.isnan(p): k += 1
    return qmap, ps

norm_cross_q, norm_cross_p = cross_qs(norm_cross)
neo_cross_q,  neo_cross_p  = cross_qs(neo_cross)   # single comparison: q == raw p
print("\n  cross-stage NORMAL (BH/2):")
for (lbl, _, _), p in zip(norm_cross, norm_cross_p):
    print(f"    {lbl}: p={p:.4g}, q={norm_cross_q[lbl]:.4g}")
print("  cross-stage NEOPLASTIC (single, no BH):")
for (lbl, _, _), p in zip(neo_cross, neo_cross_p):
    print(f"    {lbl}: p={p:.4g}")

print(f"JT(normal across NBM -> SMM):     z={z_n:.3f}, p={p_n:.4g}  [primary figure statistic]")
print("JT(neoplastic): not computed (only 2 groups after NDMM removal; JT requires >=3)")
for (stage, _, _), p in zip(pair_idx, pair_p):
    print(f"  caption: {stage} normal vs neoplastic MWU p={p:.4g}, BH-q={qmap[stage]:.4g}")
print("counts per box (NBM_n, MGUS_n, MGUS_neo, SMM_n, SMM_neo):",
      [len(d) for d in data])

# Plot ----
# Positions: NBM at 1; MGUS pair at (3, 4.6); SMM pair at (6.8, 8.4). Within-pair gap 1.6,
# between-stage gap ~2.2 -- enough room to keep "MGUS normal" / "MGUS malignant" tick labels
# from colliding at the small panel size.
pos = [1, 3, 4.6, 6.8, 8.4]
tick_labels = [f"{lab}\nn={len(d)}" for lab, d in zip(labels, data)]

# Compact figure tuned for use as a single panel in a multi-panel figure: small physical size,
# but fonts/strokes sized to survive scale-down. Within-stage gray brackets removed (all ns at
# the available n; their absence reduces clutter without losing information that the JT and
# cross-stage brackets already convey).
fig, ax = plt.subplots(figsize=(5.4, 4.1))
bp = ax.boxplot(data, positions=pos, widths=0.78, patch_artist=True,
                medianprops=dict(color='black', linewidth=1.4))
for patch, c, h in zip(bp['boxes'], colors, hatches):
    patch.set_facecolor(c); patch.set_alpha(0.6)
    if h: patch.set_hatch(h)
rng = np.random.default_rng(42)
for i, (vals, c) in enumerate(zip(data, colors)):
    if len(vals) == 0: continue
    ax.scatter(rng.normal(pos[i], 0.10, size=len(vals)), vals, s=14, color=c,
               edgecolor='black', linewidth=0.4, alpha=0.92)
ax.set_xticks(pos); ax.set_xticklabels(tick_labels, fontsize=9)
ax.set_ylabel('APRIL module score\n(per-sample mean)', fontsize=10)
ax.tick_params(axis='y', labelsize=9)

ymax_data = max(np.max(d) if len(d) else 0 for d in data)
ymin_data = min(np.min(d) if len(d) else 0 for d in data)
yr = ymax_data - ymin_data
# Bracket levels above the boxes (staggered so spanning brackets don't collide).
# Normal NBM-MGUS (1-3) and neoplastic MGUS-SMM (4-7) have no horizontal overlap -> same level.
# NBM-SMM (1-6) overlaps both -> one level higher.
norm_lvl  = [ymax_data + 0.10 * yr,  # NBM-MGUS  (short, level 1)
             ymax_data + 0.26 * yr]  # NBM-SMM   (long, level 2; above MGUS-malignant)
neo_lvl   = [ymax_data + 0.10 * yr]  # MGUS-malignant (single comparison, level 1)
jt_norm_y = ymax_data + 0.52 * yr   # JT(normal) two-line label, lifted a bit higher
ax.set_ylim(ymin_data - 0.05 * yr, jt_norm_y + 0.18 * yr)

def fmt_q(q):
    if np.isnan(q): return "n/a"
    return f"q={q:.2f}" if q >= 0.01 else f"q={q:.1e}"

def draw_bracket(x1, x2, y, label, color, lw=1.0, fontsize=8.5):
    h = 0.013 * yr
    ax.plot([x1, x1, x2, x2], [y - h, y, y, y - h], color=color, linewidth=lw)
    ax.text((x1 + x2) / 2, y + 0.4 * h, label,
            ha='center', va='bottom', fontsize=fontsize, color=color)

# Cross-stage normal-PC brackets (dark blue; BH-corrected within the 2 normal comparisons)
for (lbl, ia, ib), y in zip(norm_cross, norm_lvl):
    draw_bracket(pos[ia], pos[ib], y, fmt_q(norm_cross_q[lbl]),
                 color=NORMAL_BRACKET_COLOR, lw=1.0, fontsize=8.5)

# Cross-stage neoplastic-PC bracket (dark red; single comparison, q == raw p)
for (lbl, ia, ib), y in zip(neo_cross, neo_lvl):
    draw_bracket(pos[ia], pos[ib], y, fmt_q(neo_cross_q[lbl]),
                 color=NEOPLASTIC_BRACKET_COLOR, lw=1.0, fontsize=8.5)

center_x = 4.7  # midpoint of the 5-box layout (pos 1..8.4). Two-line JT label centered here.
# JT normal PC line: prefix in normal weight + stat (z, p) in BOLD on a second line (mathtext).
ax.text(center_x, jt_norm_y,
        "Jonckheere-Terpstra trend, normal PCs (NBM → SMM)\n"
        + r"$\mathbf{" + f"z={z_n:.2f},\\ p={p_n:.3g}" + r"}$",
        ha='center', va='center', fontsize=9, fontweight='normal',
        multialignment='center')

ax.set_title('BM plasma cells (Boiarsky GSE193531)',
             fontsize=10.5, pad=10, fontweight='bold')

# Single-row legend below the x-axis, no frame; gives the panel a clean external label strip
leg = [mpatches.Patch(facecolor='lightgray', edgecolor='black', alpha=0.6, label='normal PCs'),
       mpatches.Patch(facecolor='lightgray', edgecolor='black', alpha=0.6, hatch='///', label='malignant PCs')]
ax.legend(handles=leg, loc='upper center', bbox_to_anchor=(0.5, -0.28),
          ncol=2, fontsize=8.5, frameon=False, handlelength=1.6, columnspacing=2.0)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(FIGURES_DIR / "Figure3F.png", dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(FIGURES_DIR / "Figure3F.svg", bbox_inches='tight', facecolor='white')
plt.close()
print("Saved: Figure3F.png + Figure3F.svg")
