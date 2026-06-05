#!/usr/bin/env python3
"""Figure 3E (= manuscript Figure 3D): APRIL-responsive gene module score in peripheral B cells
across the precursor disease continuum (HD -> MGUS -> SMM), restricted to treatment-naive
participants. Statistics: Jonckheere-Terpstra ordered trend (HD->SMM) and pairwise Mann-Whitney
HD-vs-X with Benjamini-Hochberg correction; age and sex adjusted in the prose via rank-based
ANCOVA (reported in the response letter / Methods). Output: per-patient mean module score per
group.

The script requires per-cell obs annotations (Annotation_Level_1/2, Timepoint, Diagnosis,
TreatmentStatus). Age/Sex/TreatmentStatus are joined per patient from Supplementary Table 2 for
the adjusted statistics quoted in the manuscript.

Purpose:      Figure 3E (= manuscript Figure 3D): APRIL-responsive gene module score in peripheral B cells across the treatment-naive HD/MGUS/SMM continuum; JT ordered-trend test + BH-adjusted pairwise Mann-Whitney; age and sex adjusted via rank-based ANCOVA.

Inputs:       H5AD_ANNOTATED (B-cell subset, treatment-naive); APRIL-responsive gene list curated from the literature.

Outputs:      figures/Figure3E.png + per-patient module-score table.

Dependencies: Python + scanpy, numpy, pandas, scipy, statsmodels, matplotlib; reads config.py.
"""
import sys
import warnings
warnings.filterwarnings('ignore')
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import *

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats
# Match mathtext bold to the surrounding sans-serif font (used to make z/p bold inline)
mpl.rcParams['mathtext.fontset'] = 'dejavusans'
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf

# APRIL-responsive module (same 15 genes as published Fig 3D / SupFig3)
APRIL_TARGETS = {
    'Survival': ['MCL1', 'BCL2'],
    'NF-kB_pathway': ['NFKB2', 'RELB', 'NFKB1', 'NFKBIA'],
    'Adhesion': ['CD44', 'ICAM1'],
    'Immunomodulation': ['CCL3', 'CCL4', 'VEGFA', 'IL10', 'CD274', 'TGFB1', 'CXCL8'],
}
ALL_APRIL_GENES = [g for genes in APRIL_TARGETS.values() for g in genes]

# Load + filter (clean B cells, Post-2nd; matches the published cohort definition)
adata = sc.read_h5ad(H5AD_CELLXGENE)
print(f"Total cells: {adata.n_obs:,}")
_l2 = adata.obs['Annotation_Level_2'].astype(str)
clean_mask = (~_l2.isin(['QC_removed', 'CLL'])) & (~_l2.str.startswith('db:'))
bcell_mask = adata.obs['Annotation_Level_1'] == 'B'
postvx_mask = adata.obs['Timepoint'] == 'Post-2nd'
sub = adata[bcell_mask & postvx_mask & clean_mask].copy()
print(f"B/Post-2nd/clean: {sub.n_obs:,} cells")

available = [g for g in ALL_APRIL_GENES if g in sub.var_names]
print(f"APRIL genes present: {len(available)}/{len(ALL_APRIL_GENES)}")
sc.tl.score_genes(sub, gene_list=available, score_name='APRIL_score', random_state=42)

# Aggregate per patient
pat = (sub.obs.groupby('Deidentified_Patient_ID')
            .agg(APRIL_score=('APRIL_score','mean'),
                 n_cells=('APRIL_score','size'),
                 Diagnosis=('Diagnosis','first'),
                 TreatmentStatus=('TreatmentStatus','first'))
            .reset_index().rename(columns={'Deidentified_Patient_ID':'Patient_ID'}))
pat['Diagnosis'] = pat['Diagnosis'].replace({'IgM-MGUS': 'MGUS'})

# Treatment-naive cohort: HD (by definition naive) + MGUS/SMM with TreatmentStatus == 'Never_treated'
nt = pat[(pat['Diagnosis'] == 'HD') |
         ((pat['Diagnosis'].isin(['MGUS', 'SMM'])) & (pat['TreatmentStatus'] == 'Never_treated'))].copy()
print("Treatment-naive cohort by group:", nt['Diagnosis'].value_counts().to_dict())

order = ['HD', 'MGUS', 'SMM']
present = [g for g in order if g in set(nt['Diagnosis'])]
data = [nt.loc[nt['Diagnosis'] == g, 'APRIL_score'].values for g in present]
counts = [len(d) for d in data]

# Statistics ----
# Jonckheere-Terpstra (ordered alternative across HD -> MGUS -> SMM)
def jonckheere_terpstra(groups_in_order):
    g = [np.asarray(x) for x in groups_in_order]
    k = len(g)
    J = 0.0
    for i in range(k):
        for j in range(i + 1, k):
            xi, xj = g[i], g[j]
            for v in xi:
                J += float(np.sum(v < xj)) + 0.5 * float(np.sum(v == xj))
    n = [len(x) for x in g]; N = sum(n)
    mu = (N * N - sum(ni * ni for ni in n)) / 4.0
    var = (N * N * (2 * N + 3) - sum(ni * ni * (2 * ni + 3) for ni in n)) / 72.0
    z = (J - mu) / np.sqrt(var)
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    return z, p

z_jt, p_jt = jonckheere_terpstra(data)
print(f"JT (HD -> {present[-1]}): z={z_jt:.3f}, p={p_jt:.4g}")

# Pairwise HD vs X (Mann-Whitney; BH across the 2 contrasts)
hd = nt.loc[nt['Diagnosis'] == 'HD', 'APRIL_score'].values
pairs = []
for g in ['MGUS', 'SMM']:
    if g in present:
        other = nt.loc[nt['Diagnosis'] == g, 'APRIL_score'].values
        pairs.append((g, stats.mannwhitneyu(hd, other, alternative='two-sided').pvalue))
pair_q = multipletests([p for _, p in pairs], method='fdr_bh')[1]
qmap = dict(zip([g for g, _ in pairs], pair_q))
for g, p in pairs:
    print(f"  HD vs {g}: MWU p={p:.4g}, q(BH)={qmap[g]:.4g}")

# Age + sex adjusted rank-ANCOVA (Supp Table 2). Quoted in Methods + response letter.
try:
    _meta = pd.read_csv(REPO_DIR.parent / "Supplementary_Tables" / "Supplementary_Table_2_scRNAseq_sample_list.csv")
    _as = (_meta.groupby("Deidentified_Patient_ID")
                .agg(Age=("Age", "median"),
                     Sex=("Sex", lambda s: s.dropna().mode().iloc[0] if s.dropna().size else None))
                .reset_index().rename(columns={"Deidentified_Patient_ID": "Patient_ID"}))
    fit = nt.merge(_as, on='Patient_ID', how='left').dropna(subset=['Age', 'Sex'])
    fit['rk'] = fit['APRIL_score'].rank()
    fit['Diagnosis'] = pd.Categorical(fit['Diagnosis'], categories=order)
    _m = smf.ols('rk ~ C(Diagnosis, Treatment(reference="HD")) + Age + C(Sex)', data=fit).fit()
    print("  rank-ANCOVA age+sex-adjusted (n=%d):" % len(fit))
    for term in _m.pvalues.index:
        if 'Diagnosis' in term:
            lab = term.replace('C(Diagnosis, Treatment(reference="HD"))[T.', '').replace(']', '')
            print(f"    {lab} vs HD: p={_m.pvalues[term]:.4g}")
except Exception as e:
    print(f"age+sex adjustment skipped: {e}")

# Plot ----
colors_map = {'HD': '#4DBBD5', 'MGUS': '#F39B7F', 'SMM': '#E64B35'}
pos = list(range(1, len(present) + 1))
tick_labels = [f"{g}\nn={n}" for g, n in zip(present, counts)]

fig, ax = plt.subplots(figsize=(5.0, 5.0))
bp = ax.boxplot(data, positions=pos, widths=0.6, patch_artist=True,
                medianprops=dict(color='black', linewidth=1.8))
for patch, g in zip(bp['boxes'], present):
    patch.set_facecolor(colors_map[g]); patch.set_alpha(0.55)
rng = np.random.default_rng(42)
for i, (vals, g) in enumerate(zip(data, present)):
    ax.scatter(rng.normal(pos[i], 0.08, size=len(vals)), vals, s=28,
               color=colors_map[g], edgecolor='black', linewidth=0.5, alpha=0.92)

ax.set_xticks(pos); ax.set_xticklabels(tick_labels, fontsize=11)
ax.set_ylabel('APRIL-responsive module score\n(per-patient mean, B cells)', fontsize=12)
ax.tick_params(axis='y', labelsize=11)

# Headroom for q-brackets + JT
ymax_data = max(np.max(d) for d in data)
ymin_data = min(np.min(d) for d in data)
yr = ymax_data - ymin_data
bracket1_y = ymax_data + 0.10 * yr   # HD vs MGUS
bracket2_y = ymax_data + 0.22 * yr   # HD vs SMM
jt_y       = ymax_data + 0.42 * yr   # two-line JT label centered here
tick_h = 0.018 * yr
ax.set_ylim(ymin_data - 0.05 * yr, jt_y + 0.18 * yr)

def fmt_q(q):
    return f"q={q:.2f}" if q >= 0.01 else f"q={q:.1e}"

def bracket(x1, x2, y, label):
    ax.plot([x1, x1, x2, x2], [y - tick_h, y, y, y - tick_h], color='black', linewidth=1.0)
    ax.text((x1 + x2) / 2, y + 0.4 * tick_h, label, ha='center', va='bottom', fontsize=11)

hd_pos = pos[present.index('HD')]
for grp, lvl in [('MGUS', bracket1_y), ('SMM', bracket2_y)]:
    if grp in present:
        bracket(hd_pos, pos[present.index(grp)], lvl, fmt_q(qmap[grp]))

# JT trend — two lines stacked. Prefix in normal weight; the z/p statistics in bold (significant
# trend). Uses mathtext $\mathbf{...}$ for inline bold on the second line so the prefix above stays
# normal weight; mathtext.fontset='dejavusans' matches the surrounding font.
ax.text(sum(pos) / len(pos), jt_y,
        f"Jonckheere-Terpstra trend (HD → {present[-1]})\n"
        + r"$\mathbf{" + f"z={z_jt:.2f},\\ p={p_jt:.3g}" + r"}$",
        ha='center', va='center', fontsize=11.5, fontweight='normal',
        multialignment='center')

ax.set_title('APRIL-responsive module in B cells', fontsize=12, pad=8, fontweight='bold')
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(FIGURES_DIR / "Figure3E.png", dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(FIGURES_DIR / "Figure3E.svg", bbox_inches='tight', facecolor='white')
plt.close()
print("Saved: Figure3E.png + Figure3E.svg")
