#!/usr/bin/env python3
"""Figure 3E: APRIL-responsive gene module score (sample-level boxplot)."""

import sys
import warnings
warnings.filterwarnings('ignore')
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import *

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# APRIL-responsive genes (literature-validated)
APRIL_TARGETS = {
    'Survival': ['MCL1', 'BCL2'],
    'NF-kB_pathway': ['NFKB2', 'RELB', 'NFKB1', 'NFKBIA'],
    'Adhesion': ['CD44', 'ICAM1'],
    'Immunomodulation': ['CCL3', 'CCL4', 'VEGFA', 'IL10', 'CD274', 'TGFB1', 'CXCL8']
}
ALL_APRIL_GENES = [g for genes in APRIL_TARGETS.values() for g in genes]

# Load data
adata = sc.read_h5ad(H5AD_CELLXGENE)
print(f"Total cells: {adata.n_obs:,}")

available_genes = [g for g in ALL_APRIL_GENES if g in adata.var_names]
print(f"APRIL genes available: {len(available_genes)}/{len(ALL_APRIL_GENES)}")

# Filter to B cells, post-vax, HD + untreated SMM
bcell_mask = adata.obs['Annotation_Level_1'] == 'B'
postvx_mask = adata.obs['Timepoint'] == 'Post-2nd'
hd_mask = adata.obs['Diagnosis'] == 'HD'
smm_untreated_mask = (adata.obs['Diagnosis'] == 'SMM') & (adata.obs['TreatmentStatus'] == 'Never_treated')
group_mask = hd_mask | smm_untreated_mask

adata_sub = adata[bcell_mask & postvx_mask & group_mask].copy()
adata_sub.obs['Group'] = adata_sub.obs['Diagnosis'].astype(str)
adata_sub.obs['Group'] = pd.Categorical(adata_sub.obs['Group'], categories=['HD', 'SMM'], ordered=True)
print(f"Subset: {adata_sub.n_obs:,} cells (HD={(adata_sub.obs['Group']=='HD').sum():,}, SMM={(adata_sub.obs['Group']=='SMM').sum():,})")

# Calculate module score
sc.tl.score_genes(adata_sub, gene_list=available_genes, score_name='APRIL_score')

# Aggregate per patient
sample_scores = adata_sub.obs.groupby('Deidentified_Patient_ID').agg({
    'APRIL_score': 'mean', 'Group': 'first'
}).reset_index()
sample_scores.columns = ['Patient_ID', 'APRIL_score', 'Group']
sample_scores = sample_scores.dropna(subset=['APRIL_score'])

n_hd = (sample_scores['Group'] == 'HD').sum()
n_smm = (sample_scores['Group'] == 'SMM').sum()
print(f"Patients: HD={n_hd}, SMM={n_smm}")

# Sample-level statistics
hd_sample_scores = sample_scores[sample_scores['Group'] == 'HD']['APRIL_score'].values
smm_sample_scores = sample_scores[sample_scores['Group'] == 'SMM']['APRIL_score'].values
_, pval_sample = stats.mannwhitneyu(hd_sample_scores, smm_sample_scores, alternative='two-sided')
print(f"Mann-Whitney p = {pval_sample:.4e}")

# Plot (sample-level boxplot)
fig, ax = plt.subplots(figsize=(4.375, 4))

sample_plot = sample_scores.copy()
sample_plot['Group'] = pd.Categorical(sample_plot['Group'], categories=['HD', 'SMM'], ordered=True)

sns.boxplot(data=sample_plot, x='Group', y='APRIL_score',
            hue='Group', palette={'HD': '#4DBBD5', 'SMM': '#E64B35'}, ax=ax,
            width=0.6, linewidth=1.5)
sns.stripplot(data=sample_plot, x='Group', y='APRIL_score',
              color='black', alpha=0.6, size=6, ax=ax)

ax.set_title('APRIL-Responsive Genes', fontsize=14, fontweight='bold', pad=18)
ax.text(0.5, 1.01, 'Expression Score (B Cells) - PBMCs',
        transform=ax.transAxes, fontsize=10, ha='center', va='bottom')
ax.set_xlabel('')
ax.set_ylabel('Module Score', fontsize=12)
ax.tick_params(axis='both', labelsize=11)
ax.set_xticklabels([f'HD\n(n={n_hd})', f'SMM\n(n={n_smm})'], fontsize=10)

if ax.get_legend():
    ax.get_legend().remove()

# p-value bracket
ymax = sample_plot['APRIL_score'].max()
ymin = sample_plot['APRIL_score'].min()
yrange = ymax - ymin
bracket_y = ymax + 0.10 * yrange
text_y = bracket_y + 0.03 * yrange

ax.plot([0, 0, 1, 1], [bracket_y - 0.02*yrange, bracket_y, bracket_y, bracket_y - 0.02*yrange],
        'k-', linewidth=1)
pval_str = f"p = {pval_sample:.2e}" if pval_sample < 0.001 else f"p = {pval_sample:.3f}"
ax.text(0.5, text_y, pval_str, ha='center', va='bottom', fontsize=10)

ax.set_ylim(ymin - 0.1*yrange, text_y + 0.20*yrange)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(FIGURES_DIR / "Figure3E.png", dpi=300, bbox_inches='tight')
plt.close()
print("Saved: Figure3E.png")
