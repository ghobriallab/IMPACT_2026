#!/usr/bin/env python3
"""Figure 3F: External validation of APRIL-responsive genes (GSE193531, normal PCs only)."""

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
import matplotlib.pyplot as plt
import seaborn as sns

np.random.seed(42)

# APRIL-responsive genes (same as IMPACT analysis)
APRIL_TARGETS = {
    'Survival': ['MCL1', 'BCL2'],
    'NF-kB_pathway': ['NFKB2', 'RELB', 'NFKB1', 'NFKBIA'],
    'Adhesion': ['CD44', 'ICAM1'],
    'Immunomodulation': ['CCL3', 'CCL4', 'VEGFA', 'IL10', 'CD274', 'TGFB1', 'CXCL8']
}
ALL_APRIL_GENES = [g for genes in APRIL_TARGETS.values() for g in genes]

EXT_DATA_DIR = DATA_DIR / "external"

# Load data
print("Loading GSE193531 data...")
metadata = pd.read_csv(EXT_DATA_DIR / "GSE193531_cell-level-metadata.csv", index_col=0)
print(f"Metadata: {len(metadata):,} cells")

counts = pd.read_csv(EXT_DATA_DIR / "GSE193531_umi-count-matrix.csv.gz", index_col=0)
print(f"Count matrix: {counts.shape[0]:,} genes x {counts.shape[1]:,} cells")

# Create AnnData
counts_T = counts.T
common_cells = counts_T.index.intersection(metadata.index)
counts_T = counts_T.loc[common_cells]
metadata = metadata.loc[common_cells]

adata = sc.AnnData(X=csr_matrix(counts_T.values), obs=metadata,
                   var=pd.DataFrame(index=counts_T.columns))
print(f"AnnData: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

# Normalize
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Check gene availability
available_genes = [g for g in ALL_APRIL_GENES if g in adata.var_names]
print(f"APRIL genes available: {len(available_genes)}/{len(ALL_APRIL_GENES)}")

# Score APRIL targets
sc.tl.score_genes(adata, gene_list=available_genes, score_name='APRIL_score', random_state=42)

# Filter to NBM vs SMM (normal PCs only via inferCNV classification)
adata_sub = adata[adata.obs['disease_stage'].isin(['NBM', 'SMM'])].copy()
smm_neoplastic_mask = ((adata_sub.obs['disease_stage'] == 'SMM') &
                       (adata_sub.obs['normal_or_neoplastic'] != 'normal'))
adata_sub = adata_sub[~smm_neoplastic_mask].copy()

print(f"NBM cells: {(adata_sub.obs['disease_stage'] == 'NBM').sum():,}")
print(f"SMM normal PCs: {(adata_sub.obs['disease_stage'] == 'SMM').sum():,}")

# Aggregate to sample level
sample_scores = adata_sub.obs.groupby(['sample_ID', 'disease_stage'])['APRIL_score'].mean().reset_index()
sample_scores.columns = ['Sample', 'Disease', 'APRIL_score']

nbm_sample = sample_scores[sample_scores['Disease'] == 'NBM']['APRIL_score']
smm_sample = sample_scores[sample_scores['Disease'] == 'SMM']['APRIL_score']
_, pval_sample = stats.mannwhitneyu(nbm_sample, smm_sample, alternative='two-sided')

n_nbm = len(nbm_sample)
n_smm = len(smm_sample)
print(f"Sample-level: NBM n={n_nbm}, SMM n={n_smm}, p={pval_sample:.4f}")

# Plot (sample-level boxplot)
fig, ax = plt.subplots(figsize=(4.375, 4.5))

sample_plot = sample_scores.copy()
sample_plot['Disease'] = pd.Categorical(sample_plot['Disease'], categories=['NBM', 'SMM'], ordered=True)

colors_module = {'NBM': '#4DBBD5', 'SMM': '#E64B35'}

sns.boxplot(data=sample_plot, x='Disease', y='APRIL_score',
            hue='Disease', palette=colors_module, ax=ax, width=0.6, linewidth=1.5)
sns.stripplot(data=sample_plot, x='Disease', y='APRIL_score',
              color='black', alpha=0.6, size=6, ax=ax)

ax.set_title('APRIL-Responsive Genes', fontsize=14, fontweight='bold', pad=32)
ax.text(0.5, 1.065, 'Expression Score (Non-Malignant Plasma Cells) - BM',
        transform=ax.transAxes, fontsize=10, ha='center', va='bottom')
ax.text(0.5, 1.01, 'External validation dataset: GSE193531',
        transform=ax.transAxes, fontsize=9, ha='center', va='bottom', style='italic')

ax.set_xlabel('')
ax.set_ylabel('Module Score', fontsize=12)
ax.tick_params(axis='both', labelsize=11)
ax.set_xticklabels([f'NBM\n(n={n_nbm})', f'SMM\n(n={n_smm})'], fontsize=10)

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
plt.savefig(FIGURES_DIR / "Figure3F.png", dpi=300, bbox_inches='tight')
plt.close()
print("Saved: Figure3F.png")
