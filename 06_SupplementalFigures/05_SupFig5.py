#!/usr/bin/env python3
"""Supplementary Figure 5: Individual APRIL-responsive gene violin plots (HD vs SMM).

Purpose:      Supplementary Figure 5: individual APRIL-responsive gene violins (15 genes, grouped by functional module) in peripheral B cells, HD vs SMM, post-vaccination.

Inputs:       H5AD_ANNOTATED (B-cell subset post-vaccination); 15-gene APRIL-target list.

Outputs:      figures/SupFig5.png (per-gene violin grid).

Dependencies: Python + scanpy, pandas, numpy, matplotlib, seaborn, scipy, statsmodels; reads config.py.
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
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
from matplotlib.patches import Patch

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

# Filter to B cells, post-vax, HD + untreated SMM.
# exclude doublet/QC cells retained (labeled) in the comprehensive deposit object.
_l2 = adata.obs['Annotation_Level_2'].astype(str)
clean_mask = (~_l2.isin(['QC_removed', 'CLL'])) & (~_l2.str.startswith('db:'))
bcell_mask = adata.obs['Annotation_Level_1'] == 'B'
postvx_mask = adata.obs['Timepoint'] == 'Post-2nd'
hd_mask = adata.obs['Diagnosis'] == 'HD'
smm_untreated_mask = (adata.obs['Diagnosis'] == 'SMM') & (adata.obs['TreatmentStatus'] == 'Never_treated')
group_mask = hd_mask | smm_untreated_mask

adata_sub = adata[bcell_mask & postvx_mask & group_mask & clean_mask].copy()
adata_sub.obs['Group'] = adata_sub.obs['Diagnosis'].astype(str)
adata_sub.obs['Group'] = pd.Categorical(adata_sub.obs['Group'], categories=['HD', 'SMM'], ordered=True)

n_hd = (adata_sub.obs['Group'] == 'HD').sum()
n_smm = (adata_sub.obs['Group'] == 'SMM').sum()
print(f"Cells: HD={n_hd:,}, SMM={n_smm:,}")

# Individual gene analysis
gene_results = []
for gene in available_genes:
    if gene in adata_sub.var_names:
        gene_idx = adata_sub.var_names.get_loc(gene)
        expr = adata_sub.X[:, gene_idx]
        if hasattr(expr, 'toarray'):
            expr = expr.toarray().flatten()
        else:
            expr = np.array(expr).flatten()

        hd_expr = expr[adata_sub.obs['Group'].values == 'HD']
        smm_expr = expr[adata_sub.obs['Group'].values == 'SMM']
        stat, pval = stats.mannwhitneyu(hd_expr, smm_expr, alternative='two-sided')

        gene_results.append({
            'Gene': gene,
            'HD_mean': np.mean(hd_expr), 'SMM_mean': np.mean(smm_expr),
            'Log2FC_SMM_vs_HD': np.log2((np.mean(smm_expr) + 0.01) / (np.mean(hd_expr) + 0.01)),
            'pvalue': pval
        })

df_genes = pd.DataFrame(gene_results)
_, df_genes['padj'], _, _ = multipletests(df_genes['pvalue'], method='fdr_bh')
df_genes = df_genes.sort_values('pvalue')
_ndown = int(((df_genes['padj'] < 0.05) & (df_genes['SMM_mean'] < df_genes['HD_mean'])).sum())
print(f"SupFig3: {_ndown}/{len(df_genes)} genes significantly downregulated in SMM (q<0.05)")

# Build plot data
plot_data = []
for gene in available_genes:
    if gene in adata_sub.var_names:
        gene_idx = adata_sub.var_names.get_loc(gene)
        expr = adata_sub.X[:, gene_idx]
        if hasattr(expr, 'toarray'):
            expr = expr.toarray().flatten()
        else:
            expr = np.array(expr).flatten()
        for group, e in zip(adata_sub.obs['Group'].values, expr):
            plot_data.append({'Gene': gene, 'Expression': e, 'Group': group})

df_plot = pd.DataFrame(plot_data)
df_plot['Group'] = pd.Categorical(df_plot['Group'], categories=['HD', 'SMM'], ordered=True)

pval_dict = dict(zip(df_genes['Gene'], df_genes['pvalue']))
padj_dict = dict(zip(df_genes['Gene'], df_genes['padj']))
available_sorted = df_genes['Gene'].tolist()

# Violin plots
n_genes = len(available_sorted)
n_cols = 4
n_rows = int(np.ceil(n_genes / n_cols))

fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.5*n_cols, 4*n_rows))
axes = axes.flatten()

colors = {'HD': '#4DBBD5', 'SMM': '#E64B35'}

for i, gene in enumerate(available_sorted):
    ax = axes[i]
    gene_data = df_plot[df_plot['Gene'] == gene]

    sns.violinplot(data=gene_data, x='Group', y='Expression', hue='Group',
                   palette=colors, ax=ax, inner='quartile', linewidth=1.5,
                   cut=0, density_norm='width')

    padj = padj_dict.get(gene, 1)
    gene_row = df_genes[df_genes['Gene'] == gene]
    if len(gene_row) > 0:
        direction = "HD > SMM" if gene_row['HD_mean'].values[0] > gene_row['SMM_mean'].values[0] else "SMM > HD"
    else:
        direction = ""

    if padj < 0.001:
        qval_str = f"q = {padj:.1e}\n{direction}"
        sig = "***"
    elif padj < 0.01:
        qval_str = f"q = {padj:.3f}\n{direction}"
        sig = "**"
    elif padj < 0.05:
        qval_str = f"q = {padj:.3f}\n{direction}"
        sig = "*"
    else:
        qval_str = f"q = {padj:.3f}"
        sig = "ns"

    ax.set_title(f"{gene} ({sig})", fontsize=16, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('Log1p Expression', fontsize=14)
    ax.tick_params(axis='x', labelsize=13)
    ax.tick_params(axis='y', labelsize=12)
    ax.text(0.5, 0.88, qval_str, transform=ax.transAxes,
            ha='center', va='top', fontsize=13, color='#333333', linespacing=1.2)
    if ax.get_legend():
        ax.get_legend().remove()

for j in range(i+1, len(axes)):
    fig.delaxes(axes[j])

fig.suptitle('APRIL-Responsive Gene Expression (Log1p Normalized)\nHD vs SMM (Post-Vaccination B Cells)',
             fontsize=18, fontweight='bold', y=1.01)

legend_elements = [Patch(facecolor=colors['HD'], label=f'HD (n={n_hd:,})'),
                   Patch(facecolor=colors['SMM'], label=f'SMM (n={n_smm:,})')]
fig.legend(handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, -0.02),
           ncol=2, fontsize=16, frameon=False)

plt.tight_layout()
plt.subplots_adjust(top=0.92, bottom=0.06)
plt.savefig(FIGURES_DIR / "SupFig5.png", dpi=300, bbox_inches='tight')
plt.close()
print("Saved: SupFig5.png")
