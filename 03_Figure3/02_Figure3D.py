#!/usr/bin/env python3
"""Figure 3D (= manuscript Figure 3C): APRIL expression in myeloid cells (Pre-Vx vs Post-Vx boxplots).

Purpose:      Figure 3D (= manuscript Figure 3C): TNFSF13 (APRIL) expression in peripheral myeloid cells across HD/MGUS/SMM, pre and post vaccination, with age- and sex-adjusted rank-based ANCOVA + BH correction per timepoint.

Inputs:       H5AD_ANNOTATED (filtered to clean cells: drops QC_removed, doublets, Platelets, CLL).

Outputs:      figures/Figure3D.png and the per-timepoint adjusted q-values.

Dependencies: Python + scanpy, numpy, pandas, scipy, statsmodels, matplotlib; reads config.py.
"""
import sys
import warnings
warnings.filterwarnings('ignore')
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import *

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

APRIL_GENE = 'TNFSF13'
DISEASE_GROUPS = ['HD', 'MGUS', 'SMM']

def main():
    # consolidated single de-identified shared object. Read the per-cell
    # group labels straight from obs (Diagnosis, Timepoint, TreatmentStatus, Deidentified_Patient_ID)
    # instead of joining the original-ID metadata table. Untreated SMM only, matching the prior
    # SMM_treated exclusion. APRIL (TNFSF13) expression comes from the same file's normalized X.
    adata = sc.read_h5ad(H5AD_NORM, backed='r')
    april_idx = list(adata.var_names).index(APRIL_GENE)

    cells = pd.DataFrame({
        'Annotation_Level_1': adata.obs['Annotation_Level_1'].astype(str).values,
        'l2': adata.obs['Annotation_Level_2'].astype(str).values,
        'diagnosis': adata.obs['Diagnosis'].astype(str).values,
        'tp_raw': adata.obs['Timepoint'].astype(str).values,
        'treat': adata.obs['TreatmentStatus'].astype(str).values,
        'patient': adata.obs['Deidentified_Patient_ID'].astype(str).values,
    })
    cells['pos'] = np.arange(len(cells))

    tp_map = {'Pre-Vx': 'Pre-Vx', 'Post-2nd': 'Post-Vx'}  # obs timepoints -> figure timepoints
    # exclude doublet/QC cells retained (labeled) in the comprehensive deposit object.
    clean = (~cells['l2'].isin(['QC_removed', 'CLL'])) & (~cells['l2'].str.startswith('db:'))
    keep = (
        clean &
        cells['Annotation_Level_1'].isin(['Mono', 'DC']) &
        cells['tp_raw'].isin(tp_map.keys()) &
        cells['diagnosis'].isin(DISEASE_GROUPS)
    )
    # exclude treated SMM (keep untreated SMM only); HD/MGUS unaffected
    keep &= ~((cells['diagnosis'] == 'SMM') & (cells['treat'] != 'Never_treated'))
    myeloid_df = cells[keep].copy()
    myeloid_df['Timepoint'] = myeloid_df['tp_raw'].map(tp_map)
    print(f"Myeloid cells: {len(myeloid_df):,}")

    # APRIL (TNFSF13) expression from the same file, chunked (positions are sorted ascending)
    positions = myeloid_df['pos'].to_numpy()
    chunk_size = 100000
    expr = []
    for i in range(0, len(positions), chunk_size):
        sl = list(positions[i:i+chunk_size])
        chunk = adata.X[sl, april_idx]
        chunk = chunk.toarray().flatten() if hasattr(chunk, 'toarray') else np.asarray(chunk).flatten()
        expr.extend(chunk)
    myeloid_df['APRIL_expr'] = expr

    # Aggregate by patient
    patient_data = myeloid_df.groupby(['patient', 'diagnosis', 'Timepoint']).agg({
        'APRIL_expr': 'mean'
    }).reset_index()

    # join per-patient Age + Sex from Supplementary Table 2 (scRNA cohort metadata) for the
    # covariate-adjusted analysis. Supp Table 2 has 100% coverage of Age/Sex for the 118 scRNA-cohort
    # patients used here. No PHI is loaded: Supp Table 2 carries only the deid Patient ID + capped Age.
    META = REPO_DIR.parent / "Supplementary_Tables" / "Supplementary_Table_2_scRNAseq_sample_list.csv"
    _meta = pd.read_csv(META)
    _as = (_meta.groupby('Deidentified_Patient_ID')
                 .agg(Age=('Age', 'median'),
                      Sex=('Sex', lambda s: s.dropna().mode().iloc[0] if s.dropna().size else None))
                 .reset_index().rename(columns={'Deidentified_Patient_ID': 'patient'}))
    patient_data = patient_data.merge(_as, on='patient', how='left')

    # Statistical tests per timepoint: rank(expr) ~ Diagnosis + Age + Sex (rank-based ANCOVA),
    # BH-corrected across the 4 disease-vs-HD contrasts. (The unadjusted Mann-Whitney U on raw means
    # is also computed and printed for comparison; the figure brackets show the adjusted q-values.)
    stats_results = []
    pvals_all = []

    for tp in ['Pre-Vx', 'Post-Vx']:
        tp_data = patient_data[patient_data['Timepoint'] == tp].dropna(subset=['Age', 'Sex']).copy()
        tp_data['rk'] = tp_data['APRIL_expr'].rank()
        tp_data['diagnosis'] = pd.Categorical(tp_data['diagnosis'], categories=DISEASE_GROUPS)
        m = smf.ols('rk ~ C(diagnosis, Treatment(reference="HD")) + Age + C(Sex)', data=tp_data).fit()
        # raw MWU for comparison logging
        hd_vals = tp_data[tp_data['diagnosis'] == 'HD']['APRIL_expr'].values
        for comp_name, group in [('HD vs MGUS', 'MGUS'), ('HD vs SMM', 'SMM')]:
            comp_vals = tp_data[tp_data['diagnosis'] == group]['APRIL_expr'].values
            term = f'C(diagnosis, Treatment(reference="HD"))[T.{group}]'
            p_adj = float(m.pvalues[term])
            p_raw = float(stats.mannwhitneyu(hd_vals, comp_vals, alternative='two-sided').pvalue) if len(comp_vals) else np.nan
            pvals_all.append(p_adj)
            stats_results.append({
                'Timepoint': tp, 'Comparison': comp_name,
                'n1': len(hd_vals), 'n2': len(comp_vals),
                'mean1': hd_vals.mean(), 'mean2': comp_vals.mean() if len(comp_vals) else np.nan,
                'pval_adj': p_adj, 'pval_raw': p_raw,
            })
        print(f"{tp}: HD vs MGUS p_adj={stats_results[-2]['pval_adj']:.4f} (raw {stats_results[-2]['pval_raw']:.4f}), "
              f"HD vs SMM p_adj={stats_results[-1]['pval_adj']:.4f} (raw {stats_results[-1]['pval_raw']:.4f})")

    # BH correction across all 4 adjusted contrasts
    _, qvals_corrected, _, _ = multipletests(pvals_all, method='fdr_bh')
    for i, qval in enumerate(qvals_corrected):
        stats_results[i]['qval'] = qval

    stats_df = pd.DataFrame(stats_results)
    print("\n=== age+sex-adjusted, BH-corrected q-values (displayed on figure) ===")
    for _, r in stats_df.iterrows():
        print(f"  {r['Timepoint']} {r['Comparison']}: p_adj={r['pval_adj']:.4g}, q_adj(BH)={r['qval']:.4g}")

    # Two-panel boxplot figure
    diag_colors = {'HD': '#3498db', 'MGUS': '#f1c40f', 'SMM': '#e74c3c'}
    fig, axes = plt.subplots(1, 2, figsize=(6.75, 3.2), sharey=True)
    plt.subplots_adjust(wspace=0.12)

    for ax_idx, (tp, ax) in enumerate(zip(['Pre-Vx', 'Post-Vx'], axes)):
        tp_data = patient_data[patient_data['Timepoint'] == tp]
        box_data = [tp_data[tp_data['diagnosis'] == d]['APRIL_expr'].values for d in DISEASE_GROUPS]

        bp = ax.boxplot(box_data, positions=range(len(DISEASE_GROUPS)), widths=0.65,
                        patch_artist=True, showfliers=False)
        for patch, diag in zip(bp['boxes'], DISEASE_GROUPS):
            patch.set_facecolor(diag_colors[diag])
            patch.set_alpha(0.7)
        for element in ['whiskers', 'caps', 'medians']:
            plt.setp(bp[element], color='black', linewidth=1.2)

        for i, (data, diag) in enumerate(zip(box_data, DISEASE_GROUPS)):
            x = np.random.normal(i, 0.1, size=len(data))
            ax.scatter(x, data, alpha=0.6, color=diag_colors[diag], edgecolor='white',
                       s=45, zorder=3, linewidth=0.5)

        ymax = max([d.max() for d in box_data if len(d) > 0])

        # Adjusted q-value brackets (always show the value; matches the published Fig 3 panel style
        # and is more informative than a binary NS marker)
        def _fmt_q(q):
            return f'q={q:.2f}' if q >= 0.01 else f'q={q:.3f}'

        row_mgus = stats_df[(stats_df['Timepoint'] == tp) & (stats_df['Comparison'] == 'HD vs MGUS')].iloc[0]
        ax.plot([0, 1], [ymax + 0.08, ymax + 0.08], 'k-', linewidth=1)
        ax.text(0.5, ymax + 0.1, _fmt_q(row_mgus['qval']), ha='center', va='bottom', fontsize=11)

        row_smm = stats_df[(stats_df['Timepoint'] == tp) & (stats_df['Comparison'] == 'HD vs SMM')].iloc[0]
        ax.plot([0, 2], [ymax + 0.22, ymax + 0.22], 'k-', linewidth=1)
        ax.text(1, ymax + 0.24, _fmt_q(row_smm['qval']), ha='center', va='bottom', fontsize=11)

        ax.set_title(tp, fontsize=14, fontweight='bold', pad=6)
        ax.set_xticks(range(len(DISEASE_GROUPS)))
        xlabels = [f'{d}\n(n={len(tp_data[tp_data["diagnosis"]==d])})' for d in DISEASE_GROUPS]
        ax.set_xticklabels(xlabels, fontsize=11)
        ax.tick_params(axis='y', labelsize=11)

        if ax_idx == 0:
            ax.set_ylabel('APRIL Expression\n[log$_2$(CPM+1)]', fontsize=12, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlim(-0.5, len(DISEASE_GROUPS) - 0.5)

    all_data = patient_data['APRIL_expr'].values
    for ax in axes:
        ax.set_ylim(0, all_data.max() + 0.45)

    plt.suptitle('APRIL Expression (Myeloid Cells)', fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'Figure3D.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print("Saved: Figure3D.png")


if __name__ == "__main__":
    main()
