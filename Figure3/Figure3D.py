#!/usr/bin/env python3
"""Figure 3D: APRIL expression in myeloid cells (Pre-Vx vs Post-Vx boxplots)."""

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
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

APRIL_GENE = 'TNFSF13'
DISEASE_GROUPS = ['HD', 'MGUS', 'SMM']

def main():
    # Load metadata
    meta = pd.read_csv(METADATA_CSV)
    meta = meta.rename(columns={'Demux_Sample_ID': 'sample_id'})

    def get_diagnosis(row):
        disease = row['DiseaseStatus']
        treatment = row['TreatmentStatus']
        if disease == 'PROMISE_NEG':
            return 'HD'
        elif disease == 'IgM-MGUS':
            return 'MGUS'
        elif disease == 'SMM':
            if pd.isna(treatment) or treatment == 0:
                return 'SMM'
            else:
                return 'SMM_treated'
        else:
            return disease

    meta['DiagnosisBinary'] = meta.apply(get_diagnosis, axis=1)

    sample_to_diag = dict(zip(meta['sample_id'], meta['DiagnosisBinary']))
    sample_to_tp = dict(zip(meta['sample_id'], meta['VaccineTimepointFinal']))
    sample_to_patient = dict(zip(meta['sample_id'], meta['PatientID']))

    # Load annotated h5ad for cell type labels
    adata_annot = sc.read_h5ad(H5AD_ANNOTATED, backed='r')

    annot_df = pd.DataFrame({
        'Annotation_Level_1': adata_annot.obs['Annotation_Level_1'].values,
        'sample_id': adata_annot.obs['sample_id'].values
    }, index=adata_annot.obs_names)

    annot_df['diagnosis'] = annot_df['sample_id'].map(sample_to_diag)
    annot_df['timepoint'] = annot_df['sample_id'].map(sample_to_tp)
    annot_df['patient'] = annot_df['sample_id'].map(sample_to_patient)

    # Filter to myeloid cells (Mono + DC), HD/MGUS/SMM, Pre/Post-Vx
    myeloid_df = annot_df[
        (annot_df['Annotation_Level_1'].isin(['Mono', 'DC'])) &
        (annot_df['timepoint'].isin(['Before_1st', 'After_2nd'])) &
        (annot_df['diagnosis'].isin(DISEASE_GROUPS))
    ].copy()

    myeloid_df['Timepoint'] = myeloid_df['timepoint'].map({
        'Before_1st': 'Pre-Vx', 'After_2nd': 'Post-Vx'
    })
    print(f"Myeloid cells: {len(myeloid_df):,}")

    # Load full h5ad for APRIL expression
    adata_full = sc.read_h5ad(H5AD_NORM, backed='r')
    april_idx = list(adata_full.var_names).index(APRIL_GENE)

    # Extract expression in chunks
    full_bc_to_idx = {bc: i for i, bc in enumerate(adata_full.obs_names)}
    barcodes = myeloid_df.index.tolist()
    indices = [full_bc_to_idx[bc] for bc in barcodes if bc in full_bc_to_idx]
    matched = [bc for bc in barcodes if bc in full_bc_to_idx]

    chunk_size = 100000
    expr = []
    for i in range(0, len(indices), chunk_size):
        chunk = adata_full.X[indices[i:i+chunk_size], april_idx].toarray().flatten()
        expr.extend(chunk)

    myeloid_df = myeloid_df.loc[matched].copy()
    myeloid_df['APRIL_expr'] = expr

    # Aggregate by patient
    patient_data = myeloid_df.groupby(['patient', 'diagnosis', 'Timepoint']).agg({
        'APRIL_expr': 'mean'
    }).reset_index()

    # Statistical tests: HD vs MGUS, HD vs SMM within each timepoint
    stats_results = []
    pvals_all = []

    for tp in ['Pre-Vx', 'Post-Vx']:
        tp_data = patient_data[patient_data['Timepoint'] == tp]
        hd_vals = tp_data[tp_data['diagnosis'] == 'HD']['APRIL_expr'].values
        mgus_vals = tp_data[tp_data['diagnosis'] == 'MGUS']['APRIL_expr'].values
        smm_vals = tp_data[tp_data['diagnosis'] == 'SMM']['APRIL_expr'].values

        for comp_name, comp_vals in [('HD vs MGUS', mgus_vals), ('HD vs SMM', smm_vals)]:
            stat, pval = stats.mannwhitneyu(hd_vals, comp_vals, alternative='two-sided')
            pvals_all.append(pval)
            stats_results.append({
                'Timepoint': tp, 'Comparison': comp_name,
                'n1': len(hd_vals), 'n2': len(comp_vals),
                'mean1': hd_vals.mean(), 'mean2': comp_vals.mean(), 'pval': pval
            })
        print(f"{tp}: HD vs MGUS p={stats_results[-2]['pval']:.4f}, HD vs SMM p={stats_results[-1]['pval']:.4f}")

    # FDR correction
    _, qvals_corrected, _, _ = multipletests(pvals_all, method='fdr_bh')
    for i, qval in enumerate(qvals_corrected):
        stats_results[i]['qval'] = qval

    stats_df = pd.DataFrame(stats_results)

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

        # q-value brackets
        row_mgus = stats_df[(stats_df['Timepoint'] == tp) & (stats_df['Comparison'] == 'HD vs MGUS')].iloc[0]
        qval_mgus = row_mgus['qval']
        qtext_mgus = 'NS' if qval_mgus > 0.05 else (f'q={qval_mgus:.2f}' if qval_mgus >= 0.01 else f'q={qval_mgus:.3f}')
        ax.plot([0, 1], [ymax + 0.08, ymax + 0.08], 'k-', linewidth=1)
        ax.text(0.5, ymax + 0.1, qtext_mgus, ha='center', va='bottom', fontsize=11)

        row_smm = stats_df[(stats_df['Timepoint'] == tp) & (stats_df['Comparison'] == 'HD vs SMM')].iloc[0]
        qval_smm = row_smm['qval']
        qtext_smm = 'NS' if qval_smm > 0.05 else (f'q={qval_smm:.2f}' if qval_smm >= 0.01 else f'q={qval_smm:.3f}')
        ax.plot([0, 2], [ymax + 0.22, ymax + 0.22], 'k-', linewidth=1)
        ax.text(1, ymax + 0.24, qtext_smm, ha='center', va='bottom', fontsize=11)

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
