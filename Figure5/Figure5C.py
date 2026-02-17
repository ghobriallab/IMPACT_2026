#!/usr/bin/env python3
"""Figure 5C: IL-1B response signature analysis (pre vs post-vx, overall + per-celltype)."""

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
import gc

# Cell types to exclude
EXCLUDE_CELLTYPES = ['QC_removed', 'MK', 'CLL']
EXCLUDE_PATTERNS = ['db:']

# Plot colors
COLORS_PREPOST = {'Pre-Vx': '#4878A8', 'Post-Vx': '#E05555'}

# IL-1B response gene list
IL1B_GENES_PATH = DATA_DIR / "il1b_response_genes_human.csv"


def filter_untreated_smm(df):
    """Keep HD, MGUS as-is; for SMM keep only untreated. Exclude MM."""
    def keep_sample(row):
        if row['Disease'] in ['HD', 'MGUS', 'PROMISE_NEG', 'IgM-MGUS']:
            return True
        elif row['Disease'] == 'SMM':
            return row['Treatment_Status'] == 0.0 or pd.isna(row['Treatment_Status'])
        return False
    return df[df.apply(keep_sample, axis=1)].copy()


def should_exclude_celltype(ct):
    if ct in EXCLUDE_CELLTYPES:
        return True
    for pattern in EXCLUDE_PATTERNS:
        if ct.startswith(pattern):
            return True
    return False


def compute_statistics(paired_df, disease_list):
    """Compute paired Wilcoxon statistics for each disease."""
    stats_results = []
    for disease in disease_list:
        subset = paired_df[paired_df['Disease'] == disease]
        n = len(subset)
        if n >= 3:
            pre_vals = subset['Pre-Vx']
            post_vals = subset['Post-Vx']
            try:
                stat, p_value = stats.wilcoxon(post_vals, pre_vals)
            except Exception:
                stat, p_value = np.nan, np.nan
            stats_results.append({
                'Disease': disease, 'n': n,
                'Mean_PreVx': pre_vals.mean(), 'Mean_PostVx': post_vals.mean(),
                'Mean_Diff': post_vals.mean() - pre_vals.mean(),
                'Wilcoxon_stat': stat, 'p_value': p_value
            })
        else:
            stats_results.append({
                'Disease': disease, 'n': n,
                'Mean_PreVx': np.nan, 'Mean_PostVx': np.nan,
                'Mean_Diff': np.nan, 'Wilcoxon_stat': np.nan, 'p_value': np.nan
            })
    return pd.DataFrame(stats_results)


def generate_boxplot(plot_df, stats_df, output_path, disease_order=None):
    """Generate boxplot with significance brackets."""
    if disease_order is None:
        disease_order = ['HD', 'MGUS', 'SMM']
    disease_order = [d for d in disease_order if d in plot_df['Disease'].unique()]

    fig, ax = plt.subplots(figsize=(4.5, 4.2))

    sns.boxplot(data=plot_df, x='Disease', y='IL1B_score', hue='Timepoint_Clean',
                order=disease_order, hue_order=['Pre-Vx', 'Post-Vx'],
                palette=COLORS_PREPOST, width=0.65, linewidth=1.5, ax=ax)
    sns.stripplot(data=plot_df, x='Disease', y='IL1B_score', hue='Timepoint_Clean',
                  order=disease_order, hue_order=['Pre-Vx', 'Post-Vx'],
                  palette=COLORS_PREPOST, dodge=True, alpha=0.6, size=5, ax=ax, legend=False)

    y_min = plot_df['IL1B_score'].min()
    y_max = plot_df['IL1B_score'].max()
    y_range = y_max - y_min
    bracket_y = y_max + 0.02 * y_range

    for i, disease in enumerate(disease_order):
        if disease not in stats_df['Disease'].values:
            continue
        row = stats_df[stats_df['Disease'] == disease].iloc[0]
        p_val = row['p_value']
        n = int(row['n'])

        if pd.isna(p_val):
            p_str, stars = "n/a", ""
        elif p_val < 0.001:
            p_str, stars = "p<0.001", "***"
        elif p_val < 0.01:
            p_str, stars = f"p={p_val:.3f}", "**"
        elif p_val < 0.05:
            p_str, stars = f"p={p_val:.2f}", "*"
        else:
            p_str, stars = f"p={p_val:.2f}", "ns"

        label = f"{stars}\n{p_str}"
        x_left, x_right = i - 0.17, i + 0.17
        bracket_height = 0.008 * y_range
        ax.plot([x_left, x_left, x_right, x_right],
                [bracket_y, bracket_y + bracket_height, bracket_y + bracket_height, bracket_y],
                'k-', linewidth=1.2)
        ax.text(i, bracket_y + bracket_height + 0.005 * y_range, label,
                ha='center', va='bottom', fontsize=11)
        ax.text(i, -0.12, f'n={n}', ha='center', va='top', fontsize=11,
                transform=ax.get_xaxis_transform())

    ax.set_xlabel('')
    ax.set_ylabel('IL-1\u03b2 Response Score', fontsize=13)
    ax.tick_params(axis='both', labelsize=12)
    ax.tick_params(axis='x', labelsize=11)
    ax.legend(title='', fontsize=11, frameon=False, loc='upper center',
              bbox_to_anchor=(0.5, 1.18), ncol=2, columnspacing=1.5)
    ax.set_ylim(y_min - 0.08 * y_range, bracket_y + 0.22 * y_range)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def main():
    # Load IL1B response genes
    il1b_genes = pd.read_csv(IL1B_GENES_PATH)['gene'].tolist()
    print(f"IL1B response genes: {len(il1b_genes)}")

    # Load metadata
    metadata = pd.read_csv(METADATA_CSV)
    metadata = metadata.rename(columns={
        'Demux_Sample_ID': 'sample_id',
        'VaccineTimepointFinal': 'Timepoint',
        'DiseaseStatus': 'Diagnosis',
        'TreatmentStatus': 'Treatment_Status'
    })

    # Load h5ad
    print("Loading h5ad...")
    adata = sc.read_h5ad(H5AD_ANNOTATED)
    print(f"Loaded {adata.n_obs:,} cells")

    # Normalize if counts layer exists
    if 'counts' in adata.layers:
        adata.X = adata.layers['counts'].copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    available_genes = [g for g in il1b_genes if g in adata.var_names]
    print(f"IL1B genes in data: {len(available_genes)}/{len(il1b_genes)}")

    if len(available_genes) < 50:
        print("ERROR: Too few genes found")
        return

    # Score IL1B response
    sc.tl.score_genes(adata, gene_list=available_genes, score_name='IL1B_response',
                      ctrl_size=min(100, len(available_genes)))

    # Prepare cell dataframe
    cell_df = pd.DataFrame({
        'sample_id': adata.obs['sample_id'].values,
        'PatientID': adata.obs['PatientID'].values,
        'Annotation_Level_2': adata.obs['Annotation_Level_2'].values,
        'IL1B_response': adata.obs['IL1B_response'].values
    })

    cell_df['sample_id'] = cell_df['sample_id'].astype(str)
    metadata['sample_id'] = metadata['sample_id'].astype(str)
    cell_df = cell_df.merge(
        metadata[['sample_id', 'Timepoint', 'Diagnosis', 'Treatment_Status']],
        on='sample_id', how='left'
    )
    cell_df = cell_df[cell_df['Timepoint'].notna()].copy()

    timepoint_map = {'Before_1st': 'Pre-Vx', 'After_2nd': 'Post-Vx', 'After_1st': 'Post-1st'}
    cell_df['Timepoint_Clean'] = cell_df['Timepoint'].map(timepoint_map).fillna(cell_df['Timepoint'])
    cell_df.rename(columns={'Diagnosis': 'Disease'}, inplace=True)
    cell_df['Disease'] = cell_df['Disease'].replace({'PROMISE_NEG': 'HD'})

    del adata
    gc.collect()

    # Overall analysis: aggregate per patient/timepoint
    overall_scores = cell_df.groupby(
        ['PatientID', 'Timepoint_Clean', 'Disease', 'Treatment_Status']
    )['IL1B_response'].mean().reset_index()
    overall_scores.columns = ['Patient_ID', 'Timepoint_Clean', 'Disease', 'Treatment_Status', 'IL1B_score']

    overall_filtered = filter_untreated_smm(overall_scores)
    analysis_df = overall_filtered[
        (overall_filtered['Disease'].isin(['HD', 'MGUS', 'SMM'])) &
        (overall_filtered['Timepoint_Clean'].isin(['Pre-Vx', 'Post-Vx']))
    ].copy()
    # SMM label (data already filtered to untreated only)

    # Paired analysis
    paired_df = analysis_df.pivot_table(
        index=['Patient_ID', 'Disease'], columns='Timepoint_Clean', values='IL1B_score'
    ).reset_index()
    paired_df = paired_df.dropna(subset=['Pre-Vx', 'Post-Vx'])
    print(f"Paired patients: {len(paired_df)}")

    stats_df = compute_statistics(paired_df, ['HD', 'MGUS', 'SMM'])

    for _, row in stats_df.iterrows():
        if pd.notna(row['p_value']):
            print(f"  {row['Disease']} (n={int(row['n'])}): p={row['p_value']:.4f}")

    # Generate overall plot
    plot_df = analysis_df[analysis_df['Timepoint_Clean'].isin(['Pre-Vx', 'Post-Vx'])].copy()
    generate_boxplot(plot_df, stats_df, FIGURES_DIR / "Figure5C.png")
    print("Saved: Figure5C.png")


if __name__ == "__main__":
    main()
