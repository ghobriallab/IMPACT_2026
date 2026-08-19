#!/usr/bin/env python3
"""Figure 3C (manuscript Figure 3C): APRIL expression in myeloid cells (Pre-Vx vs Post-Vx boxplots).

Purpose:      Figure 3C (manuscript Figure 3C): TNFSF13 (APRIL) expression in peripheral myeloid cells across HD/MGUS/SMM, pre and post vaccination, with age- and sex-adjusted rank-based ANCOVA + BH correction per timepoint.

Inputs:       H5AD_ANNOTATED (filtered to clean cells: drops QC_removed, doublets, Platelets, CLL).

Outputs:      figures/Figure3C.png/.pdf/.svg and the per-timepoint adjusted q-values.

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


APRIL_GENE = 'TNFSF13'
DISEASE_GROUPS = ['HD', 'MGUS', 'SMM']   # raw Diagnosis values used to filter the deposit
# SMM is displayed as two columns, untreated and previously/actively treated, so that
# disease stage is separated from prior therapy (matching Figure 1B/1C/1F). Contrasts remain each
# group versus HD; BH now runs across 6 contrasts (3 per timepoint) rather than 4.
PLOT_GROUPS = ['HD', 'MGUS', 'SMM (Untreated)', 'SMM (Treated)']

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
    myeloid_df = cells[keep].copy()
    myeloid_df['Timepoint'] = myeloid_df['tp_raw'].map(tp_map)
    # Treated SMM are retained as their own group.
    myeloid_df['group'] = np.where(
        myeloid_df['diagnosis'] != 'SMM', myeloid_df['diagnosis'],
        np.where(myeloid_df['treat'] == 'Never_treated', 'SMM (Untreated)', 'SMM (Treated)'))
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
    patient_data = myeloid_df.groupby(['patient', 'group', 'Timepoint']).agg({
        'APRIL_expr': 'mean'
    }).reset_index()

    # join per-patient Age + Sex from Supplementary Table 4 (scRNA cohort metadata) for the
    # covariate-adjusted analysis. Supp Table 4 has 100% coverage of Age/Sex for the 118 scRNA-cohort
    # patients used here. No PHI is loaded: Supp Table 4 carries only the deid Patient ID + capped Age.
    # Shipped in the Zenodo deposit as metadata/Supplementary_Table_4_scRNAseq_sample_list.csv.
    # The repo-relative fallback is the in-house working copy.
    META = DATA_DIR / "metadata" / "Supplementary_Table_4_scRNAseq_sample_list.csv"
    if not META.exists():
        META = REPO_DIR.parent / "Supplementary_Tables" / "Supplementary_Table_4_scRNAseq_sample_list.csv"
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
        tp_data['group'] = pd.Categorical(tp_data['group'], categories=PLOT_GROUPS)
        m = smf.ols('rk ~ C(group, Treatment(reference="HD")) + Age + C(Sex)', data=tp_data).fit()
        # raw MWU for comparison logging
        hd_vals = tp_data[tp_data['group'] == 'HD']['APRIL_expr'].values
        for group in PLOT_GROUPS[1:]:
            comp_name = f'HD vs {group}'
            comp_vals = tp_data[tp_data['group'] == group]['APRIL_expr'].values
            term = f'C(group, Treatment(reference="HD"))[T.{group}]'
            p_adj = float(m.pvalues[term])
            p_raw = float(stats.mannwhitneyu(hd_vals, comp_vals, alternative='two-sided').pvalue) if len(comp_vals) else np.nan
            pvals_all.append(p_adj)
            stats_results.append({
                'Timepoint': tp, 'Comparison': comp_name,
                'n1': len(hd_vals), 'n2': len(comp_vals),
                'mean1': hd_vals.mean(), 'mean2': comp_vals.mean() if len(comp_vals) else np.nan,
                'pval_adj': p_adj, 'pval_raw': p_raw,
            })
        for _r in stats_results[-(len(PLOT_GROUPS) - 1):]:
            print(f"{tp}: {_r['Comparison']} p_adj={_r['pval_adj']:.4f} (raw {_r['pval_raw']:.4f})")

    # BH correction across all 4 adjusted contrasts
    _, qvals_corrected, _, _ = multipletests(pvals_all, method='fdr_bh')
    for i, qval in enumerate(qvals_corrected):
        stats_results[i]['qval'] = qval

    stats_df = pd.DataFrame(stats_results)
    print("\n=== age+sex-adjusted, BH-corrected q-values (displayed on figure) ===")
    for _, r in stats_df.iterrows():
        print(f"  {r['Timepoint']} {r['Comparison']}: p_adj={r['pval_adj']:.4g}, q_adj(BH)={r['qval']:.4g}")

    # Two-panel boxplot figure
    # Treated SMM in a darker red than untreated SMM.
    diag_colors = {'HD': '#3498db', 'MGUS': '#f1c40f',
                   'SMM (Untreated)': '#e74c3c', 'SMM (Treated)': '#A93226'}
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
    plt.subplots_adjust(wspace=0.12)

    for ax_idx, (tp, ax) in enumerate(zip(['Pre-Vx', 'Post-Vx'], axes)):
        tp_data = patient_data[patient_data['Timepoint'] == tp]
        box_data = [tp_data[tp_data['group'] == d]['APRIL_expr'].values for d in PLOT_GROUPS]

        bp = ax.boxplot(box_data, positions=range(len(PLOT_GROUPS)), widths=0.65,
                        patch_artist=True, showfliers=False)
        for patch, diag in zip(bp['boxes'], PLOT_GROUPS):
            patch.set_facecolor(diag_colors[diag])
            patch.set_alpha(0.7)
        for element in ['whiskers', 'caps', 'medians']:
            plt.setp(bp[element], color='black', linewidth=1.2)

        for i, (data, diag) in enumerate(zip(box_data, PLOT_GROUPS)):
            x = np.random.normal(i, 0.1, size=len(data))
            ax.scatter(x, data, alpha=0.6, color=diag_colors[diag], edgecolor='white',
                       s=45, zorder=3, linewidth=0.5)

        ymax = max([d.max() for d in box_data if len(d) > 0])

        # Adjusted q-value brackets (always show the value; matches the published Fig 3 panel style
        # and is more informative than a binary NS marker)
        def _fmt_q(q):
            return f'q={q:.2f}' if q >= 0.01 else f'q={q:.3f}'

        for _j, _grp in enumerate(PLOT_GROUPS[1:], start=1):
            _row = stats_df[(stats_df['Timepoint'] == tp) &
                            (stats_df['Comparison'] == f'HD vs {_grp}')].iloc[0]
            _y = ymax + 0.08 + 0.14 * (_j - 1)
            ax.plot([0, _j], [_y, _y], 'k-', linewidth=1)
            ax.text(_j / 2, _y + 0.02, _fmt_q(_row['qval']), ha='center', va='bottom', fontsize=10)

        ax.set_title(tp, fontsize=14, fontweight='bold', pad=6)
        ax.set_xticks(range(len(PLOT_GROUPS)))
        # Stack the two SMM labels so four groups fit the original panel width without colliding.
        _disp = {'SMM (Untreated)': 'SMM\nUntreated', 'SMM (Treated)': 'SMM\nTreated'}
        xlabels = [f'{_disp.get(d, d)}\n(n={len(tp_data[tp_data["group"]==d])})' for d in PLOT_GROUPS]
        ax.set_xticklabels(xlabels, fontsize=9.5)
        ax.tick_params(axis='y', labelsize=11)

        if ax_idx == 0:
            ax.set_ylabel('APRIL Expression\n[log$_2$(CPM+1)]', fontsize=12, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlim(-0.5, len(PLOT_GROUPS) - 0.5)

    all_data = patient_data['APRIL_expr'].values
    for ax in axes:
        ax.set_ylim(0, all_data.max() + 0.65)

    plt.suptitle('APRIL Expression (Myeloid Cells)', fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    save_figure("Figure3C")
    plt.close()


if __name__ == "__main__":
    main()
