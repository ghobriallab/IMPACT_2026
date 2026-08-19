#!/usr/bin/env python3
"""Supplementary Figure 3A: sample shipment and myeloid TNFSF13 (APRIL) expression.

Purpose:      Supplementary Figure 3A: shipment control for Figure 3C. Per-participant mean
              TNFSF13 expression in myeloid cells, split by whether the sample was shipped,
              within each of the four Figure 3C groups and at both timepoints. Shipment shifts
              the absolute level in treatment-naive SMM at the post-vaccination timepoint, but
              the disease comparison that Figure 3C makes is unaffected: restricting Figure 3C
              to shipped samples leaves all six contrasts non-significant (all q>=0.99, against
              all q>=0.34 on the full cohort), both being null.

Inputs:       H5AD_NORM (scRNAseq_IMPACT_Zenodo.h5ad) for TNFSF13 expression and cell labels;
              data/metadata/Supplementary_Table_4_scRNAseq_sample_list.csv for Shipment_FedEx
              and the age/sex covariates.

Outputs:      figures/SupFig3A.png (and PDF + SVG); tables/SupFig3A_stats.csv and
              tables/SupFig3A_shipped_only_Figure3C.csv.

Dependencies: Python + scanpy, pandas, numpy, scipy, statsmodels, matplotlib; reads config.py.
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
import matplotlib
matplotlib.use('Agg')
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
DISEASE_GROUPS = ['HD', 'MGUS', 'SMM']
PLOT_GROUPS = ['HD', 'MGUS', 'SMM (Untreated)', 'SMM (Treated)']
ARMS = ['Shipped', 'Not shipped']
MIN_PER_ARM = 3          # below this a group cannot be compared and is annotated instead
RNG_SEED = 2026          # jitter is seeded so the panel is reproducible

# Shipment is read from Supplementary Table 4. 1 = shipped, 0 = not shipped. A blank field is
# EXCLUDED rather than assumed shipped, so a missing record never counts as evidence either way.
SHIP_MAP = {1.0: 'Shipped', 0.0: 'Not shipped'}


def load_metadata():
    # Shipped in the Zenodo deposit as metadata/Supplementary_Table_4_scRNAseq_sample_list.csv.
    meta_path = DATA_DIR / "metadata" / "Supplementary_Table_4_scRNAseq_sample_list.csv"
    if not meta_path.exists():
        meta_path = REPO_DIR.parent / "Supplementary_Tables" / "Supplementary_Table_4_scRNAseq_sample_list.csv"
    return pd.read_csv(meta_path)


def main():
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

    # Identical cell selection to 03_Figure3/02_Figure3C.py so the two panels describe the
    # same cells: myeloid only, doublet/QC/CLL labels dropped, both vaccination timepoints.
    tp_map = {'Pre-Vx': 'Pre-Vx', 'Post-2nd': 'Post-Vx'}
    clean = (~cells['l2'].isin(['QC_removed', 'CLL'])) & (~cells['l2'].str.startswith('db:'))
    keep = (clean
            & cells['Annotation_Level_1'].isin(['Mono', 'DC'])
            & cells['tp_raw'].isin(tp_map.keys())
            & cells['diagnosis'].isin(DISEASE_GROUPS))
    myeloid = cells[keep].copy()
    myeloid['Timepoint'] = myeloid['tp_raw'].map(tp_map)
    myeloid['group'] = np.where(
        myeloid['diagnosis'] != 'SMM', myeloid['diagnosis'],
        np.where(myeloid['treat'] == 'Never_treated', 'SMM (Untreated)', 'SMM (Treated)'))

    positions = myeloid['pos'].to_numpy()
    expr = []
    for i in range(0, len(positions), 100000):
        chunk = adata.X[list(positions[i:i + 100000]), april_idx]
        chunk = chunk.toarray().flatten() if hasattr(chunk, 'toarray') else np.asarray(chunk).flatten()
        expr.extend(chunk)
    myeloid['APRIL_expr'] = expr

    pat = (myeloid.groupby(['patient', 'group', 'Timepoint'])['APRIL_expr']
                  .mean().reset_index())

    meta = load_metadata()
    meta['tp'] = meta['Timepoint'].map(tp_map)
    ship = (meta[meta['tp'].notna()]
            .groupby(['Deidentified_Patient_ID', 'tp'])['Shipment_FedEx']
            .agg(lambda s: s.dropna().unique().tolist()).reset_index())
    # A participant-timepoint with two conflicting records would be ambiguous; assert there is none.
    assert (ship['Shipment_FedEx'].map(len) <= 1).all(), "conflicting shipment records"
    ship['Shipping'] = ship['Shipment_FedEx'].map(lambda v: SHIP_MAP.get(v[0]) if len(v) == 1 else None)
    ship = ship.rename(columns={'Deidentified_Patient_ID': 'patient', 'tp': 'Timepoint'})
    pat = pat.merge(ship[['patient', 'Timepoint', 'Shipping']], on=['patient', 'Timepoint'], how='left')

    cov = (meta.groupby('Deidentified_Patient_ID')
                .agg(Age=('Age', 'median'),
                     Sex=('Sex', lambda s: s.dropna().mode().iloc[0] if s.dropna().size else None))
                .reset_index().rename(columns={'Deidentified_Patient_ID': 'patient'}))
    pat = pat.merge(cov, on='patient', how='left')

    n_blank = int(pat['Shipping'].isna().sum())
    print(f"Participant-timepoints: {len(pat)}; shipment not recorded for {n_blank} (excluded)")

    # --- statistics: shipped vs not shipped within each group and timepoint ------------------
    rows = []
    for tp in ['Pre-Vx', 'Post-Vx']:
        for g in PLOT_GROUPS:
            s = pat[(pat['Timepoint'] == tp) & (pat['group'] == g)]
            a = s.loc[s['Shipping'] == 'Shipped', 'APRIL_expr'].values
            b = s.loc[s['Shipping'] == 'Not shipped', 'APRIL_expr'].values
            rec = dict(Timepoint=tp, group=g, n_shipped=len(a), n_not_shipped=len(b),
                       median_shipped=np.median(a) if len(a) else np.nan,
                       median_not_shipped=np.median(b) if len(b) else np.nan,
                       p=np.nan, r=np.nan)
            if len(a) >= MIN_PER_ARM and len(b) >= MIN_PER_ARM:
                u = stats.mannwhitneyu(a, b, alternative='two-sided')
                n1, n2 = len(a), len(b)
                z = (u.statistic - n1 * n2 / 2) / np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
                rec['p'] = float(u.pvalue)
                rec['r'] = abs(z) / np.sqrt(n1 + n2)   # r = |Z|/sqrt(N), as elsewhere in the paper
            rows.append(rec)
    st = pd.DataFrame(rows)
    ok = st['p'].notna()
    st.loc[ok, 'q'] = multipletests(st.loc[ok, 'p'], method='fdr_bh')[1]
    print("\nShipped vs not shipped, within group:")
    print(st.round(4).to_string(index=False))

    # --- the claim that matters: does Figure 3C's comparison change on shipped samples only? --
    ship_rows = []
    for tp in ['Pre-Vx', 'Post-Vx']:
        t = pat[(pat['Timepoint'] == tp) & (pat['Shipping'] == 'Shipped')].dropna(subset=['Age', 'Sex']).copy()
        t['rk'] = t['APRIL_expr'].rank()
        t['group'] = pd.Categorical(t['group'], categories=PLOT_GROUPS)
        fit = smf.ols('rk ~ C(group, Treatment(reference="HD")) + Age + C(Sex)', data=t).fit()
        for g in PLOT_GROUPS[1:]:
            ship_rows.append(dict(Timepoint=tp, Comparison=f'HD vs {g}',
                                  n=int((t['group'] == g).sum()), n_HD=int((t['group'] == 'HD').sum()),
                                  p=float(fit.pvalues[f'C(group, Treatment(reference="HD"))[T.{g}]'])))
    sh = pd.DataFrame(ship_rows)
    sh['q'] = multipletests(sh['p'], method='fdr_bh')[1]
    print("\nFigure 3C contrasts recomputed on shipped samples only:")
    print(sh.round(4).to_string(index=False))
    print(f"minimum q on shipped samples only: {sh['q'].min():.4f}")

    # config.py exposes FIGURES_DIR but no TABLES_DIR (config.R does); derive it the same way.
    tables_dir = REPO_DIR / "tables"
    tables_dir.mkdir(exist_ok=True)
    st.to_csv(tables_dir / "SupFig3A_stats.csv", index=False)
    sh.to_csv(tables_dir / "SupFig3A_shipped_only_Figure3C.csv", index=False)

    # --- plot -------------------------------------------------------------------------------
    ARM_COLORS = {'Shipped': '#4682B4', 'Not shipped': '#EE5C42'}
    rng = np.random.default_rng(RNG_SEED)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.0), sharey=True)
    plt.subplots_adjust(wspace=0.08)
    ymax = pat['APRIL_expr'].max()

    for ax_idx, (ax, tp) in enumerate(zip(axes, ['Pre-Vx', 'Post-Vx'])):
        xt, xl = [], []
        for gi, g in enumerate(PLOT_GROUPS):
            base = gi * 2.6
            for ai, arm in enumerate(ARMS):
                x = base + ai
                s = pat[(pat['Timepoint'] == tp) & (pat['group'] == g) & (pat['Shipping'] == arm)]
                vals = s['APRIL_expr'].values
                xt.append(x)
                xl.append(f"{arm.split()[0] if arm == 'Shipped' else 'Not'}\n(n={len(vals)})")
                if len(vals) == 0:
                    continue
                bp = ax.boxplot([vals], positions=[x], widths=0.65, patch_artist=True, showfliers=False)
                bp['boxes'][0].set_facecolor(ARM_COLORS[arm])
                bp['boxes'][0].set_alpha(0.7)
                for el in ['whiskers', 'caps', 'medians']:
                    plt.setp(bp[el], color='black', linewidth=1.2)
                ax.scatter(rng.normal(x, 0.1, size=len(vals)), vals, alpha=0.6,
                           color=ARM_COLORS[arm], edgecolor='white', s=42, zorder=3, linewidth=0.5)

            row = st[(st['Timepoint'] == tp) & (st['group'] == g)].iloc[0]
            if pd.notna(row['q']):
                y = ymax + 0.10
                ax.plot([base, base + 1], [y, y], 'k-', linewidth=1)
                lab = f"q={row['q']:.2f}" if row['q'] >= 0.01 else f"q={row['q']:.3f}"
                if row['q'] < 0.1:
                    lab += f", r={row['r']:.2f}"
                ax.text(base + 0.5, y + 0.02, lab, ha='center', va='bottom', fontsize=9.5)
            else:
                # HD and MGUS have no, or a single, non-shipped participant, so no test is
                # possible. State which on the panel rather than leaving a silent gap, so the
                # asymmetry reads as a property of the cohort and not as selective display.
                n_not = int(row['n_not_shipped'])
                msg = 'no non-shipped\nsamples' if n_not == 0 else f'only {n_not} non-shipped\nsample'
                ax.text(base + 0.5, ymax + 0.10, msg, ha='center', va='bottom',
                        fontsize=7.5, style='italic', color='#555555')

        ax.set_xticks(xt)
        ax.set_xticklabels(xl, fontsize=8.5)
        for gi, g in enumerate(PLOT_GROUPS):
            ax.text(gi * 2.6 + 0.5, -0.155 * (ymax + 0.55), g.replace(' (', '\n('),
                    ha='center', va='top', fontsize=10, transform=ax.transData)
        ax.set_title(tp, fontsize=14, fontweight='bold', pad=6)
        ax.tick_params(axis='y', labelsize=11)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlim(-0.8, (len(PLOT_GROUPS) - 1) * 2.6 + 1.8)
        ax.set_ylim(0, ymax + 0.55)
        if ax_idx == 0:
            ax.set_ylabel('APRIL Expression\n[log$_2$(CPM+1)]', fontsize=12, fontweight='bold')

    plt.suptitle('APRIL Expression (Myeloid Cells) by Sample Shipment', fontsize=14,
                 fontweight='bold', y=1.00)
    save_figure("SupFig3A")


if __name__ == "__main__":
    main()
