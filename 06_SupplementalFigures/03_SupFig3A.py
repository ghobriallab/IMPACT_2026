#!/usr/bin/env python3
"""Supplementary Figure 3A: myeloid TNFSF13 (APRIL) expression by sample shipment status.

Purpose:      Shipment control for Figure 3C, part one. Per-participant mean TNFSF13 expression
              in myeloid cells, shipped versus not shipped, at each timepoint, in participants
              with SMM. SMM is the only Figure 3C group in which shipment varies: every healthy
              donor and all but one participant with MGUS were shipped, so no comparison is
              possible there. Both treatment arms are pooled, since the question is about sample
              handling rather than disease.

              Shipment does shift the absolute level at the post-vaccination timepoint. Panel B
              asks the question that matters for the paper, whether the disease comparison
              Figure 3C makes survives when shipment is held constant, and it does.

Statistics:   Two-sided Wilcoxon rank-sum, raw p, with the effect size r = |Z|/sqrt(N) appended
              when p < 0.1, matching the convention of the main figures. Benjamini-Hochberg
              across the two timepoints is written to the stats table.

Inputs:       H5AD_NORM (scRNAseq_IMPACT_Zenodo.h5ad) for TNFSF13 and the cell labels, and
              data/metadata/Supplementary_Table_4_scRNAseq_sample_list.csv for Shipment_FedEx.
              The per-participant aggregation and the shipment join are imported from
              03_SupFig3B.py so the two panels cannot describe different data.

Outputs:      figures/SupFig3A.png (and PDF + SVG); tables/SupFig3A_stats.csv.

Dependencies: Python + scanpy, pandas, numpy, scipy, statsmodels, matplotlib; reads config.py.
"""
import sys, glob, importlib.util
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import *

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as _fm

# Every figure in this paper is set in Arial; fall back to the metrically identical Liberation
# Sans and rewrite the SVG font-family. svg.fonttype="none" keeps text editable.
FONT = "Arial"
_av = {f.name for f in _fm.fontManager.ttflist}
if FONT not in _av:
    for _p in glob.glob("/usr/share/fonts/**/LiberationSans-*.ttf", recursive=True):
        _fm.fontManager.addfont(_p)
    _av = {f.name for f in _fm.fontManager.ttflist}
PLOT_FONT = FONT if FONT in _av else ("Liberation Sans" if "Liberation Sans" in _av else "sans-serif")
plt.rcParams["font.family"] = PLOT_FONT
plt.rcParams["svg.fonttype"] = "none"

# Panel B computes the per-participant values and the shipment join; import them rather than
# reimplement, so the two panels cannot drift apart.
_spec = importlib.util.spec_from_file_location(
    "supfig3b", Path(__file__).parent / "03_SupFig3B.py")
sf3b = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sf3b)

ARMS = ['Shipped', 'Not shipped']
ARM_COLORS = {'Shipped': 'steelblue', 'Not shipped': 'tomato'}
TIMEPOINTS = ['Pre-Vx', 'Post-Vx']
RNG_SEED = 2026


def save_figure(basename, dpi=300):
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


def main():
    pat = sf3b.per_patient_expression()
    pat = sf3b.attach_metadata(pat)
    smm = pat[pat['group'].str.startswith('SMM') & pat['Shipping'].notna()].copy()
    print("Participants with SMM by timepoint and shipment status:")
    print(smm.groupby(['Timepoint', 'Shipping']).size().unstack(fill_value=0).to_string())

    rows = []
    for tp in TIMEPOINTS:
        s = smm[smm['Timepoint'] == tp]
        a = s.loc[s['Shipping'] == 'Shipped', 'APRIL_expr'].values
        b = s.loc[s['Shipping'] == 'Not shipped', 'APRIL_expr'].values
        u = stats.mannwhitneyu(a, b, alternative='two-sided')
        n1, n2 = len(a), len(b)
        z = (u.statistic - n1 * n2 / 2) / np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
        rows.append(dict(Timepoint=tp, n_shipped=n1, n_not_shipped=n2,
                         median_shipped=float(np.median(a)), median_not_shipped=float(np.median(b)),
                         p=float(u.pvalue), r=abs(z) / np.sqrt(n1 + n2)))
    st = pd.DataFrame(rows)
    st['q'] = multipletests(st['p'], method='fdr_bh')[1]
    print("\nShipped versus not shipped, SMM:")
    print(st.round(4).to_string(index=False))

    tables_dir = REPO_DIR / "tables"
    tables_dir.mkdir(exist_ok=True)
    st.to_csv(tables_dir / "SupFig3A_stats.csv", index=False)

    rng = np.random.default_rng(RNG_SEED)
    fig, axes = plt.subplots(1, 2, figsize=(5.6, 4.2), sharey=True)
    plt.subplots_adjust(wspace=0.08)
    ymin = smm['APRIL_expr'].min()
    ymax = smm['APRIL_expr'].max()
    yr = ymax - ymin

    for ax, tp in zip(axes, TIMEPOINTS):
        s = smm[smm['Timepoint'] == tp]
        data = [s.loc[s['Shipping'] == arm, 'APRIL_expr'].values for arm in ARMS]
        pos = [1, 2]
        parts = ax.violinplot(data, positions=pos, widths=0.8, showextrema=False)
        for body, arm in zip(parts['bodies'], ARMS):
            body.set_facecolor(ARM_COLORS[arm]); body.set_alpha(0.35); body.set_edgecolor('none')
        bp = ax.boxplot(data, positions=pos, widths=0.4, patch_artist=True, showfliers=False,
                        medianprops=dict(color='black', linewidth=1.4))
        for patch, arm in zip(bp['boxes'], ARMS):
            patch.set_facecolor(ARM_COLORS[arm]); patch.set_alpha(0.55)
        for i, (vals, arm) in enumerate(zip(data, ARMS)):
            ax.scatter(rng.normal(pos[i], 0.07, size=len(vals)), vals, s=22,
                       color=ARM_COLORS[arm], edgecolor='black', linewidth=0.4, alpha=0.9, zorder=3)

        row = st[st['Timepoint'] == tp].iloc[0]
        y = ymax + 0.06 * yr
        h = 0.018 * yr
        ax.plot([pos[0], pos[0], pos[1], pos[1]], [y - h, y, y, y - h], color='black', linewidth=1.0)
        lab = f"p={row['p']:.1e}" if row['p'] < 0.01 else f"p={row['p']:.2f}"
        if row['p'] < 0.1:
            lab += f", r={row['r']:.2f}"
        ax.text(np.mean(pos), y + 0.4 * h, lab, ha='center', va='bottom', fontsize=9)

        ax.set_xticks(pos)
        ax.set_xticklabels([f"{arm}\n(n={len(v)})" for arm, v in zip(ARMS, data)], fontsize=9)
        ax.set_title(tp, fontsize=11)
        ax.tick_params(axis='y', labelsize=9)
        ax.set_ylim(ymin - 0.06 * yr, ymax + 0.26 * yr)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)

    axes[0].set_ylabel('Myeloid TNFSF13 (APRIL) expression\n(per-participant mean)', fontsize=9)
    save_figure("SupFig3A")


if __name__ == "__main__":
    main()
