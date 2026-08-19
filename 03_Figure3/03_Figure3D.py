#!/usr/bin/env python3
"""Figure 3D (manuscript Figure 3D): APRIL-responsive module score in non-malignant plasma cells
across disease stage, SWIFT-seq cohort (Lightbody et al., Nat Cancer 2025).

Purpose:      Cross-sectional decline of the 15-gene APRIL-responsive module in NON-MALIGNANT
              plasma cells along NBM -> MGUS -> SMM -> NDMM, shown separately for bone marrow and
              peripheral blood.

Inputs:       data/external/swiftseq_april_persample_deid.csv -- de-identified per-sample summary
              (sample and patient keys, compartment, grouped disease stage, serial timepoint index,
              per-population cell counts and mean module scores). Cell-level scoring was performed
              upstream with scanpy sc.tl.score_genes() using the gene list in Supplementary Table 5;
              source data are in dbGaP phs003855.v1.p1 (controlled access).

Outputs:      figures/Figure3D.png (and PDF + SVG).

Dependencies: Python + pandas, numpy, scipy, statsmodels, matplotlib; reads config.py.
"""
import sys, glob, itertools
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import *

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib as mpl, matplotlib.pyplot as plt
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
mpl.rcParams['mathtext.fontset'] = 'dejavusans'

MIN_CELLS = 10                     # non-malignant plasma cells per sample, as in Figure 3E
ORDER = ['NBM', 'MGUS', 'SMM', 'NDMM']
# Palette matched to Figure 3C (its diag_colors) and to Supplementary Figure 4B, so the same
# disease stage is the same colour wherever it appears. NBM takes the HD blue, since it is
# the healthy reference of this cohort. NDMM is brown, one step darker along the progression
# axis and distinct from the SMM red.
SHADES = {'NBM': '#3498db', 'MGUS': '#f1c40f', 'SMM': '#e74c3c', 'NDMM': '#8B4513'}
COMPS = [('BM', 'Bone marrow'), ('PB', 'Peripheral blood')]
RNG_SEED = 42


def jt(groups):
    """Jonckheere-Terpstra ordered-trend statistic, normal approximation, two-sided."""
    J = 0
    for i, j in itertools.combinations(range(len(groups)), 2):
        a, b = groups[i], groups[j]
        J += sum((x < y) + 0.5 * (x == y) for x in a for y in b)
    n = [len(g) for g in groups]
    N = sum(n)
    mu = (N * N - sum(v * v for v in n)) / 4.0
    var = (N * N * (2 * N + 3) - sum(v * v * (2 * v + 3) for v in n)) / 72.0
    z = (J - mu) / np.sqrt(var)
    return z, 2 * (1 - stats.norm.cdf(abs(z)))


df = pd.read_csv(DATA_DIR / "external" / "swiftseq_april_persample_deid.csv")
df = df[df.n_normal_PC.fillna(0) >= MIN_CELLS]

# ONE SAMPLE PER PARTICIPANT PER COMPARTMENT, the earliest serial timepoint, matching Figure 3E.
# Several participants were sampled repeatedly and would otherwise contribute more than once to
# a cross-sectional comparison. Ties are technical replicates of the same timepoint and are broken
# on Sample_ID so the choice is deterministic.
_before = len(df)
df = (df.sort_values(['Patient_ID', 'Compartment', 'Timepoint_Index', 'Sample_ID'])
        .drop_duplicates(['Patient_ID', 'Compartment'], keep='first'))
print(f"One sample per participant per compartment: kept {len(df)} of {_before} qualifying samples")
_nostage = int(df.Disease_Stage.isna().sum())
df = df[df.Disease_Stage.notna()]
print(f"Excluded {_nostage} samples whose disease stage is not recorded in the source metadata")

fig, axes = plt.subplots(1, 2, figsize=(6.0, 4.6), sharey=True)
for ax, (key, nice) in zip(axes, COMPS):
    s = df[df.Compartment == key]
    data = [s.loc[s.Disease_Stage == st, 'APRIL_normal'].values for st in ORDER]
    z, pj = jt(data)
    ps = [stats.mannwhitneyu(data[0], v, alternative='two-sided').pvalue for v in data[1:]]
    qs = multipletests(ps, method='fdr_bh')[1]
    print(f"{nice}: n={[len(v) for v in data]}, medians="
          f"{[round(float(np.median(v)), 4) for v in data]}, JT z={z:.3f} p={pj:.3g}, "
          f"q(NBM vs MGUS/SMM/NDMM)={[float(f'{q:.3g}') for q in qs]}")

    pos = list(range(1, len(ORDER) + 1))
    bp = ax.boxplot(data, positions=pos, widths=0.6, patch_artist=True, showfliers=False,
                    medianprops=dict(color='black', linewidth=1.8))
    for patch, st in zip(bp['boxes'], ORDER):
        patch.set_facecolor(SHADES[st]); patch.set_alpha(0.55)
    rng = np.random.default_rng(RNG_SEED)
    for i, (vals, st) in enumerate(zip(data, ORDER)):
        ax.scatter(rng.normal(pos[i], 0.08, size=len(vals)), vals, s=18, color=SHADES[st],
                   edgecolor='black', linewidth=0.4, alpha=0.92)
    ax.set_xticks(pos)
    ax.set_xticklabels([f"{st}\nn={len(v)}" for st, v in zip(ORDER, data)], fontsize=8.5)
    ax.tick_params(axis='y', labelsize=9)

    ymax = max(np.max(v) for v in data); ymin = min(np.min(v) for v in data); yr = ymax - ymin
    lvls = [ymax + 0.08 * yr, ymax + 0.20 * yr, ymax + 0.32 * yr]

    def bracket(x1, x2, y, lab):
        h = 0.018 * yr
        ax.plot([x1, x1, x2, x2], [y - h, y, y, y - h], color='black', linewidth=1.0)
        ax.text((x1 + x2) / 2, y + 0.4 * h, lab, ha='center', va='bottom', fontsize=8)
    for j, q in enumerate(qs):
        bracket(pos[0], pos[j + 1], lvls[j], f"q={q:.2f}" if q >= 0.01 else f"q={q:.1e}")
    # JT annotation in the same style as Figure 2E: italic, with a colon after the test name.
    ax.text(np.mean(pos), ymax + 0.50 * yr,
            f"Jonckheere-Terpstra:\nz={z:.2f}, p={pj:.1e}", ha='center', va='bottom',
            fontsize=8, fontstyle='italic')
    ax.set_ylim(ymin - 0.05 * yr, ymax + 0.74 * yr)
    ax.set_title(nice, fontsize=10, pad=5)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)

axes[0].set_ylabel('APRIL-responsive module score\n(non-malignant plasma cells)', fontsize=9)
plt.tight_layout()
for ext in ("png", "pdf", "svg"):
    out = FIGURES_DIR / f"Figure3D.{ext}"
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    if ext == "svg" and PLOT_FONT != FONT:
        t = Path(out).read_text(encoding="utf-8")
        for qm in ('"', "'"):
            t = t.replace(f"font-family: {qm}{PLOT_FONT}{qm}", f"font-family: {FONT}")
        t = t.replace(f"font-family: {PLOT_FONT}", f"font-family: {FONT}")
        Path(out).write_text(t, encoding="utf-8")
    print(f"Saved: {out.name}")
