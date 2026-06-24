#!/usr/bin/env python3
# ============================================================================
# Purpose:      Figure 2F: paired pre/post-teclistamab plasma Olink for APRIL,
#               BAFF, and sBCMA in 10 patients enrolled on the teclistamab arm
#               of Immuno-PRISM (NCT05469893, DFCI 22-154). The 3-protein
#               composite is cited as a single Figure 2F in the manuscript.
#               Three paired-design boxplots with patient-level connector lines;
#               paired Wilcoxon signed-rank tests with Benjamini-Hochberg
#               correction across the three proteins; median delta-NPX
#               annotation.
# Inputs:       data/external/PrePostTEC_olink_deid.csv (deid: P01..P10).
# Outputs:      Figure2F_PrePostTEC_3panels.{png,pdf,svg} in figures/.
# Dependencies: Python + pandas, scipy, matplotlib. Sources ../config.py.
# ============================================================================
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
try:
    from config import DATA_DIR, FIGURES_DIR
except Exception:
    DATA_DIR = HERE.parent / "data"
    FIGURES_DIR = HERE.parent / "figures"

IN = Path(DATA_DIR) / "external" / "PrePostTEC_olink_deid.csv"
FIGURES_DIR = Path(FIGURES_DIR)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

PROTEINS = ["APRIL", "BAFF", "sBCMA"]
GENE_LABELS = {"APRIL": "TNFSF13 (APRIL)",
               "BAFF": "TNFSF13B (BAFF)",
               "sBCMA": "TNFRSF17 (sBCMA)"}
BOX_COLORS = ["#A6CEE3", "#FBB4AE"]  # light blue baseline, light red post-TEC

# ---- BH adjustment helper (no external dep) ----
def bh_qvalues(p_array):
    """Benjamini-Hochberg adjusted q-values (FDR)."""
    p_array = np.asarray(p_array, dtype=float)
    n = len(p_array)
    order = np.argsort(p_array)
    ranks_0idx = np.empty(n, dtype=int)
    ranks_0idx[order] = np.arange(n)
    raw_q = p_array * n / (ranks_0idx + 1)
    q = np.empty(n)
    min_so_far = 1.0
    for r in range(n - 1, -1, -1):
        idx = order[r]
        min_so_far = min(min_so_far, raw_q[idx])
        q[idx] = min(min_so_far, 1.0)
    return q

# ---- Load ----
df = pd.read_csv(IN)
df["NPX"] = pd.to_numeric(df["NPX"], errors="coerce")
print(f"Loaded {len(df)} rows; {df['Patient_ID'].nunique()} patients; proteins: {df['Protein'].unique()}")

# ---- Compute paired Wilcoxon stats for all three proteins UP FRONT (so we can BH-adjust together) ----
stats_data = {}
for protein in PROTEINS:
    sub = df[df["Protein"] == protein].copy()
    wide = sub.pivot(index="Patient_ID", columns="Timepoint", values="NPX").dropna()
    baseline = wide["Baseline"].values
    post = wide["PostTEC"].values
    delta = post - baseline
    n = len(wide)
    w_two = stats.wilcoxon(post, baseline, alternative="two-sided", zero_method="wilcox")
    median_delta = float(np.median(delta))
    median_fc_linear = float(2 ** median_delta)
    n_up = int(np.sum(delta > 0))
    stats_data[protein] = dict(
        baseline=baseline, post=post, delta=delta, n=n,
        p=float(w_two.pvalue), median_delta=median_delta,
        median_fc_linear=median_fc_linear, n_up=n_up,
    )
# Benjamini-Hochberg across the 3 proteins (family-wise correction for the 3 hypothesis tests
# on this targeted APRIL/BAFF/sBCMA panel; matches q<0.1 threshold convention in IMPACT).
q_array = bh_qvalues([stats_data[p]["p"] for p in PROTEINS])
for i, protein in enumerate(PROTEINS):
    stats_data[protein]["q"] = float(q_array[i])

# ---- One panel per protein ----
fig, axes = plt.subplots(1, 3, figsize=(8.5, 3.6), dpi=200)
positions = [1, 2]
np.random.seed(42)

for ax, protein in zip(axes, PROTEINS):
    s = stats_data[protein]
    baseline = s["baseline"]; post = s["post"]; delta = s["delta"]; n = s["n"]
    median_delta = s["median_delta"]; median_fc_linear = s["median_fc_linear"]; n_up = s["n_up"]

    # Boxplot
    bp = ax.boxplot([baseline, post], positions=positions, widths=0.55, showfliers=False,
                    patch_artist=True,
                    medianprops=dict(color="black", linewidth=1.4),
                    whiskerprops=dict(color="black", linewidth=1),
                    capprops=dict(color="black", linewidth=1))
    for patch, color in zip(bp["boxes"], BOX_COLORS):
        patch.set_facecolor(color)
        patch.set_edgecolor("black")
        patch.set_linewidth(1)
        patch.set_zorder(1)

    # Connector lines + points (lines OVER boxes)
    jitter = np.random.uniform(-0.06, 0.06, size=n)
    for i in range(n):
        ax.plot([positions[0] + jitter[i], positions[1] + jitter[i]],
                [baseline[i], post[i]],
                color="grey", linewidth=0.6, alpha=0.7, zorder=2.5)
    for i in range(n):
        ax.scatter([positions[0] + jitter[i], positions[1] + jitter[i]],
                   [baseline[i], post[i]],
                   s=22, color="#333333", edgecolor="black", linewidth=0.4, zorder=3)

    # Significance annotation bar
    ymax = max(post.max(), baseline.max())
    ymin = min(post.min(), baseline.min())
    y_span = ymax - ymin
    bar_y = ymax + 0.10 * y_span
    ax.plot([positions[0], positions[1]], [bar_y, bar_y], color="black", linewidth=0.9)
    ax.plot([positions[0], positions[0]], [bar_y - 0.025*y_span, bar_y], color="black", linewidth=0.9)
    ax.plot([positions[1], positions[1]], [bar_y - 0.025*y_span, bar_y], color="black", linewidth=0.9)
    q_val = s["q"]
    q_str = f"q = {q_val:.3f}" if q_val >= 0.001 else f"q = {q_val:.1e}"
    # Magnitude: median ΔNPX. NPX is log2-scale, so ΔNPX is the directly-measured paired
    # difference. Natively signed (no arrow needed); same units as the y-axis.
    delta_str = f"ΔNPX = {median_delta:+.1f}"
    ax.text((positions[0] + positions[1]) / 2, bar_y + 0.05 * y_span,
            f"{q_str}\n{delta_str}",
            ha="center", va="bottom", fontsize=8)

    # Cosmetics
    ax.set_title(GENE_LABELS[protein], fontsize=10, loc="left", fontweight="bold")
    ax.set_xticks(positions)
    ax.set_xticklabels(["Pre-TEC", "Post-TEC"], fontsize=9)
    ax.set_ylabel("Plasma NPX (Olink)", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=8)
    ax.set_ylim(ymin - 0.10*y_span, ymax + 0.40*y_span)

    print(f"{protein:>6}: median delta NPX = {median_delta:+.3f}, "
          f"fold-change = {median_fc_linear:.2f}, paired Wilcoxon p = {s['p']:.3g}, "
          f"BH q = {s['q']:.3g}, {n_up}/{n} directional")

fig.suptitle("Pre vs post teclistamab plasma proteomics (n=10 SMM patients, Immuno-PRISM)",
             fontsize=11, y=1.02)
plt.tight_layout()

for ext in ("png", "pdf", "svg"):
    out = FIGURES_DIR / f"Figure2F_PrePostTEC_3panels.{ext}"
    plt.savefig(out, dpi=200, bbox_inches="tight")
    print(f"Saved: {out}")
plt.close()
