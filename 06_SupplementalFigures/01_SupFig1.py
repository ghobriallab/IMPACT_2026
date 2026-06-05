#!/usr/bin/env python3
"""Supplementary Figure 1: Vaccine response by SMM 20/2/20 risk tier.

Design:
  - Panels A-C: tier reproduction of Fig 1B / 1C / 1F. Columns are HD / MGUS / LR-SMM /
    IR-SMM / HR-SMM / MM. HD, MGUS, MM use the FULL eligible ELISA cohort. SMM is filtered
    to the 20/2/20 risk-classed subset (the metadata union with the older scRNA shipping-status
    file was rejected after audit -- 10/11 patients DISAGREED across the two curators).
    EACH SMM TIER COLUMN'S n LABEL IS ANNOTATED WITH THE NUMBER OF TREATED PATIENTS in
    parentheses, so the treatment confound is visible. Treatment confounding is severe in
    HR-SMM (4/4 treated at peak post-2nd dose; this cohort enrolled HR-SMM predominantly
    in treatment-arm contexts) and notable in IR/LR-SMM, which is why we cannot use
    SupFig1 to detect a clean tier-stratified vaccine-response gradient at peak.
  - Continuous tumor-burden panels were dropped after the user judged the
    M-spike / BM PC% / FLC ratio scatters within treatment-naive SMM to be
    underpowered and inconclusive across metrics (M-spike rho=-0.21 p=0.28,
    BM PC% rho=+0.15 p=0.48, FLC ratio rho=-0.31 p=0.081); the numerical
    results are reported in the manuscript text and in the response letter only.
  - Cytogenetics panels are NOT included: only 4 untreated risk-classed SMM
    at peak; n_pos <= 2 per CNA, not testable.
  - Fig 1E (waning slope) NOT reproduced -- per-tier serial-sample n = 2-6.
  - Age + sex adjustment via rank-based ANCOVA (Fig 1 main convention). JT ordered-trend
    test on age+sex residuals across HD < MGUS < LR-SMM < IR-SMM < HR-SMM; "Jonckheere-Terpstra"
    spelled out in full per manuscript convention. Statistics drawn IN-PANEL as brackets.

Purpose:      Supplementary Figure 1: vaccine response by SMM 20/2/20 risk tier. Three panels (Fig 1B/C/F equivalents) split SMM into LR/IR/HR with per-tier n labels annotating treated/treatment-naive counts; statistics use age+sex-adjusted rank-based ANCOVA + JT ordered-trend on HD<MGUS<LR<IR<HR.

Inputs:       data/elisa/elisa_spike_post2nd.csv, data/elisa/elisa_spike_post3rd.csv, data/smm_risk_strat.csv (de-identified 20/2/20 tier table).

Outputs:      figures/SupFig1.png.

Dependencies: Python + pandas, numpy, matplotlib, scipy, statsmodels; reads config.py.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import *

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

RNG_SEED = 2026

ELISA_2ND = DATA_DIR / 'elisa' / 'elisa_spike_post2nd.csv'
ELISA_3RD = DATA_DIR / 'elisa' / 'elisa_spike_post3rd.csv'
RISK      = DATA_DIR / 'smm_risk_strat.csv'

TIER_ORDER  = ['HD', 'MGUS', 'LRSMM', 'IRSMM', 'HRSMM', 'MM']
TIER_LABELS = {'HD':'HD', 'MGUS':'MGUS', 'LRSMM':'LR-SMM', 'IRSMM':'IR-SMM',
               'HRSMM':'HR-SMM', 'MM':'MM'}
TIER_COLORS = ['steelblue','orange','#FCD0A1','#EE5C42','#8B0000','#458B00']

risk = pd.read_csv(RISK)

# Helpers ---
def make_panel_df(elisa_file, t_lo, t_hi, day_col):
    """Per-patient first-sample DataFrame.

    NOTE: All risk-classed SMM are KEPT regardless of treatment status
    for the tier panels; the treatment count per tier is annotated in the x-label so
    the confound is transparent.
    """
    d = pd.read_csv(elisa_file)
    d = d[(d[day_col] >= t_lo) & (d[day_col] <= t_hi)].copy()
    d = d.sort_values(['Common_ID', day_col]).drop_duplicates('Common_ID', keep='first')
    d['Group'] = d['Disease'].map({'Healthy':'HD','MGUS':'MGUS','IgM-MGUS':'MGUS','MM':'MM','SMM':None})
    smm_rows = d['Disease'] == 'SMM'
    d.loc[smm_rows, 'Group'] = d.loc[smm_rows, 'Common_ID'].map(risk.set_index('Common_ID')['SMM_risk'])
    return d[d['Group'].isin(TIER_ORDER)]

def jt_trend(df, value_col, group_col, order, adjust_with):
    sub = df.dropna(subset=[value_col, group_col] + list(adjust_with)).copy()
    sub['rk'] = sub[value_col].rank()
    formula = 'rk ~ ' + ' + '.join([f'C({c})' if sub[c].dtype == object else c for c in adjust_with])
    sub['resid'] = smf.ols(formula, data=sub).fit().resid
    sub = sub[sub[group_col].isin(order)]
    groups = [sub[sub[group_col]==g]['resid'].values for g in order]
    if any(len(g) < 1 for g in groups) or sum(len(g) for g in groups) < 3:
        return float('nan'), float('nan'), 0
    U = 0
    for i in range(len(order)-1):
        for j in range(i+1, len(order)):
            for a in groups[i]:
                for b in groups[j]:
                    if   b > a: U += 1
                    elif b == a: U += 0.5
    n = [len(g) for g in groups]; N = sum(n)
    mu = (N*N - sum(ni*ni for ni in n)) / 4.0
    sigma2 = (N*N*(2*N+3) - sum(ni*ni*(2*ni+3) for ni in n)) / 72.0
    z = (U - mu) / np.sqrt(sigma2) if sigma2 > 0 else float('nan')
    p = 2 * stats.norm.sf(abs(z))
    return z, p, N

def adjusted_q_vs_hd(df, value_col, group_col, ref='HD'):
    sub = df.dropna(subset=[value_col, group_col, 'Age', 'Sex']).copy()
    sub['rk'] = sub[value_col].rank()
    sub[group_col] = pd.Categorical(sub[group_col], categories=[ref] + [g for g in TIER_ORDER if g != ref])
    m = smf.ols(f"rk ~ C({group_col}) + Age + C(Sex)", data=sub).fit()
    ps, names = [], []
    for g in TIER_ORDER:
        if g == ref: continue
        key = f"C({group_col})[T.{g}]"
        if key in m.pvalues.index:
            ps.append(m.pvalues[key]); names.append(g)
    if not ps: return {}
    q = multipletests(ps, method='fdr_bh')[1]
    return dict(zip(names, q))

def fmt_p(p):
    if pd.isna(p): return 'NA'
    return f"{p:.1e}" if p < 0.001 else f"{p:.3f}"

def draw_bracket(ax, x1, x2, y, label, tip_frac=0.012, fontsize=7.2, linewidth=0.6,
                  text_x=None):
    """Draw a square bracket between x1 and x2 at height y with a label above.
    If text_x is provided, the label is positioned at that x-coordinate (default = midpoint).
    """
    ylim = ax.get_ylim(); tip = (ylim[1] - ylim[0]) * tip_frac
    ax.plot([x1, x1, x2, x2], [y - tip, y, y, y - tip],
            color='black', linewidth=linewidth, clip_on=False, solid_capstyle='butt')
    tx = (x1 + x2) / 2 if text_x is None else text_x
    ax.text(tx, y + tip * 0.6, label, ha='center', va='bottom',
            fontsize=fontsize, clip_on=False)

# ============================================================
# Tier panels A/B/C — Fig 1B / 1C / 1F (ALL risk-classed SMM, treatment annotated)
# ============================================================
b = make_panel_df(ELISA_2ND, 14, 60,    'Days_post2nd')
c = make_panel_df(ELISA_2ND, 60.0001, 120, 'Days_post2nd')
f_ = make_panel_df(ELISA_3RD, 0, 100000, 'Days_post3rd')

def tier_panel_data(df):
    d = df.dropna(subset=['ELISA_Titer'])
    counts = {t: (d['Group']==t).sum() for t in TIER_ORDER}
    treated = {t: ((d['Group']==t) & (d.get('Ever_treated','') == 'Yes')).sum()
               for t in TIER_ORDER}
    q = adjusted_q_vs_hd(d, 'ELISA_Titer', 'Group')
    z, p_jt, _ = jt_trend(d, 'ELISA_Titer', 'Group',
                          order=['HD','MGUS','LRSMM','IRSMM','HRSMM'],
                          adjust_with=['Age','Sex'])
    return {'df': d, 'counts': counts, 'treated': treated, 'q': q, 'jt_z': z, 'jt_p': p_jt}

tier_results = [
    ('A', tier_panel_data(b),  '2 weeks – 2 months\npost-2nd dose (Fig 1B)'),
    ('B', tier_panel_data(c),  '2 – 4 months\npost-2nd dose (Fig 1C)'),
    ('C', tier_panel_data(f_), 'Post-3rd dose (Fig 1F)'),
]

# ============================================================
# FIGURE — 1 row × 3 cols
# ============================================================
fig = plt.figure(figsize=(15.0, 6.5))
gs = fig.add_gridspec(1, 3, hspace=0.30, wspace=0.28,
                       left=0.06, right=0.985, top=0.85, bottom=0.20)

for col, (lbl, r, title) in enumerate(tier_results):
    ax = fig.add_subplot(gs[0, col])
    d = r['df']
    data = [d[d['Group']==t]['ELISA_Titer'].values for t in TIER_ORDER]
    bp = ax.boxplot(data, positions=range(len(TIER_ORDER)), widths=0.65,
                    patch_artist=True, showfliers=False)
    for patch, c0 in zip(bp['boxes'], TIER_COLORS):
        patch.set_facecolor(c0); patch.set_alpha(0.55)
    for el in ['whiskers','caps','medians']:
        plt.setp(bp[el], color='black', linewidth=1.0)
    for i, (v, c0) in enumerate(zip(data, TIER_COLORS)):
        if len(v):
            ax.scatter(np.random.default_rng(RNG_SEED+i).normal(i, 0.08, len(v)), v,
                       s=18, color=c0, edgecolor='white', linewidth=0.4, alpha=0.85, zorder=3)
    ax.set_xticks(range(len(TIER_ORDER)))
    xlab = []
    for t in TIER_ORDER:
        ntot = r['counts'][t]; ntx = r['treated'][t]
        if t in ('LRSMM','IRSMM','HRSMM') and ntot > 0:
            xlab.append(f"{TIER_LABELS[t]}\n(n={ntot}; {ntx} tx)")
        else:
            xlab.append(f"{TIER_LABELS[t]}\n(n={ntot})")
    ax.set_xticklabels(xlab, fontsize=8.5, rotation=30, ha='right')
    if col == 0: ax.set_ylabel('Spike IgG titer (OD$_{450-570}$)', fontsize=10)
    ax.set_title(title, fontsize=10.5, fontweight='bold', pad=6)
    ax.spines[['top','right']].set_visible(False)

    # In-panel brackets: q vs HD (5 brackets, stacked) + JT trend (top bracket)
    flat_all = np.concatenate([v for v in data if len(v)])
    y_top = float(np.max(flat_all)); y_min = float(np.min(flat_all))
    step  = (y_top - y_min) * 0.16
    for k, g in enumerate(['MGUS','LRSMM','IRSMM','HRSMM','MM']):
        if g not in r['q']: continue
        y = y_top + (k + 1) * step
        draw_bracket(ax, 0, k + 1, y, f"q={fmt_p(r['q'][g])}", fontsize=9.0)
    jt_y = y_top + 6.4 * step
    draw_bracket(ax, 0, 4, jt_y,
                 f"Jonckheere-Terpstra (HD→HR-SMM, age+sex-adj):\nz={r['jt_z']:.2f}, p={fmt_p(r['jt_p'])}",
                 fontsize=9.0, text_x=2.5)
    ax.set_ylim(top=jt_y + 3.2 * step)

# Panel-letter labels
for lbl, spec in [('A', gs[0,0]),('B', gs[0,1]),('C', gs[0,2])]:
    pos = spec.get_position(fig)
    fig.text(pos.x0 - 0.012, pos.y1 + 0.005, lbl, fontsize=13, fontweight='bold')

fig.suptitle('Vaccine response by SMM 20/2/20 risk tier (treated counts annotated under n)',
             fontsize=11.5, fontweight='bold', y=0.96)

OUT_PNG = FIGURES_DIR / 'SupFig1_RiskGenetics_VaccineResponse.png'
plt.savefig(OUT_PNG, dpi=300, bbox_inches='tight')
plt.savefig(str(OUT_PNG).replace('.png','.pdf'), bbox_inches='tight')
plt.savefig(str(OUT_PNG).replace('.png','.svg'), bbox_inches='tight')
print(f"\nSaved: {OUT_PNG}")

# Console dump ---
print("\n=== TIER PANELS (all risk-classed SMM, treated count annotated) ===")
for lbl, r, title in tier_results:
    print(f"  {lbl} ({title.splitlines()[0]}):")
    print(f"    counts (treated/total):")
    for t in TIER_ORDER:
        nt = r['counts'][t]; tr = r['treated'][t]
        if t in ('LRSMM','IRSMM','HRSMM'):
            print(f"      {TIER_LABELS[t]:8s} n={nt:3d}  treated={tr}/{nt}")
        else:
            print(f"      {TIER_LABELS[t]:8s} n={nt:3d}")
    q_items = ', '.join([f"{k}={v:.3g}" for k, v in r['q'].items()])
    print(f"    q vs HD (age+sex adj, BH): {{ {q_items} }}")
    print(f"    JT (HD→HR-SMM): z={r['jt_z']:.2f}, p={r['jt_p']:.4g}")
plt.close(fig)
