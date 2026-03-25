#!/usr/bin/env python3
"""
plot_cooccurrence_report.py — Rapport HTML co-occurrence ARG × MGE
Pipeline MetagenAMR

Figures :
  1. Bubble plot  — Top ARGs × MGE category (tous samples agrégés)
  2. Network graph — ARG ↔ MGE category (edge list pondéré)

Inputs  : results/r_analysis/cooccurrence/01_cooccurrence_pairs.tsv
          results/r_analysis/cooccurrence/03_cooccurrence_network.tsv
Output  : results/figures/cooccurrence/cooccurrence_report.html
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import base64, io, os, sys, textwrap
from datetime import datetime

# ── Snakemake I/O ─────────────────────────────────────────────────────────────
if 'snakemake' in dir():
    input_pairs   = snakemake.input.pairs
    input_network = snakemake.input.network
    output_html   = snakemake.output.report
    DPI           = snakemake.params.get('dpi', 150)
else:
    input_pairs   = "results/r_analysis/cooccurrence/01_cooccurrence_pairs.tsv"
    input_network = "results/r_analysis/cooccurrence/03_cooccurrence_network.tsv"
    output_html   = "results/figures/cooccurrence/cooccurrence_report.html"
    DPI           = 150

os.makedirs(os.path.dirname(output_html), exist_ok=True)

# ── Palette MGE ───────────────────────────────────────────────────────────────
MGE_COLORS = {
    "integration/excision":             "#4E79A7",
    "transfer":                         "#F28E2B",
    "phage":                            "#59A14F",
    "replication/recombination/repair": "#E15759",
    "stability/transfer/defense":       "#B07AA1",
}
MGE_SHORT = {
    "integration/excision":             "Int/Exc",
    "transfer":                         "Transfer",
    "phage":                            "Phage",
    "replication/recombination/repair": "RRR",
    "stability/transfer/defense":       "STD",
}

# ── Chargement données ────────────────────────────────────────────────────────
pairs   = pd.read_csv(input_pairs,   sep='\t')
network = pd.read_csv(input_network, sep='\t')

print(f"[OK] Pairs   : {len(pairs):,} lignes")
print(f"[OK] Network : {len(network):,} edges")

# ── Helper : fig → base64 ─────────────────────────────────────────────────────
def fig_to_b64(fig, dpi=DPI):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=dpi,
                bbox_inches='tight', facecolor='white')
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Bubble plot Top ARGs × MGE category
# ══════════════════════════════════════════════════════════════════════════════
print("\n[Fig 1] Bubble plot ARGs × MGE category...")

# Agréger tous samples — top 20 ARGs par nombre total de paires
arg_totals = pairs.groupby('best_hit_aro').size().sort_values(ascending=False)
top_args   = arg_totals.head(20).index.tolist()

# Matrice : ARG × MGE category — nombre de contigs co-porteurs uniques
bubble_df = (
    pairs[pairs['best_hit_aro'].isin(top_args)]
    .groupby(['best_hit_aro', 'major_category'])
    .agg(n_contigs=('contig_id', 'nunique'),
         n_samples=('sample_id', 'nunique'))
    .reset_index()
)

mge_cats  = list(MGE_COLORS.keys())
arg_order = (bubble_df.groupby('best_hit_aro')['n_contigs']
             .sum().sort_values(ascending=True).index.tolist())

fig1, ax = plt.subplots(figsize=(10, 7))

for _, row in bubble_df.iterrows():
    x = mge_cats.index(row['major_category']) if row['major_category'] in mge_cats else None
    y = arg_order.index(row['best_hit_aro'])
    if x is None:
        continue
    size  = row['n_contigs'] * 18
    color = MGE_COLORS.get(row['major_category'], '#999999')
    ax.scatter(x, y, s=size, color=color, alpha=0.80,
               edgecolors='white', linewidths=0.6, zorder=3)
    if row['n_contigs'] >= 3:
        ax.text(x, y, str(int(row['n_contigs'])),
                ha='center', va='center',
                fontsize=6.5, color='white', fontweight='bold')

# Grid
ax.set_xlim(-0.6, len(mge_cats) - 0.4)
ax.set_ylim(-0.6, len(arg_order) - 0.4)
ax.set_xticks(range(len(mge_cats)))
ax.set_xticklabels([MGE_SHORT[c] for c in mge_cats],
                   fontsize=9, fontweight='bold')
ax.set_yticks(range(len(arg_order)))
ax.set_yticklabels([a for a in arg_order], fontsize=8, style='italic')
ax.yaxis.grid(True, linestyle='--', alpha=0.3, color='grey', linewidth=0.5)
ax.xaxis.grid(True, linestyle='--', alpha=0.3, color='grey', linewidth=0.5)
ax.set_axisbelow(True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel("MGE Functional Category", fontsize=10, labelpad=10)
ax.set_ylabel("Antibiotic Resistance Gene", fontsize=10, labelpad=10)
ax.set_title("Top 20 ARGs Co-localized with MGEs\n(bubble size = number of co-carrying contigs)",
             fontsize=11, fontweight='bold', pad=12)

# Légende taille
for s, label in [(5*18, "5"), (15*18, "15"), (30*18, "30")]:
    ax.scatter([], [], s=s, color='#888888', alpha=0.7,
               label=f'{label} contigs', edgecolors='white', linewidths=0.5)
ax.legend(title="Co-carrying contigs", fontsize=8, title_fontsize=8,
          loc='lower right', framealpha=0.9, edgecolor='#cccccc')

plt.tight_layout()
b64_fig1 = fig_to_b64(fig1)
plt.close(fig1)
print("  [OK] Bubble plot généré")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Network graph ARG ↔ MGE category
# ══════════════════════════════════════════════════════════════════════════════
print("\n[Fig 2] Network graph ARG ↔ MGE...")

# Filtrer top ARGs pour lisibilité (top 25 par poids total)
top_net = (network.groupby('source')['weight']
           .sum().sort_values(ascending=False).head(25).index.tolist())
net_df  = network[network['source'].isin(top_net)].copy()

# Layout manuel : MGE categories à gauche, ARGs à droite
mge_nodes  = list(MGE_COLORS.keys())
arg_nodes  = (net_df.groupby('source')['weight']
              .sum().sort_values(ascending=False).index.tolist())

n_mge = len(mge_nodes)
n_arg = len(arg_nodes)

mge_pos = {m: (0.0, 1 - i / max(n_mge - 1, 1)) for i, m in enumerate(mge_nodes)}
arg_pos = {a: (1.0, 1 - i / max(n_arg - 1, 1)) for i, a in enumerate(arg_nodes)}

fig2, ax = plt.subplots(figsize=(12, 9))
ax.set_xlim(-0.15, 1.25)
ax.set_ylim(-0.08, 1.08)
ax.axis('off')

# Normaliser épaisseur edges
max_w = net_df['weight'].max()
min_w = net_df['weight'].min()

for _, row in net_df.iterrows():
    src = row['source']
    tgt = row['target']
    if src not in arg_pos or tgt not in mge_pos:
        continue
    x1, y1 = arg_pos[src]
    x2, y2 = mge_pos[tgt]
    lw    = 0.3 + 3.5 * (row['weight'] - min_w) / max(max_w - min_w, 1)
    color = MGE_COLORS.get(tgt, '#aaaaaa')
    ax.plot([x1, x2], [y1, y2], color=color,
            alpha=0.35, linewidth=lw, zorder=1)

# Noeuds MGE (gauche)
for m, (x, y) in mge_pos.items():
    total_w = net_df[net_df['target'] == m]['weight'].sum()
    size    = 800 + total_w * 0.08
    ax.scatter(x, y, s=min(size, 4000), color=MGE_COLORS[m],
               edgecolors='white', linewidths=1.5, zorder=4)
    ax.text(x - 0.04, y, MGE_SHORT[m],
            ha='right', va='center', fontsize=9,
            fontweight='bold', color=MGE_COLORS[m])

# Noeuds ARG (droite)
for a, (x, y) in arg_pos.items():
    total_w = net_df[net_df['source'] == a]['weight'].sum()
    size    = 200 + total_w * 0.04
    ax.scatter(x, y, s=min(size, 1500), color='#475569',
               edgecolors='white', linewidths=1, zorder=4, alpha=0.85)
    ax.text(x + 0.04, y, a, ha='left', va='center',
            fontsize=7.5, style='italic', color='#1e293b')

# Labels colonnes
ax.text(0.0,  1.06, "MGE Category", ha='center', fontsize=10,
        fontweight='bold', color='#334155')
ax.text(1.0,  1.06, "ARG (Top 25)", ha='center', fontsize=10,
        fontweight='bold', color='#334155')
ax.axvline(0.5, color='#e2e8f0', linewidth=1, linestyle='--', zorder=0)

ax.set_title("ARG–MGE Co-occurrence Network\n(edge width ∝ number of co-carrying contigs, all samples)",
             fontsize=11, fontweight='bold', pad=16)

# Note bas de page
ax.text(0.5, -0.06,
        "RRR: Replication/Recombination/Repair  |  STD: Stability/Transfer/Defense",
        ha='center', fontsize=7.5, color='#94a3b8',
        style='italic', transform=ax.transAxes)

plt.tight_layout()
b64_fig2 = fig_to_b64(fig2)
plt.close(fig2)
print("  [OK] Network graph généré")

# ══════════════════════════════════════════════════════════════════════════════
# STATISTIQUES RÉSUMÉ
# ══════════════════════════════════════════════════════════════════════════════
n_contigs   = pairs.groupby(['sample_id', 'contig_id']).ngroups
n_args      = pairs['best_hit_aro'].nunique()
n_mge_types = pairs['mobileog_id'].nunique()
n_samples   = pairs['sample_id'].nunique()
top_arg     = pairs['best_hit_aro'].value_counts().index[0]
top_mge     = pairs['major_category'].value_counts().index[0]

# Top associations
top_assoc = (pairs.groupby(['best_hit_aro', 'major_category'])
             .agg(n_contigs_=('contig_id', 'nunique'))
             .reset_index()
             .sort_values('n_contigs_', ascending=False)
             .head(10))

assoc_rows = ""
for _, r in top_assoc.iterrows():
    assoc_rows += f"""
        <tr>
          <td style="font-style:italic">{r['best_hit_aro']}</td>
          <td><span class="badge" style="background:{MGE_COLORS.get(r['major_category'],'#888')};color:white">
              {MGE_SHORT.get(r['major_category'], r['major_category'])}</span></td>
          <td style="text-align:right;font-weight:600">{int(r['n_contigs_'])}</td>
        </tr>"""

# ══════════════════════════════════════════════════════════════════════════════
# HTML
# ══════════════════════════════════════════════════════════════════════════════
print("\n[HTML] Génération du rapport...")

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Co-occurrence ARG × MGE — MetagenAMR</title>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: 'Segoe UI', Arial, sans-serif; background: #f8fafc;
          color: #1e293b; line-height: 1.6; }}
  .header {{ background: linear-gradient(135deg, #1e293b 0%, #334155 100%);
             color: white; padding: 36px 48px; }}
  .header h1 {{ font-size: 1.8em; font-weight: 700; letter-spacing: -0.01em; }}
  .header p  {{ opacity: 0.7; font-size: 0.9em; margin-top: 6px; }}
  .container {{ max-width: 1200px; margin: 0 auto; padding: 40px 32px; }}
  .stats-grid {{ display: grid; grid-template-columns: repeat(4, 1fr);
                 gap: 16px; margin-bottom: 40px; }}
  .stat-card {{ background: white; border-radius: 10px; padding: 20px 24px;
                border: 1px solid #e2e8f0;
                box-shadow: 0 1px 4px rgba(0,0,0,0.06); }}
  .stat-val  {{ font-size: 2em; font-weight: 700; color: #0f172a; }}
  .stat-lbl  {{ font-size: 0.8em; color: #64748b; margin-top: 4px;
                text-transform: uppercase; letter-spacing: 0.05em; }}
  .section   {{ background: white; border-radius: 12px; padding: 32px;
                border: 1px solid #e2e8f0; margin-bottom: 32px;
                box-shadow: 0 1px 4px rgba(0,0,0,0.06); }}
  .section h2 {{ font-size: 1.15em; font-weight: 700; color: #0f172a;
                 margin-bottom: 8px; }}
  .section .caption {{ font-size: 0.82em; color: #64748b; margin-bottom: 20px;
                       font-style: italic; }}
  .fig-wrap  {{ text-align: center; }}
  .fig-wrap img {{ max-width: 100%; border-radius: 8px;
                   border: 1px solid #f1f5f9; }}
  table      {{ width: 100%; border-collapse: collapse; font-size: 0.88em; }}
  th         {{ background: #f1f5f9; padding: 10px 14px; text-align: left;
                font-weight: 600; font-size: 0.8em; text-transform: uppercase;
                letter-spacing: 0.05em; color: #475569; }}
  td         {{ padding: 9px 14px; border-bottom: 1px solid #f1f5f9; }}
  tr:hover td {{ background: #f8fafc; }}
  .badge     {{ display: inline-block; padding: 2px 10px; border-radius: 4px;
                font-size: 0.78em; font-weight: 600; }}
  .footer    {{ text-align: center; padding: 24px; font-size: 0.75em;
                color: #94a3b8; border-top: 1px solid #e2e8f0; }}
</style>
</head>
<body>

<div class="header">
  <h1>ARG × MGE Co-occurrence Report</h1>
  <p>MetagenAMR Pipeline · Contig-level co-localization · mobileOG-db × CARD · {datetime.now().strftime("%Y-%m-%d")}</p>
</div>

<div class="container">

  <!-- Stats -->
  <div class="stats-grid">
    <div class="stat-card">
      <div class="stat-val">{n_contigs}</div>
      <div class="stat-lbl">Co-carrying contigs</div>
    </div>
    <div class="stat-card">
      <div class="stat-val">{n_args}</div>
      <div class="stat-lbl">Unique ARGs</div>
    </div>
    <div class="stat-card">
      <div class="stat-val">{n_mge_types}</div>
      <div class="stat-lbl">Unique MGE hits</div>
    </div>
    <div class="stat-card">
      <div class="stat-val">{n_samples}</div>
      <div class="stat-lbl">Samples</div>
    </div>
  </div>

  <!-- Fig 1 : Bubble plot -->
  <div class="section">
    <h2>Figure 1 — Bubble Plot: Top ARGs × MGE Functional Category</h2>
    <p class="caption">Top 20 ARGs by co-occurrence frequency across all samples.
    Bubble size and label indicate the number of contigs co-carrying both an ARG and an MGE.
    All four hospital wastewater samples aggregated.</p>
    <div class="fig-wrap">
      <img src="data:image/png;base64,{b64_fig1}" alt="Bubble plot ARG x MGE">
    </div>
  </div>

  <!-- Fig 2 : Network -->
  <div class="section">
    <h2>Figure 2 — Co-occurrence Network: ARG ↔ MGE Category</h2>
    <p class="caption">Bipartite network linking the top 25 ARGs (right) to MGE functional
    categories (left). Edge width is proportional to the number of contigs on which both
    elements co-occur. Node size reflects total co-occurrence weight.</p>
    <div class="fig-wrap">
      <img src="data:image/png;base64,{b64_fig2}" alt="Network ARG MGE">
    </div>
  </div>

  <!-- Top associations table -->
  <div class="section">
    <h2>Top 10 ARG–MGE Associations</h2>
    <p class="caption">Ranked by number of unique contigs co-carrying both elements (all samples).</p>
    <table>
      <thead>
        <tr>
          <th>ARG</th>
          <th>MGE Category</th>
          <th style="text-align:right">Co-carrying contigs</th>
        </tr>
      </thead>
      <tbody>{assoc_rows}</tbody>
    </table>
  </div>

</div>

<div class="footer">
  MetagenAMR Pipeline · RGI main (CARD) × mobileOG-db beatrix-1.6 ·
  Contig-level co-localization · {datetime.now().strftime("%Y-%m-%d %H:%M")}
</div>

</body>
</html>"""

with open(output_html, 'w') as f:
    f.write(html)

print(f"\n[OK] Rapport : {output_html}")
print(f"     Taille  : {os.path.getsize(output_html)/1024:.0f} KB")

