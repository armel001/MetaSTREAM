#!/usr/bin/env python3
"""
plot_mge_report.py — MGE HTML Report
Inputs  : 07_mge_normalized.tsv + 01_mge_hits_raw.tsv (sorties mge_analysis.R)
Figures : 4 (composition relative, counts normalisés, heatmap contigs, overview)
Nomenclature : catégories natives mobileOG-db (integration/excision, phage, etc.)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
import base64, os, sys
from io import BytesIO
from datetime import datetime

# ── Snakemake params ──────────────────────────────────────────────────────────
normalized_path = snakemake.input.normalized       # 07_mge_normalized.tsv
hits_raw_path   = snakemake.input.hits_raw         # 01_mge_hits_raw.tsv
out_report      = snakemake.output.report
dpi             = snakemake.params.get("dpi", 150)

log = open(snakemake.log[0], 'w')
sys.stderr = sys.stdout = log

report_date = datetime.now().strftime("%Y-%m-%d %H:%M")
print("=" * 70)
print(f"MGE HTML Report — {report_date}")
print("=" * 70)

os.makedirs(os.path.dirname(out_report), exist_ok=True)

# ── Catégories natives mobileOG-db ────────────────────────────────────────────
# Ordre fixé par abondance décroissante typique (eaux usées hospitalières)
MGE_CATEGORIES = [
    "integration/excision",
    "transfer",
    "phage",
    "replication/recombination/repair",
    "stability/transfer/defense",
]

MGE_LABELS = {
    "integration/excision":             "Integration/Excision",
    "transfer":                         "Transfer",
    "phage":                            "Phage",
    "replication/recombination/repair": "Replication/Recombination/Repair",
    "stability/transfer/defense":       "Stability/Transfer/Defense",
}

MGE_COLORS = {
    "integration/excision":             "#4E79A7",
    "transfer":                         "#F28E2B",
    "phage":                            "#59A14F",
    "replication/recombination/repair": "#E15759",
    "stability/transfer/defense":       "#B07AA1",
}

SAMPLE_COLORS = ["#4E79A7", "#F28E2B", "#59A14F", "#E15759",
                 "#B07AA1", "#76B7B2", "#EDC948", "#9C755F"]

# ── Helpers ───────────────────────────────────────────────────────────────────

def fig_to_b64(fig):
    buf = BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight', facecolor='white')
    buf.seek(0)
    enc = base64.b64encode(buf.read()).decode('utf-8')
    plt.close(fig)
    return f"data:image/png;base64,{enc}"

def style_ax(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.grid(True, linestyle='--', alpha=0.35, color='grey')
    ax.set_axisbelow(True)

# ── Chargement des données ────────────────────────────────────────────────────

print("\n[1/5] Loading data...")

norm = pd.read_csv(normalized_path, sep='\t')
hits = pd.read_csv(hits_raw_path,   sep='\t')

samples   = sorted(norm['sample_id'].unique().tolist())
n_samples = len(samples)

print(f"  [OK] normalized  : {norm.shape[0]} lignes")
print(f"  [OK] hits_raw    : {hits.shape[0]} hits")
print(f"  [OK] Samples     : {', '.join(samples)}")

# Vérifier catégories présentes
cats_present = norm['major_category'].unique().tolist()
MGE_CATEGORIES = [c for c in MGE_CATEGORIES if c in cats_present]
# Ajouter catégories non prévues
for c in cats_present:
    if c not in MGE_CATEGORIES:
        MGE_CATEGORIES.append(c)
        MGE_LABELS[c]  = c.title()
        MGE_COLORS[c]  = "#999999"

# Stats globales
total_hits_all  = int(norm['total_hits'].sum())
total_contigs   = int(norm['n_contigs'].sum())
total_orfs      = int(hits['orf_global_id'].nunique()) if 'orf_global_id' in hits.columns else 0
dominant_cat    = norm.groupby('major_category')['total_hits'].sum().idxmax()
dominant_pct    = norm.groupby('major_category')['total_hits'].sum().max() / total_hits_all * 100

print(f"  [OK] Total hits    : {total_hits_all:,}")
print(f"  [OK] Total contigs : {total_contigs:,}")
print(f"  [OK] Dominant cat  : {dominant_cat} ({dominant_pct:.1f}%)")

# ── Matrices pour figures ─────────────────────────────────────────────────────

# Matrice relative (%) : samples × catégories
pct_matrix = norm.pivot_table(
    index='sample_id', columns='major_category',
    values='relative_abundance', fill_value=0
).reindex(columns=MGE_CATEGORIES, fill_value=0)

# Matrice normalisée (copies/Gb) : samples × catégories
norm_matrix = norm.pivot_table(
    index='sample_id', columns='major_category',
    values='normalized_copies', fill_value=0
).reindex(columns=MGE_CATEGORIES, fill_value=0)

# ── FIGURE 1 : Composition relative ──────────────────────────────────────────

print("\n[2/5] Generating figures...")
print("  -> Fig 1: MGE relative composition (%)")

fig1, ax = plt.subplots(figsize=(10, 6))
bottom = np.zeros(n_samples)
sample_list = pct_matrix.index.tolist()

for cat in MGE_CATEGORIES:
    vals = pct_matrix[cat].values
    ax.bar(sample_list, vals, bottom=bottom,
           color=MGE_COLORS.get(cat, "#999"),
           label=MGE_LABELS.get(cat, cat),
           width=0.55, edgecolor='white', linewidth=0.5)
    bottom += vals

ax.set_xlabel("Sample",                    fontsize=12, labelpad=10)
ax.set_ylabel("Relative composition (%)",  fontsize=12, labelpad=10)
ax.set_title("MGE Composition per Sample",
             fontsize=14, fontweight='bold', pad=15)
ax.set_ylim(0, 115)
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:.0f}%"))
handles = [mpatches.Patch(color=MGE_COLORS.get(c,"#999"),
           label=MGE_LABELS.get(c,c)) for c in MGE_CATEGORIES]
ax.legend(handles=handles, title="MGE category", title_fontsize=10,
          fontsize=9, loc='upper left', bbox_to_anchor=(1.01, 1),
          frameon=True, framealpha=0.9, edgecolor='grey')
style_ax(ax)
plt.tight_layout()
img_f1 = fig_to_b64(fig1)

# ── FIGURE 2 : Copies / Gb par catégorie ─────────────────────────────────────

print("  -> Fig 2: Normalized counts (copies/Gb)")

n_cats    = len(MGE_CATEGORIES)
bar_width = 0.75 / n_samples
offsets   = np.linspace(-(n_samples-1)/2, (n_samples-1)/2, n_samples) * bar_width
x         = np.arange(n_cats)

fig2, ax = plt.subplots(figsize=(12, 6))

for i, (sid, color) in enumerate(zip(samples, SAMPLE_COLORS[:n_samples])):
    if sid not in norm_matrix.index:
        continue
    vals = [norm_matrix.loc[sid, c] if c in norm_matrix.columns else 0
            for c in MGE_CATEGORIES]
    bars = ax.bar(x + offsets[i], vals, width=bar_width,
                  color=color, alpha=0.85, label=sid,
                  edgecolor='white', linewidth=0.5)
    for bar, val in zip(bars, vals):
        if val > 100:
            ax.text(bar.get_x() + bar.get_width()/2,
                    bar.get_height() + 20,
                    f"{val:,.0f}", ha='center', va='bottom',
                    fontsize=7, rotation=90)

ax.set_xlabel("MGE category",        fontsize=12, labelpad=10)
ax.set_ylabel("Normalized copies/Gb", fontsize=12, labelpad=10)
ax.set_title("MGE Abundance — Normalized (copies/Gb)",
             fontsize=14, fontweight='bold', pad=15)
ax.set_xticks(x)
xlabels = [MGE_LABELS.get(c, c) for c in MGE_CATEGORIES]
ax.set_xticklabels(xlabels, fontsize=9, rotation=20, ha='right')
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x):,}"))
ax.legend(title="Sample", fontsize=10, title_fontsize=10,
          loc='upper right', framealpha=0.9)
style_ax(ax)
plt.tight_layout()
img_f2 = fig_to_b64(fig2)

# ── FIGURE 3 : Heatmap contigs porteurs ──────────────────────────────────────

print("  -> Fig 3: Per-contig heatmap (top contigs by total hits)")

TOP_N = 40

# Agréger hits par contig (tous samples) × catégorie
contig_cat = (
    hits.groupby(['contig', 'major_category'])
        .size()
        .reset_index(name='hits')
)
contig_pivot = contig_cat.pivot_table(
    index='contig', columns='major_category', values='hits', fill_value=0
)
contig_pivot = contig_pivot.reindex(
    columns=[c for c in MGE_CATEGORIES if c in contig_pivot.columns],
    fill_value=0
)
contig_pivot['_total'] = contig_pivot.sum(axis=1)
top = contig_pivot.nlargest(TOP_N, '_total').drop(columns='_total')

if len(top) >= 2:
    log_mat   = np.log10(top.values + 0.1)
    row_link  = linkage(pdist(log_mat, metric='euclidean'), method='ward')
    row_order = leaves_list(row_link)
    top = top.iloc[row_order]

fig_h = max(8, len(top) * 0.30)
fig3, ax = plt.subplots(figsize=(8, fig_h))

img = ax.imshow(np.log10(top.values + 0.1),
                aspect='auto', cmap='YlOrRd', interpolation='nearest')
cbar = fig3.colorbar(img, ax=ax, shrink=0.35, pad=0.02)
cbar.set_label('log₁₀(hits + 0.1)', fontsize=9)
cbar.ax.tick_params(labelsize=8)

col_labels = [MGE_LABELS.get(c, c) for c in top.columns]
ax.set_xticks(range(len(top.columns)))
ax.set_xticklabels(col_labels, fontsize=9, rotation=25, ha='left')
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')

ax.set_yticks(range(len(top)))
ax.set_yticklabels(top.index, fontsize=7)
ax.set_xticks(np.arange(-0.5, len(top.columns), 1), minor=True)
ax.set_yticks(np.arange(-0.5, len(top), 1), minor=True)
ax.grid(which='minor', color='white', linewidth=0.5)
ax.tick_params(which='minor', bottom=False, left=False)
ax.set_title(f"Top {len(top)} Contigs — MGE Hit Distribution",
             fontsize=12, fontweight='bold', pad=45)
plt.tight_layout()
img_f3 = fig_to_b64(fig3)

# ── FIGURE 4 : Overview — hits totaux + ORFs uniques ─────────────────────────

print("  -> Fig 4: Per-sample overview")

total_per_sample = norm.groupby('sample_id')['total_hits'].sum().reindex(samples)
orfs_per_sample  = (
    hits.groupby('sample_id')['orf_global_id'].nunique().reindex(samples)
    if 'orf_global_id' in hits.columns
    else pd.Series(0, index=samples)
)

fig4, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, (values, title, ylabel, fmt) in zip(axes, [
    (total_per_sample, "Total MGE Hits per Sample",    "Total hits",     "{:,.0f}"),
    (orfs_per_sample,  "Unique MGE ORFs per Sample",   "Unique ORFs",    "{:,.0f}"),
]):
    colors = SAMPLE_COLORS[:n_samples]
    bars   = ax.bar(samples, values.values, color=colors,
                    width=0.55, edgecolor='white', linewidth=0.5)
    for bar, val in zip(bars, values.values):
        ax.text(bar.get_x() + bar.get_width()/2,
                bar.get_height() * 1.02,
                fmt.format(val),
                ha='center', va='bottom', fontsize=9, fontweight='600')
    ax.set_title(title, fontsize=13, fontweight='bold', pad=12)
    ax.set_xlabel("Sample",  fontsize=11)
    ax.set_ylabel(ylabel,    fontsize=11)
    ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, _: f"{int(x):,}"))
    style_ax(ax)

plt.tight_layout()
img_f4 = fig_to_b64(fig4)

# ── HTML helpers ──────────────────────────────────────────────────────────────

def figure_block(b64, caption, fid=""):
    return (f'<div class="fig-block" id="{fid}">'
            f'<img src="{b64}" style="max-width:100%;border-radius:8px;" alt="{caption}">'
            f'<div class="fig-caption">{caption}</div></div>')

def section(title, content, icon="📊"):
    return (f'<div class="section" id="sec-{title.lower().replace(" ","-")}">'
            f'<h2>{icon} {title}</h2>{content}</div>')

sections_nav = ["Summary","Composition","Counts","Heatmap","Overview","Statistics"]
nav_html = "".join(
    f'<a href="#sec-{s.lower()}">{s}</a>' for s in sections_nav)

# ── Section Summary ───────────────────────────────────────────────────────────

top5_cat = (norm.groupby('major_category')['total_hits']
                .sum().sort_values(ascending=False).head(5))

stat_cards = "".join(f"""
  <div style="background:#f0f4f8;border-radius:10px;padding:18px 22px;
              min-width:140px;flex:1;text-align:center">
    <div style="font-size:1.8em;font-weight:700;color:#1a2634">{v}</div>
    <div style="font-size:0.82em;color:#666;margin-top:4px">{l}</div>
  </div>"""
for v, l in [
    (f"{n_samples}",           "Samples analysés"),
    (f"{total_hits_all:,}",    "Total hits MGE"),
    (f"{total_contigs:,}",     "Contigs avec MGE"),
    (f"{total_orfs:,}",        "ORFs uniques"),
    (MGE_LABELS.get(dominant_cat, dominant_cat), f"Catégorie dominante<br>Moy. {dominant_pct:.1f}%"),
])

top5_html = "".join(
    f'<li style="margin:4px 0"><strong>{MGE_LABELS.get(c,c)}</strong>'
    f' — {v:,} hits ({v/total_hits_all*100:.1f}%)</li>'
    for c, v in top5_cat.items()
)

sec_summary = section("Summary", f"""
  <div style="display:flex;gap:16px;flex-wrap:wrap;margin-bottom:20px">
    {stat_cards}
  </div>
  <div style="display:grid;grid-template-columns:1fr 1fr;gap:20px;margin-top:8px">
    <div>
      <div style="font-weight:700;margin-bottom:8px;color:#1a2634">
        Top catégories MGE</div>
      <ul style="list-style:none;padding:0;font-size:0.9em">{top5_html}</ul>
    </div>
    <div style="background:#f8f9ff;border-left:4px solid #4E79A7;
                padding:12px 16px;border-radius:6px;font-size:0.88em">
      <strong>Note méthodologique.</strong><br>
      Détection par <em>mobileOG-db beatrix-1.6</em> (alignement DIAMOND) sur
      contigs assemblés par Flye à partir de lectures Nanopore.
      Catégories fonctionnelles natives mobileOG-db :
      Integration/Excision, Transfer, Phage,
      Replication/Recombination/Repair, Stability/Transfer/Defense.
      Filtres appliqués : Pident ≥ 90%, couverture query ≥ 90%.
    </div>
  </div>
""", "📊")

# ── Section Composition ───────────────────────────────────────────────────────

sec_composition = section("Composition", f"""
  <div style="font-size:0.9em;color:#555;margin-bottom:16px">
    Proportion relative (%) de chaque catégorie MGE par échantillon,
    calculée sur le total des hits filtrés (Pident ≥ 90%, cov ≥ 90%).
  </div>
  {figure_block(img_f1,
    "Fig. 1 — Composition relative du mobilome par échantillon. "
    "Chaque barre représente 100% des hits MGE, décomposés par catégorie fonctionnelle mobileOG-db.",
    "fig1")}
""", "🧬")

# ── Section Counts normalisés ─────────────────────────────────────────────────

sec_counts = section("Counts", f"""
  <div style="font-size:0.9em;color:#555;margin-bottom:16px">
    Abondance normalisée (copies/Gb) par catégorie et par échantillon.
    La normalisation par gigabase séquencée permet la comparaison inter-samples
    indépendamment de la profondeur de séquençage.
  </div>
  {figure_block(img_f2,
    "Fig. 2 — Abondance normalisée des MGE (copies/Gb) par catégorie fonctionnelle. "
    "Normalisation : nombre de hits / taille du dataset en Gb.",
    "fig2")}
""", "📊")

# ── Section Heatmap ───────────────────────────────────────────────────────────

sec_heatmap = section("Heatmap", f"""
  <div style="background:#f8f9ff;border-left:4px solid #F28E2B;
              padding:12px 16px;border-radius:6px;margin-bottom:18px;font-size:0.92em">
    Clustering hiérarchique (méthode Ward, distance euclidienne sur log₁₀(hits + 0.1)).
    Seuls les <strong>top {TOP_N} contigs</strong> par nombre total de hits sont affichés.
    Les contigs porteurs de plusieurs catégories (éléments composites) sont visibles
    par leur distribution multi-colonnes.
  </div>
  {figure_block(img_f3,
    f"Fig. 3 — Heatmap des {len(top)} contigs présentant le plus grand nombre de hits MGE "
    "(agrégés sur l'ensemble des échantillons, échelle log₁₀). "
    "Les contigs avec signal dans plusieurs catégories représentent des éléments composites potentiels.",
    "fig3")}
""", "🗺️")

# ── Section Overview ──────────────────────────────────────────────────────────

sec_overview = section("Overview", f"""
  {figure_block(img_f4,
    "Fig. 4 — Vue d'ensemble par échantillon : hits MGE totaux (gauche) "
    "et ORFs uniques porteurs de MGE (droite). "
    "Les valeurs absolues reflètent à la fois la charge MGE et la profondeur d'assemblage.",
    "fig4")}
""", "📈")

# ── Section Statistics table ──────────────────────────────────────────────────

# Construire la table depuis norm (07_mge_normalized.tsv)
norm_wide = norm.pivot_table(
    index='sample_id', columns='major_category',
    values='total_hits', fill_value=0
).reindex(columns=MGE_CATEGORIES, fill_value=0)

pct_wide = norm.pivot_table(
    index='sample_id', columns='major_category',
    values='relative_abundance', fill_value=0
).reindex(columns=MGE_CATEGORIES, fill_value=0)

contig_wide = norm.pivot_table(
    index='sample_id', columns='major_category',
    values='n_contigs', fill_value=0
).reindex(columns=MGE_CATEGORIES, fill_value=0)

hits_per_sample  = norm.groupby('sample_id')['total_hits'].sum()
contig_per_sample = norm.groupby('sample_id')['n_contigs'].sum()
orfs_per_sample2 = (
    hits.groupby('sample_id')['orf_global_id'].nunique()
    if 'orf_global_id' in hits.columns
    else pd.Series(0, index=samples)
)

header = ("<tr><th>Sample</th>"
          + "".join(f"<th>{MGE_LABELS.get(c,c)}</th>" for c in MGE_CATEGORIES)
          + "<th>Total hits</th><th>Contigs MGE+</th><th>ORFs uniques</th></tr>")

rows = ""
for sid in samples:
    cells = "".join(
        f'<td>{int(norm_wide.loc[sid, c]):,}'
        f'<span style="color:#888;font-size:0.85em"> '
        f'({pct_wide.loc[sid, c]:.1f}%)</span></td>'
        for c in MGE_CATEGORIES
    )
    rows += (f'<tr><td><strong>{sid}</strong></td>{cells}'
             f'<td><strong>{int(hits_per_sample[sid]):,}</strong></td>'
             f'<td>{int(contig_per_sample[sid]):,}</td>'
             f'<td>{int(orfs_per_sample2.get(sid, 0)):,}</td></tr>')

# Ligne totaux/moyennes
sum_row = "<tr style='background:#f0f4f8;font-weight:600'><td>Total / Mean</td>"
for c in MGE_CATEGORIES:
    t = int(norm_wide[c].sum())
    m = pct_wide[c].mean()
    sum_row += f'<td>{t:,} <span style="color:#888;font-size:0.85em">({m:.1f}%)</span></td>'
sum_row += (f'<td>{total_hits_all:,}</td>'
            f'<td>{total_contigs:,}</td>'
            f'<td>{total_orfs:,}</td></tr>')

sec_statistics = section("Statistics", f"""
  <div style="overflow-x:auto">
    <table class="data-table">
      <thead>{header}</thead>
      <tbody>{rows}{sum_row}</tbody>
    </table>
  </div>
  <div style="margin-top:14px;font-size:0.83em;color:#666;font-style:italic">
    Les pourcentages indiquent la proportion relative de chaque catégorie
    sur le total des hits de l'échantillon.
    Abondances normalisées disponibles dans 07_mge_normalized.tsv.
  </div>
""", "📋")

# ── Assemble HTML ─────────────────────────────────────────────────────────────

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>MGE Report — MetagenAMR</title>
  <style>
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
      background: #f4f6f9; color: #2c3e50; line-height: 1.6;
    }}
    .navbar {{
      position: sticky; top: 0; z-index: 1000;
      background: #1a2634; padding: 0 32px;
      display: flex; align-items: center; justify-content: space-between;
      box-shadow: 0 2px 10px rgba(0,0,0,0.2); height: 56px;
    }}
    .navbar .brand {{ color: white; font-size: 1.1em; font-weight: 700; }}
    .navbar .brand span {{ color: #F28E2B; }}
    .navbar nav {{ display: flex; gap: 4px; }}
    .navbar nav a {{
      color: #cdd9e5; text-decoration: none; padding: 6px 14px;
      border-radius: 6px; font-size: 0.88em; font-weight: 500;
      transition: background 0.2s, color 0.2s;
    }}
    .navbar nav a:hover {{ background: #2d3f52; color: white; }}
    .container {{ max-width: 1200px; margin: 0 auto; padding: 32px 24px; }}
    .report-header {{
      background: linear-gradient(135deg, #1a2634 0%, #2d4a1a 100%);
      color: white; padding: 36px 40px; border-radius: 14px;
      margin-bottom: 28px; box-shadow: 0 4px 20px rgba(0,0,0,0.15);
    }}
    .report-header h1 {{ font-size: 1.9em; font-weight: 700; margin-bottom: 6px; }}
    .report-header .subtitle {{ color: #b8d4b0; font-size: 0.95em; margin-bottom: 16px; }}
    .report-header .meta {{ display: flex; gap: 24px; flex-wrap: wrap; margin-top: 12px; }}
    .report-header .meta span {{
      background: rgba(255,255,255,0.1); padding: 4px 14px;
      border-radius: 20px; font-size: 0.85em;
    }}
    .section {{
      background: white; border-radius: 12px;
      padding: 28px 32px; margin-bottom: 24px;
      box-shadow: 0 2px 10px rgba(0,0,0,0.06);
    }}
    .section h2 {{
      font-size: 1.25em; font-weight: 700; color: #1a2634;
      margin-bottom: 20px; padding-bottom: 10px;
      border-bottom: 2px solid #e8ecf0;
    }}
    .data-table {{ width: 100%; border-collapse: collapse; font-size: 0.87em; }}
    .data-table th {{
      background: #1a2634; color: white;
      padding: 9px 14px; text-align: left;
      font-weight: 600; font-size: 0.84em;
    }}
    .data-table td {{ padding: 8px 14px; border-bottom: 1px solid #e8ecf0; color: #444; }}
    .data-table tr:hover td {{ background: #f8fafc; }}
    .data-table tr:last-child td {{ border-bottom: none; }}
    .fig-block {{ margin: 16px 0; text-align: center; }}
    .fig-caption {{
      font-size: 0.83em; color: #666; margin-top: 10px;
      font-style: italic; max-width: 900px;
      margin-left: auto; margin-right: auto;
    }}
    .footer {{ text-align: center; padding: 24px; color: #888; font-size: 0.82em; }}
    @media (max-width: 768px) {{
      .navbar nav {{ display: none; }}
      .container {{ padding: 16px; }}
      .report-header {{ padding: 24px; }}
    }}
  </style>
</head>
<body>
<div class="navbar">
  <div class="brand">MetagenAMR <span>| MGE</span></div>
  <nav>{nav_html}</nav>
</div>
<div class="container">
  <div class="report-header">
    <h1>🧬 Mobile Genetic Elements Report</h1>
    <div class="subtitle">Hospital Wastewater Metagenomics — mobileOG-db beatrix-1.6</div>
    <div class="meta">
      <span>📅 {report_date}</span>
      <span>🧪 {n_samples} samples</span>
      <span>🎯 {total_hits_all:,} MGE hits</span>
      <span>🔗 {total_contigs:,} contigs porteurs</span>
      <span>🧩 {total_orfs:,} ORFs uniques</span>
    </div>
  </div>
  {sec_summary}
  {sec_composition}
  {sec_counts}
  {sec_heatmap}
  {sec_overview}
  {sec_statistics}
  <div class="footer">
    Generated by MetagenAMR pipeline — {report_date}<br>
    mobileOG-db beatrix-1.6 | DIAMOND alignment | Flye assembly | Nanopore long reads
  </div>
</div>
</body>
</html>"""

print("\n[5/5] Writing HTML report...")
with open(out_report, 'w', encoding='utf-8') as f:
    f.write(html)

size_kb = os.path.getsize(out_report) / 1024
print(f"  [OK] {out_report}  ({size_kb:.0f} KB)")
print("\n" + "=" * 70)
print("MGE HTML REPORT COMPLETED")
print("=" * 70)
log.close()
