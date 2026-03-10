#!/usr/bin/env python3
"""
Generate taxonomic analysis HTML report.
Sections:
  1. Summary
  2. Composition — Genus level (stacked barplot)
  3. Composition — Species level (stacked barplot)
  4. Alpha Diversity (multi-panel barplot)
  5. Heatmap — Top N species (clustered)
  6. Priority Pathogens (grouped barplot)
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
import base64
import os
import sys
from io import BytesIO
from datetime import datetime

# ── Snakemake params ──────────────────────────────────────────────────────────
genus_matrix_path   = snakemake.input.genus_matrix
species_matrix_path = snakemake.input.species_matrix
alpha_div_path      = snakemake.input.alpha_div
out_report          = snakemake.output.report

top_genera  = snakemake.params.get("top_genera",  15)
top_species = snakemake.params.get("top_species", 20)
top_heatmap = snakemake.params.get("top_heatmap", 30)
dpi         = snakemake.params.get("dpi",         150)   # lower for HTML embedding
pathogens   = snakemake.params.get("pathogens",   [])

# ── Logging ───────────────────────────────────────────────────────────────────
log = open(snakemake.log[0], 'w')
sys.stderr = sys.stdout = log

report_date = datetime.now().strftime("%Y-%m-%d %H:%M")
print("=" * 70)
print(f"Taxonomy HTML Report — {report_date}")
print("=" * 70)

os.makedirs(os.path.dirname(out_report), exist_ok=True)

# ── Color palettes ────────────────────────────────────────────────────────────
PALETTE = [
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
    "#8CD17D", "#B6992D", "#499894", "#86BCB6", "#D37295",
    "#FABFD2", "#A0CBE8", "#FFBE7D", "#F1CE63", "#D4E6F1",
    "#C7C7C7"
]
SAMPLE_COLORS = ["#4E79A7", "#F28E2B", "#59A14F", "#E15759",
                 "#B07AA1", "#76B7B2", "#EDC948", "#9C755F"]

PRIORITY_PATHOGENS = [
    "Escherichia coli", "Klebsiella pneumoniae",
    "Pseudomonas aeruginosa", "Clostridioides difficile",
    "Clostridium perfringens", "Acinetobacter baumannii",
    "Enterococcus faecium", "Bacteroides fragilis"
]
BENEFICIAL_SPECIES = [
    "Faecalibacterium prausnitzii",
    "Agathobacter rectalis",
    "Roseburia intestinalis"
]

# ── Helper functions ──────────────────────────────────────────────────────────

def fig_to_base64(fig):
    buf = BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight',
                facecolor='white')
    buf.seek(0)
    encoded = base64.b64encode(buf.read()).decode('utf-8')
    plt.close(fig)
    return f"data:image/png;base64,{encoded}"


def load_matrix(path):
    df = pd.read_csv(path, sep='\t', index_col=0)
    if 'taxonomy_id' in df.columns:
        df = df.drop(columns=['taxonomy_id'])
    return df


def build_top_n_df(df, top_n):
    mean_abund = df.mean(axis=1)
    top_taxa   = mean_abund.nlargest(top_n).index.tolist()
    other_taxa = [t for t in df.index if t not in top_taxa]
    plot_df    = df.loc[top_taxa].copy()
    if other_taxa:
        others_row      = df.loc[other_taxa].sum(axis=0)
        others_row.name = "Others"
        plot_df = pd.concat([plot_df, others_row.to_frame().T])
    return plot_df.T   # samples × taxa


def style_ax(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.grid(True, linestyle='--', alpha=0.35, color='grey')
    ax.set_axisbelow(True)


# ── FIGURE 1 & 2 : Stacked barplots ──────────────────────────────────────────

def make_stacked_barplot(df_matrix, top_n, title):
    plot_df = build_top_n_df(df_matrix, top_n)
    taxa    = plot_df.columns.tolist()
    colors  = PALETTE[:len(taxa) - 1] + ["#C7C7C7"]

    fig, ax = plt.subplots(figsize=(10, 7))
    bottom  = np.zeros(len(plot_df))

    for i, taxon in enumerate(taxa):
        values = plot_df[taxon].values
        ax.bar(plot_df.index, values, bottom=bottom,
               color=colors[i], label=taxon,
               width=0.60, edgecolor='white', linewidth=0.5)
        bottom += values

    ax.set_xlabel("Sample",                 fontsize=12, labelpad=10)
    ax.set_ylabel("Relative Abundance (%)", fontsize=12, labelpad=10)
    ax.set_title(title,                     fontsize=14, fontweight='bold', pad=15)
    ax.set_ylim(0, 108)
    ax.set_xticks(range(len(plot_df)))
    ax.set_xticklabels(plot_df.index, fontsize=11)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:.0f}%"))
    style_ax(ax)

    handles = [mpatches.Patch(color=colors[i], label=t)
               for i, t in enumerate(taxa)]
    legend = ax.legend(
        handles=handles, title="Taxa", title_fontsize=10,
        fontsize=8, loc='upper left', bbox_to_anchor=(1.01, 1),
        borderaxespad=0, frameon=True, framealpha=0.9, edgecolor='grey'
    )
    for text in legend.get_texts():
        if not text.get_text().startswith("Others"):
            text.set_style('italic')

    plt.tight_layout()
    return fig


# ── FIGURE 3 : Alpha diversity ────────────────────────────────────────────────

def make_alpha_diversity(alpha_df):
    metrics = [
        {"col": "shannon_diversity",  "title": "Shannon Diversity (H')",
         "ylabel": "H'",                 "color": "#4E79A7", "fmt": ".3f"},
        {"col": "simpson_diversity",  "title": "Simpson Diversity (1-D)",
         "ylabel": "1-D",                "color": "#F28E2B", "fmt": ".3f"},
        {"col": "observed_richness",  "title": "Observed Richness (S)",
         "ylabel": "Number of species",  "color": "#59A14F", "fmt": ".0f"},
        {"col": "pielou_evenness",    "title": "Pielou's Evenness (J')",
         "ylabel": "J'",                 "color": "#E15759", "fmt": ".3f"},
    ]

    fig = plt.figure(figsize=(13, 9))
    fig.suptitle("Alpha Diversity Metrics — Hospital Wastewater",
                 fontsize=15, fontweight='bold', y=1.01)
    gs = gridspec.GridSpec(2, 2, hspace=0.50, wspace=0.38)

    for idx, m in enumerate(metrics):
        ax      = fig.add_subplot(gs[idx // 2, idx % 2])
        values  = alpha_df[m["col"]].values
        samples = alpha_df['sample'].values
        ylim_max = max(values) * 1.25

        bars = ax.bar(samples, values, color=m["color"], alpha=0.85,
                      width=0.55, edgecolor='white', linewidth=0.8)

        for bar, val in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + ylim_max * 0.02,
                    format(val, m["fmt"]),
                    ha='center', va='bottom', fontsize=9, fontweight='bold')

        ax.set_title(m["title"],   fontsize=12, fontweight='bold', pad=10)
        ax.set_ylabel(m["ylabel"], fontsize=10)
        ax.set_xlabel("Sample",    fontsize=10)
        ax.set_ylim(0, ylim_max)
        ax.set_xticks(range(len(samples)))
        ax.set_xticklabels(samples, fontsize=9, rotation=15, ha='right')
        style_ax(ax)

    plt.tight_layout()
    return fig


# ── FIGURE 4 : Heatmap ────────────────────────────────────────────────────────

def make_heatmap(species_df, top_n):
    mean_abund  = species_df.mean(axis=1)
    top_taxa    = mean_abund.nlargest(top_n).index.tolist()
    heatmap_df  = species_df.loc[top_taxa]
    log_df      = np.log10(heatmap_df + 0.01)

    row_link    = linkage(pdist(log_df.values, metric='euclidean'), method='ward')
    row_order   = leaves_list(row_link)
    log_ordered = log_df.iloc[row_order]
    taxa_ordered = [top_taxa[i] for i in row_order]

    fig_h = max(10, top_n * 0.38)
    fig, ax = plt.subplots(figsize=(8, fig_h))

    img = ax.imshow(log_ordered.values, aspect='auto',
                    cmap='YlOrRd', interpolation='nearest')
    cbar = fig.colorbar(img, ax=ax, shrink=0.35, pad=0.02)
    cbar.set_label('log₁₀(Relative abundance % + 0.01)', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    ax.set_xticks(range(len(log_ordered.columns)))
    ax.set_xticklabels(log_ordered.columns, fontsize=11)
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')

    ax.set_yticks(range(len(taxa_ordered)))
    ytick_labels = ax.set_yticklabels(taxa_ordered, fontsize=8, style='italic')
    for label, species in zip(ytick_labels, taxa_ordered):
        if any(p.lower() in species.lower() for p in PRIORITY_PATHOGENS):
            label.set_color('#CC0000')
            label.set_fontweight('bold')
        elif any(b.lower() in species.lower() for b in BENEFICIAL_SPECIES):
            label.set_color('#1a7a1a')
            label.set_fontweight('bold')

    ax.set_xticks(np.arange(-0.5, len(log_ordered.columns), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(taxa_ordered), 1),        minor=True)
    ax.grid(which='minor', color='white', linewidth=0.6)
    ax.tick_params(which='minor', bottom=False, left=False)
    ax.set_title("Top Species Abundance Heatmap",
                 fontsize=13, fontweight='bold', pad=40)

    legend_elements = [
        Line2D([0], [0], color='#CC0000', lw=2.5, label='Priority pathogen'),
        Line2D([0], [0], color='#1a7a1a', lw=2.5, label='Beneficial species'),
        Line2D([0], [0], color='black',   lw=2.5, label='Commensal'),
    ]
    ax.legend(handles=legend_elements, loc='lower right',
              bbox_to_anchor=(1.0, -0.06), fontsize=8,
              frameon=True, framealpha=0.9)

    plt.tight_layout()
    return fig


# ── FIGURE 5 : Priority pathogens ─────────────────────────────────────────────

def make_pathogens(species_df, pathogens):
    found   = [p for p in pathogens if p in species_df.index]
    missing = [p for p in pathogens if p not in species_df.index]

    if missing:
        print(f"  [!!] Pathogens below filter threshold (not plotted):")
        for p in missing: print(f"       - {p}")

    if not found:
        print("  [!!] No priority pathogens found in matrix — skipping figure")
        return None, []

    path_df   = species_df.loc[found]
    samples   = species_df.columns.tolist()
    n_samples = len(samples)
    n_path    = len(found)
    bar_width = 0.75 / n_samples
    offsets   = np.linspace(
        -(n_samples - 1) / 2,
         (n_samples - 1) / 2,
        n_samples
    ) * bar_width
    x = np.arange(n_path)

    fig, ax = plt.subplots(figsize=(max(10, n_path * 1.6), 7))

    for i, (sample, color) in enumerate(zip(samples, SAMPLE_COLORS[:n_samples])):
        values = path_df[sample].values
        bars   = ax.bar(x + offsets[i], values, width=bar_width,
                        color=color, alpha=0.85, label=sample,
                        edgecolor='white', linewidth=0.5)
        for bar, val in zip(bars, values):
            if val > 0.1:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.08,
                        f"{val:.1f}%", ha='center', va='bottom',
                        fontsize=7, rotation=90)

    max_val = path_df.values.max()
    ax.axhline(y=1.0, color='red', linestyle='--', alpha=0.5, linewidth=1.2)
    ax.text(n_path - 0.5, 1.05, 'Alert threshold (1%)',
            color='red', fontsize=8, alpha=0.7, ha='right')

    ax.set_xlabel("Priority Pathogens",                    fontsize=12, labelpad=10)
    ax.set_ylabel("Relative Abundance (%)\n(Gb-normalized)", fontsize=12, labelpad=10)
    ax.set_title("Priority Pathogen Abundance",            fontsize=14,
                 fontweight='bold', pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(found, fontsize=9, style='italic', rotation=30, ha='right')
    ax.set_ylim(0, max(max_val * 1.3, 2.0))
    ax.legend(title="Sample", fontsize=10, title_fontsize=10,
              loc='upper right', framealpha=0.9)
    style_ax(ax)
    plt.tight_layout()
    return fig, found


# ── Load data ─────────────────────────────────────────────────────────────────

print("\n[1/5] Loading data...")
genus_df   = load_matrix(genus_matrix_path)
species_df = load_matrix(species_matrix_path)
alpha_df   = pd.read_csv(alpha_div_path, sep='\t')

samples    = list(genus_df.columns)
n_samples  = len(samples)
n_genera   = genus_df.shape[0]
n_species  = species_df.shape[0]

print(f"  [OK] Genus matrix  : {n_genera} genera  x {n_samples} samples")
print(f"  [OK] Species matrix: {n_species} species x {n_samples} samples")
print(f"  [OK] Alpha div     : {len(alpha_df)} samples")
print(f"  [OK] Samples       : {', '.join(samples)}")

# ── Generate figures ──────────────────────────────────────────────────────────

print("\n[2/5] Generating figures...")

print("  -> Fig 1: Genus composition")
fig1    = make_stacked_barplot(genus_df, top_genera,
                               "Bacterial Composition — Genus Level")
img_f1  = fig_to_base64(fig1)

print("  -> Fig 2: Species composition")
fig2    = make_stacked_barplot(species_df, top_species,
                               "Bacterial Composition — Species Level")
img_f2  = fig_to_base64(fig2)

print("  -> Fig 3: Alpha diversity")
required_cols = ['sample', 'shannon_diversity', 'simpson_diversity',
                 'observed_richness', 'pielou_evenness']
missing_cols  = [c for c in required_cols if c not in alpha_df.columns]
if missing_cols:
    print(f"  [!!] Missing alpha columns: {missing_cols}")
    sys.exit(1)
fig3   = make_alpha_diversity(alpha_df)
img_f3 = fig_to_base64(fig3)

print("  -> Fig 4: Species heatmap")
fig4   = make_heatmap(species_df, top_heatmap)
img_f4 = fig_to_base64(fig4)

print("  -> Fig 5: Priority pathogens")
fig5_result = make_pathogens(species_df, pathogens)
fig5, found_pathogens = fig5_result if fig5_result[0] is not None else (None, [])
img_f5 = fig_to_base64(fig5) if fig5 is not None else None


# ── Summary stats ─────────────────────────────────────────────────────────────

print("\n[3/5] Computing summary stats...")

# Top 5 genera across all samples
top5_genera = genus_df.mean(axis=1).nlargest(5)

# Top 5 species across all samples
top5_species = species_df.mean(axis=1).nlargest(5)

# Alpha diversity summary
shannon_mean = alpha_df['shannon_diversity'].mean()
shannon_std  = alpha_df['shannon_diversity'].std()
richness_mean = alpha_df['observed_richness'].mean()

# Pathogen detection summary
pathogen_detection = {}
for p in pathogens:
    if p in species_df.index:
        row = species_df.loc[p]
        pathogen_detection[p] = {
            'mean': row.mean(),
            'max':  row.max(),
            'detected_in': int((row > 0).sum())
        }


# ── Build HTML ────────────────────────────────────────────────────────────────

print("\n[4/5] Building HTML report...")

def stat_card(value, label, color="#4E79A7", sub=None):
    sub_html = f'<div style="font-size:0.78em;color:#888;margin-top:3px">{sub}</div>' if sub else ""
    return f"""
    <div style="background:#fff;border-radius:10px;padding:18px 22px;
                box-shadow:0 2px 8px rgba(0,0,0,0.08);text-align:center;
                border-left:4px solid {color};min-width:130px">
      <div style="font-size:1.9em;font-weight:700;color:{color}">{value}</div>
      <div style="font-size:0.85em;color:#555;margin-top:4px">{label}</div>
      {sub_html}
    </div>"""


def section(title, content, icon=""):
    return f"""
    <div class="section" id="sec-{title.lower().replace(' ','-')}">
      <h2>{icon} {title}</h2>
      {content}
    </div>"""


def figure_block(img_b64, caption, fig_id=""):
    return f"""
    <div class="fig-block" id="{fig_id}">
      <img src="{img_b64}" alt="{caption}" style="max-width:100%;border-radius:8px;
           box-shadow:0 2px 12px rgba(0,0,0,0.10)">
      <div class="fig-caption">{caption}</div>
    </div>"""


# ── Navigation ────────────────────────────────────────────────────────────────

nav_items = [
    ("Summary",           "#sec-summary"),
    ("Genus Composition", "#sec-genus-composition"),
    ("Species Composition","#sec-species-composition"),
    ("Alpha Diversity",   "#sec-alpha-diversity"),
    ("Heatmap",           "#sec-heatmap"),
    ("Pathogens",         "#sec-pathogens"),
]
nav_html = "\n".join(
    f'<a href="{href}">{label}</a>'
    for label, href in nav_items
)


# ── Section: Summary ──────────────────────────────────────────────────────────

sample_badges = " ".join(
    f'<span style="background:{SAMPLE_COLORS[i % len(SAMPLE_COLORS)]};'
    f'color:white;padding:4px 12px;border-radius:20px;'
    f'font-size:0.9em;font-weight:600">{s}</span>'
    for i, s in enumerate(samples)
)

cards_html = f"""
<div style="display:flex;flex-wrap:wrap;gap:16px;margin:20px 0">
  {stat_card(n_samples,  "Samples",           "#4E79A7")}
  {stat_card(n_genera,   "Genera detected",   "#F28E2B")}
  {stat_card(n_species,  "Species detected",  "#59A14F")}
  {stat_card(f"{shannon_mean:.2f} ± {shannon_std:.2f}", "Shannon H' (mean ± sd)", "#E15759")}
  {stat_card(f"{richness_mean:.0f}", "Mean richness (S)",  "#B07AA1")}
  {stat_card(len(found_pathogens), "Pathogens detected", "#CC0000")}
</div>
<div style="margin:12px 0">
  <strong>Samples:</strong>&nbsp;&nbsp;{sample_badges}
</div>
<div style="margin:18px 0">
  <strong>Top 5 genera (mean relative abundance):</strong>
  <table class="data-table" style="margin-top:10px">
    <thead><tr><th>Genus</th>{"".join(f"<th>{s}</th>" for s in samples)}<th>Mean (%)</th></tr></thead>
    <tbody>
      {"".join(
          f'<tr><td style="font-style:italic">{g}</td>'
          + "".join(f'<td>{genus_df.loc[g, s]:.2f}%</td>'
                    if g in genus_df.index and s in genus_df.columns else '<td>—</td>'
                    for s in samples)
          + f'<td><strong>{v:.2f}%</strong></td></tr>'
          for g, v in top5_genera.items()
      )}
    </tbody>
  </table>
</div>
<div style="margin:18px 0">
  <strong>Top 5 species (mean relative abundance):</strong>
  <table class="data-table" style="margin-top:10px">
    <thead><tr><th>Species</th>{"".join(f"<th>{s}</th>" for s in samples)}<th>Mean (%)</th></tr></thead>
    <tbody>
      {"".join(
          f'<tr><td style="font-style:italic">{sp}</td>'
          + "".join(f'<td>{species_df.loc[sp, s]:.2f}%</td>'
                    if sp in species_df.index and s in species_df.columns else '<td>—</td>'
                    for s in samples)
          + f'<td><strong>{v:.2f}%</strong></td></tr>'
          for sp, v in top5_species.items()
      )}
    </tbody>
  </table>
</div>
"""

sec_summary = section("Summary", cards_html, "📊")


# ── Section: Figures ──────────────────────────────────────────────────────────

sec_genus = section(
    "Genus Composition",
    figure_block(img_f1,
        f"Fig. 1 — Stacked barplot of relative abundances at genus level. "
        f"Top {top_genera} genera shown; remaining taxa aggregated as 'Others'.",
        "fig1"),
    "🦠"
)

sec_species = section(
    "Species Composition",
    figure_block(img_f2,
        f"Fig. 2 — Stacked barplot of relative abundances at species level. "
        f"Top {top_species} species shown.",
        "fig2"),
    "🔬"
)

sec_alpha = section(
    "Alpha Diversity",
    f"""
    <div style="background:#f8f9ff;border-left:4px solid #4E79A7;
                padding:12px 16px;border-radius:6px;margin-bottom:18px;font-size:0.93em">
      <strong>Metrics:</strong>
      Shannon H' (information-theoretic diversity),
      Simpson 1-D (dominance-corrected),
      Observed richness S (number of detected species),
      Pielou's J' (evenness relative to maximum possible diversity).
    </div>
    """ + figure_block(img_f3,
        "Fig. 3 — Alpha diversity metrics across samples. "
        "Higher Shannon and Pielou values indicate more diverse and evenly distributed communities.",
        "fig3"),
    "📈"
)

sec_heatmap = section(
    "Heatmap",
    f"""
    <div style="background:#f8f9ff;border-left:4px solid #F28E2B;
                padding:12px 16px;border-radius:6px;margin-bottom:18px;font-size:0.93em">
      Hierarchical clustering (Ward's method, Euclidean distance on log₁₀-transformed abundances).
      <span style="color:#CC0000;font-weight:600">Red</span> = priority pathogen,
      <span style="color:#1a7a1a;font-weight:600">Green</span> = beneficial species,
      Black = commensal.
    </div>
    """ + figure_block(img_f4,
        f"Fig. 4 — Clustered heatmap of top {top_heatmap} species by mean abundance. "
        f"Values are log₁₀-transformed relative abundances.",
        "fig4"),
    "🗺️"
)

if img_f5:
    pathogen_table_rows = ""
    for p, stats in pathogen_detection.items():
        alert = "⚠️" if stats['max'] > 1.0 else ""
        pathogen_table_rows += (
            f'<tr><td style="font-style:italic">{p}</td>'
            f'<td>{stats["mean"]:.3f}%</td>'
            f'<td>{stats["max"]:.3f}% {alert}</td>'
            f'<td>{stats["detected_in"]}/{n_samples}</td></tr>'
        )
    sec_pathogens = section(
        "Pathogens",
        f"""
        <div style="background:#fff5f5;border-left:4px solid #CC0000;
                    padding:12px 16px;border-radius:6px;margin-bottom:18px;font-size:0.93em">
          ⚠️ Alert threshold set at <strong>1%</strong> relative abundance.
          Species are from the WHO priority pathogen list and local hospital surveillance targets.
        </div>
        <table class="data-table" style="margin-bottom:24px">
          <thead>
            <tr><th>Pathogen</th><th>Mean abundance</th>
                <th>Max abundance</th><th>Detected in</th></tr>
          </thead>
          <tbody>{pathogen_table_rows}</tbody>
        </table>
        """ + figure_block(img_f5,
            "Fig. 5 — Priority pathogen abundances per sample. "
            "Red dashed line = alert threshold (1%). "
            "Values above threshold warrant clinical review.",
            "fig5"),
        "⚠️"
    )
else:
    sec_pathogens = section(
        "Pathogens",
        '<div style="color:#888;padding:20px">No priority pathogens detected above filter threshold.</div>',
        "⚠️"
    )


# ── Assemble HTML ─────────────────────────────────────────────────────────────

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Taxonomy Report — MetagenAMR</title>
  <style>
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
      background: #f4f6f9;
      color: #2c3e50;
      line-height: 1.6;
    }}
    /* ── Navbar ── */
    .navbar {{
      position: sticky; top: 0; z-index: 1000;
      background: #1a2634;
      padding: 0 32px;
      display: flex; align-items: center; justify-content: space-between;
      box-shadow: 0 2px 10px rgba(0,0,0,0.2);
      height: 56px;
    }}
    .navbar .brand {{
      color: white; font-size: 1.1em; font-weight: 700;
      letter-spacing: 0.5px;
    }}
    .navbar .brand span {{ color: #76B7B2; }}
    .navbar nav {{ display: flex; gap: 4px; }}
    .navbar nav a {{
      color: #cdd9e5; text-decoration: none;
      padding: 6px 14px; border-radius: 6px;
      font-size: 0.88em; font-weight: 500;
      transition: background 0.2s, color 0.2s;
    }}
    .navbar nav a:hover {{ background: #2d3f52; color: white; }}
    /* ── Main layout ── */
    .container {{
      max-width: 1200px; margin: 0 auto;
      padding: 32px 24px;
    }}
    /* ── Header ── */
    .report-header {{
      background: linear-gradient(135deg, #1a2634 0%, #2d4a6e 100%);
      color: white; padding: 36px 40px; border-radius: 14px;
      margin-bottom: 28px;
      box-shadow: 0 4px 20px rgba(0,0,0,0.15);
    }}
    .report-header h1 {{ font-size: 1.9em; font-weight: 700; margin-bottom: 6px; }}
    .report-header .subtitle {{
      color: #a8c4d8; font-size: 0.95em; margin-bottom: 16px;
    }}
    .report-header .meta {{
      display: flex; gap: 24px; flex-wrap: wrap; margin-top: 12px;
    }}
    .report-header .meta span {{
      background: rgba(255,255,255,0.1); padding: 4px 14px;
      border-radius: 20px; font-size: 0.85em;
    }}
    /* ── Sections ── */
    .section {{
      background: white; border-radius: 12px;
      padding: 28px 32px; margin-bottom: 24px;
      box-shadow: 0 2px 10px rgba(0,0,0,0.06);
    }}
    .section h2 {{
      font-size: 1.25em; font-weight: 700;
      color: #1a2634; margin-bottom: 20px;
      padding-bottom: 10px;
      border-bottom: 2px solid #e8ecf0;
    }}
    /* ── Tables ── */
    .data-table {{
      width: 100%; border-collapse: collapse;
      font-size: 0.88em;
    }}
    .data-table th {{
      background: #1a2634; color: white;
      padding: 9px 14px; text-align: left;
      font-weight: 600; font-size: 0.85em;
    }}
    .data-table td {{
      padding: 8px 14px;
      border-bottom: 1px solid #e8ecf0;
      color: #444;
    }}
    .data-table tr:hover td {{ background: #f8fafc; }}
    .data-table tr:last-child td {{ border-bottom: none; }}
    /* ── Figures ── */
    .fig-block {{
      margin: 16px 0;
      text-align: center;
    }}
    .fig-caption {{
      font-size: 0.83em; color: #666;
      margin-top: 10px; font-style: italic;
      max-width: 800px; margin-left: auto; margin-right: auto;
    }}
    /* ── Footer ── */
    .footer {{
      text-align: center; padding: 24px;
      color: #888; font-size: 0.82em;
    }}
    @media (max-width: 768px) {{
      .navbar nav {{ display: none; }}
      .container {{ padding: 16px; }}
      .report-header {{ padding: 24px; }}
    }}
  </style>
</head>
<body>

<div class="navbar">
  <div class="brand">MetagenAMR <span>| Taxonomy</span></div>
  <nav>{nav_html}</nav>
</div>

<div class="container">

  <div class="report-header">
    <h1>🦠 Taxonomic Analysis Report</h1>
    <div class="subtitle">Hospital Wastewater Metagenomics — Nanopore Long Reads</div>
    <div class="meta">
      <span>📅 {report_date}</span>
      <span>🧪 {n_samples} samples</span>
      <span>🔬 {n_species} species detected</span>
      <span>🌿 {n_genera} genera detected</span>
    </div>
  </div>

  {sec_summary}
  {sec_genus}
  {sec_species}
  {sec_alpha}
  {sec_heatmap}
  {sec_pathogens}

  <div class="footer">
    Generated by MetagenAMR pipeline — {report_date}<br>
    Kraken2 / Bracken taxonomic profiling | Nanopore long reads
  </div>

</div>
</body>
</html>"""

# ── Write output ──────────────────────────────────────────────────────────────

print("\n[5/5] Writing HTML report...")
with open(out_report, 'w', encoding='utf-8') as f:
    f.write(html)

size_kb = os.path.getsize(out_report) / 1024
print(f"  [OK] Report saved : {out_report}  ({size_kb:.0f} KB)")
print(f"  [OK] Figures      : 5 (embedded as base64 PNG)")
print(f"  [OK] Sections     : Summary | Genus | Species | Alpha | Heatmap | Pathogens")

print("\n" + "=" * 70)
print("TAXONOMY HTML REPORT COMPLETED")
print("=" * 70)

sink_close = log.close()
