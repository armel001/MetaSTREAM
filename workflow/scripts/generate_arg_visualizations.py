#!/usr/bin/env python3
"""
Generate ARG visualization report (HTML) for rgi main or rgi_bwt
Sections: Summary | Drug Classes | Diversity | Heatmap | Mechanisms | Rarefaction | PCoA | Top ARGs | Samples
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import MDS
from sklearn.metrics import pairwise_distances
from scipy.stats import entropy
from scipy.cluster.hierarchy import dendrogram, linkage
import base64
import sys
from io import BytesIO
from pathlib import Path
from datetime import datetime

# ─── Logging ─────────────────────────────────────────────────────────────────

log = open(snakemake.log[0], 'w')
sys.stderr = sys.stdout = log

print("=" * 70)
print("ARG Visualization Report")
print("=" * 70)

# ─── Parameters ──────────────────────────────────────────────────────────────

tool         = snakemake.params.tool
min_identity = snakemake.params.min_identity
min_coverage = snakemake.params.min_coverage
tool_label   = "RGI Main (Assembly)" if tool == "rgi" else "RGI BWT (Reads)"
report_date  = datetime.now().strftime("%Y-%m-%d %H:%M")

print(f"\n  Tool    : {tool_label}")
print(f"  Date    : {report_date}")

# ─── Style ───────────────────────────────────────────────────────────────────

sns.set_style("ticks")
master_palette = sns.color_palette("tab20", 25)
plt.rcParams.update({
    'figure.dpi':        300,
    'savefig.dpi':       300,
    'font.size':         12,
    'axes.labelweight':  'bold',
    'axes.titleweight':  'bold',
    'figure.facecolor':  'white',
    'axes.facecolor':    'white',
})


# ─── Charger le mapping manuel si disponible ─────────────────────────────────
ARO_MAPPING_FILE = Path("config/aro_short_names.tsv")
arg_norm_full = pd.read_csv(snakemake.input.arg_normalized, sep='\t')
if ARO_MAPPING_FILE.exists():
    aro_mapping = pd.read_csv(ARO_MAPPING_FILE, sep='\t',
                              index_col='full_name')['short_name'].to_dict()
    print(f"  ✓ ARO mapping loaded: {len(aro_mapping)} entries")
else:
    aro_mapping = {}
    print(f"  ⚠ No ARO mapping file found — using auto short_name()")

def short_name(full_name):
    # 1. Mapping manuel en priorité
    if full_name in aro_mapping:
        return aro_mapping[full_name]

    # 2. Nom déjà court (pas d'espace)
    tokens = full_name.split()
    if len(tokens) == 1:
        return full_name

    # 3. rRNA genes
    if "rRNA" in tokens:
        rrna_idx = tokens.index("rRNA")
        rrna_type = tokens[rrna_idx - 1]
        genus   = tokens[0][0]
        species = tokens[1][:8]
        return f"{genus}.{species}_{rrna_type}"

    # 4. Gène en position 3
    if len(tokens) >= 3:
        candidate = tokens[2]
        if candidate not in SKIP_TOKENS and not candidate[0].islower():
            return candidate

    # 5. Avant "conferring" ou "with"
    for keyword in ["conferring", "with"]:
        if keyword in tokens:
            idx = tokens.index(keyword)
            if idx > 0 and tokens[idx - 1] not in SKIP_TOKENS:
                return tokens[idx - 1]

    # Fallback
    if len(tokens) >= 2:
        return f"{tokens[0][0]}.{tokens[1][:8]}"
    return full_name[:20]

# ─── Helpers ─────────────────────────────────────────────────────────────────

def fig_to_b64(fig):
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()


def img_tag(b64, alt="figure"):
    return (f'<div class="fig-container">'
            f'<img src="data:image/png;base64,{b64}" alt="{alt}">'
            f'<a class="dl-btn" href="data:image/png;base64,{b64}" '
            f'download="{alt}.png">⬇ Download PNG</a>'
            f'</div>')


def section_header(title, icon=""):
    return f'<div class="section-header">{icon} {title}</div>'


SKIP_TOKENS = {"with", "conferring", "mutation", "rRNA", "16S", "23S",
               "resistance", "to", "gene", "variant", "and"}


def short_name(full_name):
    tokens = full_name.split()
    if len(tokens) == 1:
        return full_name
    if len(tokens) >= 3:
        candidate = tokens[2]
        if candidate not in SKIP_TOKENS and not candidate[0].islower():
            return candidate
    for keyword in ["conferring", "with"]:
        if keyword in tokens:
            idx = tokens.index(keyword)
            if idx > 0 and tokens[idx - 1] not in SKIP_TOKENS:
                return tokens[idx - 1]
    return full_name[:20]


def make_unique_labels(columns):
    seen = {}
    result = {}
    for col in columns:
        label = short_name(col)
        if label in seen:
            seen[label] += 1
            result[col] = f"{label}_{seen[label]}"
        else:
            seen[label] = 0
            result[col] = label
    return result

# ─── Load data ───────────────────────────────────────────────────────────────

print("\n[1/8] Loading data...")

arg_matrix_norm  = pd.read_csv(snakemake.input.arg_matrix_normalized, sep='\t', index_col=0)
arg_matrix_counts = pd.read_csv(snakemake.input.arg_matrix_counts,    sep='\t', index_col=0)
arg_matrix_pres  = pd.read_csv(snakemake.input.arg_matrix_presence,   sep='\t', index_col=0)
drug_class_df    = pd.read_csv(snakemake.input.drug_class,             sep='\t')
mechanism_df     = pd.read_csv(snakemake.input.mechanism,              sep='\t')
family_df        = pd.read_csv(snakemake.input.family,                 sep='\t')
arg_counts_df    = pd.read_csv(snakemake.input.arg_counts,             sep='\t')

n_samples = len(arg_matrix_norm)
n_args    = len(arg_matrix_norm.columns)
n_classes = drug_class_df['drug_class'].nunique() if 'drug_class' in drug_class_df.columns else 0

print(f"  ✓ Matrix: {n_samples} samples × {n_args} ARGs")
print(f"  ✓ Drug classes: {n_classes}")

# ─── Figure 1: Drug Classes Distribution ─────────────────────────────────────

print("\n[2/8] Figure 1: Drug Classes Distribution...")

# Séparer les drug classes combinées si pas déjà fait en amont
drug_class_clean = (drug_class_df
                    .assign(drug_class=drug_class_df['drug_class'].str.split('; '))
                    .explode('drug_class')
                    .assign(drug_class=lambda x: x['drug_class'].str.strip()))

# Réagréger après split
drug_class_clean = (drug_class_clean
                    .groupby(['sample_id', 'drug_class'])['count']
                    .sum()
                    .reset_index())
drug_class_clean['total'] = drug_class_clean.groupby('sample_id')['count'].transform('sum')
drug_class_clean['relative_abundance'] = (drug_class_clean['count'] /
                                           drug_class_clean['total'] * 100)

# Garder uniquement les top N classes pour lisibilité
TOP_N = 15
top_classes = (drug_class_clean
               .groupby('drug_class')['count']
               .sum()
               .nlargest(TOP_N)
               .index)

drug_class_plot = drug_class_clean[drug_class_clean['drug_class'].isin(top_classes)].copy()

# Recalculer relative_abundance sur les top N uniquement
other = (drug_class_clean[~drug_class_clean['drug_class'].isin(top_classes)]
         .groupby('sample_id')['count']
         .sum()
         .reset_index()
         .assign(drug_class='Other'))
other['total'] = drug_class_clean.groupby('sample_id')['count'].sum().values[0]
other['relative_abundance'] = other['count'] / other['total'] * 100

drug_class_plot = pd.concat([drug_class_plot, other], ignore_index=True)

drug_pivot = (drug_class_plot
              .pivot(index='sample_id', columns='drug_class', values='relative_abundance')
              .fillna(0))

# Trier les colonnes par abondance moyenne décroissante
col_order = drug_pivot.mean().sort_values(ascending=False).index
drug_pivot = drug_pivot[col_order]

fig, ax = plt.subplots(figsize=(12, 6))
drug_pivot.plot(kind='bar', stacked=True, ax=ax,
                color=master_palette[:len(drug_pivot.columns)])

ax.set_title(f"Antibiotic Resistance Classes Distribution — {tool_label}",
             fontsize=14, pad=15, weight='bold')
ax.set_ylabel("Relative Abundance (%)", fontsize=12, weight='bold')
ax.set_xlabel("")
ax.set_ylim(0, 105)
ax.tick_params(axis='x', rotation=0, labelsize=11)
ax.tick_params(axis='y', labelsize=11)
ax.margins(x=0.01)
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left',
          title='Drug Class', frameon=False, fontsize=9,
          ncol=1)
sns.despine()
plt.tight_layout()
b64_fig1 = fig_to_b64(fig)
plt.close(fig)
print("  ✓ Done")


# ─── Figure 2: Alpha Diversity ───────────────────────────────────────────────

print("\n[3/8] Figure 2: Alpha Diversity...")

def shannon_diversity(row):
    row = row[row > 0]
    if len(row) == 0:
        return 0
    p = row / row.sum()
    return entropy(p, base=np.e)

shannon = arg_matrix_norm.apply(shannon_diversity, axis=1)
richness = (arg_matrix_norm > 0).sum(axis=1)

diversity_df = pd.DataFrame({
    'sample_id': shannon.index,
    'shannon':   shannon.values,
    'richness':  richness.values
})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

for ax, col, title, ylabel, color in [
    (ax1, 'shannon',  "Alpha Diversity: Shannon Index", "Shannon Index",        "#2e6da4"),
    (ax2, 'richness', "ARG Richness",                   "Number of Unique ARGs", "#e07b39")
]:
    values = diversity_df[col].values
    samples = diversity_df['sample_id'].values
    mean_val = values.mean()

    bars = ax.bar(samples, values,
                  color=color, edgecolor='black',
                  linewidth=0.7, alpha=0.85, width=0.5)

    # Valeur sur chaque barre
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(values) * 0.01,
                f'{val:.2f}' if col == 'shannon' else f'{int(val)}',
                ha='center', va='bottom', fontsize=10, fontweight='bold')

    # Ligne moyenne
    ax.axhline(mean_val, color='gray', linestyle='--',
               linewidth=1.5, alpha=0.6,
               label=f'Mean = {mean_val:.2f}' if col == 'shannon'
                     else f'Mean = {int(mean_val)}')

    ax.set_title(title, fontsize=13, pad=12, weight='bold')
    ax.set_ylabel(ylabel, fontsize=11, weight='bold')
    ax.set_xlabel("")
    ax.tick_params(axis='x', rotation=30, labelsize=10)
    ax.tick_params(axis='y', labelsize=10)
    ax.set_ylim(0, max(values) * 1.15)
    ax.legend(frameon=False, fontsize=9)
    ax.margins(x=0.1)

sns.despine(fig=fig)
plt.suptitle(f"ARG Alpha Diversity — {tool_label}",
             fontsize=14, weight='bold', y=1.02)
plt.tight_layout()
b64_fig2 = fig_to_b64(fig)
plt.close(fig)
print("  ✓ Done")


# ─── Figure 3: Heatmap Top 30 ────────────────────────────────────────────────

print("\n[4/8] Figure 3: Heatmap Top 30 ARGs...")

# Top 30 par somme d'abondance
top_args = arg_matrix_norm.sum(axis=0).nlargest(30).index
arg_top  = arg_matrix_norm[top_args].copy()

# Noms courts
unique_labels = make_unique_labels(arg_top.columns)
arg_top = arg_top.rename(columns=unique_labels)

# Transpose : ARGs en lignes, samples en colonnes
df_heatmap = arg_top.T

# Z-score par ligne (row-wise) — met en évidence quel sample
# a le plus d'un ARG donné, indépendamment des valeurs absolues
g = sns.clustermap(
    df_heatmap,
    cmap="mako",
    z_score=0,                          # normalisation par ligne
    figsize=(max(8, len(arg_top) * 1.2), 14),
    linewidths=0.4,
    linecolor='black',
    cbar_pos=(0.02, 0.85, 0.03, 0.1),
    dendrogram_ratio=(0.1, 0.15),
    xticklabels=True,
    yticklabels=True,
)

g.ax_heatmap.set_yticklabels(
    g.ax_heatmap.get_ymajorticklabels(),
    fontsize=9, rotation=0
)
g.ax_heatmap.set_xticklabels(
    g.ax_heatmap.get_xmajorticklabels(),
    fontsize=11, weight='bold', rotation=45
)
g.ax_heatmap.set_xlabel("Samples", fontsize=12, weight='bold', labelpad=10)
g.ax_heatmap.set_ylabel("Top 30 ARGs", fontsize=12, weight='bold', labelpad=10)

g.fig.suptitle(f"Top 30 ARGs Heatmap (Row Z-score) — {tool_label}",
               fontsize=14, weight='bold', y=1.02)

b64_fig3 = fig_to_b64(g.fig)
plt.close(g.fig)
print("  ✓ Done")


# ─── Figure 4: Resistance Mechanisms ─────────────────────────────────────────

print("\n[5/8] Figure 4: Resistance Mechanisms...")

# Séparer les mécanismes combinés (même logique que drug_class)
mech_expanded = (mechanism_df
                 .assign(resistance_mechanism=mechanism_df['resistance_mechanism']
                         .str.split('; '))
                 .explode('resistance_mechanism')
                 .assign(resistance_mechanism=lambda x: x['resistance_mechanism'].str.strip()))

mech_total = (mech_expanded
              .groupby('resistance_mechanism')['count']
              .sum()
              .sort_values(ascending=True))

# Wrapper les labels longs
import textwrap
wrapped_labels = [textwrap.fill(label, width=35) for label in mech_total.index]

fig, ax = plt.subplots(figsize=(12, max(5, len(mech_total) * 0.55 + 1)))

# Palette avec contraste suffisant
palette = sns.color_palette("Set2", len(mech_total))
bars = ax.barh(range(len(mech_total)), mech_total.values,
               color=palette, edgecolor='black', linewidth=0.7, height=0.6)

# Valeur sur chaque barre
for bar, val in zip(bars, mech_total.values):
    ax.text(bar.get_width() + mech_total.max() * 0.01,
            bar.get_y() + bar.get_height() / 2,
            f'{int(val):,}', va='center', ha='left', fontsize=9)

ax.set_yticks(range(len(mech_total)))
ax.set_yticklabels(wrapped_labels, fontsize=10)
ax.set_title(f"Resistance Mechanisms — {tool_label}",
             fontsize=14, pad=15, weight='bold')
ax.set_xlabel("Number of ARG Observations", fontsize=12, weight='bold')
ax.margins(y=0.02)
ax.set_xlim(0, mech_total.max() * 1.12)  # espace pour les labels
sns.despine()
plt.tight_layout()
b64_fig4 = fig_to_b64(fig)
plt.close(fig)
print("  ✓ Done")

# ─── Figure 5: Rarefaction Curves ────────────────────────────────────────────

print("\n[6/8] Figure 5: Rarefaction Curves...")

def rarefaction_curve(row, step=10):
    total = int(row.sum())
    if total == 0:
        return [0], [0]
    if total < step:
        return [total], [int((row > 0).sum())]
    counts = row[row > 0].values
    observations = np.repeat(np.arange(len(counts)), counts.astype(int))
    steps = list(range(step, total, step))
    if not steps or steps[-1] != total:
        steps.append(total)
    richness_vals = []
    for n in steps:
        np.random.seed(42)
        sub = np.random.choice(observations, size=n, replace=False)
        richness_vals.append(len(np.unique(sub)))
    return steps, richness_vals

fig, ax = plt.subplots(figsize=(10, 7))
for idx, (sample, row) in enumerate(arg_matrix_counts.iterrows()):
    try:
        steps, rich = rarefaction_curve(row, step=10)
        if len(steps) > 1:
            ax.plot(steps, rich, label=sample,
                    color=master_palette[idx % len(master_palette)],
                    linewidth=2.5, alpha=0.8)
        else:
            ax.scatter(steps, rich, label=sample,
                       color=master_palette[idx % len(master_palette)],
                       s=100, zorder=3)
    except Exception as e:
        print(f"  ⚠ Rarefaction skipped for {sample}: {e}")

ax.set_title(f"ARG Rarefaction Curves — {tool_label}",
             fontsize=14, pad=15, weight='bold')
ax.set_xlabel("Number of ARG Observations", fontsize=12, weight='bold')
ax.set_ylabel("Number of Unique ARGs",       fontsize=12, weight='bold')
ax.tick_params(axis='both', labelsize=11)
ax.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False, fontsize=9)
sns.despine()
plt.tight_layout()
b64_fig5 = fig_to_b64(fig)
plt.close(fig)
print("  ✓ Done")

# ─── Figure 6: PCoA Beta Diversity ───────────────────────────────────────────

print("\n[7/8] Figure 6: PCoA Beta Diversity...")

def bray_curtis(u, v):
    denom = np.sum(u + v)
    return np.sum(np.abs(u - v)) / denom if denom > 0 else 0

distances = pairwise_distances(arg_matrix_norm.values, metric=bray_curtis)
mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
coords = mds.fit_transform(distances)

fig, ax = plt.subplots(figsize=(9, 7))
for idx, sample in enumerate(arg_matrix_norm.index):
    ax.scatter(coords[idx, 0], coords[idx, 1], s=150,
               c=[master_palette[idx % len(master_palette)]],
               edgecolors='black', linewidth=1.5, alpha=0.8, label=sample)
    ax.annotate(sample, (coords[idx, 0], coords[idx, 1]),
                xytext=(8, 8), textcoords='offset points',
                fontsize=9, weight='bold')
ax.set_title(f"PCoA — Bray-Curtis Distance\n{tool_label}",
             fontsize=14, pad=15, weight='bold')
ax.set_xlabel("PCo1", fontsize=12, weight='bold')
ax.set_ylabel("PCo2", fontsize=12, weight='bold')
ax.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', frameon=False, fontsize=9)
sns.despine()
plt.tight_layout()
b64_fig6 = fig_to_b64(fig)
plt.close(fig)
print("  ✓ Done")

# ─── Figure 7: Top 20 ARGs ────────────────────────────────────────────────────

print("\n[8/8] Figure 7: Top 20 ARGs...")

top20 = (arg_counts_df
         .groupby('best_hit_aro')['count']
         .sum()
         .nlargest(20)
         .sort_values(ascending=True))

fig, ax = plt.subplots(figsize=(10, 8))
colors = [master_palette[i % len(master_palette)] for i in range(len(top20))]
ax.barh(top20.index, top20.values, color=colors, edgecolor='black', linewidth=0.7)
ax.set_title(f"Top 20 Most Abundant ARGs — {tool_label}",
             fontsize=14, pad=15, weight='bold')
ax.set_xlabel("Number of Observations", fontsize=12, weight='bold')
ax.tick_params(axis='both', labelsize=10)
ax.margins(y=0.01)
sns.despine()
plt.tight_layout()
b64_fig7 = fig_to_b64(fig)
plt.close(fig)
print("  ✓ Done")



# ─── Figure 8: Bubble plot abondance × prévalence ────────────────────────────

print("\n[8/8] Figure 8: Bubble plot Family Abundance × Prevalence...")

family_prefix_df = pd.read_csv(snakemake.input.family_prefix, sep='\t')

# Agréger par famille
bubble_df = (family_prefix_df
             .groupby('family_prefix')
             .agg(
                 mean_abundance  = ('normalized_abundance', 'mean'),
                 total_abundance = ('normalized_abundance', 'sum'),
                 prevalence      = ('sample_id', 'nunique'),
                 n_variants      = ('n_variants', 'sum')
             )
             .reset_index())

n_samples_total = family_prefix_df['sample_id'].nunique()

# Normaliser taille des bulles
size_scale = (bubble_df['total_abundance'] /
              bubble_df['total_abundance'].max() * 800 + 50)

fig, ax = plt.subplots(figsize=(12, 8))

scatter = ax.scatter(
    bubble_df['mean_abundance'],
    bubble_df['prevalence'],
    s=size_scale,
    c=bubble_df['mean_abundance'],
    cmap='YlOrRd',
    alpha=0.7,
    edgecolors='black',
    linewidths=0.6
)

# Colorbar
cbar = plt.colorbar(scatter, ax=ax, shrink=0.6)
cbar.set_label('Mean Normalized Abundance (reads/Gb)', fontsize=10)

# Annoter les familles prioritaires
# Priorité = abondant (top 30%) ET prévalent (> 50% des samples)
abund_thresh = bubble_df['mean_abundance'].quantile(0.70)
prev_thresh  = n_samples_total * 0.5

priority = bubble_df[
    (bubble_df['mean_abundance'] >= abund_thresh) |
    (bubble_df['prevalence'] >= prev_thresh)
]

for _, row in priority.iterrows():
    ax.annotate(
        row['family_prefix'],
        (row['mean_abundance'], row['prevalence']),
        xytext=(8, 4), textcoords='offset points',
        fontsize=8, fontweight='bold',
        arrowprops=dict(arrowstyle='-', color='gray', lw=0.5)
    )

# Ligne de prévalence 100%
ax.axhline(n_samples_total, color='red', linestyle='--',
           linewidth=1, alpha=0.5, label=f'All samples (n={n_samples_total})')

# Quadrants
ax.axvline(bubble_df['mean_abundance'].median(),
           color='gray', linestyle=':', linewidth=1, alpha=0.4)
ax.axhline(n_samples_total / 2,
           color='gray', linestyle=':', linewidth=1, alpha=0.4)

ax.set_xlabel("Mean Normalized Abundance (reads/Gb)", fontsize=12, weight='bold')
ax.set_ylabel("Prevalence (number of samples)", fontsize=12, weight='bold')
ax.set_title(f"ARG Family Abundance × Prevalence\n{tool_label}",
             fontsize=14, pad=15, weight='bold')
ax.set_yticks(range(0, n_samples_total + 1))
ax.legend(frameon=False, fontsize=9)
sns.despine()
plt.tight_layout()
b64_fig8 = fig_to_b64(fig)
plt.close(fig)
print("  ✓ Done")

# ─── ARG Table HTML ───────────────────────────────────────────────────────────

arg_table = arg_counts_df.copy()
arg_table['prevalence'] = arg_table.groupby('best_hit_aro')['sample_id'] \
                                    .transform('nunique')
arg_table = arg_table.rename(columns={
    'best_hit_aro': 'ARG',
    'sample_id':    'Sample',
    'count':        'Count',
    'mapped_reads': 'Mapped Reads',
    'drug_class':   'Drug Class',
})
arg_table = arg_table.sort_values(['ARG', 'Sample']).reset_index(drop=True)

# Index de la colonne Drug Class dans le tableau HTML (0-based)
col_order = list(arg_table.columns)
drug_col_idx = col_order.index('Drug Class') if 'Drug Class' in col_order else -1

# Split drug_class combinées pour le dropdown — valeurs atomiques uniquement
if 'Drug Class' in arg_table.columns:
    drug_classes = sorted(set(
        dc.strip()
        for val in arg_table['Drug Class'].dropna()
        for dc in str(val).split(';')
        if dc.strip()
    ))
else:
    drug_classes = []

drug_options = ''.join(f'<option value="{dc}">{dc}</option>'
                       for dc in drug_classes)

arg_table_html = f"""
<div style="margin-bottom:1em; display:flex; gap:1em; align-items:center; flex-wrap:wrap;">
  <label style="font-weight:bold;">Filter by Drug Class:</label>
  <select id="drug-filter" onchange="filterTable()"
          style="padding:6px 12px; border-radius:6px; border:1px solid #ccc; min-width:200px;">
    <option value="">— All —</option>
    {drug_options}
  </select>
  <label style="font-weight:bold;">Search ARG:</label>
  <input id="arg-search" type="text" oninput="filterTable()"
         placeholder="e.g. CTX-M..."
         style="padding:6px 12px; border-radius:6px; border:1px solid #ccc; width:180px;">
  <button onclick="exportCSV()"
          style="padding:6px 16px; background:#2e6da4; color:white;
                 border:none; border-radius:20px; cursor:pointer; font-size:0.85em;">
    ⬇ Export CSV
  </button>
</div>

<div style="max-height:500px; overflow-y:auto;">
{arg_table.to_html(index=False, border=0, classes="stats-table", table_id="arg-table")}
</div>

<script>
function filterTable() {{
  const drugVal   = document.getElementById('drug-filter').value.toLowerCase();
  const searchVal = document.getElementById('arg-search').value.toLowerCase();
  const drugIdx   = {drug_col_idx};
  const rows = document.querySelectorAll('#arg-table tbody tr');
  rows.forEach(row => {{
    const cells    = row.querySelectorAll('td');
    const argText  = cells[1] ? cells[1].innerText.toLowerCase() : '';
    const drugText = (drugIdx >= 0 && cells[drugIdx])
                     ? cells[drugIdx].innerText.toLowerCase() : '';
    const show = (!drugVal  || drugText.includes(drugVal)) &&
                 (!searchVal || argText.includes(searchVal));
    row.style.display = show ? '' : 'none';
  }});
}}

function exportCSV() {{
  const rows = document.querySelectorAll('#arg-table tr');
  const csv  = Array.from(rows)
    .filter(r => r.style.display !== 'none')
    .map(r => Array.from(r.querySelectorAll('th,td'))
                   .map(c => '"' + c.innerText.replace(/"/g, '""') + '"')
                   .join(','))
    .join('\\n');
  const blob = new Blob([csv], {{type: 'text/csv'}});
  const a    = document.createElement('a');
  a.href     = URL.createObjectURL(blob);
  a.download = '{tool}_arg_table.csv';
  a.click();
}}
</script>
"""


# ─── Summary table ───────────────────────────────────────────────────────────

sample_stats = []
for sample in arg_matrix_norm.index:
    row = arg_matrix_norm.loc[sample]
    sample_stats.append({
        'Sample':        sample,
        'Unique ARGs':   int((row > 0).sum()),
        'Shannon Index': f"{shannon_diversity(row):.3f}",
        'Top Drug Class': (drug_class_df[drug_class_df['sample_id'] == sample]
                           .sort_values('count', ascending=False)
                           ['drug_class'].iloc[0]
                           if len(drug_class_df[drug_class_df['sample_id'] == sample]) > 0
                           else 'N/A'),
        'Top Mechanism':  (mechanism_df[mechanism_df['sample_id'] == sample]
                           .sort_values('count', ascending=False)
                           ['resistance_mechanism'].iloc[0]
                           if len(mechanism_df[mechanism_df['sample_id'] == sample]) > 0
                           else 'N/A'),
    })

stats_table_html = pd.DataFrame(sample_stats).to_html(index=False, border=0,
                                                        classes="stats-table")

# ─── HTML Report ─────────────────────────────────────────────────────────────

CSS = """
<style>
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body {
    font-family: 'Segoe UI', Arial, sans-serif;
    background: #f8f9fa;
    color: #222;
  }
  header {
    background: linear-gradient(135deg, #1a3a5c, #2e6da4);
    color: white;
    padding: 2em 3em;
  }
  header h1 { font-size: 1.8em; margin-bottom: 0.3em; }
  header p  { opacity: 0.85; font-size: 0.95em; }
  .badge {
    display: inline-block;
    background: rgba(255,255,255,0.2);
    border-radius: 12px;
    padding: 3px 12px;
    font-size: 0.85em;
    margin-top: 0.5em;
    margin-right: 0.4em;
  }
  .container { max-width: 1100px; margin: 0 auto; padding: 2em; }

  /* Summary cards */
  .cards { display: flex; gap: 1em; flex-wrap: wrap; margin: 1.5em 0; }
  .card {
    background: white;
    border-radius: 10px;
    padding: 1.2em 1.8em;
    box-shadow: 0 2px 8px rgba(0,0,0,0.08);
    min-width: 140px;
  }
  .card b { display: block; font-size: 2em; color: #2e6da4; }
  .card span { font-size: 0.85em; color: #666; }

  /* Params box */
  .params {
    background: white;
    border-left: 4px solid #2e6da4;
    border-radius: 0 8px 8px 0;
    padding: 1em 1.5em;
    margin: 1.5em 0;
    box-shadow: 0 2px 8px rgba(0,0,0,0.06);
  }
  .params h3 { color: #1a3a5c; margin-bottom: 0.6em; }
  .params table { font-size: 0.9em; }
  .params td { padding: 2px 12px 2px 0; }
  .params td:first-child { font-weight: bold; color: #555; }

  /* Tabs */
  .tabs { display: flex; gap: 4px; flex-wrap: wrap; margin: 2em 0 0; }
  .tab-btn {
    background: #dde3ea;
    border: none;
    border-radius: 8px 8px 0 0;
    padding: 0.6em 1.2em;
    cursor: pointer;
    font-size: 0.9em;
    font-weight: bold;
    color: #444;
    transition: background 0.2s;
  }
  .tab-btn:hover   { background: #c5d0db; }
  .tab-btn.active  { background: white; color: #1a3a5c; border-bottom: 2px solid white; }
  .tab-content {
    display: none;
    background: white;
    border-radius: 0 8px 8px 8px;
    padding: 2em;
    box-shadow: 0 2px 8px rgba(0,0,0,0.07);
    margin-bottom: 2em;
  }
  .tab-content.active { display: block; }

  /* Figures */
  .fig-container { text-align: center; margin: 1em 0; }
  .fig-container img { max-width: 100%; border-radius: 6px; }
  .dl-btn {
    display: inline-block;
    margin-top: 0.6em;
    background: #2e6da4;
    color: white;
    text-decoration: none;
    padding: 5px 14px;
    border-radius: 20px;
    font-size: 0.82em;
    transition: background 0.2s;
  }
  .dl-btn:hover { background: #1a3a5c; }

  /* Stats table */
  .stats-table { width: 100%; border-collapse: collapse; font-size: 0.88em; }
  .stats-table th {
    background: #1a3a5c; color: white;
    padding: 8px 12px; text-align: left;
  }
  .stats-table td { border-bottom: 1px solid #eee; padding: 7px 12px; }
  .stats-table tr:hover td { background: #f0f5fa; }

  .section-desc { color: #555; font-size: 0.92em; margin-bottom: 1.2em; }
</style>
"""

JS = """
<script>
function showTab(tabId) {
  document.querySelectorAll('.tab-content').forEach(t => t.classList.remove('active'));
  document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
  document.getElementById(tabId).classList.add('active');
  event.target.classList.add('active');
}
</script>
"""

# Résumé exécutif
exec_summary = f"""
<div class="cards">
  <div class="card"><b>{n_samples}</b><span>Samples</span></div>
  <div class="card"><b>{n_args}</b><span>Unique ARGs</span></div>
  <div class="card"><b>{n_classes}</b><span>Drug Classes</span></div>
  <div class="card"><b>{diversity_df['shannon'].mean():.2f}</b><span>Mean Shannon</span></div>
  <div class="card"><b>{int(diversity_df['richness'].mean())}</b><span>Mean Richness</span></div>
</div>
"""

params_box = f"""
<div class="params">
  <h3>Pipeline Parameters</h3>
  <table>
    <tr><td>Tool</td><td>{tool_label}</td></tr>
    <tr><td>Min Identity</td><td>{min_identity}%</td></tr>
    <tr><td>Min Coverage</td><td>{min_coverage}%</td></tr>
    <tr><td>Report date</td><td>{report_date}</td></tr>
  </table>
</div>
"""

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>ARG Report — {tool_label}</title>
  {CSS}
</head>
<body>

<header>
  <h1>ARG Analysis Report</h1>
  <p>MetagenAMR Pipeline — Antibiotic Resistance Genes</p>
  <span class="badge">{tool_label}</span>
  <span class="badge">{report_date}</span>
  <span class="badge">{n_samples} samples</span>
</header>

<div class="container">

  <!-- Summary tab always visible -->
  <div style="margin-top:1.5em">
    <h2 style="color:#1a3a5c; margin-bottom:0.8em">Executive Summary</h2>
    {exec_summary}
    {params_box}
  </div>

  <!-- Tabs -->
  <div class="tabs">
    <button class="tab-btn active" onclick="showTab('tab-drug')">Drug Classes</button>
    <button class="tab-btn" onclick="showTab('tab-diversity')">Diversity</button>
    <button class="tab-btn" onclick="showTab('tab-heatmap')">Heatmap</button>
    <button class="tab-btn" onclick="showTab('tab-mech')">Mechanisms</button>
    <button class="tab-btn" onclick="showTab('tab-rarefaction')">Rarefaction</button>
    <button class="tab-btn" onclick="showTab('tab-pcoa')">Beta Diversity</button>
    <button class="tab-btn" onclick="showTab('tab-top20')">Top ARGs</button>
    <button class="tab-btn" onclick="showTab('tab-samples')">Samples</button>
    <button class="tab-btn" onclick="showTab('tab-bubble')">ARG Families</button>
    <button class="tab-btn" onclick="showTab('tab-argtable')">ARG Table</button>
  </div>

  <div id="tab-drug" class="tab-content active">
    <p class="section-desc">Relative abundance of antibiotic resistance classes per sample.</p>
    {img_tag(b64_fig1, f"{tool}_drug_classes")}
  </div>
  
  <div id="tab-bubble" class="tab-content">
  <p class="section-desc">ARG family abundance vs prevalence across samples. 
  Bubble size reflects total abundance. Annotated families are high-priority 
  (abundant and/or ubiquitous).</p>
  {img_tag(b64_fig8, f"{tool}_bubble_plot")}
  </div>

  <div id="tab-argtable" class="tab-content">
  <p class="section-desc">Per-sample ARG detection table. 
  Filter by drug class or search by ARG name. Export visible rows as CSV.</p>
  {arg_table_html}
  </div>

  <div id="tab-diversity" class="tab-content">
    <p class="section-desc">Alpha diversity (Shannon index) and ARG richness per sample.</p>
    {img_tag(b64_fig2, f"{tool}_alpha_diversity")}
  </div>

  <div id="tab-heatmap" class="tab-content">
    <p class="section-desc">Hierarchical clustering of top 30 ARGs across samples (log₁₀ normalized copy number).</p>
    {img_tag(b64_fig3, f"{tool}_heatmap_top30")}
  </div>

  <div id="tab-mech" class="tab-content">
    <p class="section-desc">Distribution of resistance mechanisms detected across all samples.</p>
    {img_tag(b64_fig4, f"{tool}_mechanisms")}
  </div>

  <div id="tab-rarefaction" class="tab-content">
    <p class="section-desc">Rarefaction curves showing ARG discovery saturation per sample.</p>
    {img_tag(b64_fig5, f"{tool}_rarefaction")}
  </div>

  <div id="tab-pcoa" class="tab-content">
    <p class="section-desc">Principal Coordinates Analysis (PCoA) based on Bray-Curtis dissimilarity.</p>
    {img_tag(b64_fig6, f"{tool}_pcoa")}
  </div>

  <div id="tab-top20" class="tab-content">
    <p class="section-desc">Top 20 most frequently observed ARGs across all samples.</p>
    {img_tag(b64_fig7, f"{tool}_top20_args")}
  </div>

  <div id="tab-samples" class="tab-content">
    <p class="section-desc">Per-sample statistics summary.</p>
    {stats_table_html}
  </div>

</div>

{JS}
</body>
</html>"""

Path(snakemake.output.report).parent.mkdir(parents=True, exist_ok=True)
with open(snakemake.output.report, 'w') as f:
    f.write(html)

print(f"\n✓ Report saved → {snakemake.output.report}")
print("\n" + "=" * 70)
print(f"REPORT COMPLETED  [{tool_label}]")
print("=" * 70)
print(f"  Samples : {n_samples}")
print(f"  ARGs    : {n_args}")
print(f"  Classes : {n_classes}")
print("=" * 70 + "\n")

log.close()
