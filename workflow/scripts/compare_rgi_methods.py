#!/usr/bin/env python3
"""
Compare les ARGs détectés par rgi main vs rgi bwt
Outputs : TSV de comparaison + rapport HTML (venn + corrélations)
"""

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from scipy import stats
import base64
from io import BytesIO
import sys
from pathlib import Path


def fig_to_b64(fig):
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()


def img_tag(b64):
    return f'<img src="data:image/png;base64,{b64}" style="max-width:600px; display:block; margin:1em 0;">'


def main():
    input_rgi     = snakemake.input.rgi
    input_rgi_bwt = snakemake.input.rgi_bwt
    output_tsv    = snakemake.output.tsv
    output_html   = snakemake.output.html
    log_file      = snakemake.log[0]

    log = open(log_file, 'w')
    sys.stderr = sys.stdout = log

    print("=" * 70)
    print("RGI Comparison — rgi main vs rgi bwt")
    print("=" * 70)

    # ─── Load ────────────────────────────────────────────────────────────────

    print(f"\nLoading rgi main     → {input_rgi}")
    rgi = pd.read_csv(input_rgi, sep='\t', low_memory=False)
    print(f"  ✓ {len(rgi):,} rows")

    print(f"Loading rgi bwt      → {input_rgi_bwt}")
    bwt = pd.read_csv(input_rgi_bwt, sep='\t', low_memory=False)
    print(f"  ✓ {len(bwt):,} rows")

    # Colonnes normalisées post-filter (lowercase + _)
    # rgi main  : best_hit_aro / best_identities / percentage_length_of_reference_sequence
    # rgi bwt   : aro_term     / average_percent_identity / average_percent_coverage

    RGI_ARO  = "best_hit_aro"
    RGI_ID   = "best_identities"
    RGI_COV  = "percentage_length_of_reference_sequence"

    BWT_ARO  = "aro_term"
    BWT_ID   = "average_percent_identity"
    BWT_COV  = "average_percent_coverage"

    # ─── Sets pour Venn ──────────────────────────────────────────────────────

    set_rgi = set(rgi[RGI_ARO].dropna().unique()) if RGI_ARO in rgi.columns else set()
    set_bwt = set(bwt[BWT_ARO].dropna().unique()) if BWT_ARO in bwt.columns else set()

    only_rgi = set_rgi - set_bwt
    only_bwt = set_bwt - set_rgi
    shared   = set_rgi & set_bwt

    print(f"\nOverlap:")
    print(f"  - rgi main only : {len(only_rgi)}")
    print(f"  - rgi bwt only  : {len(only_bwt)}")
    print(f"  - Shared        : {len(shared)}")

    # ─── TSV de comparaison ──────────────────────────────────────────────────

    rows = []
    for aro in set_rgi | set_bwt:
        rows.append({
            "aro_term":     aro,
            "in_rgi_main":  aro in set_rgi,
            "in_rgi_bwt":   aro in set_bwt,
            "exclusive_to": "rgi_main" if aro in only_rgi
                            else "rgi_bwt" if aro in only_bwt
                            else "shared"
        })

    comparison = pd.DataFrame(rows).sort_values("exclusive_to").reset_index(drop=True)

    # Ajout identité/coverage moyens pour ARGs partagés
    if shared and RGI_ID in rgi.columns and RGI_COV in rgi.columns:
        rgi_shared = (rgi[rgi[RGI_ARO].isin(shared)]
                      .groupby(RGI_ARO)[[RGI_ID, RGI_COV]]
                      .mean()
                      .rename(columns={RGI_ID: "identity_rgi", RGI_COV: "coverage_rgi"}))
        comparison = comparison.merge(rgi_shared, left_on="aro_term",
                                      right_index=True, how="left")

    if shared and BWT_ID in bwt.columns and BWT_COV in bwt.columns:
        bwt_shared = (bwt[bwt[BWT_ARO].isin(shared)]
                      .groupby(BWT_ARO)[[BWT_ID, BWT_COV]]
                      .mean()
                      .rename(columns={BWT_ID: "identity_rgi_bwt", BWT_COV: "coverage_rgi_bwt"}))
        comparison = comparison.merge(bwt_shared, left_on="aro_term",
                                      right_index=True, how="left")

    Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(output_tsv, sep='\t', index=False)
    print(f"\n✓ TSV saved → {output_tsv}")

    # ─── Figure 1 : Venn ─────────────────────────────────────────────────────

    fig1, ax1 = plt.subplots(figsize=(5, 4))
    venn2([set_rgi, set_bwt], set_labels=("rgi main", "rgi bwt"), ax=ax1)
    ax1.set_title("ARGs détectés — overlap")
    b64_venn = fig_to_b64(fig1)
    plt.close(fig1)

    # ─── Figures 2-3 : Corrélations (ARGs partagés) ──────────────────────────

    corr_blocks = []

    shared_df = comparison[comparison["exclusive_to"] == "shared"].copy()

    for metric, label in [("identity", "% Identity"), ("coverage", "% Coverage")]:
        col_r = f"{metric}_rgi"
        col_b = f"{metric}_rgi_bwt"

        if col_r not in shared_df.columns or col_b not in shared_df.columns:
            continue

        sub = shared_df[[col_r, col_b]].dropna()
        if len(sub) < 3:
            corr_blocks.append(f"<h3>{label}</h3><p>Pas assez de données partagées ({len(sub)} ARGs).</p>")
            continue

        x, y = sub[col_r], sub[col_b]
        r, p = stats.pearsonr(x, y)

        fig, ax = plt.subplots(figsize=(5, 4))
        ax.scatter(x, y, alpha=0.6, edgecolors="k", linewidths=0.4)
        lims = [min(x.min(), y.min()) - 2, max(x.max(), y.max()) + 2]
        ax.plot(lims, lims, "r--", lw=1, label="y = x")
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_xlabel(f"rgi main — {label}")
        ax.set_ylabel(f"rgi bwt — {label}")
        ax.set_title(f"{label}  (r = {r:.2f}, p = {p:.2e}, n = {len(sub)})")
        ax.legend(fontsize=8)
        corr_blocks.append(f"<h3>{label}</h3>{img_tag(fig_to_b64(fig))}")
        plt.close(fig)

    corr_html = "\n".join(corr_blocks) if corr_blocks else \
        "<p>Pas d'ARGs partagés avec données suffisantes.</p>"

    # ─── Tableau HTML ARGs exclusifs ─────────────────────────────────────────

    def subset_table(df, tool):
        sub = df[df["exclusive_to"] == tool][["aro_term"]].copy()
        if sub.empty:
            return "<p>Aucun.</p>"
        return sub.to_html(index=False, border=0)

    # ─── Rapport HTML ────────────────────────────────────────────────────────

    html = f"""<!DOCTYPE html>
<html lang="fr">
<head>
  <meta charset="UTF-8">
  <title>Comparaison RGI main vs RGI BWT</title>
  <style>
    body  {{ font-family: Arial, sans-serif; max-width: 960px; margin: 2em auto; color: #222; }}
    h1    {{ border-bottom: 3px solid #333; padding-bottom: 6px; }}
    h2    {{ border-bottom: 1px solid #aaa; padding-bottom: 4px; margin-top: 2em; }}
    table {{ border-collapse: collapse; width: 100%; font-size: 0.83em; margin-top: 0.5em; }}
    th, td {{ border: 1px solid #ccc; padding: 5px 10px; text-align: left; }}
    th    {{ background: #f0f0f0; }}
    tr:nth-child(even) {{ background: #fafafa; }}
    .stat {{ display: inline-block; background: #f5f5f5; border: 1px solid #ddd;
             border-radius: 6px; padding: 0.6em 1.2em; margin: 0.3em; font-size: 1.1em; }}
    .stat b {{ display: block; font-size: 1.6em; color: #2a6ebb; }}
  </style>
</head>
<body>

<h1>Comparaison rgi main vs rgi bwt</h1>

<h2>1. Overlap des ARGs détectés</h2>

<div>
  <span class="stat"><b>{len(set_rgi)}</b> ARGs — rgi main</span>
  <span class="stat"><b>{len(set_bwt)}</b> ARGs — rgi bwt</span>
  <span class="stat"><b>{len(shared)}</b> Partagés</span>
  <span class="stat"><b>{len(only_rgi)}</b> rgi main uniquement</span>
  <span class="stat"><b>{len(only_bwt)}</b> rgi bwt uniquement</span>
</div>

{img_tag(b64_venn)}

<h2>2. ARGs exclusifs à rgi main</h2>
{subset_table(comparison, "rgi_main")}

<h2>3. ARGs exclusifs à rgi bwt</h2>
{subset_table(comparison, "rgi_bwt")}

<h2>4. Corrélation identité / coverage (ARGs partagés, n={len(shared)})</h2>
{corr_html}

<h2>5. Tableau complet</h2>
{comparison.to_html(index=False, border=0)}

</body>
</html>"""

    with open(output_html, 'w') as fh:
        fh.write(html)
    print(f"✓ HTML report → {output_html}")

    print("\n" + "=" * 70)
    print("COMPARISON COMPLETED")
    print("=" * 70)
    print(f"  rgi main only : {len(only_rgi)}")
    print(f"  rgi bwt only  : {len(only_bwt)}")
    print(f"  Shared        : {len(shared)}")
    print("=" * 70 + "\n")

    log.close()


if __name__ == "__main__":
    main()
