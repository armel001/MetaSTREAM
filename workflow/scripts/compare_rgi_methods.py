#!/usr/bin/env python3
"""
compare_rgi_methods.py — ARG overlap: rgi main vs rgi bwt
Outputs: TSV comparison + HTML report (Venn + correlations)
"""

import sys, base64
from io import BytesIO
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from scipy import stats


# ── Helpers ───────────────────────────────────────────────────────────────────
def fig_to_b64(fig):
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()

def img_tag(b64):
    return f'<img src="data:image/png;base64,{b64}" style="max-width:600px;display:block;margin:1em 0;">'

def df_to_html(df):
    return df.to_html(index=False, border=0) if not df.empty else "<p>None.</p>"


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    input_rgi     = snakemake.input.rgi
    input_rgi_bwt = snakemake.input.rgi_bwt
    output_tsv    = snakemake.output.tsv
    output_html   = snakemake.output.html

    log = open(snakemake.log[0], 'w')
    sys.stderr = sys.stdout = log

    print("=" * 70)
    print("RGI Comparison — rgi main vs rgi bwt")
    print("=" * 70)

    # Column names (post-normalization: lowercase + underscore)
    RGI_ARO, RGI_ID, RGI_COV = "best_hit_aro", "best_identities", "percentage_length_of_reference_sequence"
    BWT_ARO, BWT_ID, BWT_COV = "aro_term", "average_percent_identity", "average_percent_coverage"

    # ── Load ──────────────────────────────────────────────────────────────────
    rgi = pd.read_csv(input_rgi,     sep='\t', low_memory=False)
    bwt = pd.read_csv(input_rgi_bwt, sep='\t', low_memory=False)
    print(f"\n  rgi main : {len(rgi):,} rows | rgi bwt : {len(bwt):,} rows")

    # ── ARG sets ──────────────────────────────────────────────────────────────
    set_rgi  = set(rgi[RGI_ARO].dropna().unique()) if RGI_ARO in rgi.columns else set()
    set_bwt  = set(bwt[BWT_ARO].dropna().unique()) if BWT_ARO in bwt.columns else set()
    only_rgi = set_rgi - set_bwt
    only_bwt = set_bwt - set_rgi
    shared   = set_rgi & set_bwt

    print(f"\n  rgi main only: {len(only_rgi)} | rgi bwt only: {len(only_bwt)} | Shared: {len(shared)}")

    # ── Comparison TSV ────────────────────────────────────────────────────────
    comparison = pd.DataFrame([{
        "aro_term":     aro,
        "in_rgi_main":  aro in set_rgi,
        "in_rgi_bwt":   aro in set_bwt,
        "exclusive_to": "rgi_main" if aro in only_rgi else "rgi_bwt" if aro in only_bwt else "shared"
    } for aro in set_rgi | set_bwt]).sort_values("exclusive_to").reset_index(drop=True)

    # Join mean identity/coverage for shared ARGs
    for aro_col, id_col, cov_col, df, suffix in [
        (RGI_ARO, RGI_ID, RGI_COV, rgi, "rgi"),
        (BWT_ARO, BWT_ID, BWT_COV, bwt, "rgi_bwt"),
    ]:
        if shared and id_col in df.columns and cov_col in df.columns:
            agg = (df[df[aro_col].isin(shared)]
                   .groupby(aro_col)[[id_col, cov_col]].mean()
                   .rename(columns={id_col: f"identity_{suffix}", cov_col: f"coverage_{suffix}"}))
            comparison = comparison.merge(agg, left_on="aro_term", right_index=True, how="left")

    Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(output_tsv, sep='\t', index=False)
    print(f"\n  [OK] TSV → {output_tsv}")

    # ── Figure 1 : Venn ───────────────────────────────────────────────────────
    fig1, ax1 = plt.subplots(figsize=(5, 4))
    venn2([set_rgi, set_bwt], set_labels=("rgi main", "rgi bwt"), ax=ax1)
    ax1.set_title("ARG overlap — rgi main vs rgi bwt")
    b64_venn = fig_to_b64(fig1)
    plt.close(fig1)

    # ── Figures 2–3 : Correlations (shared ARGs) ──────────────────────────────
    shared_df   = comparison[comparison["exclusive_to"] == "shared"]
    corr_blocks = []

    for metric, label in [("identity", "% Identity"), ("coverage", "% Coverage")]:
        col_r, col_b = f"{metric}_rgi", f"{metric}_rgi_bwt"
        if col_r not in shared_df.columns or col_b not in shared_df.columns:
            continue
        sub = shared_df[[col_r, col_b]].dropna()
        if len(sub) < 3:
            corr_blocks.append(f"<h3>{label}</h3><p>Insufficient shared data ({len(sub)} ARGs).</p>")
            continue

        x, y   = sub[col_r], sub[col_b]
        r, p   = stats.pearsonr(x, y)
        lims   = [min(x.min(), y.min()) - 2, max(x.max(), y.max()) + 2]

        fig, ax = plt.subplots(figsize=(5, 4))
        ax.scatter(x, y, alpha=0.6, edgecolors="k", linewidths=0.4)
        ax.plot(lims, lims, "r--", lw=1, label="y = x")
        ax.set_xlim(lims); ax.set_ylim(lims)
        ax.set_xlabel(f"rgi main — {label}"); ax.set_ylabel(f"rgi bwt — {label}")
        ax.set_title(f"{label}  (r={r:.2f}, p={p:.2e}, n={len(sub)})"); ax.legend(fontsize=8)
        corr_blocks.append(f"<h3>{label}</h3>{img_tag(fig_to_b64(fig))}")
        plt.close(fig)

    corr_html = "\n".join(corr_blocks) or "<p>No shared ARGs with sufficient data.</p>"

    # ── HTML report ───────────────────────────────────────────────────────────
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>RGI main vs RGI BWT Comparison</title>
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
<h1>RGI main vs RGI BWT — ARG Comparison</h1>

<h2>1. ARG Overlap</h2>
<div>
  <span class="stat"><b>{len(set_rgi)}</b> rgi main</span>
  <span class="stat"><b>{len(set_bwt)}</b> rgi bwt</span>
  <span class="stat"><b>{len(shared)}</b> Shared</span>
  <span class="stat"><b>{len(only_rgi)}</b> rgi main only</span>
  <span class="stat"><b>{len(only_bwt)}</b> rgi bwt only</span>
</div>
{img_tag(b64_venn)}

<h2>2. ARGs exclusive to rgi main</h2>
{df_to_html(comparison[comparison["exclusive_to"]=="rgi_main"][["aro_term"]])}

<h2>3. ARGs exclusive to rgi bwt</h2>
{df_to_html(comparison[comparison["exclusive_to"]=="rgi_bwt"][["aro_term"]])}

<h2>4. Identity / Coverage correlation (shared ARGs, n={len(shared)})</h2>
{corr_html}

<h2>5. Full comparison table</h2>
{comparison.to_html(index=False, border=0)}
</body>
</html>"""

    with open(output_html, 'w') as fh:
        fh.write(html)
    print(f"  [OK] HTML → {output_html}")
    print(f"\n══ COMPARISON COMPLETED ══\n  rgi main only: {len(only_rgi)} | rgi bwt only: {len(only_bwt)} | Shared: {len(shared)}\n")
    log.close()


if __name__ == "__main__":
    main()
