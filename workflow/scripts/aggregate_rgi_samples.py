#!/usr/bin/env python3
"""
aggregate_rgi_samples.py — Aggregate RGI results (main or bwt) across all samples.
"""

import sys
from pathlib import Path
import pandas as pd

# Column names after normalize step (lowercase + underscore)
COLUMN_MAP = {
    "rgi":     {"aro": "best_hit_aro", "drug": "drug_class",
                "mechanism": "resistance_mechanism", "family": "amr_gene_family"},
    "rgi_bwt": {"aro": "aro_term",     "drug": "drug_class",
                "mechanism": "resistance_mechanism", "family": "amr_gene_family"},
}


def main():
    input_files  = snakemake.input.selected
    output_file  = snakemake.output.aggregated
    summary_file = snakemake.output.summary
    sample_names = snakemake.params.samples
    tool         = snakemake.wildcards.tool

    log = open(snakemake.log[0], 'w')
    sys.stderr = sys.stdout = log

    cols = COLUMN_MAP.get(tool)
    if cols is None:
        sys.exit(f"ERROR: Unknown tool '{tool}'")
    col_aro, col_drug, col_mech = cols["aro"], cols["drug"], cols["mechanism"]

    print(f"══ RGI Aggregation [{tool}] ══  {len(sample_names)} samples → {output_file}\n")

    # ── Load samples ──────────────────────────────────────────────────────────
    all_data, sample_stats = [], []

    for i, (path, sid) in enumerate(zip(input_files, sample_names), 1):
        print(f"  [{i}/{len(sample_names)}] {sid}", end=" ")
        if not Path(path).exists():
            print("— SKIPPED (file not found)")
            continue
        try:
            df = pd.read_csv(path, sep='\t', low_memory=False)
            df['sample_id'] = sid
            n_unique = df[col_aro].nunique() if col_aro in df.columns else 0
            print(f"— {len(df):,} rows | {n_unique} unique ARGs")
            sample_stats.append({"sample_id": sid, "n_args": len(df), "n_unique_args": n_unique})
            all_data.append(df)
        except Exception as e:
            print(f"— ERROR: {e}")

    if not all_data:
        sys.exit("ERROR: No valid sample files loaded")

    # ── Combine ───────────────────────────────────────────────────────────────
    df = pd.concat(all_data, ignore_index=True)
    if 'sample_id' in df.columns:
        df = df[['sample_id'] + [c for c in df.columns if c != 'sample_id']]

    print(f"\n  Combined: {len(df):,} rows × {len(df.columns)} columns")

    # ── Global stats ──────────────────────────────────────────────────────────
    if col_aro in df.columns:
        print(f"\n  Unique ARGs: {df[col_aro].nunique():,}")
        print("  Top 10 ARGs:")
        for j, (arg, n) in enumerate(df[col_aro].value_counts().head(10).items(), 1):
            print(f"    {j:2d}. {arg}: {n:,} ({100*n/len(df):.1f}%)")

    if col_drug in df.columns:
        expanded = df[col_drug].str.split('; ').explode()
        print(f"\n  Drug classes ({expanded.nunique()}) — Top 10:")
        for j, (dc, n) in enumerate(expanded.value_counts().head(10).items(), 1):
            print(f"    {j:2d}. {dc}: {n:,} ({100*n/len(expanded):.1f}%)")

    if col_mech in df.columns:
        print("\n  Mechanisms — Top 10:")
        for mech, n in df[col_mech].value_counts().head(10).items():
            print(f"    • {mech}: {n:,} ({100*n/len(df):.1f}%)")

    # ── Save aggregated TSV ───────────────────────────────────────────────────
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    try:
        df.to_csv(output_file, sep='\t', index=False)
        size_mb = Path(output_file).stat().st_size / 1024 / 1024
        print(f"\n  [OK] {output_file}  ({size_mb:.2f} MB)")
    except Exception as e:
        sys.exit(f"ERROR: Failed to save aggregated file — {e}")

    # ── Summary file ──────────────────────────────────────────────────────────
    try:
        lines = [
            f"══ RGI Aggregation Summary [{tool}] ══",
            f"Samples  : {len(all_data)}",
            f"Total ARGs: {len(df):,}",
        ]
        if col_aro in df.columns:
            lines.append(f"Unique ARGs: {df[col_aro].nunique():,}")
        lines += ["", f"{'Sample':<20} {'N_ARGs':>10} {'Unique_ARGs':>12}",
                  "-" * 44]
        for s in sample_stats:
            lines.append(f"{s['sample_id']:<20} {s['n_args']:>10,} {s['n_unique_args']:>12}")
        lines += ["", "── Per-sample detail ──"]
        for sid in sample_names:
            sd = df[df['sample_id'] == sid]
            lines.append(f"\n{sid}: {len(sd):,} ARGs")
            if col_aro in sd.columns:
                lines.append(f"  Unique ARGs: {sd[col_aro].nunique()}")
            if col_drug in sd.columns:
                top3 = sd[col_drug].str.split('; ').explode().value_counts().head(3)
                lines += [f"  Top 3 drug classes:"] + \
                         [f"    • {dc}: {n}" for dc, n in top3.items()]

        Path(summary_file).write_text("\n".join(lines) + "\n")
        print(f"  [OK] {summary_file}")
    except Exception as e:
        print(f"  [WARN] Summary not saved — {e}")

    print(f"\n══ AGGREGATION COMPLETED ══  {len(all_data)} samples | {len(df):,} rows | {size_mb:.2f} MB\n")
    log.close()


if __name__ == "__main__":
    main()
