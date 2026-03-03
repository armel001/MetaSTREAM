#!/usr/bin/env python3
"""
Agrège les résultats RGI (main ou bwt) de tous les échantillons en un seul fichier
"""

import pandas as pd
import sys
from pathlib import Path


# ─── Mapping colonnes selon tool ─────────────────────────────────────────────
# NB: après filter_rgi_results.py, les colonnes sont normalisées (lowercase, _ )

COLUMN_MAP = {
    "rgi": {
        "aro":       "best_hit_aro",
        "drug":      "drug_class",
        "mechanism": "resistance_mechanism",
        "family":    "amr_gene_family",
    },
    "rgi_bwt": {
        "aro":       "aro_term",
        "drug":      "drug_class",
        "mechanism": "resistance_mechanism",
        "family":    "amr_gene_family",
    }
}


def main():
    input_files  = snakemake.input.selected
    output_file  = snakemake.output.aggregated
    summary_file = snakemake.output.summary
    sample_names = snakemake.params.samples
    tool         = snakemake.wildcards.tool
    log_file     = snakemake.log[0]

    log = open(log_file, 'w')
    sys.stderr = sys.stdout = log

    cols = COLUMN_MAP.get(tool)
    if cols is None:
        print(f"ERROR: Unknown tool '{tool}'.")
        sys.exit(1)

    col_aro  = cols["aro"]
    col_drug = cols["drug"]
    col_mech = cols["mechanism"]

    print("=" * 70)
    print(f"RGI Results Aggregation  [{tool}]")
    print("=" * 70)
    print(f"\nNumber of samples: {len(sample_names)}")
    print(f"Output: {output_file}")
    print("\n" + "-" * 70)

    all_data     = []
    sample_stats = []

    print(f"\nLoading individual sample files...")
    for i, (file_path, sample_name) in enumerate(zip(input_files, sample_names), 1):
        print(f"\n[{i}/{len(sample_names)}] Processing: {sample_name}")

        if not Path(file_path).exists():
            print(f"  ⚠ WARNING: File not found, skipping: {file_path}")
            continue

        try:
            df = pd.read_csv(file_path, sep='\t', low_memory=False)
            df['sample_id'] = sample_name

            n_rows      = len(df)
            n_unique    = df[col_aro].nunique() if col_aro in df.columns else 0

            print(f"  ✓ Loaded: {n_rows:,} rows, {n_unique} unique ARGs")

            sample_stats.append({
                'sample_id':     sample_name,
                'n_args':        n_rows,
                'n_unique_args': n_unique
            })
            all_data.append(df)

        except Exception as e:
            print(f"  ✗ ERROR: {str(e)}")
            continue

    if not all_data:
        print(f"\nERROR: No valid sample files found!")
        sys.exit(1)

    print(f"\n" + "-" * 70)
    print(f"Successfully loaded {len(all_data)} sample(s)")

    print(f"\nCombining all samples...")
    df_combined = pd.concat(all_data, ignore_index=True)

    # sample_id en premier
    cols_order = df_combined.columns.tolist()
    if 'sample_id' in cols_order:
        cols_order.remove('sample_id')
        df_combined = df_combined[['sample_id'] + cols_order]

    print(f"  ✓ Total rows:    {len(df_combined):,}")
    print(f"  ✓ Total columns: {len(df_combined.columns)}")

    # ─── Statistiques globales ───────────────────────────────────────────────

    print(f"\n" + "-" * 70)
    print(f"Global Statistics:")
    print(f"  - Total ARG observations: {len(df_combined):,}")
    print(f"  - Unique samples:         {df_combined['sample_id'].nunique()}")

    if col_aro in df_combined.columns:
        print(f"  - Unique ARGs detected:   {df_combined[col_aro].nunique():,}")
        print(f"\n  Top 10 most frequent ARGs:")
        for j, (arg, count) in enumerate(df_combined[col_aro].value_counts().head(10).items(), 1):
            print(f"    {j:2d}. {arg}: {count:,} ({100*count/len(df_combined):.1f}%)")

    if col_drug in df_combined.columns:
        drug_expanded = df_combined[col_drug].str.split('; ').explode()
        print(f"\n  Unique drug classes: {drug_expanded.nunique()}")
        print(f"\n  Top 10 drug classes:")
        for j, (dc, count) in enumerate(drug_expanded.value_counts().head(10).items(), 1):
            print(f"    {j:2d}. {dc}: {count:,} ({100*count/len(drug_expanded):.1f}%)")

    if col_mech in df_combined.columns:
        print(f"\n  Resistance mechanisms distribution:")
        for mech, count in df_combined[col_mech].value_counts().head(10).items():
            print(f"    • {mech}: {count:,} ({100*count/len(df_combined):.1f}%)")

    # ─── Sauvegarde ──────────────────────────────────────────────────────────

    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"\n" + "-" * 70)
    print(f"Saving aggregated data...")
    try:
        df_combined.to_csv(output_file, sep='\t', index=False)
        file_size = output_path.stat().st_size / 1024 / 1024
        print(f"  ✓ Saved → {output_file} ({file_size:.2f} MB)")
    except Exception as e:
        print(f"\nERROR: Failed to save: {str(e)}")
        sys.exit(1)

    # ─── Fichier summary ─────────────────────────────────────────────────────

    print(f"\nCreating summary file...")
    try:
        with open(summary_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write(f"RGI Aggregation Summary  [{tool}]\n")
            f.write("=" * 70 + "\n\n")
            f.write(f"Total samples processed:  {len(all_data)}\n")
            f.write(f"Total ARG observations:   {len(df_combined):,}\n")
            if col_aro in df_combined.columns:
                f.write(f"Unique ARGs detected:     {df_combined[col_aro].nunique():,}\n\n")

            f.write("-" * 70 + "\n")
            f.write("Per-sample statistics:\n")
            f.write("-" * 70 + "\n")
            f.write(f"{'Sample':<20} {'N_ARGs':>10} {'Unique_ARGs':>15}\n")
            f.write("-" * 70 + "\n")
            for stat in sample_stats:
                f.write(f"{stat['sample_id']:<20} {stat['n_args']:>10,} {stat['n_unique_args']:>15}\n")

            f.write("\n" + "=" * 70 + "\n")
            f.write("\nDetailed statistics by sample:\n")
            f.write("=" * 70 + "\n")
            for sample in sample_names:
                sd = df_combined[df_combined['sample_id'] == sample]
                f.write(f"\n{sample}:\n")
                f.write(f"  - Total ARGs:  {len(sd):,}\n")
                if col_aro in sd.columns:
                    f.write(f"  - Unique ARGs: {sd[col_aro].nunique()}\n")
                if col_drug in sd.columns:
                    top3 = sd[col_drug].str.split('; ').explode().value_counts().head(3)
                    f.write(f"  - Top 3 drug classes:\n")
                    for dc, count in top3.items():
                        f.write(f"      • {dc}: {count}\n")

        print(f"  ✓ Summary → {summary_file}")

    except Exception as e:
        print(f"\nWARNING: Failed to save summary: {str(e)}")

    # ─── Résumé final ────────────────────────────────────────────────────────

    print("\n" + "=" * 70)
    print("AGGREGATION COMPLETED SUCCESSFULLY")
    print("=" * 70)
    print(f"  Tool:        {tool}")
    print(f"  Samples:     {len(all_data)}")
    print(f"  Total ARGs:  {len(df_combined):,}")
    if col_aro in df_combined.columns:
        print(f"  Unique ARGs: {df_combined[col_aro].nunique():,}")
    print(f"  Output size: {file_size:.2f} MB")
    print("=" * 70 + "\n")

    log.close()


if __name__ == "__main__":
    main()
