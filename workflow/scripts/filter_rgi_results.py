#!/usr/bin/env python3

import pandas as pd
import sys
from pathlib import Path


# ─── Mapping colonnes selon tool ─────────────────────────────────────────────

COLUMN_MAP = {
    "rgi": {
        "aro":        "Best_Hit_ARO",
        "identity":   "Best_Identities",
        "coverage":   "Percentage Length of Reference Sequence",
        "model_type": "Model_type",
    },
    "rgi_bwt": {
        "aro":        "ARO Term",
        "identity":   "Average Percent Identity",
        "coverage":   "Average Percent Coverage",
        "model_type": "Model Type",
    }
}


def main():
    input_file  = snakemake.input.rgi
    output_file = snakemake.output.filtered
    min_identity = snakemake.params.min_identity
    min_coverage  = snakemake.params.min_coverage
    model_types   = snakemake.params.model_types
    tool          = snakemake.wildcards.tool
    log_file      = snakemake.log[0]

    log = open(log_file, 'w')
    sys.stderr = sys.stdout = log

    cols = COLUMN_MAP.get(tool)
    if cols is None:
        print(f"ERROR: Unknown tool '{tool}'. Expected 'rgi' or 'rgi_bwt'.")
        sys.exit(1)

    print("=" * 70)
    print(f"RGI Results Quality Filtering  [{tool}]")
    print("=" * 70)
    print(f"\nInput:  {input_file}")
    print(f"Output: {output_file}")
    print(f"\nParameters:")
    print(f"  - Min identity:  {min_identity}%")
    print(f"  - Min coverage:  {min_coverage}%")
    print(f"  - Model types:   {', '.join(model_types)}")
    print("\n" + "-" * 70)

    if not Path(input_file).exists():
        print(f"\nERROR: Input file not found: {input_file}")
        sys.exit(1)

    print(f"\nLoading RGI results...")
    df = pd.read_csv(input_file, sep='\t', low_memory=False)
    print(f"✓ Loaded {len(df):,} rows")
    print(f"✓ Found {len(df.columns)} columns")

    print(f"\nInitial statistics:")
    print(f"  - Total hits: {len(df):,}")

    col_aro   = cols["aro"]
    col_id    = cols["identity"]
    col_cov   = cols["coverage"]
    col_model = cols["model_type"]

    if col_aro in df.columns:
        print(f"  - Unique ARGs: {df[col_aro].nunique():,}")

    if col_id in df.columns:
        print(f"  - Identity range: {df[col_id].min():.1f}% - {df[col_id].max():.1f}%")
        print(f"  - Identity mean:  {df[col_id].mean():.1f}%")

    if col_cov in df.columns:
        print(f"  - Coverage range: {df[col_cov].min():.1f}% - {df[col_cov].max():.1f}%")
        print(f"  - Coverage mean:  {df[col_cov].mean():.1f}%")

    if col_model in df.columns:
        print(f"\n  Model types found:")
        for model, count in df[col_model].value_counts().items():
            print(f"    • {model}: {count:,} ({100*count/len(df):.1f}%)")

    print(f"\nApplying filters:")
    df_filtered = df.copy()

    # Filtre 1 : Identité
    if col_id in df.columns:
        before = len(df_filtered)
        df_filtered = df_filtered[df_filtered[col_id] >= min_identity]
        removed = before - len(df_filtered)
        print(f"  1. Identity ≥ {min_identity}%")
        print(f"     → Removed: {removed:,} | Retained: {len(df_filtered):,}")

    # Filtre 2 : Couverture
    if col_cov in df.columns:
        before = len(df_filtered)
        if before > 0:
            df_filtered = df_filtered[df_filtered[col_cov] >= min_coverage]
            removed = before - len(df_filtered)
            print(f"  2. Coverage ≥ {min_coverage}%")
            print(f"     → Removed: {removed:,} | Retained: {len(df_filtered):,}")
        else:
            print(f"  2. Coverage ≥ {min_coverage}% → Skipped (no hits remaining)")

    # Filtre 3 : Model type
    if col_model in df.columns and model_types:
        before = len(df_filtered)
        if before > 0:
            df_filtered = df_filtered[df_filtered[col_model].isin(model_types)]
            removed = before - len(df_filtered)
            print(f"  3. Model type filter")
            print(f"     → Removed: {removed:,} | Retained: {len(df_filtered):,}")
        else:
            print(f"  3. Model type filter → Skipped (no hits remaining)")

    # Normalisation des noms de colonnes
    print(f"\nCleaning column names...")
    df_filtered.columns = (df_filtered.columns
                           .str.replace(' ', '_')
                           .str.replace('(', '')
                           .str.replace(')', '')
                           .str.replace('%', 'percent')
                           .str.lower())
    print(f"  ✓ Standardized {len(df_filtered.columns)} column names")

    # Sauvegarde
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df_filtered.to_csv(output_file, sep='\t', index=False)
    print(f"  ✓ Saved → {output_file} ({output_path.stat().st_size/1024:.2f} KB)")

    # Warning résultat vide
    if len(df_filtered) == 0:
        print(f"\n{'!' * 70}")
        print(f"WARNING: No ARGs passed quality filters!")
        print(f"{'!' * 70}")
        if col_id in df.columns:
            print(f"\n  Identity range in raw data: "
                  f"{df[col_id].min():.1f}% - {df[col_id].max():.1f}%")
        print(f"  Consider lowering thresholds in config.yaml [{tool}]")

    print("\n" + "=" * 70)
    print("FILTERING COMPLETED")
    print("=" * 70)
    if len(df) > 0:
        print(f"  Input:  {len(df):,} hits")
        print(f"  Output: {len(df_filtered):,} hits ({100*len(df_filtered)/len(df):.1f}% retained)")
    print("=" * 70 + "\n")

    log.close()


if __name__ == "__main__":
    main()
