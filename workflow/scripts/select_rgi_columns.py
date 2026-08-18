#!/usr/bin/env python3
"""
select_rgi_columns.py — Retain only useful columns from filtered RGI results.
"""

import sys
from pathlib import Path
import pandas as pd


def main():
    tool            = getattr(snakemake.wildcards, 'tool', 'rgi')
    input_file      = snakemake.input.filtered
    output_file     = snakemake.output.selected
    columns_to_keep = snakemake.params.columns

    log = open(snakemake.log[0], 'w')
    sys.stderr = sys.stdout = log

    print(f"══ RGI Column Selection [{tool}] ══\n  Input : {input_file}\n  Output: {output_file}")
    print(f"  Requested columns ({len(columns_to_keep)}): {', '.join(columns_to_keep)}\n")

    # ── Load ──────────────────────────────────────────────────────────────────
    if not Path(input_file).exists():
        sys.exit(f"ERROR: Input file not found: {input_file}")

    try:
        df = pd.read_csv(input_file, sep='\t', low_memory=False)
        print(f"  Loaded {len(df):,} rows × {len(df.columns)} columns")
    except Exception as e:
        sys.exit(f"ERROR: Failed to load file — {e}")

    # ── Column validation ─────────────────────────────────────────────────────
    available = [c for c in columns_to_keep if c in df.columns]
    missing   = [c for c in columns_to_keep if c not in df.columns]

    if missing:
        print(f"  [WARN] {len(missing)} column(s) not found: {', '.join(missing)}")
    if not available:
        sys.exit("ERROR: None of the requested columns exist — check config.yaml")

    # ── Select & save ─────────────────────────────────────────────────────────
    df_sel = df[available]
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    try:
        df_sel.to_csv(output_file, sep='\t', index=False)
        file_kb = Path(output_file).stat().st_size / 1024
    except Exception as e:
        sys.exit(f"ERROR: Failed to save — {e}")

    # ── Summary ───────────────────────────────────────────────────────────────
    n_removed = len(df.columns) - len(available)
    print(f"\n  Columns : {len(df.columns)} → {len(available)} selected ({n_removed} removed)")
    print(f"  Rows    : {len(df_sel):,}")
    print(f"  Size    : {file_kb:.2f} KB")
    print(f"\n  Preview :\n{df_sel.head(3).to_string()}")
    print(f"\n══ COLUMN SELECTION COMPLETED ══\n")

    log.close()


if __name__ == "__main__":
    main()
