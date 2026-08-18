#!/usr/bin/env python3
"""
combine_bracken.py — Combine Bracken outputs across samples into an abundance matrix.
"""

import sys
from pathlib import Path
import pandas as pd

input_files     = snakemake.input
output_matrix   = snakemake.output.abundance_matrix
output_metadata = snakemake.output.metadata
level           = snakemake.wildcards.level

if not input_files:
    sys.exit("ERROR: No input files provided")

REQUIRED_COLS = ['name', 'taxonomy_id', 'new_est_reads', 'fraction_total_reads']

print(f"══ Bracken combine [{level}] ══  {len(input_files)} files\n")

# ── Load samples ──────────────────────────────────────────────────────────────
all_data, sample_stats = [], []

for filepath in input_files:
    sid = Path(filepath).parent.parent.name
    try:
        df = pd.read_csv(filepath, sep='\t')
        missing = [c for c in REQUIRED_COLS if c not in df.columns]
        if missing:
            print(f"  [WARN] {sid} — missing columns: {missing}", file=sys.stderr)
            continue
        sub = df[REQUIRED_COLS].copy()
        sub['sample'] = sid
        all_data.append(sub)
        sample_stats.append({'sample': sid, 'n_taxa': len(df),
                              'total_reads': int(df['new_est_reads'].sum()),
                              'taxonomic_level': level})
        print(f"  {sid}: {len(df)} taxa | {int(df['new_est_reads'].sum()):,} reads")
    except Exception as e:
        print(f"  [ERROR] {sid} — {e}", file=sys.stderr)

if not all_data:
    sys.exit("ERROR: No data loaded")

# ── Build abundance matrix (taxa × samples) ───────────────────────────────────
combined = pd.concat(all_data, ignore_index=True)

abundance = (combined
    .pivot_table(index='name', columns='sample',
                 values='new_est_reads', fill_value=0)
    .join(combined[['name', 'taxonomy_id']].drop_duplicates().set_index('name')))

# taxonomy_id first column, sorted by total abundance descending
abundance = abundance[['taxonomy_id'] + [c for c in abundance.columns if c != 'taxonomy_id']]
abundance = abundance.loc[abundance.drop('taxonomy_id', axis=1).sum(axis=1)
                                   .sort_values(ascending=False).index]

# ── Save ──────────────────────────────────────────────────────────────────────
abundance.to_csv(output_matrix, sep='\t')
print(f"\n  [OK] Matrix  : {output_matrix}  ({len(abundance)} taxa × {len(abundance.columns)-1} samples)")

metadata = pd.DataFrame(sample_stats).sort_values('sample')
metadata.to_csv(output_metadata, sep='\t', index=False)
print(f"  [OK] Metadata: {output_metadata}")

# ── Summary ───────────────────────────────────────────────────────────────────
print(f"\n══ SUMMARY [{level}] ══")
print(f"  Taxa     : {len(abundance)}")
print(f"  Samples  : {len(sample_stats)}")
print(f"  Total reads  : {metadata['total_reads'].sum():,}")
print(f"  Mean reads/sample : {metadata['total_reads'].mean():.0f}")
print(f"  Mean taxa/sample  : {metadata['n_taxa'].mean():.1f}\n")
