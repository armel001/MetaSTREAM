#!/usr/bin/env python3
"""
Selectionne uniquement les colonnes utiles des resultats RGI filtres
"""

import pandas as pd
import sys
from pathlib import Path

def main():
    tool = getattr(snakemake.wildcards, 'tool', 'rgi')  # fallback safe
    # Recuperer les parametres
    input_file = snakemake.input.filtered
    print(f"RGI Column Selection  [{tool}]")
    output_file = snakemake.output.selected
    columns_to_keep = snakemake.params.columns
    log_file = snakemake.log[0]
    
    # Configurer le logging
    log = open(log_file, 'w')
    sys.stderr = sys.stdout = log
    
    print("=" * 70)
    print("RGI Column Selection")
    print("=" * 70)
    print(f"\nInput:  {input_file}")
    print(f"Output: {output_file}")
    print(f"\nColumns to keep: {len(columns_to_keep)}")
    for col in columns_to_keep:
        print(f"  • {col}")
    print("\n" + "-" * 70)
    
    # Verifier que le fichier existe
    if not Path(input_file).exists():
        print(f"\nERROR: Input file not found: {input_file}")
        sys.exit(1)
    
    # Charger les donnees
    print(f"\nLoading filtered RGI results...")
    try:
        df = pd.read_csv(input_file, sep='\t', low_memory=False)
        print(f"✓ Loaded {len(df):,} rows")
        print(f"✓ Found {len(df.columns)} columns")
    except Exception as e:
        print(f"\nERROR: Failed to load file")
        print(f"  {str(e)}")
        sys.exit(1)
    
    # Afficher les colonnes disponibles
    print(f"\nAvailable columns in input file:")
    for i, col in enumerate(df.columns, 1):
        print(f"  {i:2d}. {col}")
    
    # Verifier que toutes les colonnes demandees existent
    print(f"\nChecking requested columns...")
    missing_columns = []
    available_columns = []
    
    for col in columns_to_keep:
        if col in df.columns:
            available_columns.append(col)
            print(f"  ✓ {col}")
        else:
            missing_columns.append(col)
            print(f"  ✗ {col} (NOT FOUND)")
    
    if missing_columns:
        print(f"\nWARNING: {len(missing_columns)} requested columns not found:")
        for col in missing_columns:
            print(f"  • {col}")
        print(f"\nContinuing with {len(available_columns)} available columns...")
    
    # Si aucune colonne valide, erreur
    if len(available_columns) == 0:
        print(f"\nERROR: None of the requested columns exist in the input file!")
        print(f"\nPlease check your column names in config.yaml")
        sys.exit(1)
    
    # Select colonnes
    print(f"\nSelecting columns...")
    df_selected = df[available_columns].copy()
    
    print(f"  ✓ Selected {len(df_selected.columns)} columns")
    print(f"  ✓ Retained {len(df_selected):,} rows")
    
    # Statistiques
    print(f"\nColumn statistics:")
    print(f"  - Original columns: {len(df.columns)}")
    print(f"  - Selected columns: {len(df_selected.columns)}")
    print(f"  - Reduction: {len(df.columns) - len(df_selected.columns)} columns removed")
    
    # Calcul
    import io
    original_size = len(df.to_csv(sep='\t', index=False))
    selected_size = len(df_selected.to_csv(sep='\t', index=False))
    reduction_percent = 100 * (1 - selected_size / original_size)
    
    print(f"\nFile size reduction:")
    print(f"  - Original size: {original_size / 1024:.2f} KB")
    print(f"  - Selected size: {selected_size / 1024:.2f} KB")
    print(f"  - Reduction: {reduction_percent:.1f}%")
    
    # Data preview
    print(f"\nData preview:")
    print(df_selected.head(3).to_string())
    
    # Creer le repertoire de sortie
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Sauvegarder
    print(f"\nSaving selected columns...")
    try:
        df_selected.to_csv(output_file, sep='\t', index=False)
        file_size = output_path.stat().st_size / 1024
        print(f"  ✓ Saved to: {output_file}")
        print(f"  ✓ File size: {file_size:.2f} KB")
    except Exception as e:
        print(f"\nERROR: Failed to save file")
        print(f"  {str(e)}")
        sys.exit(1)
    
    # Resume
    print("\n" + "=" * 70)
    print("COLUMN SELECTION COMPLETED SUCCESSFULLY")
    print("=" * 70)
    print(f"  Rows:    {len(df_selected):,}")
    print(f"  Columns: {len(df_selected.columns)}")
    print(f"  Size:    {file_size:.2f} KB")
    print("=" * 70 + "\n")
    
    log.close()


if __name__ == "__main__":
    main()
