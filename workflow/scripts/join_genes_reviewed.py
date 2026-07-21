import pandas as pd

df = pd.read_csv('results/summary/rgi_bwt_all_samples.tsv', sep='\t')
mapping = pd.read_excel('results/summary/unique_genes_REVIEWED.xlsx')

# Vérifier les clean_name laissés vides
empty_mask = mapping['clean_name'].isna() | (mapping['clean_name'].astype(str).str.strip() == '')
if empty_mask.any():
    print(f"ATTENTION : {empty_mask.sum()} gènes sans clean_name assigné :")
    print(mapping.loc[empty_mask, ['aro_accession', 'aro_term']].to_string(index=False))
    print("Ces lignes garderont leur aro_term d'origine.\n")
    mapping.loc[empty_mask, 'clean_name'] = mapping.loc[empty_mask, 'aro_term']

# Informer (pas bloquant) des fusions volontaires de plusieurs accessions
# vers un même nom — agrégées correctement en aval par le pipeline R
dup_counts = mapping['clean_name'].value_counts()
collisions = dup_counts[dup_counts > 1]
if len(collisions) > 0:
    print(f"Note : {len(collisions)} noms regroupent plusieurs ARO accessions (OK, seront sommés) :")
    for name, n in collisions.items():
        accs = mapping.loc[mapping['clean_name'] == name, 'aro_accession'].tolist()
        print(f"  '{name}' <- {n} accessions : {accs}")

# Jointure par aro_accession (ID CARD stable — pas de risque de mismatch texte)
mapping_lookup = mapping.set_index('aro_accession')['clean_name'].to_dict()

df['aro_term_original'] = df['aro_term']          # traçabilité conservée
df['aro_term'] = df['aro_accession'].map(mapping_lookup).fillna(df['aro_term'])

output_path = 'results/summary/rgi_bwt_all_samples_CLEAN.tsv'
df.to_csv(output_path, sep='\t', index=False)
print(f"\nSauvegardé : {output_path}")
print(f"  {len(df)} lignes, {df['aro_term'].nunique()} noms de gènes uniques après nettoyage")
