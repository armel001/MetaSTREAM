import pandas as pd

df = pd.read_csv('results/summary/rgi_bwt_all_samples.tsv', sep='\t')

unique_genes = (df[['aro_accession', 'aro_term', 'amr_gene_family', 'drug_class']]
                .drop_duplicates(subset='aro_accession')
                .sort_values('aro_term'))
unique_genes['clean_name'] = ''

output_path = 'results/summary/unique_genes_FOR_REVIEW.xlsx'
unique_genes.to_excel(output_path, index=False)

print(f"Saved: {output_path}")
print(f"{len(unique_genes)} gènes uniques exportés pour révision")
