import requests
import pandas as pd

bed_file1 = 'new_intersect_1y.bed'
bed_file2 = 'closest_no_overlap_1y.bed'
columns = ['chr_peak', 'start_peak', 'end_peak', 'chr_gene', 'start_gene', 'end_gene',
           'ensembl_gene_id', 'score', 'strand']

# Read the BED file
df1 = pd.read_csv(bed_file1, sep='\t', names=columns)
df2 = pd.read_csv(bed_file2, sep='\t', names=columns)

df1['ensembl_gene_id'] = df1['ensembl_gene_id'].str.replace(';', '')
df2['ensembl_gene_id'] = df2['ensembl_gene_id'].str.replace(';', '')


# Annotating ensemble names of files with actual gene names
def fetch_gene_names(ensembl_ids):
    url = "https://rest.ensembl.org/lookup/id"
    headers = {"Content-Type": "application/json"}
    results = {}
    for ensembl_id in ensembl_ids:
        response = requests.get(f"{url}/{ensembl_id}", headers=headers)
        if response.ok:
            data = response.json()
            results[ensembl_id] = data.get('display_name', 'Unknown')
        else:
            results[ensembl_id] = 'Not Found'
    return results


ids = list(set(list(df1['ensembl_gene_id']) + list(df2['ensembl_gene_id'])))
gene_names = fetch_gene_names(ids)

# Apply the replacement to the specific column
df1['ensembl_gene_id'] = df1['ensembl_gene_id'].apply(lambda x: gene_names.get(x, x))
df2['ensembl_gene_id'] = df2['ensembl_gene_id'].apply(lambda x: gene_names.get(x, x))


print("Saving files 1")
df1.to_csv('intersection_annotated_with_genes_on_peaks.bed', sep='\t', header=False, index=False)
print("Saving files 2")
df2.to_csv('Nearest_promoter_to_peaks_with_genes.bed', sep='\t', header=False, index=False)
print("Saving files done")

columns = ['chr_peak', 'start_peak', 'end_peak', 'chr_gene', 'start_gene', 'end_gene',
           'gene_name', 'score', 'strand']
df1 = pd.read_csv("Nearest_promoter_to_peaks_with_genes.bed", sep='\t', names=columns)
df2 = pd.read_csv("intersection_annotated_with_genes_on_peaks.bed", sep='\t', names=columns)

# Merge DataFrames on the first three columns (chrom, start, end)
merged_df = pd.merge(df1, df2, on=['chr_peak', 'start_peak', 'end_peak'], how='inner', suffixes=('_left', '_right'))

# Show the merged DataFrame
tf_prom = merged_df[['gene_name_left', 'gene_name_right']]
tf_prom = tf_prom.drop_duplicates()
tf_prom = tf_prom[(tf_prom['gene_name_left'] != 'Not Found') & (tf_prom['gene_name_right'] != 'Not Found')]
tf_prom.columns = ['tf', 'nearestProm']
# print(tf_prom)

tf_prom.to_csv('atac_processed.csv', sep='\t', header=False, index=False)
