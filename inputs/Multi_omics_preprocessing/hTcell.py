import pandas as pd
import scanpy as sc
import scipy
from scipy.stats import spearmanr
import requests

# Input files obtained from https://zenodo.org/records/6383269#.ZCWsCOzMJqs
# The preprocessing steps where tested on cd34 cells, similar steps can be followed for the bone marrow

rna_data = sc.read_h5ad('cd34_multiome_rna.h5ad')
atac_data = sc.read_h5ad('cd34_multiome_atac.h5ad')

# filter counts
sc.pp.filter_cells(atac_data, min_counts=1000)
sc.pp.filter_cells(rna_data, min_counts=1000)

sc.pp.filter_genes(atac_data, min_counts=10)
sc.pp.filter_genes(rna_data, min_counts=10)

atac_data.X = sc.pp.normalize_total(atac_data, target_sum=1e4, inplace=False)['X']
atac_data.X = scipy.sparse.csr_matrix(atac_data.X)

sc.pp.normalize_total(rna_data, target_sum=1e4)
sc.pp.log1p(rna_data)

matching_genes = rna_data.var.index.intersection(atac_data.var['nearestGene'])
#gene_peak_correlations = []
gene_peak_corr_df = pd.DataFrame(columns=['gene', 'corr'])

common_cells = atac_data.obs.index.intersection(rna_data.obs.index)
atac_data = atac_data[common_cells, :]
rna_data = rna_data[common_cells, :]

# Generating regulators file
# finding correlations with respective peaks for genes in matching genes

for gene in matching_genes:
    atac_indices = atac_data.var['nearestGene'] == gene
    rna_expression = rna_data[:, gene].X.toarray()
    atac_signal = atac_data[:, atac_indices].X.toarray()
    # Extract RNA-seq expression data for the gene
    rna_expression = rna_data[:, gene].X.toarray()
    # Flatten arrays to 1D for correlation calculation
    atac_signal_flat = atac_signal.flatten()
    rna_expression_flat = rna_expression.flatten()
    # Ensure both arrays have the same size after flattening
    if atac_signal_flat.size == rna_expression_flat.size:
        correlation, p_value = spearmanr(atac_signal_flat, rna_expression_flat)
        gene_peak_corr_df.append({gene: correlation}, ignore_index=True)

gene_peak_corr_sorted = gene_peak_corr_df.sort_values(by='corr', ascending=False)
gene_peak_corr_sorted.to_csv("regulators.csv", index=False)
print("Regulators file has been created. ")


# Generating tf_binding ot processed ATAC file
def convert_ucsc_to_ensembl(ucsc_ids):
    server = "https://rest.ensembl.org"
    endpoint = "/xrefs/symbol/human/{}"

    gene_names = {}

    for ucsc_id in ucsc_ids:
        # Make a GET request to Ensembl REST API
        response = requests.get(server + endpoint.format(ucsc_id), headers={"Content-Type": "application/json"})

        if response.ok:
            data = response.json()
            # Extract the gene name
            if data:
                gene_names[ucsc_id] = data[0]["id"]
            else:
                gene_names[ucsc_id] = "Not found"
        else:
            gene_names[ucsc_id] = "Error"

    return gene_names


# list unique nearest TSS column

ucsc = list(set(atac_data.var['nearestTSS']))
ensembls = convert_ucsc_to_ensembl(ucsc)

# Unique ensemble values obtained

ensembl_results = list(set(ensembls.values()))


def ensembl_to_gene_names(ensembl_ids):
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


gene_names_dict = ensembl_to_gene_names(ensembl_results)

# Map UCSC to gene names :
ucsc_to_gene_name = {ucsc_id: gene_names_dict.get(ensembl_id, ensembl_id) for ucsc_id, ensembl_id in ensembls.items()}
atac_data.var['nearestTSS'] = atac_data.var['nearestTSS'].map(ucsc_to_gene_name)

atac_tfs = atac_data.var[['nearestGene', 'nearestTSS', 'score']]
atac_tf_filtered = atac_tfs[(atac_tfs['nearestTSS'] != 'Unknown') & (atac_tfs['nearestTSS'] != 'Not Found')]
atac_tf_sorted = atac_tf_filtered.sort_values(by='score', ascending=False)
atac_tf_sorted.to_csv("atac.csv", index=False)

# Writing the RNA seq file
rna_matrix = rna_data.X

# Convert to a DataFrame for easier manipulation and export
rna_df = pd.DataFrame(rna_matrix.toarray(), index=rna_data.obs_names, columns=rna_data.var_names)

# Save the RNA-seq matrix to a CSV file
rna_df.to_csv('rna_seq_matrix.csv')

