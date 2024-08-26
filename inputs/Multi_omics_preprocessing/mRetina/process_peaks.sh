#!/bin/bash

# Download and unzip peak calls file
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6062732&format=file&file=GSM6062732%5Fmultiomics%5Fmouse%5Fwt%5FGeneExpression%5FPeaks%5Flists%2ERDS%2Egz
gunzip GSM6062732_multiomics_mouse_wt_GeneExpression_Peaks_lists.RDS.gz

# Download the GTF file
wget ftp://ftp.ensembl.org/pub/current_gtf/mus_musculus/Mus_musculus.GRCm39.104.gtf.gz
# Unzip the GTF file
gunzip Mus_musculus.GRCm39.104.gtf.gz

# Input RDS file
RDS_FILE=GSM6062732_multiomics_mouse_wt_GeneExpression_Peaks_lists.RDS

# Run the R script with the input RDS file
Rscript process_peaks.R "$RDS_FILE"

echo "Processing completed. Output files: gene_expr_matrix.csv and peaks.bed"

# Install necessary tools
sudo apt-get update
sudo apt-get install -y bedops bedtools

# Convert GTF to BED format
gff2bed < Mus_musculus.GRCm39.104.gtf > genes.bed

# Filter peaks and process BED files
less -S peaks.bed | grep -v "GL" | grep -v "JL" | grep -v "JH" > peaks_filtered.bed

bedtools sort -i peaks_filtered.bed > peaks_sorted.bed
bedtools sort -i genes.bed > genes_sorted.bed

# Find the closest gene to each peak
# bedtools closest -a peaks_sorted.bed -b genes_sorted.bed -wa -wb > closest_genes.bed

# Find the closest gene to each peak with no overlap
bedtools closest -io -a peaks_sorted.bed -b genes_sorted.bed > closest_no_overlap.bed

# Intersect peaks with genes
bedtools intersect -a peaks_filtered.bed -b genes.bed -wa -wb > intersect_genes.bed

python final_process.py

# Clear intermediate files
rm *peaks.bed *genes.bed *closest_1y.bed *intersection_annotated_with_genes_on_peaks.bed  *intersect_genes.bed
rm *genes_sorted.bed *peaks_filtered.bed *Nearest_promoter_to_peaks_with_genes.bed *closest_no_overlap.bed