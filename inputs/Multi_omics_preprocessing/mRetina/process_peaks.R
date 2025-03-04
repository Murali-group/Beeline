
# The dataset was downloaded from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6062732&format=file&file=GSM6062732%5Fmultiomics%5Fmouse%5Fwt%5FGeneExpression%5FPeaks%5Flists%2ERDS%2Egz
# The multiomics file is used for the preprocessing and GRN construction
# The Mus musculus gtf file is downloaded from
# This is used for annotation of peaks and genes.


# library(Seurat)
# atac_data <- readRDS("GSM6062731_retina_10x.atac.RDS")
#
# # Load RNA-seq data
# rna_data <- readRDS("GSM6062732_retina_10x.rna.RDS")

# Load multiomics data
multi_data <- readRDS("GSM6062732_multiomics_mouse_wt_GeneExpression_Peaks_lists.RDS")
gene_expr_matrix <- multi_data[["Gene Expression"]]

write.csv(gene_expr_matrix, "gene_expr_matrix.csv")

# ATAC data
peaks_matrix <- multi_data[["Peaks"]]
# Extracting the row names which are in the format "chr:start-end"
peak_coords <- rownames(peaks_matrix)

# Split the coordinates into chr, start, and end
peak_coords_split <- strsplit(peak_coords, "[:-]")

# Create a data frame with chr, start, and end columns
peaks_bed <- do.call(rbind, peak_coords_split)
colnames(peaks_bed) <- c("chr", "start", "end")
peaks_bed <- data.frame(
  chr = peaks_bed[, 1],
  start = as.integer(peaks_bed[, 2]),
  end = as.integer(peaks_bed[, 3])
)

# Write the BED file
write.table(peaks_bed, "peaks.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



