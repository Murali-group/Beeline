    # test package library and install
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE))
  BiocManager::install("EnsDb.Hsapiens.v86")
if (!requireNamespace("EnsDb.Mmusculus.v79", quietly = TRUE))
  BiocManager::install("EnsDb.Mmusculus.v79")
if (!requireNamespace("scater", quietly = TRUE))
  BiocManager::install("scater")
if (!requireNamespace("bluster", quietly = TRUE))
  BiocManager::install("bluster")
if (!requireNamespace("GenomeInfoDb", quietly = TRUE))
  BiocManager::install("GenomeInfoDb")
if (!requireNamespace("GenomeInfoDb", quietly = TRUE))
  BiocManager::install("GenomeInfoDb")
if (!requireNamespace("IRanges", quietly = TRUE))
  BiocManager::install("IRanges")
if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
  BiocManager::install("ComplexHeatmap")
if (!requireNamespace("rtracklayer", quietly = TRUE))
  BiocManager::install ("rtracklayer")
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
if (!requireNamespace("RColorBrewer", quietly = TRUE))
  install.packages('RColorBrewer')
if (!requireNamespace("reticulate", quietly = TRUE))
  install.packages("reticulate")
if (!requireNamespace("plyr", quietly = TRUE))
  install.packages('plyr')
if (!requireNamespace("dsb", quietly = TRUE))
  install.packages('dsb')
if (!requireNamespace("Seurat", quietly = TRUE))
  install.packages('Seurat')
if (!requireNamespace("Signac", quietly = TRUE))
  install.packages("Signac")
if (!requireNamespace("cluster", quietly = TRUE))
  install.packages("cluster")
if (!requireNamespace("igraph", quietly = TRUE))
  install.packages ("igraph")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages ("ggplot2")
if (!requireNamespace("Matrix", quietly = TRUE))
  install.packages ("Matrix")
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
if (!requireNamespace("tinytex", quietly = TRUE))
  install.packages("tinytex")
if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages ( "tidyverse" )
if (!requireNamespace("devtools", quietly = TRUE))
  library(devtools)
if (!requireNamespace("MAESTRO", quietly = TRUE))
  install_github("liulab-dfci/MAESTRO")
if (!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table")
if (!requireNamespace("igraph", quietly = TRUE))
  install.packages ("igraph")
if (!requireNamespace("parallel", quietly = TRUE))
  install.packages("parallel")
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
if (!requireNamespace("Hmisc", quietly = TRUE))
  install.packages("Hmisc")
if (!requireNamespace("CellChat", quietly = TRUE))
  devtools::install_github("sqjin/CellChat")
if (!requireNamespace("patchwork", quietly = TRUE))
  devtools::install_github("thomasp85/patchwork")
# Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
if (!requireNamespace("Seurat", quietly = TRUE))
  install.packages('Seurat')
if (!requireNamespace("plyr", quietly = TRUE))
  install.packages('plyr')
if (!requireNamespace("dsb", quietly = TRUE))
  install.packages('dsb')
if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
  BiocManager::install("ComplexHeatmap")
if (!requireNamespace("RColorBrewer", quietly = TRUE))
  install.packages('RColorBrewer')
if (!requireNamespace("reticulate", quietly = TRUE))
  install.packages("reticulate")
if (!requireNamespace("CellChat", quietly = TRUE))
  devtools::install_github("sqjin/CellChat")
if (!requireNamespace("patchwork", quietly = TRUE))
  devtools::install_github("thomasp85/patchwork")
if (!requireNamespace("scater", quietly = TRUE))
  BiocManager::install("scater")

library(MAESTRO)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(scater)
library(Seurat)
library(Signac)
library(cluster)
library(bluster)
library(GenomeInfoDb)
library(igraph)
library(GenomicRanges)
library(IRanges)
library(ggplot2)
library(Matrix)
library(dplyr)
library(tinytex)
library(tidyverse)
library(rtracklayer)
library(reticulate)
library(dplyr)
library(parallel)
library(igraph)
library(data.table)
library(Hmisc)
library(dplyr)
library(Seurat)
library(patchwork)
library(cluster)

## Setup path and environment
#source("/Algorithms/DEEPMAPS/scRNA_scATAC1.r")
source("/scRNA_scATAC1.r")
# set python environment
Sys.setenv(RETICULATE_PYTHON = "/home/user/miniconda/bin/python")
use_python("/home/user/miniconda/bin/python")
py_config()
lisa_path <- "/lisa_output/"
jaspar_path <- "/home/user/miniconda/lib/python3.8/site-packages/lisa/data/"

option_list <- list(
              make_option(c("-e","--expressionFile"), type = 'character',
              help= "Path to the heirarchical file with rna matrix and atac files as matrices. The gene-by-cell sample
              matrix with cell names as the first row and gene names as the first column. The peaks-by-cell sample matrix
              is also inculded. Usually an h5 or RDS file."),

              make_option(c("-r","--rnaFile"), type = 'character',
              help= "Path to comma separated file or tab seperated file containing gene-by-cell
              matrix with cell names as the first row and gene names as  the first column. "),

              make_option(c("-a","--atacFile"), type = 'character',
              help= "Path to ATAC-seq dataset with Chromatin Accessibility. Or an file with peaks by cell sample formatted."),

              make_option(c("-c","--veloFile"), type = 'character',
              help= "Path to the gene by peak matrix file ATAC-seq dataset with Chromatin Accessibility. Or an file with peaks by cell sample formatted."),

              make_option(c("-o","--outPrefix"), type = 'character',
              help= "Path to write output files. Required."),

              make_option(c("","--outFile"), default = 'MICA_Out', type = 'character',
              help= "outFile name to write the output ranked edges. Required.")

              )

parser <- OptionParser(option_list = option_list)
arguments <- parse_args(parser, positional_arguments = FALSE)


if (length(arguments$outPrefix) == 0){
 stop('outPrefix is required.
      Run Rscript runMICA.R -h for more details.')
}

if (length(arguments$expressionFile) == 0 and length(arguments$rnaFile) + length(arguments$atacFile) < 2 ){
    stop("Please enter the path to expression file or matrix (single cell RNA and ATAC) using the -e flag.
          Run Rscript runMICA.R -h for more details.")
}

# Based on the input files provided reading and creating object
if (length(arguments$expressionFile) == 0 and length(arguments$rnaFile) + length(arguments$atacFile) > 1 ){
   # Matrix is read
   rna_matrix = readRDS(arguments$rnaFile)
   atac_matrix = readRDS(arguments$atacFile)
   lymph_obj <- ReadData(rna_matrix = rna_matrix, atac_matrix = atac_matrix, dataFormat = "scRNA_scATAC", min_cell = 0.001)
}

if (length(arguments$expressionFile) != 0 and length(arguments$rnaFile) + length(arguments$atacFile) > 1 ){
   rna_matrix = readRDS(arguments$rnaFile)
   atac_matrix = readRDS(arguments$atacFile)
   lymph_obj <- ReadData(rna_matrix = rna_matrix, atac_matrix = atac_matrix, dataFormat = "scRNA_scATAC",
        dataFormat = "h5", min_cell = 0.001, h5Path = arguments$expressionFile )

}

if (length(arguments$expressionFile) != 0 and length(arguments$rnaFile) + length(arguments$atacFile) < 2 ){
    lymph_obj <- ReadData(dataFormat = "scRNA_scATAC", dataFormat = "h5", min_cell = 0.001, h5Path = arguments$expressionFile)
}

ATAC_gene_peak <- CalGenePeakScore(peak_count_matrix = lymph_obj@assays$ATAC@counts,organism = "GRCh38")

# Calculate gene active score (integration)
velo <-arguments$veloFile
GAS_obj <- calculate_GAS_v1(ATAC_gene_peak = ATAC_gene_peak , obj = lymph_obj, method = "velo", veloPath = velo)
GAS <- GAS_obj[[1]]
lymph_obj <- GAS_obj[[2]]

# Perform heterogeneous graph transformer (HGT) model on GAS matrix
HGT_result <- run_HGT(GAS = as.matrix(GAS),result_dir=arguments$outPrefix, data_type='scRNA_scATAC', envPath=NULL,
 lr=0.2, epoch=30, n_hid=128, n_heads=16)

# Cluster the cells
cell_hgt_matrix <- HGT_result[['cell_hgt_matrix']]
rownames(cell_hgt_matrix) <- colnames(GAS)

lymph_obj <- lymph_obj[, colnames(GAS)]
cell_hgt_matrix <- cell_hgt_matrix[colnames(GAS),]

HGT_embedding <-
  CreateDimReducObject(embeddings = cell_hgt_matrix,
                       key = "HGT_",
                       assay = "RNA")
lymph_obj@reductions[['HGT']] <- HGT_embedding
lymph_obj <-
  FindVariableFeatures(lymph_obj, selection.method = "vst", nfeatures = 2000)
lymph_obj <- ScaleData(lymph_obj, features = VariableFeatures(lymph_obj))
lymph_obj <-
  RunUMAP(
    lymph_obj,
    reduction = 'HGT',
    dims = 1:ncol(cell_hgt_matrix),
    reduction.name = "umap.rna",
    reduction.key = "rnaUMAP_"
  )
lymph_obj <-
  FindNeighbors(lymph_obj,
                reduction = "HGT",
                dims = 1:ncol(cell_hgt_matrix))
lymph_obj <- FindClusters(lymph_obj, resolution = 1)
graph.out <- as.factor(lymph_obj$seurat_clusters)

DefaultAssay(lymph_obj) <- "RNA"
png("plot.png")
DimPlot(lymph_obj, reduction = 'umap.rna')
dev.off()

# Calculate cell cluster active gene modules and run LISA for TF infer.
co <- get_gene_module(obj = lymph_obj, GAS = GAS, att = HGT_result[['attention']],method = 'SFP' )

dir.create(lisa_path, showWarnings = F)
write_GM(co = co, lisa_path = lisa_path)

# Run LISA on generated gene modules
system(
  paste0(
    "/home/user/miniconda/bin/python /run_lisa.py --path ",
    #"/home/user/miniconda/bin/python Algorithms/DEEPMAPS/run_lisa.py --path ",
    lisa_path,
    " --species ",
    "hg38"
  )
)

# Filter gene with no accessible peak in promoter
gene_peak_pro <- AccPromoter(obj = lymph_obj, gene_peak = ATAC_gene_peak, GAS = GAS, species = 'hg38')

pre_regulon_res <- Calregulon(GAS = GAS, co = co,gene_peak_pro = gene_peak_pro, species = "hg38", jaspar_path = jaspar_path, lisa_path = lisa_path)
BA_score <- pre_regulon_res[[1]]
ct_regulon_v1 <- pre_regulon_res[[2]]
TFinGAS<- pre_regulon_res[[3]]

# Combine same TFs
peak_TF <- uni(gene_peak_pro = gene_peak_pro, BA_score = BA_score)

# Calculate regulatory Intensive (RI) score in cell level and infer cell type active regulon
RI_C <- RI_cell(obj = lymph_obj, ct_regulon = ct_regulon_v1, GAS = GAS, gene_peak_pro = gene_peak_pro, peak_TF = peak_TF, graph.out = graph.out)

# Calculate regulon active score1 (RAS1)
regulon_res <- calRAS(RI_C = RI_C, ct_regulon = ct_regulon_v1, graph.out = graph.out, TFinGAS = TFinGAS)
RAS_CT <- regulon_res[[1]]
RI_CT <- regulon_res[[2]]
ct_regulon <- regulon_res[[3]]
RAS_C1 <- regulon_res[[4]]

# Calculate RAS2 same topology in a row
RAS_C2 <- CalRAS2(ct_regulon, graph.out)

# Calculate master TF
masterTF <- masterFac(ct_regulon = ct_regulon, RI_CT = RI_CT)
TF_cen <- masterTF[[1]]
gene_cen <- masterTF[[2]]
network <- masterTF[[3]]


# Constructing GRN for ALL cell types :
for (ct_name in names(network)) {
  # Extract the adjacency matrix for the current cell type
  network_matrix <- network[[ct_name]]

  # Save the adjacency matrix to a CSV file
  write.csv(network_matrix, file = paste0(arguments$outPrefix,arguments$outFile,ct_name,"_network.csv"), row.names = TRUE)

  # Convert the adjacency matrix to an edge list
  # Get the indices of non-zero entries
  non_zero_indices <- which(network_matrix != 0, arr.ind = TRUE)

  # Create a data frame for the edge list
  edge_list <- data.frame(
    source = rownames(network_matrix)[non_zero_indices[, 1]],
    target = colnames(network_matrix)[non_zero_indices[, 2]],
    weight = network_matrix[non_zero_indices]
  )

  # Save the edge list to a CSV file
  write.csv(edge_list, file = paste0(arguments$outPrefix,arguments$outFile,ct_name,"_edge_list.csv"), row.names = FALSE)
}


# Incorporating tf and gene centrality outcomes into the network

# for (ct_name in names(network)) {
#   # Extract the adjacency matrix and centrality measures
#   network_matrix <- network[[ct_name]]
#   tf_cent <- TF_cen[[ct_name]]
#   gene_cent <- gene_cen[[ct_name]]
#
#   # Enhance the network by multiplying edges by centrality scores
#   enhanced_network_matrix <- network_matrix * outer(tf_cent, gene_cent)
#
#   # Save the enhanced adjacency matrix
#   write.csv(enhanced_network_matrix, file = paste0(ct_name, "_enhanced_network.csv"), row.names = TRUE)
#
#   # Convert to edge list
#   non_zero_indices <- which(enhanced_network_matrix != 0, arr.ind = TRUE)
#   edge_list <- data.frame(
#     source = rownames(enhanced_network_matrix)[non_zero_indices[, 1]],
#     target = colnames(enhanced_network_matrix)[non_zero_indices[, 2]],
#     weight = enhanced_network_matrix[non_zero_indices]
#   )
#
#   # Save the enhanced edge list
#   write.csv(edge_list, file = paste0(ct_name, "_enhanced_edge_list.csv"), row.names = FALSE)
# }









