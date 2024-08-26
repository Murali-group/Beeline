# package library and install
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
## ----------------------------------------------------------------------------------------------------------------
# Read data
##' Read matched scRNA + scATAC data from H5 file
#input:
# h5Path: the path of h5 file
# rna_matrix: an expression matrix with gene * cell
# atac_matrix: an accessibility matrix with peak * cell
# min_cell: the peak / gene will be removed if the value in the gene / peak with more than min_cell cell is equal to zero
# dataFormat: the format of input data, default as "h5"
#output:
#a seurat object
ReadData <- function(h5Path = NULL, rna_matrix = NULL, atac_matrix =NULL, min_cell = 0.1, dataFormat = "h5"){
  if (dataFormat == "h5") {
    tmp_obj <- Read10Xdata(h5Path = h5Path, min_cell = min_cell)
  }else{
    tmp_obj <- readmatrix(rna_matrix = rna_matrix, atac_matrix = atac_matrix, min_cell = min_cell)
  }
  tmp_obj <- filterCell(tmp_obj)
  return(tmp_obj)
}


# Creat a seuart object
# input:
#  h5Path: address of  h5 file to read. (When dataFormat is 'h5', it couldn't be NULL.)
#  rna_matrix: a RNA matrix where the rows correspond to genes and the columns correspond to cells. (When dataFormat is 'matrixs', it couldn't be NULL.)
#  atac_matrix: a matrix where the rows correspond to peaks and the columns correspond to cells. (When data_type is 'RNA_ATAC'and dataFormat is 'matrixs', it couldn't be NULL.)
#  adt_matrix: a matrix where the rows correspond to proteins and the columns correspond to cells. (When data_type is 'CITE' and dataFormat is 'matrixs', it couldn't be NULL.)
#  data_type: 'CITE', 'RNA_ATAC'
#  dataFormat: 'matrixs' or 'h5'
#  min_cell: the peak/gene will be removed if the value in the gene/peak with more than min_cell cell is equal to zero
#  nmad: a numeric scalar, specifying the minimum number of MADs away from median required for a value to be called an outlier
#  gene.filter: if do gene filtering
#  gene.filter: if do cell filtering

# output:
#  obj: a seuart object for cite-seq data after normalizing and scaling (When data_type is 'CITE'); a seuart object for scRNA-seq and scATAC-seq (When data_type is 'RNA_ATAC')


ReadData <- function(h5Path = NULL, rna_matrix = NULL, atac_matrix = NULL, adt_matrix = NULL, data_type = NULL, dataFormat = NULL, min_cell=0.1, nmad=3, gene.filter=TRUE, cell.filter=TRUE){
  if (data_type=='CITE'){
    # obtain RNA and ADT matrixs
    if (dataFormat=='matrixs'){
      rna <- rna_matrix
      adt <- adt_matrix
    }else if(dataFormat=='h5'){
      h5 <- Read10X_h5(h5Path)
      rna <- h5$`Gene Expression`
      adt <- h5$`Antibody Capture`
      
    }
    # gene filtering
    if (gene.filter==TRUE){
      binaryrna <- rna
      binaryadt <- adt
      binaryrna[binaryrna>0] <-1
      binaryadt[binaryadt>0] <-1
      rna <- rna[which(rowSums(binaryrna) > ncol(binaryrna)*min_cell),]
      adt <- adt[which(rowSums(binaryadt) > ncol(binaryadt)*min_cell),]
    }
    
    #Setup a Seurat object, add the RNA and protein data
    obj <- CreateSeuratObject(counts = rna)
    obj [["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    adt_assay <- CreateAssayObject(counts = adt)
    obj[["ADT"]] <- adt_assay
    # cell filtering
    if (cell.filter==TRUE){
      adtn<-isOutlier(
        obj$nCount_ADT,
        nmads = nmad,
        log = F,
        type = "both"
      )
      rnan<-isOutlier(
        obj$nCount_RNA,
        nmads = nmad,
        log = F,
        type = "both"
      )
      mito<-isOutlier(
        obj$percent.mt,
        nmads = nmad,
        log = F,
        type = "both"
      )
      obj <-
        AddMetaData(obj, adtn, col.name = "adtn")
      obj <-
        AddMetaData(obj, rnan, col.name = "rnan")
      obj <-
        AddMetaData(obj, mito, col.name = "mito")
      obj<-subset(
        x = obj,
        subset = adtn == F &
          rnan == F &
          mito == F
      )
    }
    # normalization and scaling 
    DefaultAssay(obj) <- 'RNA'
    obj <- NormalizeData(obj) %>% FindVariableFeatures()
    all.genes <- rownames(obj)
    obj <- ScaleData(obj, features = all.genes)
    DefaultAssay(obj) <- 'ADT'
    VariableFeatures(obj) <- rownames(obj[["ADT"]])
    obj <- NormalizeData(obj, normalization.method = 'CLR', margin = 2) %>% 
      ScaleData()
    DefaultAssay(obj) <- 'RNA'
  }else if (data_type=='scRNA_scATAC'){
    if (gene.filter==FALSE){
      min_cell = 0
    }
    if (dataFormat == "h5") {
      obj <- Read10Xdata(h5Path = h5Path, min_cell = min_cell)
    }else{
      obj <- readmatrix(rna_matrix = rna_matrix, atac_matrix = atac_matrix, min_cell = min_cell)
    }
    if (cell.filter==TRUE){
      obj <- filterCell(obj, nmad= nmad)
    }
  }
  return(obj)
}




# Read 10X data
##' Read matched scRNA + scATAC data from H5 file
#input:
#  h5Path: the path of h5 file
#  min_cell: the peak / gene will be removed if the value in the gene / peak with more than min_cell cell is equal to zero
#output:
#  a seurat object
Read10Xdata <-
  function(h5Path,
           #annoObj = NULL,
           #fragmentsPath = NULL,
           #hintPath = NULL,
           min_cell = 0.01) {
    inputdata.10x <- Read10X_h5(h5Path)
    rna_counts <- inputdata.10x$`Gene Expression`
    atac_counts <- inputdata.10x$Peaks
    grange.counts <-
      StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
    grange.use <-
      seqnames(grange.counts) %in% standardChromosomes(grange.counts)
    atac_counts <- atac_counts[as.vector(grange.use), ]
    chrom_assay <- CreateChromatinAssay(
      counts = atac_counts,
      sep = c(":", "-"),
      min.cells = ncol(atac_counts) * min_cell,
      #fragments = fragmentsPath,
      #annotation = anno
      #min.feature = 300,
    )
    tmp_obj <- CreateSeuratObject(counts = chrom_assay,
                                  assay = "ATAC")
    
    exp_assay <-
      CreateAssayObject(counts = rna_counts,
                        min.cells = ncol(rna_counts) * min_cell)
    tmp_obj[["RNA"]] <- exp_assay
    DefaultAssay(tmp_obj) <- "RNA"
    tmp_obj[["percent.mt"]] <-
      PercentageFeatureSet(tmp_obj, pattern = "^MT-")
    return (tmp_obj)
  }


# Read data with matrix format
##' Read matched scRNA+scATAC data from matrix format
#input:
#  rna_matrix: an expression matrix with gene * cell
#  atac_matrix: an accessibility matrix with peak * cell
#  min_cell: the peak/gene will be removed if the value in the gene/peak with more than min_cell cell is equal to zero
#output:
#  a seurat object
readmatrix <- function(rna_matrix, atac_matrix, min_cell = 0.1) {
  rna_matrix <-
    rna_matrix[, intersect(colnames(atac_matrix), colnames(rna_matrix))]
  atac_matrix <-
    atac_matrix[, intersect(colnames(atac_matrix), colnames(rna_matrix))]
  grange.counts <-
    StringToGRanges(rownames(atac_matrix), sep = c(":", "-"))
  grange.use <-
    seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_matrix <- atac_matrix[as.vector(grange.use), ]
  atac_matrix <-
    atac_matrix[lengths(strsplit(gsub(":", "-", rownames(atac_matrix)) , split = "-")) ==
                  3,]
  rna_matrix <- rna_matrix[unique(rownames(rna_matrix)),]
  #cell_type<-rna_cell$V7[-1][grepl("*RNA*",rna_gene$V2[-1])]
  #min_cell=0.01
  chrom_assay <- CreateChromatinAssay(
    counts = atac_matrix,
    sep = c(":", "-"),
    #genome = annota,
    #fragments = fragments,
    min.cells =  ncol(atac_matrix) * min_cell,
    #min.feature = 300,
    #annotation = annotations
  )
  obj <- CreateSeuratObject(counts = chrom_assay,
                            assay = "ATAC")
  exp_assay <-
    CreateAssayObject(counts = rna_matrix,
                      min.cells = ncol(rna_matrix) * min_cell )
  obj[["RNA"]] <- exp_assay
  DefaultAssay(obj) <- "RNA"
  obj[["percent.mt"]] <-
    PercentageFeatureSet(obj, pattern = "^mt-")
  return(obj)
}


## ----------------------------------------------------------------------------------------------------------------
# Filter abnormal cells
##' Filter abnormal cells
#input:
#  obj: a seurat object
#  nmad: a numeric scalar, specifying the minimum number of MADs away from median required for a value to be called an outlier
# output:
#  a seurat object

filterCell <- function(obj, nmad = 3, data_type = 'scRNA_scATAC') {
  if (data_type == "scRNA_scATAC"){
    
    atac <- isOutlier(obj$nCount_ATAC,
                    nmads = nmad,
                    log = F,
                    type = "both")
  rna <- isOutlier(obj$nCount_RNA,
                   nmads = nmad,
                   log = F,
                   type = "both")
  
  mito <- isOutlier(obj$percent.mt,
                    nmads = nmad,
                    log = F,
                    type = "both")
  obj <-
    AddMetaData(obj, atac, col.name = "atac")
  obj <-
    AddMetaData(obj, rna, col.name = "rna")
  obj <-
    AddMetaData(obj, mito, col.name = "mito")
  obj <- subset(x = obj,
                subset = atac == F &
                  rna == F &
                  mito == F)
  }
  
  if (data_type == "multipleRNA"){
      rna<-isOutlier(
        obj$nCount_RNA,
        nmads = nmad,
        log = F,
        type = "both"
      )
      
      obj <-
        AddMetaData(obj, rna, col.name = "rna")
      obj<-subset(
        x = obj,
        subset = rna == F
      )
      
  }
  return(obj)
  
}


## ----------------------------------------------------------------------------------------------------------------
##' Calculate gene active score matrix
# input:
#  peak_count_matrix: a peak_count matrix from scATAC-seq with peak * cell which return from filterCell function
#  organism: species type GRCh38 / GRCm38
# output:
#  a gene * peak matrix, the elements represent the regulatory potential for peak to gene

CalGenePeakScore <-
  function(peak_count_matrix, organism = "GRCh38") {
    pbmc_peak <- peak_count_matrix
    n <- nrow(pbmc_peak)
    dia <- diag(n)
    rownames(dia) <- rownames(pbmc_peak)
    colnames(dia) <- 1:ncol(dia)
    gene_peak <-
      ATACCalculateGenescore(dia,
                             organism = organism,
                             decaydistance = 10000,
                             model = "Enhanced")
    colnames(gene_peak) <- rownames(peak_count_matrix)
    return (gene_peak)
  }

## ----------------------------------------------------------------------------------------------------------------
##' Calculate gene active score matrix
#input:
#  ATAC_gene_peak: a matrix with gene * peak which return from CalGenePeakScore fucntion
#  obj: a seurat object after data preprocessing which return from filterCell function
#  method: the method to integrate scRNA-seq and scATAC-seq velo (velocity) / WNN (weighted nearest neighbor)
#  veloPath: if use velocity method, the veloPath should be provided
#  return.weight: if return.weight = T, return modality integrated weight, else return GAS
#output:
#  GAS matrix with gene * peak, the elements represent the gene activity score in each cell
#  obj: a seurat object with obj[['ATAC_active']] a gene * peak matrix, the elements represent the regulatory potential for peak to gene
calculate_GAS_v1 <-
  function(ATAC_gene_peak,
           obj,
           method = "velo",
           return.weight = F,
           veloPath = NULL) {
    peak_count <- obj@assays$ATAC@counts
    gene_count <- obj@assays$RNA@counts
    peak_count[peak_count > 0] = 1
    WA <- ATAC_gene_peak %*% peak_count
    
    colnames(WA) <- colnames(peak_count)
    rownames(WA) <- rownames(ATAC_gene_peak)
    #print(colnames(WA)[1:3])
    #print(colnames(obj)[1:3])
    
    
    
    WA <- WA[which(rowSums(as.matrix(WA)) > 0),]
    gene_count <-
      gene_count[which(rowSums(as.matrix(gene_count)) > 0),]
    commongene <-
      intersect(x = rownames(WA), y = rownames(gene_count))
    WA <- as.matrix(WA)
    WA <- WA[commongene,]
    atac_active <- CreateAssayObject(counts = WA,
                                     min.cells = 0)
    obj[['ATAC_active']] <- atac_active
    gene_count <- gene_count[commongene,]
    gene_rowsum <- rowSums(gene_count)
    peak_rowsum <- rowSums(WA)
    norm_gene_count <- gene_count / rowSums(gene_count)
    norm_WBinary <- WA / rowSums(WA)
    #norm_gene_count<-NormalizeData(CreateSeuratObject(counts = gene_count))$ RNA@data
    gene_count <- norm_gene_count
    #norm_WBinary<-NormalizeData(CreateSeuratObject(counts = WA))$RNA@data
    peak_count <- norm_WBinary
    #print(str(peak_count))
    if (method == "velo") {
      velo <- read.csv(veloPath, header = TRUE)
      velo <- as.matrix(velo)
      rownames(velo) <- velo[, 1]
      
      temp <-
        cbind(toupper(rownames(gene_count)), rownames(gene_count))
      rownames(temp) <- temp[, 1]
      temp <-
        temp[intersect(toupper(rownames(gene_count)), rownames(velo)),]
      velo <- velo[temp[, 1],]
      rownames(velo) <- temp[, 2]
      
      velo <- velo[,-1]
      colnames(velo) <- gsub("\\.", "-", colnames(velo))
      vv <- matrix(as.numeric(velo), dim(velo)[1], dim(velo)[2])
      rownames(vv) <- rownames(velo)
      colnames(vv) <- colnames(velo)
      velo <- vv
      rm(vv)
      
      rna <- gene_count
      atac <- peak_count
      velo <-
        velo[intersect(rownames(rna), rownames(velo)), intersect(colnames(rna), colnames(velo))]
      rna <-
        rna[intersect(rownames(rna), rownames(velo)), intersect(colnames(rna), colnames(velo))]
      atac <-
        atac[intersect(rownames(rna), rownames(velo)), intersect(colnames(rna), colnames(velo))]
      gene_rowsum <- gene_rowsum[rownames(rna)]
      peak_rowsum <- peak_rowsum[rownames(rna)]
      #str(velo)
      
      genes <- dim(velo)[1]
      cells <- dim(velo)[2]
      # rank matrix
      rank_cell <- velo
      rank_gene <- velo
      rank_cell <- apply(velo, 2, rank)
      rank_gene <- t(apply(velo, 1, rank))
      
      rank_cell[velo > 0] = genes - rank_cell[velo > 0]
      rank_cell[velo < 0] = rank_cell[velo < 0] - 1
      rank_cell[velo == 0] = 0
      rank_gene[velo > 0] = cells - rank_gene[velo > 0]
      rank_gene[velo < 0] = rank_gene[velo < 0] - 1
      rank_gene[velo == 0] = 0
      
      # number of positive/negative for each gene/cell
      cell_posi_num <- colSums(velo > 0)
      cell_nega_num <- colSums(velo < 0)
      gene_posi_num <- rowSums(velo > 0)
      gene_nega_num <- rowSums(velo < 0)
      
      # weights
      weights = ((rank_cell ^ 2 + rank_gene ^ 2) / ((t((t(velo > 0)) * cell_posi_num + (t(velo <
                                                                                            0)) * cell_nega_num
      )) ^ 2 + ((velo > 0) * gene_posi_num + (velo < 0) * gene_nega_num
      ) ^ 2 + (velo == 0))) ^ 0.5
      weights[velo < 0] = weights[velo < 0] * (-1)
      
      # GAS
      GAS = rna * gene_rowsum + ((1 + weights) * atac) * ((1 + weights) *
                                                            peak_rowsum)
    }
    if (method == "wnn") {
      obj <- obj[, colnames(gene_count)]
      DefaultAssay(obj) <- "RNA"
      obj <-
        FindVariableFeatures(obj,
                             selection.method = "vst",
                             nfeatures = 2000)
      obj <- ScaleData(obj, features = VariableFeatures(obj))
      obj <- RunPCA(obj)
      # ATAC analysis
      # We exclude the first dimension as this is typically correlated with sequencing depth
      DefaultAssay(obj) <- "ATAC"
      obj <- RunTFIDF(obj)
      obj <- FindTopFeatures(obj, min.cutoff = 'q0')
      obj <- RunSVD(obj)
      obj <-
        RunUMAP(
          obj,
          reduction = 'lsi',
          dims = 2:50,
          reduction.name = "umap.atac",
          reduction.key = "atacUMAP_"
        )
      obj <-
        FindMultiModalNeighbors(obj,
                                reduction.list = list("pca", "lsi"),
                                dims.list = list(1:50, 2:50))
      GAS <-
        gene_count * gene_rowsum * obj$RNA.weight + peak_count * obj$ATAC.weight *
        peak_rowsum
    }
    if (isTRUE(return.weight)) {
      return(weights)
    } else {
      m <- list()
      m[[1]] <- GAS
      m[[2]] <- obj
      return(m)
    }
  }




# required package -- reticulate
# input:
#  GAS: the spliced and normalized matrix obtained from CLR function
#  result_dir: The address for storing the models and optimization results(Type:str)
#  epoch:(Type:int)
#  lr: learning rate(Type:float)
#  n_hid: Number of hidden dimension(Type:int)
#  n_heads: Number of attention head(Type:int)
#  cuda: 0 use GPU0 else cpu(Type:int)
#  data_type: 'CITE', 'RNA_ATAC', or 'multiple RNA'
#  envPath: The address for environment to use if use.env is TRUE(Type:str)
# output:
#  HGT_result: a list containing requried results of HGT model as follows:
#  parameters: given parameters from user --epoch, lr, n_hid, n_heads, cuda
#  cell_hgt_matrix: cell embedding matrix 
#  feature_hgt_matrix : gene embedding matrix and protein embedding matrix when data_type is 'CITE';
#  attention: attention meassage for features and cells
#  data_type: 'CITE', 'RNA_ATAC', or 'multiple RNA'
#  result_dir: The address for storing the models and optimization results
#  GAS: the spliced and normalized matrix obtained from CLR function

run_HGT <- function(GAS,result_dir,data_type,envPath=NULL,lr=NULL, epoch=NULL, n_hid=NULL, n_heads=NULL,cuda=0){
  if (data_type == 'CITE') {
    if (is.null(lr)){lr = 0.1}
    if (is.null(epoch)){epoch = 50}
    if (is.null(n_hid)){n_hid = 104}
    if (is.null(n_heads)){n_heads = 13}
  }
  if (data_type == 'scRNA_scATAC') {
    if (is.null(lr)){lr = 0.1}
    if (is.null(epoch)){epoch = 100}
    if (is.null(n_hid)){n_hid = 128}
    if (is.null(n_heads)){n_heads = 16}
  }
  if (data_type == 'multipleRNA') {
    if (is.null(lr)){lr = 0.1}
    if (is.null(epoch)){epoch = 100}
    if (is.null(n_hid)){n_hid = 104}
    if (is.null(n_heads)){n_heads = 13}
  }
  print(epoch)
  cat(lr, epoch, n_hid, n_heads, cuda)
  if (!is.null(envPath)){use_condaenv(envPath)}
  list_in <- assign("list_in", list(lr=lr, epoch=epoch, n_hid=n_hid, n_heads=n_heads, result_dir=result_dir, cuda=cuda,
  data_type=data_type, cell_gene=GAS, gene_name=rownames(GAS), cell_name=colnames(GAS)), envir = .GlobalEnv)
    source_python('./arg.py')
    cell_hgt_matrix <- py$cell_matrix
    gene_hgt_matrix <- py$gene_matrix
    attention <- py$df2
    rownames(cell_hgt_matrix) <- list_in$cell_name
    rownames(gene_hgt_matrix) <- list_in$gene_name
    HGT_result <- list()
    HGT_result[['parameters']] <- data.frame(lr,epoch,n_hid,n_heads,cuda)
    HGT_result[['GAS']] <- GAS
    HGT_result[['cell_hgt_matrix']] <- cell_hgt_matrix
    HGT_result[['feature_hgt_matrix']] <- gene_hgt_matrix
    HGT_result[['attention']] <- attention
    HGT_result[['result_dir']] <- result_dir
    HGT_result[['data_type']] <- data_type
    return(HGT_result)
  
}



## ----------------------------------------------------------------------------------------------------------------
# CT active gene modules calculation
# input:
#  GAS: a gene active matrix with gene * cell which return from calculate_GAS function
#  cell_hgt_matrix: cell-embedding matrix which otains from HGT function
#  att: attention matrix with gene-cell * head which obtain from HGT function
#  gene_hgt_matrix: gene_embedding matrix with gene-cell * head which obtain from HGT function
#  cutoff: the threshold of gene module scale, the higher value the fewer number gene in the module. Default as 1.6.
#output:
#  co (variable 1): a biological gene module. a list with name CT-i and active gene list in CT-i

get_gene_module <-
  function(obj, GAS, att, cutoff = 1.6, method = NULL) {
    if (method == 'SFp'){

      `%!in%` <- Negate(`%in%`) # define the negation of %in%
      
      
      n.matching <- 10 # the number of predicting maximum matching
      m<-inp(GAS, HGT_result[['attention']], HGT_result[['cell_hgt_matrix']], HGT_result[['feature_hgt_matrix']], l=1.2)
      df1 <- m[[1]]
      df2 <- m[[2]]
      terminals <- m[[3]]
      cells <- m[[4]]
      graph.out<-m[[5]]
      
      cut.G <- global_matching_graph(df1, df2, n.matching = 10, cells = cells, terminals)
      cat ('The cutted graph contains', length(V(cut.G)), 'nodes and', length(E(cut.G)), 'edges.\n')
      steiner.ig <- set_cover_mst(G = cut.G, terminals = terminals)
      mods <- get_modules(steiner.ig = steiner.ig, terminals = terminals, 
                          out.file = out.file)
      co<-list()
      for (i in unique(mods$terminal)){
        a<-as.numeric(unlist(strsplit(as.character(unique(mods$terminal)[i]), split ='[.]'))[1])-1
        co[[paste0('ct_',a)]]<-unique(mods$steiner_node[mods$terminal==i])
      }
      
    }else{
      graph.out <- Idents(obj)
      nhead <- ncol(att)
      gene_name <- rownames(GAS)[att$gene + 1]
      cell_name <- colnames(GAS)[att$cell + 1]
      
      att$ct <- graph.out[cell_name]
      att$gene_name <- gene_name
      att$cell_name <- cell_name
      mod <- function(x) {
        return(sqrt(sum(c(x ^ 2))))
      }
      nor <- function(x) {
        return((x - min(x)) / (max(x) - min(x)))
      }
      
      att[, 3:nhead] <- nor(att[, 3:nhead])
      attention <-
        aggregate(x = as.list(att[, 3:nhead]),
                  by = list(att$ct, att$gene_name),
                  mean)
      #att[,4:nhead]<-1-att[,4:nhead]
      weight <- apply(att[, 3:nhead], 1, mod)
      df <-
        data.frame(
          'node1' = att$gene_name,
          'node2' = att$cell_name,
          'weight' = weight,
          'ct' = att$ct
        )
      attention <-
        aggregate(x = df$weight, by = list(df$ct, df$node1), mean)
      co <- list()
      for (i in (0:(length(unique(att$ct)) - 1))) {
        t <-
          mean(attention[attention$Group.1 == i, ]$x) + 1.6 * sd(attention[attention$Group.1 ==
                                                                             i, ]$x)
        co[[paste('ct', i, sep = "_")]] <-
          attention[attention$Group.1 == i, ]$Group.2[attention[attention$Group.1 ==
                                                                  i, ]$x > t]
      }
    }
    
    return (co)
    
  }


## ----------------------------------------------------------------------------------------------------------------
# gene module save
#input:
#  co: the active gene module from get_gene_module function
#  lisa_path: the path of active gene module to save
#result
#  write gene module to the lisa_path

write_GM <- function(co, lisa_path) {
  if (length(dir(path = lisa_path, pattern = ".csv")) >
      0) {
    system(paste0("rm ", lisa_path, "*.csv "))
    
  }
  if (length(dir(path = lisa_path, pattern = ".txt")) >
      0) {
    system(paste0("rm ", lisa_path, "*.txt "))
  }
  
  for (j in (1:length(co))) {
    if (length(unique(co[[j]])) < 20 |
        length(unique(co[[j]])) > 20000) {
      next
    } else{
      ct <- unlist(strsplit(names(co[j]), split = "_"))[1]
      
      write.table(
        co[[j]],
        paste0(lisa_path, names(co[j]), ".txt"),
        quote = F,
        sep = "\t",
        row.names = F,
        col.names = F
      )
    }
  }
}



## ----------------------------------------------------------------------------------------------------------------
# Filter gene with no accessible peak in promoter
# input:
#  obj: a seurat object which return from filterCell function
#  gene_peak: a matrix with gene * peak from scATAC-seq which return from filterCell function
#  GAS: the GAS matrix with gene * cell which return calculate_GAS function
#  species: human / mouse
#output:
#  a matrix with gene * peak. The gene with no accessible peak will be removed

AccPromoter <- function(obj, gene_peak, GAS, species = "hg38") {
  peak_cell <- obj@assays$ATAC@counts
  if (species == "hg38") {
    gene.ranges <- genes(EnsDb.Hsapiens.v86)
  } else{
    gene.ranges <- genes(EnsDb.Mmusculus.v79)
  }
  
  gene.use <-
    seqnames(gene.ranges) %in% standardChromosomes(gene.ranges)[standardChromosomes(gene.ranges) !=
                                                                  "MT"]
  gene.ranges <- gene.ranges[as.vector(gene.use)]
  gene.ranges <-
    gene.ranges[gene.ranges$gene_name %in% rownames(GAS)]
  genebodyandpromoter.coords <-
    Extend(x = gene.ranges,
           upstream = 2000,
           downstream = 0)
  #str(genebodyandpromoter.coords)
  x <- as.data.frame(genebodyandpromoter.coords@ranges)
  peaks <-
    GRanges(
      seqnames = paste("chr", genebodyandpromoter.coords@seqnames, sep = ""),
      ranges = IRanges(start = , x$start,
                       width = x$width)
    )
  
  peak_name <-
    colnames(gene_peak)[lengths(strsplit(gsub(":", "-", colnames(gene_peak)) , split = "-")) ==
                          3]
  peak_name <-
    do.call(what = rbind, strsplit(gsub(":", "-", peak_name) , split = "-"))
  peak_name <- as.data.frame(peak_name)
  names(peak_name) <- c("chromosome", 'start', 'end')
  peak_name <- GenomicRanges::makeGRangesFromDataFrame(peak_name)
  #str(peaks)
  over <- findOverlaps(peak_name, peaks)
  #str(over)
  promoter_gene <-
    genebodyandpromoter.coords$gene_name[unique(over@to)]
  str(promoter_gene)
  gene_peak <- gene_peak[promoter_gene,]
  
  return(gene_peak)
}

## ----------------------------------------------------------------------------------------------------------------
##' infer ct active regulons
# input:
#  GAS: the GAS matrix with gene * cell which return calculate_GAS function
#  co: a list of bio network which reture from gene_ function
#  gene_peak_pro: the matrix with gene * peak which return AccPromoter function
#  species: human / mouse (human = "hg38", mouse = "mm10" )
#  humanPath: if species == human, the TF binding RData absolute path of hg38 should be provided
#  mousePath: if species == mouse, the TF binding RData absolute path of mm10 should be provided
#output:
#  BA_score a TF binding affinity matrix with TF * peak, the elements in the matrix is the binding power of TF to peak
#  ct_regulon: candidate cell type active regulon
#  TFinGAS: TRUE / FALSE, if number of the intersection of candidate TF from LISA and gene in GAS > 50 TFinGAS will be true, else it will be false

Calregulon <-
  function(GAS,
           co,
           gene_peak_pro,
           species = "hg38",
           jaspar_path = "/scratch/deepmaps/jaspar",
           lisa_path = "/home/wan268/hgt/RNA_ATAC/lymph_14k/") {
    if (species == "hg38") {
      tfbs_df <- qs::qread(paste0(jaspar_path, "jaspar_hg38_500.qsave"))
    }
    else {
      tfbs_df <- qs::qread(paste0(jaspar_path, "jaspar_mm10_500.qsave"))
    }
    BA_score <-
      matrix(0, ncol(gene_peak_pro), length(unique(tfbs_df$V4)))
    colnames(BA_score) <- unique(tfbs_df$V4)
    rownames(BA_score) <- colnames(gene_peak_pro)
    gene_TF <-
      matrix(0, nrow(gene_peak_pro), length(unique(tfbs_df$V4)))
    colnames(gene_TF) <- unique(tfbs_df$V4)
    rownames(gene_TF) <- rownames(gene_peak_pro)
    
    peak <- tfbs_df[, 1:3]
    colnames(peak) <- c("chromosome", 'start', 'end')
    peak <- GenomicRanges::makeGRangesFromDataFrame(peak)
    
    ct_subregulon <- list()
    ct_regulon <- list()
    coexp_tf <- list()
    No <- 1
    for (i in (1:length(co))) {
      if (length(co[[i]]) > 0) {
        co[[i]] <- intersect(co[[i]], rownames(gene_peak_pro))
        a <- which(gene_peak_pro[co[[i]],] > 0, arr.ind = T)
        op <- colnames(gene_peak_pro)[unname(a[, 'col'])]
        peak_name <-
          op[lengths(strsplit(gsub(":", "-", op) , split = "-")) == 3]
        peak_name <-
          do.call(what = rbind, strsplit(gsub(":", "-", peak_name) , split = "-"))
        peak_name <- as.data.frame(peak_name)
        names(peak_name) <- c("chromosome", 'start', 'end')
        peak_name <-
          GenomicRanges::makeGRangesFromDataFrame(peak_name)
        over <- findOverlaps(peak_name, peak)
        #print(i)
        p <- op[over@from]
        pp <- tfbs_df$V5[over@to] / 100
        df <- data.frame(p, pp, tfbs_df$V4[over@to])
        hh <- df[!duplicated(df[,-2]),]
        for (k1 in (1:nrow(hh))) {
          BA_score[hh[k1,]$p, hh[k1,]$tfbs_df.V4.over.to.] <- hh[k1,]$pp
        }
        
        gene_TF <- gene_peak_pro %*% BA_score
        TF <- unique(tfbs_df$V4[over@to])
        
        if (length(co[[i]]) < 20000 & length(co[[i]]) > 20) {
          tf <-
            read.csv(paste(lisa_path,
                           names(co[i]),
                           ".txt.csv",
                           sep = ""))
          tf_pval_0.05 <- unique(tf[, 3][tf$summary_p_value < 0.05])
          TF <- intersect(unique(TF), tf_pval_0.05)
        }
        #print(length(TF))
        
        #gene_TF[co[[i]],pp]<-gene_peak_pro[co[[i]],p] %*% BA_score[p,TF]
        coexp_tf[[names(co[i])]] <- TF
        #print(length(TF))
        if (length(TF) > 0) {
          for (k in 1:length(TF)) {
            if (length(intersect(TF[k], rownames(GAS))) > 50) {
              TFinGAS <- T
              if (TF[k] %in% rownames(GAS)) {
                #print(TF[k])
                a <- unlist(strsplit(names(co[i]), "_"))
                a <- paste0(a[1], a[2])
                h <- paste(TF[k], a, sep = "_")
                ct_subregulon[[h]] <-
                  co[[i]][gene_TF[co[[i]], TF[k]] > 0]
              }
            } else{
              a <- unlist(strsplit(names(co[i]), "_"))
              TFinGAS <- F
              a <- paste0(a[1], a[2])
              h <- paste(TF[k], a, sep = "_")
              ct_subregulon[[h]] <-
                co[[i]][gene_TF[co[[i]], TF[k]] > 0]
            }
          }
        }
      }
    }
    
    m <- list()
    m[[1]] <- BA_score
    ct_subregulon <- ct_subregulon[lengths(ct_subregulon) > 10]
    m[[2]] <- ct_subregulon
    m[[3]] <- TFinGAS
    
    return(m)
  }


## ----------------------------------------------------------------------------------------------------------------
# combine same TF
#input:
#  gene_peak_pro: a matrix with gene * peak. The gene with no accessible peak will be removed which return from AccPromoter function
#  BA_score: a TF binding affinity matrix with TF * peak, the elements in the matrix is the binding power of TF to peak which returen from Calregulon function
#output:
#  peak_TF: a matrix with peak * TF without repeat TF

uni <- function(gene_peak_pro, BA_score) {
  gene_TF <- gene_peak_pro %*% BA_score
  rownames(gene_TF) <- rownames(gene_peak_pro)
  mat <-
    matrix(0, nrow = length(rownames(gene_TF)), ncol = length(unique(colnames(gene_TF))))#peak_TF
  mat1 <-
    matrix(0, nrow = length(rownames(BA_score)), ncol = length(unique(colnames(BA_score))))#gene_TF
  rownames(mat1) <- rownames(BA_score)
  colnames(mat1) <- unique(colnames(BA_score))
  #mat1: peak_TF score
  for (x in unique(colnames(BA_score))) {
    if (is.null(nrow(gene_TF[, colnames(BA_score) == x]))) {
      #mat<-rbind(mat, unlist(gene_peak_matrix[rownames(gene_peak_matrix)==x,]))
      mat1[, x] <-
        unname(unlist(BA_score[, colnames(BA_score) == x]))
    } else{
      #mat<-rbind(mat,unlist(colSums(gene_peak_matrix[rownames(gene_peak_matrix)==x,])))
      mat1[, x] <-
        unname(unlist(rowSums(BA_score[, colnames(BA_score) == x])))
    }
  }
  
  
  
  rownames(mat) <- rownames(gene_TF)
  colnames(mat) <- unique(colnames(gene_TF))
  
  for (x in unique(colnames(gene_TF))) {
    if (is.null(nrow(gene_TF[, colnames(gene_TF) == x]))) {
      #print("111")
      #mat<-rbind(mat, unlist(gene_peak_matrix[rownames(gene_peak_matrix)==x,]))
      mat[, x] <- unname(unlist(gene_TF[, colnames(gene_TF) == x]))
    } else{
      #mat<-rbind(mat,unlist(colSums(gene_peak_matrix[rownames(gene_peak_matrix)==x,])))
      mat[, x] <-
        unname(unlist(rowSums(gene_TF[, colnames(gene_TF) == x])))
    }
  }
  #m = list()
  #m[[1]] <- mat
  #m[[2]] <- mat1
  return (mat1)
}



## ----------------------------------------------------------------------------------------------------------------
# Regulatory Intensive (RI) score in cell level
# input:
# 1 - obj: a seurat object which return from filterCell function
# 2 - ct_regulon: cell type active regulon
# 3 - GAS: the GAS matrix with gene * cell which return calculate_GAS function
# 4 - gene_peak_pro: the matrix with gene * peak which return from AccPromoter function
# 5 - peak_TF: the matrix with peak * TF which return from uni function
# 6 - graph.out: a factor variable. The predict cell cluster which return from get_gene_module function
#output:
#  - RI_C: a regulatory intensive matrix with TF-gene pair * cell, the element means the intensity of TF to gene in each cell

RI_cell <-
  function(obj,
           ct_regulon,
           GAS,
           gene_peak_pro,
           peak_TF,
           graph.out) {
    v <- vector()
    for (i in (1:length(ct_regulon))) {
      v <-
        append(v, paste(unlist(strsplit(
          names(ct_regulon[i]), "_"
        ))[1], ct_regulon[[i]], sep = "_"))
    }
    peak_cell <- obj@assays$ATAC@counts
    peak_cell[peak_cell > 0] = 1
    peak_cell <- (peak_cell[, colnames(GAS)])
    
    #graph.out<-Idents(obj)
    #TG_cell <- foreach(i=1:length(ct_regulon), .packages='Matrix',.combine='c') %dopar%{
    TG_cell <- matrix(0, length(unique(v)), length(graph.out))
    rownames(TG_cell) <- unique(v)
    t <-
      unlist(strsplit(unique(v), "_"))[seq(1, length(unlist(strsplit(unique(
        v
      ), "_"))), 2)]
    g <-
      unlist(strsplit(unique(v), "_"))[seq(2, length(unlist(strsplit(unique(
        v
      ), "_"))), 2)]
    t <- unique(t)
    g <- unique(g)
    t1 <-
      unlist(strsplit(unique(v), "_"))[seq(1, length(unlist(strsplit(unique(
        v
      ), "_"))), 2)]
    g1 <-
      unlist(strsplit(unique(v), "_"))[seq(2, length(unlist(strsplit(unique(
        v
      ), "_"))), 2)]
    #TfGene_cell <- foreach(i=1:length(graph.out), .packages='Matrix',.combine='c') %dopar%
    #  {
    gene_peak <- gene_peak_pro
    for (j in (1:length(graph.out))) {
      hhh <-
        (peak_cell[, j] * gene_peak[g,]) %*% (peak_cell[, j] * peak_TF[, t])
      #hhh<-gene_peak[g,] %*%peak_TF[,t]
      bb = data.frame('gene' = g1, 'tf' = t1)
      bb$tf = as.factor(bb$tf)
      bb$gene = as.factor(bb$gene)
      levels(bb$gene) = 1:length(levels(bb$gene))
      levels(bb$tf) = 1:length(levels(bb$tf))
      bb <- as.matrix(bb)
      bb <- apply(bb, 2, as.numeric)
      TG_cell[, j] <- hhh[bb]
    }
    return(TG_cell)
  }

# calculate regulon active score in cell level / cell type level
# input:
#  RI_C: a regulatory intensive matrix with TF-gene pair * cell, the element means the intensity of TF to gene in each cell
#  ct_regulon: cell type active regulon
#  graph.out: a factor variable. The cell cluster which return from get_gene_module function
#  TFinGAS: TRUE / FALSE, if number of the intersection of candidate TF from LISA and gene in GAS > 50 TFinGAS will be true, else it will be false which return from Calregulon function
#output:
#  RAS_CT: regulon active score in cell type level
#  RI_CT: regulatory intensive score in cell type level
#  ct-regulon: cell type active regulon
#  RAS_C1: a matrix regulon-CT * cell, regulatory active score in cell level

calRAS <- function(RI_C, ct_regulon, graph.out, TFinGAS = T) {
  a <- unlist(strsplit(names(ct_regulon), "_"))
  CT <- unique(a[seq(2, length(a), 2)])
  TF <- unique(a[seq(1, length(a), 2)])
  RI_CT <- matrix(0, nrow(RI_C), length(CT))
  rownames(RI_CT) <- rownames(RI_C)
  colnames(RI_CT) <- CT
  for (i in CT) {
    a <- graph.out == as.numeric(substring(i, 3, nchar(i)))
    RI_CT[, i] <- rowSums((RI_C[, a])) / length(graph.out[a])
  }
  g <-
    unlist(strsplit(rownames(RI_C), "_"))[seq(2, length(unlist(strsplit(
      rownames(RI_C), "_"
    ))), 2)]
  RAS_E <- RI_C * as.matrix(GAS[g, ])
  colnames(RAS_E) <- colnames(RI_C)
  #RAS_C TF(regulon)*cell
  RAS_C <-
    as(matrix(0, nrow = length(TF), ncol = length(graph.out)), "sparseMatrix")
  
  rownames(RAS_C) <- TF
  colnames(RAS_C) <- names(graph.out)
  for (i in (1:length(ct_regulon))) {
    tf <- unlist(strsplit(names(ct_regulon[i]), "_"))[1]
    #ct<-unlist(strsplit(names(ct_regulon[i]),"_"))[2]
    RAS_C[tf, ] <-
      colMeans(RAS_E[paste(tf, ct_regulon[[i]], sep = "_"), ])
  }
  
  #RAS TF*CT
  RAS <-
    as(matrix(0, nrow = length(TF), ncol = length(CT)), "sparseMatrix")
  rownames(RAS) <- TF
  colnames(RAS) <- CT
  for (i in (1:length(ct_regulon))) {
    tf <- unlist(strsplit(names(ct_regulon[i]), "_"))[1]
    ct <- unlist(strsplit(names(ct_regulon[i]), "_"))[2]
    if (TFinGAS == T)
    {
      if (sum(GAS[tf, graph.out == substring(ct, 3, nchar(ct))]) == 0) {
        RAS[tf, ct] <- 0
      }
      else{
        RAS[tf, ct] <-
          mean(RAS_E[paste(tf, ct_regulon[[i]], sep = "_"), graph.out == substring(ct, 3, nchar(ct))])
      }
      
    } else{
      RAS[tf, ct] <-
        mean(RAS_E[paste(tf, ct_regulon[[i]], sep = "_"), graph.out == substring(ct, 3, nchar(ct))])
    }
    
    #print(i)
  }
  ct_re <- list()
  for (i in (1:length(ct_regulon))) {
    tf <- unlist(strsplit(names(ct_regulon[i]), "_"))[1]
    ct <- unlist(strsplit(names(ct_regulon[i]), "_"))[2]
    re <-
      names(RI_CT[paste(tf, ct_regulon[[i]], sep = "_"), ct][RI_CT[paste(tf, ct_regulon[[i]], sep =
                                                                           "_"), ct] > 0])
    if (length(re) > 0) {
      ct_re[[names(ct_regulon[i])]] <-
        unlist(strsplit(re, "_"))[seq(2, length(unlist(strsplit(re, "_"))), 2)]
    }
  }
  if (length(ct_re[lengths(ct_re) > 10]) < 10)
  {
    ct_re <- ct_re[lengths(ct_re) > 3]
  } else{
    ct_re <- ct_re[lengths(ct_re) > 10]
  }
  
  m <- list()
  m[[1]] <- RAS
  m[[2]] <- RI_CT[rowSums(RI_CT) > 0, ]
  m[[3]] <- ct_re
  m[[4]] <- RAS_C
  return (m)
}

## ----------------------------------------------------------------------------------------------------------------
##calculate master TF
#find master TF and master gene, build a cell type gene regulatory network
#input:
# ct_regulon: cell type active regulon
# RI_CT: regulatory intensive score in cell type level which return from calRAS funciton 
#output:
# TF_cen: TF centrality score in each CT
# gene_cen: gene centrality score in each CT
# network: adjcent matrix in a CT

masterFac <- function(ct_regulon, RI_CT) {
  TF_CT <- unlist(strsplit(names(ct_regulon), split = "_"))
  TF <- unique(TF_CT[seq(1, length(TF_CT), 2)])
  CT <- unique(TF_CT[seq(2, length(TF_CT), 2)])
  TF_R <- TF_CT[seq(1, length(TF_CT), 2)]
  CR_R <- TF_CT[seq(2, length(TF_CT), 2)]
  TF_cen <- list()
  network <- list()
  gene_cen <- list()
  for (i in (CT)) {
    TF_CT <- unlist(strsplit(names(ct_regulon[CR_R == i]), split = "_"))
    TF_CT <- unique(TF_CT[seq(1, length(TF_CT), 2)])
    gene_CT <- unique((union(unlist(ct_regulon[CR_R == i]), TF_CT)))
    adj <- matrix(0, nrow = length(gene_CT), ncol = length(gene_CT))
    
    rownames(adj) <- gene_CT
    colnames(adj) <- gene_CT
    for (k in (1:length(ct_regulon[CR_R == i]))) {
      j <- ct_regulon[CR_R == i][k]
      adj1 <-
        matrix(0, nrow = length(gene_CT), ncol = length(gene_CT))
      rownames(adj1) <- gene_CT
      colnames(adj1) <- gene_CT
      RI_CT[paste(unlist(strsplit(names(j), split = "_"))[1], unlist(unname(j)), sep =
                    "_"), i]
      adj1[unlist(strsplit(names(j), split = "_"))[1], unlist(unname(j))] <-
        RI_CT[paste(unlist(strsplit(names(j), split = "_"))[1], unlist(unname(j)), sep =
                      "_"), i]
      adj1[unlist(unname(j)), unlist(strsplit(names(j), split = "_"))[1]] <-
        RI_CT[paste(unlist(strsplit(names(j), split = "_"))[1], unlist(unname(j)), sep =
                      "_"), i]
      adj = adj1 + adj
    }
    network[[i]] <- adj
    g <- graph.adjacency(adj, mode = "undirected", weighted = T)
    cen <- igraph::evcent(g)$vector
    TF_cen[[i]] <- cen[TF_CT]
    gene_cen[[i]] <- cen[setdiff(names(cen), TF_CT)]
    
  }
  m <- list()
  m[[1]] <- TF_cen
  m[[2]] <- gene_cen
  m[[3]] <- network
  return (m)
}
##calculate RAS(2) same topology in a row
# calculate regulon active score in cell level / cell type level
# input:
#  ct_regulon: a matrix with gene * peak from scATAC-seq which return from get_gene_module function
#  graph.out: a factor variable. The predict cell cluster which return from get_gene_module function 
#output:
#  RAS_C2: a matrix with regulon-CT * cell. Regulon active score in cell type level

CalRAS2 <- function(ct_regulon, graph.out) {
  RAS_C2 <-
    as(matrix(0, nrow = length(ct_regulon), ncol = length(graph.out)), "sparseMatrix")
  rownames(RAS_C2) <- names(ct_regulon)
  colnames(RAS_C2) <- names(graph.out)
  for (i in (1:length(ct_regulon))) {
    if (i %% 100 == 0) {
      print(i)
    }
    tf <- unlist(strsplit(names(ct_regulon[i]), "_"))[1]
    ct <- unlist(strsplit(names(ct_regulon[i]), "_"))[2]
    g <- ct_regulon[[i]]
    RAS_C2[names(ct_regulon[i]), ] <-
      colMeans(RI_C[paste(tf, g, sep = "_"), ] * as.matrix(GAS[g, ]))
  }
  return(RAS_C2)
}



## ----------------------------------------------------------------------------------------------------------------
##calculate DR
#infer cell type specific regulon
#input:
# RAS_C2: a matrix with regulon-CT * cell. Regulon active score in cell type level
# graph.out: a factor variable. The predict cell cluster which return from get_gene_module function 
# only.pos: only print positive regulons
# lfcThre: logFoldchange threshold, default = 0.25
# pvalThres: pvalue threshold, default = 0.05
#output:
# DR: differnent regulons. 

mean.fxn <- function(x) {
  return(log(x = rowMeans(x = x) + 1, base = 2))
}


calDR <-
  function(RAS_C1,
           graph.out,
           only.pos = T,
           lfcThres = 0.25,
           pvalThres = 0.05,
           ident.1 = 1,
           ident.2 = c(2, 3, 4)) {
    ras_obj <- CreateSeuratObject(counts = RAS_C1)
    Idents(ras_obj) <- graph.out
    DR <-
      FindAllMarkers(
        ras_obj,
        only.pos = only.pos,
        min.pct = 0,
        logfc.threshold = lfcThres,
        min.cells.feature = 0,
        min.cells.group = 0,
        return.thresh = 1,
        mean.fxn = mean.fxn
      )
    
    DR1 <- DR %>%
      dplyr::filter(p_val <= pvalThres) %>%
      tidyr::separate(gene, c("gene", 'from'), "-") %>%
      dplyr::filter(from == paste0("ct",cluster)) 
    
    #m <-
    #  unlist(strsplit(DR$gene, "-"))[seq(2, length(unlist(strsplit(DR$gene, "-"))), 2)]
    #m <-
    #  unlist(strsplit(m, "ct"))[seq(2, length(unlist(strsplit(m, "ct"))), 2)]
    #DR <- DR[as.numeric(m) == DR$cluster,]
    #DR$gene <- str_remove(DR$gene, "-.*")
    return (DR1)
  }


calDR_v2 <-
  function(RAS_C1,
           graph.out,
           FindALLMarkers = T,
           lfcThres = 0,
           ident.1 = 1,
           ident.2 = c(2, 3, 4)) {
    pbmc <- CreateSeuratObject(counts = RAS_C1)
    Idents(pbmc) <- graph.out
    if (FindALLMarkers == T) {
      DR <-
        FindAllMarkers(
          pbmc,
          only.pos = T,
          logfc.threshold = lfcThres,
          min.pct = 0,
          mean.fxn = mean.fxn
        )
      DR <- DR[DR$p_val < 0.05,]
      m <-
        unlist(strsplit(DR$gene, "-"))[seq(2, length(unlist(strsplit(DR$gene, "-"))), 2)]
      m <-
        unlist(strsplit(m, "ct"))[seq(2, length(unlist(strsplit(m, "ct"))), 2)]
      DR <- DR[as.numeric(m) == DR$cluster,]
    } else{
      DR <-
        FindMarkers(pbmc,
                    ident.1 = ident.1,
                    ident.2 = ident.2,
                    min.pct = 0.25)
    }
    return (DR)
  }


# get the edge list connecting terminal nodes and steiner nodes
# steiner nodes refer to the non-terminals in the steiner forest

get_modules <- function(steiner.ig, terminals, out.file = 'examples/SFP_edges.csv') {
  
  # steiner.ig : the graph
  # terminals : the terminal list
  # steiner.es : the edge list in the steiner forest
  
  #  steiner.vs <- V(steiner.ig)[setdiff(V(steiner.ig)[unique(as.vector(ends(graph = steiner.ig, es = steiner.es)))], 
  #                             V(steiner.ig)[Reduce(union, terminals)])]
  
  steiner.vs <- setdiff(as_ids(V(steiner.ig)), Reduce(union, terminals))
  # build a data frame for the steiner edges
  
  # retain the edges connecting the terminal nodes and steiner nodes
  between.ll <- lapply(terminals, function(ter) {
    return(E(steiner.ig)[ter %--% steiner.vs]) # get the steiner edges between
    # the terminal nodes and the steiner nodes
  })
  
  # arrange the edge list into the format of four columns, i.e. 
  # - terminal - terminal id, e.g., 1, 2, or 3
  # - terminal_node - the node in a terminal
  # - steiner_node - the non-terminal node connected to a terminal node, i.e., steiner node
  # - weight - the weight of the edge connecting a terminal node and steiner node, e.g., 1 and 2
  es.info <- do.call('rbind', lapply(seq_along(between.ll), function(i) {
    x <- between.ll[[i]] # an edge list
    x.df <- data.frame(ends(graph = steiner.ig, es = x))
    # convert edge list into a data frame
    
    # rearrange the order of edge ends
    suppressMessages(library(data.table))
    xx.df <- rbindlist(apply(x.df, 1, function(y) {
      ifelse (y[1] %in% terminals[[i]], return(y), return(list(y[2], y[1])))
    }), fill = F)
    
    colnames(xx.df) <- c('terminal_node', 'steiner_node') # name the columns
    xx.df <- cbind(terminal = rep(as.numeric(names(between.ll)[i]), nrow(xx.df)), xx.df, 
                   weight = x$weight)
    # add the terminal ids and edge weights as a column, respectively
    
    return(xx.df)
  }))
  
  #write.csv(es.info, file = out.file, quote = F) # save the result to a csv file
  
  return(es.info)
}


# this algorithm is partitioned into three steps: 
# - 1 - get all the neighbors of a terminal
# - 2 - find a minimum spanning tree on the subgraph composed of the terminal nodes 
#       and their neighbors as well as the induced edges
# - 3 - remove the non-terminal leaves
# 
# return an edge sequence of the steiner forest
set_cover_mst <- function(G, terminals) {
  
  # G : the weighted/unweighted undirected graph saved in an igraph object
  # terminals : the terminal set saved in a nested list
  
  suppressMessages(library(parallel))
  mst.ll <- mclapply(terminals, function(ter) {
    ter.nbr <- Reduce(union, lapply(ter, function(v) {
      return(as_ids(neighbors(graph = G, v = v, mode = 'all')))}))
    # get neighbors of a terminal
    
    sub.G <- induced_subgraph(graph = G, v = union(ter, ter.nbr))
    # build a subgraph
    
    ter.mst <- mst(graph = sub.G, weights = E(sub.G)$weight, algorithm = 'prim')
    # minimum spanning tree for a subgraph
    
    mst.deg <- degree(graph = ter.mst, v = V(ter.mst), mode = 'all')
    # get the dgree of nodes in the minimum spanning tree
    
    leaves <- names(mst.deg)[mst.deg == 1] # get all the leaves
    return(induced_subgraph(graph = ter.mst, 
                            v = setdiff(as_ids(V(ter.mst)), 
                                        leaves[leaves %!in% ter])))
    # remove the nodes not belonging to the terminal and get the edges
  }, mc.cores = detectCores() - 1)
  
  Ugraph <- Reduce(`%u%`, mst.ll) # get the union of all the 
  # minimum spanning trees
  
  edge.ids <- get.edge.ids(graph = G, vp = as.vector(t(ends(graph = Ugraph, 
                                                            es = E(Ugraph))))) # get edge ids
  E(Ugraph)$weight <- E(G)[edge.ids]$weight
  # assign edge weights according to that on the original graph
  
  steiner.ig <- mst(graph = Ugraph, weights = E(Ugraph)$weight, algorithm = 'prim')
  # minimum spanning tree for a subgraph
  
  Uter <- Reduce(union, terminals) # the union of all terminals
  
  while (1) {
    Udeg <- degree(graph = steiner.ig, v = V(steiner.ig), mode = 'all')
    # get the dgree of nodes in the minimum spanning tree
    
    leaves <- names(Udeg)[Udeg == 1] # get all the leaves
    steiner.leaves <- leaves[leaves %!in% Uter] # the non-terminal nodes with degree one
    
    #        cat (steiner.leaves, '\n\n')
    
    if (length(steiner.leaves) < 1) {
      break
    }
    
    #  steiner.ig <- induced_subgraph(graph = Umst, 
    #                          v = setdiff(as_ids(V(Umst)), 
    #                                      leaves[leaves %!in% Uter]))
    
    steiner.ig <- delete_vertices(graph = steiner.ig, 
                                  v = steiner.leaves) # delete nodes
    #        cat (length(V(steiner.ig)), '\n\n')
  }
  
  cat(length(E(steiner.ig)), 'edges are identified for the Steiner Forest Problem.\n')
  
  return(steiner.ig)
}


# building a relatively sparse network using several times of maximum matching
global_matching_graph <- function(df1, df2, n.matching = 10, cells, terminals) {
  genes <- unique(c(df1$node1, df2$node1, df2$node2))
  cat ('There are in total', length(genes), 'genes.\n')
  
  # which terminal does each cell belong to?
  cell.ter <- Reduce(c, lapply(names(terminals), function(t) {
    t.v <- rep(t, length(terminals[[t]]))
    names(t.v) <- terminals[[t]]
    t.v
  }))
  
  tmp.cells <- cells # uncovered cells
  tmp.df <- df1[df1$node1 %in% genes & df1$node2 %in% tmp.cells, ]
  # the data frame to build the bipartite graph
  
  selected.edges.df <- data.frame(node1 = character(), node2 = character())
  # the data frame to hold selected edges
  
  nn <- 0 # the number of iterations
  
  cat('The size of the data frame to build graph is', dim(tmp.df), '.\n')
  
  # select edges iteractively
  while (1) {
    cat('There are in total', length(tmp.cells), 'uncovered cells.\n')
    if (length(tmp.cells) < 1 | nn >= n.matching) {
      break
    }
    
    nn <- nn + 1
    biG <- graph_from_data_frame(tmp.df, directed = F)
    
    if (length(V(biG)) == 0 | length(E(biG)) == 0) {
      break
    }
    
    cat ('The bipartite graph contains', length(V(biG)), 'nodes and', 
         length(E(biG)), 'edges.\n')
    V(biG)[unique(tmp.df$node1)]$type <- rep(T, length(unique(tmp.df$node1)))
    V(biG)[unique(tmp.df$node2)]$type <- rep(F, length(unique(tmp.df$node2)))
    E(biG)$weight <- tmp.df$weight
    
    max.matching <- max_bipartite_match(
      biG,
      types = V(biG)$type,
      weights = E(biG)$weight
    )
    
    matching.es <- max.matching$matching
    cat ('There are in total', length(matching.es), 'selected edges.\n')
    tmp.cells <- setdiff(names(matching.es[is.na(matching.es)]), genes)
    matching.es <- matching.es[!is.na(matching.es)]
    matching.df <- data.frame(node1 = names(matching.es), 
                              node2 = matching.es)
    
    matching.df <- rbindlist(mclapply(1:nrow(matching.df), function(i) {
      n1 <- matching.df[i, 1]
      n2 <- matching.df[i, 2]
      
      ifelse (n1 %!in% genes, return(list(node1 = n2, node2 = n1)), 
              return(list(node1 = n1, node2 = n2)))
    }, mc.cores = detectCores()), fill = T) %>% 
      distinct()
    
    rownames(matching.df) <- NULL
    selected.edges.df <- dplyr::union(selected.edges.df, matching.df) # add new edges
    cat ('After merging, there are in total', nrow(selected.edges.df), 'selectd edges and', 
         length(unique(selected.edges.df$node1)), 'genes.\n')
    tmp.df <- tmp.df[tmp.df$node2 %in% tmp.cells, ] # update the uncovered cells
    # tmp.df <- tmp.df[tmp.df$node2 %in% tmp.cells & tmp.df$node1 %!in% unique(selected.edges.df$node1), ] # update the uncovered cells
  }
  
  # global.df <- rbind(df1, df2)
  G <- graph_from_data_frame(df1, directed = F) # build graph
  E(G)$weight <- df1$weight # assign weights
  
  # add the remaining uncovered cells
  if (length(tmp.cells) > 0) {
    selected.edges.df <- dplyr::union(selected.edges.df, 
                                      rbindlist(mclapply(tmp.cells, function(u) {
                                        u.ter <- cell.ter[[u]]
                                        u.edges <- incident_edges(G, u, 'all')
                                        return(list(node1 = ends(G, u.edges[[1]][which.max(u.edges[[1]]$weight)])[1, 1], 
                                                    node2 = ends(G, u.edges[[1]][which.max(u.edges[[1]]$weight)])[1, 2]))
                                      }, mc.cores = detectCores()), fill = T))
  } # all edges linked to the uncovered cells
  
  rownames(selected.edges.df) <- paste0(selected.edges.df$node1, '_', 
                                        selected.edges.df$node2)
  rownames(df1) <- paste0(df1$node1, '_', df1$node2)
  final.df <- df1[rownames(selected.edges.df), ]
  final.df$weight <- (max(final.df$weight) - final.df$weight) / 
    (max(final.df$weight) - min(final.df$weight))
  df2$weight <- (max(df2$weight) - df2$weight) / 
    (max(df2$weight) - min(df2$weight))
  merged.df <- rbind(final.df, df2)
  
  return(graph_from_data_frame(merged.df, directed = F))
} 


inp<-function(GAS = GAS, att , cell_hgt_matrixPath, genePath, wdf1=1, thdf2=0, l=1.2){
  GAS<-read.table(GASpath)
  kidney<-CreateSeuratObject(GAS)
  cell_hgt_matrix <- read.table(cell_hgt_matrixPath)
  cell_hgt_matrix <- as.matrix(cell_hgt_matrix)
  rownames(cell_hgt_matrix) <-colnames(GAS)
  kidney<-kidney[,colnames(GAS)]
  HGT_embedding <-
    CreateDimReducObject(
      embeddings = cell_hgt_matrix,
      key = "HGT_",
      assay = "RNA"
    )
  kidney@reductions[['HGT']] <- HGT_embedding
  kidney <- FindVariableFeatures(kidney, selection.method = "vst", nfeatures = 2000)
  kidney<-ScaleData(kidney,features=VariableFeatures(kidney))
  kidney <- RunUMAP(kidney, reduction = 'HGT', dims = 1:ncol(cell_hgt_matrix), reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
  kidney <- FindNeighbors(kidney, reduction = "HGT",dims=1:ncol(cell_hgt_matrix))
  kidney <- FindClusters(kidney,resolution = 0.5)
  sil1<-silhouette(as.numeric(Idents(kidney)), dist(cell_hgt_matrix))
  print(summary(sil1)$ avg.width)
  sil2<-silhouette(as.numeric(Idents(kidney)), dist(kidney@reductions$umap.rna@cell.embeddings))
  print(summary(sil2)$ avg.width)
  #p<-DimPlot(kidney, reduction = 'umap.rna')
  #ggsave(paste0('/fs/ess/scratch/PCON0022/wxy/umap_deepmaps/',name,'_deepmaps.svg'),p, width = 5, height = 4)
  #ARI<-igraph::compare(Idents(kidney),as.factor(label),method="adjusted.rand")
  #print(ARI)
  graph.out <-Idents(kidney)
  att<-read.csv(attPath)
  nhead<-ncol(att)
  gene_name<-rownames(GAS)[att$gene+1]
  cell_name<-colnames(GAS)[att$cell+1]
  
  att$ct<-graph.out[cell_name]
  att$gene_name<-gene_name
  att$cell_name<-cell_name
  #attention<-aggregate(x=list(att$head1,att$head2,att$head3,att$head4,att$head5,att$head6,att$head7,att$head8), by=list(att$ct,att$gene_name),mean)
  mod<-function(x){
    return(sqrt(sum(c(x^2))))
  }
  nor<-function(x){
    return((x-min(x))/(max(x)-min(x)))
  }
  alpha<-function(x){
    return(x>(mean(x)+l*sd(x)))
  }
  softmax<-function(x){
    return( exp(-x)/sum(exp(-x)))
  }
  
  att[,4:nhead]<-nor(att[,4:nhead])
  att[,4:nhead]<-att[,4:nhead]
  weight<-apply(att[,4:nhead],1,mod)
  
  
  df<-data.frame('node1'=att$gene_name,'node2'=att$cell_name,'weight'=weight)
  #att sort by cell
  link<-aggregate(df,by=list(df$node2),FUN=sort)
  t<-lapply(link$weight,alpha)
  n1<-unlist(mapply(function(x,y) x[y],link$node1,t))
  n2<-unlist(mapply(function(x,y) x[y],link$node2,t))
  w<-unlist(mapply(function(x,y) x[y],link$weight,t))
  df1<-data.frame('node1'=n1,'node2'=n2,'weight'=w*wdf1)
  gene_cor<-read.table(genePath)
  rownames(gene_cor)<-rownames(GAS)
  gene_cor<-gene_cor[unique(att$gene_name),]
  
  gene_cor<-cor(t(gene_cor))
  #g_g_cor<-rcorr(t(gene_cor),type="pearson")[[1]]
  #g_g_p<-rcorr(t(gene_cor),type="pearson")[[3]]
  #g_g_p[is.na(g_g_p)]=0
  #g_g_cor[g_g_p>0.00000001]=0
  gene_cor[gene_cor<=thdf2]=0
  #g_g_cor[g_g_p==Inf]=Inf
  g_g_cor<-gene_cor
  #weight2<-1-as.vector(gene_gene)
  g_g_cor[!upper.tri(g_g_cor,diag=FALSE)]<-0
  node1<-rep(rownames(g_g_cor),nrow(gene_cor))
  node2<-rep(rownames(g_g_cor),each=nrow(gene_cor))
  #weight2<-1-as.vector(g_g_cor)
  weight2<-as.vector(g_g_cor)
  df2<-data.frame('node1'=node1,'node2'=node2,'weight'=weight2)
  df2<-df2[df2$weight!=0,]
  h<-rbind(df1,df2)
  edges.df<-h
  G <- graph_from_data_frame(d = edges.df[, 1:2], directed = F) # construct graph
  E(G)$weight <- edges.df$weight # define the weights
  a<-list()
  for (i in unique(graph.out)){
    a[[i]]<-names(graph.out[graph.out==i])
  }
  m <- list()
  m[[1]]<-df1
  m[[2]]<-df2
  m[[3]]<-a
  m[[4]]<-names(graph.out)
  m[[5]]<-graph.out
  return (m)
}

####################RNA_RNA######################
# Integration of multiple scRNA-seq
# input:
#  ifnb.list: a list with multiple seurat object
#output:
#  object : a seurat object with a integrated matrix
integratedRNA <- function(ifnb.list){
  
  ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  features <- SelectIntegrationFeatures(object.list = ifnb.list)
  
  ## Perform integration
  obj.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
  obj.combined <- IntegrateData(anchorset = obj.anchors)
  
  ## Perform an integrated analysis
  DefaultAssay(obj.combined) <- "integrated"
  return(obj.combined)
}

####################################################################
#Cite-seq 

# splice togethor the RNA matrix and ADT matrix and do normalization
# input:
# obj: a seuart object obtained from creatobject fuction
# output:
# GAS: the spliced and normalized matrix as the input of HGT model
CLR <- function(obj){
  m1 <- obj@assays$RNA@counts[obj@assays$RNA@var.features,]
  m2 <- obj@assays$ADT@counts
  m3 <- m2
  # Add 'p_' to the name of the protein 
  rownames(m3) <- paste0('p_',rownames(m2))
  m <- NormalizeData(rbind( m1, m3), normalization.method = 'CLR', margin = 2)
  return(m)
}


# Cluster cells based on HGT embedding  
# input: 
# cell_hgt_matrixPath: the path of cell_embedding matrix 
# nhid: hyper-parameter of HGT model
# resolution: resolution of cell clustering
# obj: a seurat object obtained from creatobject fuction
# output:
# obj: a seurat object containing cell clustering result 
HGT_cluster <- function(obj, cell_hgt_matrix, nhid, resolution){
  cell_hgt_matrix<-as.matrix(cell_hgt_matrix)
  rownames(cell_hgt_matrix) <-colnames(rna)
  HGT_embedding <-
    CreateDimReducObject(
      embeddings = cell_hgt_matrix,
      key = "HGT_",
      assay = "RNA"
    )
  obj<-obj[,colnames(rna)]
  obj@reductions[['HGT']] <- HGT_embedding
  for(j in 1:ncol(obj@reductions$HGT@cell.embeddings)){
    obj[[paste0('HGT_',j)]] <- obj@reductions$HGT@cell.embeddings[,paste0('HGT_',j)]
  }
  obj <- FindNeighbors(obj, reduction = "HGT",dims=1:nhid)
  obj <- FindClusters(obj, resolution = resolution)
  
  return(obj)
}


# Draw a heatmap of marker proteins or genes
# required packages: dsb, ComplexHeatmap
# input: 
# obj: a seurat object obtained from HGT_cluster fuction 
# marker: character of markers(the order of markers corresponding to ctorder) 
# ctorder: the order of celltypes to display
# assays: RNA or ADT 
# output:
# a heatmap of marker proteins or genes

MarkerHeatmap <- function(obj,marker,ctorder,assays){
  sortmarker <- 1:length(marker)
  names(sortmarker) <- marker
  if (assays == 'RNA'){
    rnamarker <- intersect(rownames(obj@assays$RNA@counts),marker)
    rnamarker <- rnamarker[order(sortmarker[rnamarker])]
    rnam <- AverageExpression(obj,assays = 'RNA', slot='counts',group.by = 'cell_type',features = rnamarker)
    rnam <- rnam[[1]][unique(rnamarker),ctorder]
    t <- DSBNormalizeProtein(rnam,rnam,use.isotype.control =FALSE,denoise.counts =FALSE)
  }else if(assays == 'ADT'){
    adtmarker <- intersect(rownames(obj@assays$ADT@counts),marker)
    adtmarker <- adtmarker[order(sortmarker[adtmarker])]
    adtm <- AverageExpression(obj,assays = 'ADT', slot='counts',group.by = 'cell_type',features = adtmarker)
    adtm <- adtm[[1]][unique(adtmarker),ctorder]
    t <- DSBNormalizeProtein(adtm,adtm,use.isotype.control =FALSE,denoise.counts =FALSE)
  }
  Heatmap(t(t),c("white","blue"),cluster_rows = F,cluster_columns = F,
          heatmap_legend_param = list( title = "normalized.expression",title_position = "leftcenter-rot" ))
}


# take a subset of obj and preprocess data again
# input: 
# obj: a seurat object obtained from HGT_cluster fuction 
# I: cell names subset
# output:
# obj: a subset of obj after preprocessing

subobject <- function(obj,I){
  obj <- obj[,I]
  DefaultAssay(obj) <- 'RNA'
  obj <- NormalizeData(obj) %>% FindVariableFeatures()
  all.genes <- rownames(obj)
  obj <- ScaleData(obj, features = all.genes)
  DefaultAssay(obj) <- 'ADT'
  VariableFeatures(obj) <- rownames(obj[["ADT"]])
  obj <- NormalizeData(obj, normalization.method = 'CLR', margin = 2) %>% 
    ScaleData()
  return(obj)
}


# find features mostly associated with a given feature based on the expression of a type of features and two types of cells. 
# input:
# feature: a given feature (gene or protein)
# ident: an expression matrix of a type of features and two types of cells. 
# output:
# co: features mostly associated with the given feature
Findcofeatures <- function(feature,ident){
  y <- ident[feature,]
  co <- c()
  for (f in rownames(ident)){
    co <- rbind(co, c(f,cor(ident[f,],y,method = "spearman")))
    if (is.na(cor(ident[f,],y,method = "spearman"))){
      print(f)
    }
  }
  co <- co[order(co[,2],decreasing =TRUE),]
  print(which(co[,1]==feature))
  co <- co[which(co[,1]!=feature),]
  co <- co[which(co[,2]>0),]
  return(co)
}


# find all features mostly associated with top DE features based on the expression of two types of cells. 
# input:
# rna.markers: the result of DE genes analysis between two cell types from FindMarkers function
# rna_ident: an RNA expression matrix of two types of cells after discarding zero expression genes. 
# adt.markers: the result of DE proteins analysis between two cell types from FindMarkers function
# adt_ident: an ADT abundance matrix of two types of cells after discarding zero abundance proteins.
# output:
# allcofeatures: features mostly associated with the given features
Findallcofeatures <- function(rna.markers,rna_ident,adt.markers,adt_ident){
  allcofeatures <- list()
  options(scipen = 100)
  corna1 <- list()
  corna2 <- list()
  coadt1 <- list()
  coadt2 <- list()
  for (i in 1:10){
    feature <- rownames(rna.markers)[i]
    print(feature)
    co <- Findcofeatures(feature,rna_ident)
    co <- co[which(co[,1] %in% rownames(rna.markers)[1:10]==FALSE),]
    if (rna.markers[feature,"avg_log2FC"]>0){
      corna1[[feature]] <- co[1:5,1]
    }else{
      corna2[[feature]] <- co[1:5,1]
    }
  }
  for (i in 1:3){
    feature <- rownames(adt.markers)[i]
    print(feature)
    co <- Findcofeatures(feature,adt_ident)
    co <- co[which(co[,1] %in% rownames(adt.markers)[1:3]==FALSE),]
    if (adt.markers[feature,"avg_log2FC"]>0){
      coadt1[[feature]] <- co[1:5,1]
    }else{
      coadt2[[feature]] <- co[1:5,1]
    }
  }
  allcofeatures[['corna1']] <- corna1
  allcofeatures[['corna2']] <- corna2
  allcofeatures[['coadt1']] <- coadt1
  allcofeatures[['coadt2']] <- coadt2
  
  return(allcofeatures)
}


# calculate subset of obj under two cell types, the DE features and top associated features
# input: 
# obj: a seurat object obtained from HGT_cluster fuction 
# ident.1, ident.2: two cell types
# output:
# obj_co: a list of the subset of obj and selected features
subobjco <- function(obj, ident.1, ident.2){
  obj_co <- list()
  # DEG analysis between ident.1 and ident.2 
  DefaultAssay(obj) <- 'RNA'
  rna.markers <-  FindMarkers(obj, ident.1 = ident.1 ,ident.2 = ident.2)
  DefaultAssay(obj) <- 'ADT'
  adt.markers <-  FindMarkers(obj, ident.1 = ident.1,ident.2 = ident.2)
  
  I0 <- colnames(obj)[which(obj@meta.data$cell_type %in% c(ident.1,ident.2))]
  obj0 <- subobject(obj,I0)
  
  # expression data for rnas and adts of two cell types
  rna_ident  <- obj0@assays$RNA@data
  adt_ident  <- obj0@assays$ADT@data
  rna_ident <- rna_ident[which(rowSums(rna_ident)!=0),]
  adt_ident <- adt_ident[which(rowSums(adt_ident)!=0),]
  
  cofeatures <- Findallcofeatures(rna.markers = rna.markers, rna_ident = rna_ident, adt.markers = adt.markers, adt_ident = adt_ident)
  corna1 <- cofeatures[['corna1']]
  corna2 <- cofeatures[['corna2']]
  coadt1 <- cofeatures[['coadt1']]
  coadt2 <- cofeatures[['coadt2']]
  obj_co[['cofeatures']] <- cofeatures
  obj_co[['obj']] <- obj0
  obj_co[['I']] <- I0
  obj_co[['rna.markers']] <- rna.markers
  obj_co[['adt.markers']] <- adt.markers
  return(obj_co)
}


# draw heatmap on correlation matrix between selected features 
# required packages: RColorBrewer
# input: 
# obj_co: a list obtained from subbmco function
# output:
# heatmap: heatmap on correlation matrix between selected features 
#library(RColorBrewer)
DEfeaturesheatmap <- function(obj_co){
  obj0 <- obj_co$obj
  corna1 <- obj_co$cofeatures$corna1
  corna2 <- obj_co$cofeatures$corna2
  coadt1 <- obj_co$cofeatures$coadt1 
  coadt2 <- obj_co$cofeatures$coadt2
  corna <- c(corna1,corna2)
  coadt <- c(coadt1,coadt2)
  # mm - correlations matrix as input of Heatmap
  m <-rbind(obj0@assays$RNA@scale.data[unique(c(unlist(corna), names(corna))),], obj0@assays$ADT@scale.data[unique(c(unlist(coadt), names(coadt))),])
  rownames(m) <- c(unique(c(unlist(corna), names(corna))),unique(c(unlist(coadt), names(coadt))))
  mm <- cor(t(m))
  # lab - the coressponding cel types of features
  lab <- list()
  for (f in rownames(m)){
    if (f %in% unique(c(unlist(corna1), names(corna1),unlist(coadt1), names(coadt1)))){
      lab[[f]] <- color1
    }else{
      lab[[f]] <- color2
    }
  }
  p <- heatmap(mm,RowSideColors=unlist(lab),col=colorRampPalette(c(rep("royalblue",1),'white','#F6DBDB',rep("firebrick3",2)))(56))
  return(p)
}


# draw a line representing a given feature
# input锛? 
# f: a given feature
# dff: a dataframe containing expression of selected features and a dimension values of HGT embedding.
# color : color of a line 
# output:
# draw a line representing a given feature
line <- function(f,dff,color){
  df <- data.frame(cbind(dff[,f],dff$index))
  model1=loess(X1 ~ X2,data=df,span=1 )
  model <- stats::predict(model1)
  lines(model, x=df$X2, col=color)
}


# draw the lines by a loess smoothing function based on the corresponding embedding and scaled gene expressions in cells
# input: 
# obj: a seurat object obtained from HGT_cluster fuction 
# obj_co: a list obtained from subbmco function
# HGT.dim: a dimension of HGT embedding
# color1, color2: two colors coressponding to two cell types
# output:
# a plot show the relationship between HGT embedding and feature expression 
DEfeaturesloess <- function(obj, obj_co, HGT.dim = n, color1, color2){
  obj0 <- obj_co$obj
  I0 <- obj_co$I
  corna1 <- obj_co$cofeatures$corna1
  corna2 <- obj_co$cofeatures$corna2
  coadt1 <- obj_co$cofeatures$coadt1 
  coadt2 <- obj_co$cofeatures$coadt2
  corna <- c(corna1,corna2)
  coadt <- c(coadt1,coadt2)
  
  m <-obj0@assays$RNA@scale.data[unique(c(unlist(corna), names(corna))),]
  df0 <- data.frame(t(m))
  colnames(df0) <- rownames(m)
  df0$index <-unlist(obj[[paste0('HGT_', HGT.dim)]][I0,])
  df0  <- df0[order(df0$index),]
  m <-obj0@assays$ADT@scale.data[unique(c(unlist(coadt), names(coadt))),]
  df1 <- data.frame(t(m))
  colnames(df1) <- rownames(m)
  df1$index <-unlist(obj[[paste0('HGT_', HGT.dim)]][I0,])
  df1  <- df1[order(df1$index),]
  
  rna.markers <- obj_co$rna.markers
  adt.markers <- obj_co$adt.markers
  par(pin = c(3,3))
  y=seq(-2,2,0.1)
  x=seq(min(df0$index),max(df0$index),length.out=length(y))
  #x=seq(-0.1,0.05,length.out=length(y))
  p <- plot(x,y,col="white",xlab = paste0('HGT_', HGT.dim), ylab="scaled experssion",type="l",main="Loess Smoothing and Prediction")
  for (j in 1:10){
    feature <- rownames(rna.markers)[j]
    if (rna.markers[feature,'avg_log2FC']>0){
      color <- color1
    }else{
      color <- color2
    }
    line (feature,df0,color)
    for (cof in unlist(corna[[feature]])[1:4]){
      line (cof,df0,color)
    }
  }
  for (j in 1:3){
    feature <- rownames(adt.markers)[j]
    if (adt.markers[feature,'avg_log2FC']>0){
      color <- color1
    }else{
      color <- color2
    }
    line (feature,df1,color)
    for (cof in unlist(coadt[[feature]])[1:4]){
      line (cof,df1,color)
    }
  }
  return(p)
}




