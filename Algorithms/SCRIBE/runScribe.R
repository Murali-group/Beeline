suppressPackageStartupMessages(library(monocle, warn.conflicts = FALSE , quietly = TRUE))
suppressPackageStartupMessages(library(Scribe, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library (optparse, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library (igraph, warn.conflicts = FALSE, quietly = TRUE))

cal_ncenter <- function(ncells, ncells_limit = 100){
  round(2 * ncells_limit * log(ncells)/ (log(ncells) + log(ncells_limit)))
}

option_list <- list (
              make_option(c("-e","--expressionFile"), type = 'character',
              help= "Path to comma separated file containing gene-by-cell matrix with
              cell names as the first row and gene names as 
              the first column. Required if -n flag is not used."),
              
              make_option(c("-c","--cellFile"), type = 'character',
              help= "Path to comma separated file containing data on cells.
              First column is Cell name (without a header). Subsequent columns
              contain information on cells (experiment, time, etc).
              Required if -n flag is not used."),
              
              make_option(c("-g","--geneFile"), type = 'character',
              help= "Path to comma separated file containing data on genes.
              First column is gene name (without a header). Subsequent columns
              contain information on genes such as gene_short_name (required). 
              Required if -n flag is not used."),
              
              make_option(c("-n","--newCellDataSet"), type = 'character',
              help= "Path to .RDS file containing an object of type
              newCellDataSet."),
              
              make_option(c("-l","--lowerDetectionLimit"), default = 0.0, type = 'double',
              help= "Single float value to pass as an argument for newCellDataSet function.
              default = %default ."),
              
              make_option(c("-x","--expressionFamily"), default = 'uninormal', type = 'character',
              help= "VGAM family function name (without parantheses) to be used for 
              expression response variables. 
              See more here: http://cole-trapnell-lab.github.io/monocle-release/docs
              default = %default ."),
              
              make_option(c("-m","--method"), default = 'ucRDI', type = 'character',
              help= "Method name for Scribe. Can be any one of 'RDI', 'cRDI', 'uRDI', or 'ucRDI'.
              default = %default ."),
              
              make_option(c("-d","--delay"), default = '1', type = 'character',
              help= "Comma separated list of delay values for Scribe. Maximum delay
              value should be the total number of cells in the dataset. default = %default ."),
              
              make_option(c("--log"), action = 'store_true', default = FALSE, type = 'character',
              help= "Log transform expression values. default = %default ."),
              
              make_option(c("-o","--outPrefix"), , type = 'character',
              help= "Path to write output files. Required."),
              
              make_option(c("","--outFile"), , type = 'character',
              help= "outFile name to write the output ranked edges. Required."),
    
              make_option(c("-i","--ignorePT"), action = 'store_true', default = FALSE, 
              type = 'character',
              help= "Ignores pseudotime computed using monocle and uses experiment time.")
              )

parser <- OptionParser(option_list = option_list)
arguments <- parse_args(parser, positional_arguments = FALSE)
print(arguments)


if (length(arguments$outPrefix) == 0){
 stop('outPrefix is required.
      Run Rscript runScribe.R -h for more details.') 
}

# if the new cell dataset file is not provided create one
if (length(arguments$newCellDataSet) == 0){
  
if (length(arguments$expressionFile) == 0){
    stop("Please enter the path to expression file using the -e flag. 
  Run Rscript runScribe.R -h for more details.")
  }
  
if (length(arguments$cellFile) == 0){
    stop("Please enter the path to cell file using the -c flag.
       Run Rscript runScribe.R -h for more details.")
  }
  
  if (length(arguments$geneFile) == 0){
    stop("Please enter the path to gene file using the -g flag.
       Run Rscript runScribe.R -h for more details.")
  }
  
  # Read data
  # Set check.names to False to avoid R adding an 'X' to the beginning of columns that start with an integer
  exprMatrix <- read.delim(arguments$expressionFile, row.names = 1, sep = ',', check.names=FALSE)
  cellData <- read.delim(arguments$cellFile, row.names = 1, sep = ',')
  geneData <- read.delim(arguments$geneFile, row.names = 1, sep = ',')
  cd <- new("AnnotatedDataFrame", data = cellData)
  gd <- new("AnnotatedDataFrame", data = geneData)

# Use uninormal if it is simulated data
if (arguments$expressionFamily == 'uninormal'){
  cat("Using uninormal() as expression family.\n")
  CDS <- newCellDataSet(as(as.matrix(exprMatrix), "sparseMatrix"),
                       phenoData = cd,
                       featureData = gd,
                       lowerDetectionLimit = arguments$lowerDetectionLimit,
                       expressionFamily = uninormal())

sizeFactors(CDS) <- 1 # Same as the one in neuronal_sim_cCDS
} else{
  # For scRNA-Seq data (counts, RPKM/FPKM)
  cat("Using negbinomial.size() as expression family.\n")
CDS <- newCellDataSet(as(as.matrix(exprMatrix), "sparseMatrix"),
                       phenoData = cd,
                       featureData = gd,
                       lowerDetectionLimit = arguments$lowerDetectionLimit,
                       expressionFamily = negbinomial.size())
  
CDS <- estimateSizeFactors(CDS)
CDS <- estimateDispersions(CDS)
disp_table <- dispersionTable(CDS)
CDS <- setOrderingFilter(CDS, ordering_genes)
}

if (arguments$ignorePT == TRUE){
  cat("Using experimental time instead of PseudoTime computed using
      monocle.\n")
  CDS$Pseudotime <- CDS$Time
  CDS@phenoData@data$Pseudotime <- CDS@phenoData@data$Time
  CDS@phenoData@data$State <- 1
  CDS$State <- 1
  
}else{
  
cat("Computing pseudotime.\n")
CDS <- reduceDimension(CDS, norm_method ="none")
CDS <- orderCells(CDS)

saveRDS(CDS, file= paste0(arguments$outPrefix,'dataset.RDS'))
write.csv(CDS@phenoData@data, file= paste0(arguments$outPrefix,'PseudoTime.csv'), quote = FALSE)
}
}

# If newcelldataset RDS is available already
if (length(arguments$newCellDataSet) != 0){
  CDS <- readRDS(arguments$newCellDataSet)
}

### Run Scribe

cat("Computing",arguments$method,"\n")


delay <- as.numeric(strsplit(arguments$delay, ",")[[1]])

if (arguments$method == 'uRDI'){
  net <- calculate_rdi(CDS, delays = delay, method = 2, uniformalize = TRUE, log = arguments$log)
  netOut <- net$max_rdi_value
  # computes CLR if we use uRDI
  # TODO: Make this an optional
  netOut <- clr(netOut)
} else if (arguments$method == 'ucRDI'){
  net <- calculate_rdi(CDS, delays = delay, method = 2, uniformalize = TRUE, log = arguments$log)
  netOut <- calculate_conditioned_rdi(CDS, rdi_list = net, uniformalize = TRUE, log = arguments$log)
} else if (arguments$method == 'RDI'){
  net <- calculate_rdi(CDS, delays = delay, method = 2, uniformalize = FALSE, log = arguments$log)
  netOut <- net$max_rdi_value
  # computes CLR if we use RDI
  # TODO: Make this optional
  netOut <- clr(netOut)
} else if (arguments$method == 'cRDI'){
  net <- calculate_rdi(CDS, delays = delay, method = 2, uniformalize = FALSE, log = arguments$log)
  netOut <- calculate_conditioned_rdi(CDS, rdi_list = net, uniformalize = FALSE, log = arguments$log)
} else{
  stop("Method must be one of RDI, cRDI, uRDI, or ucRDI. 
       Run Rscript runScribe.R -h for more details.")
}
outGraph <- graph_from_adjacency_matrix(netOut, mode = 'directed', weighted=T)
write.graph(outGraph, paste0(arguments$outPrefix,arguments$outFile),"ncol")
cat("Done.\n")
