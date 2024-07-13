options(repos = c(CRAN = "https://cran.r-project.org"))
if (!requireNamespace("optparse", quietly = TRUE))
  install.packages("optparse", quiet = TRUE)
if (!requireNamespace("igraph", quietly = TRUE))
  install.packages("igraph", quiet = TRUE)
if (!requireNamespace("rsample", quietly = TRUE))
  install.packages("rsample", quiet = TRUE)
if (!requireNamespace("purrr", quietly = TRUE))
  install.packages("purrr", quiet = TRUE)
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr", quiet = TRUE)

suppressPackageStartupMessages(library(optparse, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(igraph, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(rsample,  warn.conflicts = FALSE, quietly = TRUE ))
suppressPackageStartupMessages(library(purrr,  warn.conflicts = FALSE, quietly = TRUE ))
suppressPackageStartupMessages(library(dplyr,  warn.conflicts = FALSE, quietly = TRUE ))

option_list <- list(
              make_option(c("-e","--expressionFile"), type = 'character',
              help= "Path to comma separated file or tab seperated file containing gene-by-cell
              matrix with cell names as the first row and gene names as  the first column. "),

              make_option(c("-a","--atacFile"), type = 'character',
              help= "Path to ATAC-seq dataset with Chromatin Accessibility."),

              make_option(c("-r", "--regFile"), type = 'character',
              help = "Path to file containing information on regulatory genes among all
              genes in expression file"),

              make_option(c("-m","--method"), default = 'mica', type = 'character',
              help= "Method name for MICA. Can be any one of 'mica': Mutual Information and Chromatin Accessibility,
              'l0l2' : Sparse Regression, 'spearman': Spearman correlation , or 'genie3' : GENIE3 R module.
              default = %default ."),

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

if (length(arguments$expressionFile) == 0){
    stop("Please enter the path to expression file using the -e flag.
          Run Rscript runMICA.R -h for more details.")
}

exprMatr <- read.delim(arguments$expressionFile,  row.names = 1, sep = ',', check.names=FALSE)
exprMatr <- exprMatr[rowSums(exprMatr) > 0, ]

if (length(arguments$atacFile) != 0){
  atac <- read.delim(arguments$atacFile, row.names = 1, sep = ',', check.names = FALSE)
  atac <- atac %>%
    select(tf, nearest_prom) %>%
    rename(regulatoryGene = tf, targetGene = nearest_prom)

  exprMatr <- exprMatr[intersect(rownames(exprMatr),
                                 union(atac$regulatoryGene, atac$targetGene)), ]

}
regulators = NULL
if (length(arguments$regFile) != 0){
  regulators <- read.delim(arguments$regFile, row.names = 1, sep = ',', check.names = FALSE)
  regulators <- regulators[regulators %in% rownames(exprMatr)]

}

cat("Computing",arguments$method,"\n")

if (arguments$method == 'mica'){
    source('mutual_info.R')
    net = mutual_info(exprMatr, regulators, 1)

} else if (arguments$method == 'l0l2'){
    source('l0l2.R')
    net = l0l2(exprMatr, regulators, 1)


} else if (arguments$method == 'spearman'){
    source('spearman.R')
    net = spearman(exprMatr, regulators, 1)

} else if (arguments$method == 'genie3'){
    source('genie3.R')
    net = refined_genie3(exprMatr)

} else{
  stop("Method must be one of mica, l0l2, spearman or genie3.
       Run Rscript runMICA.R -h for more details.")
}


write.csv(net, paste0(arguments$outPrefix,arguments$outFile))

cat("Done.\n")