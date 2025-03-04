options(repos = c(CRAN = "https://cran.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", quiet = TRUE)
if (!requireNamespace("GENIE3", quietly = TRUE))
  BiocManager::install("GENIE3", quiet = TRUE)
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr", quiet = TRUE)

suppressPackageStartupMessages(library(GENIE3, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE))

#' GENIE3 regulatory network inference
#'
#' Given an n x m gene expression (n genes and m samples) and a list of
#' regulators assign likelihood scores to all possible regulator x gene
#' interactions using GENIE3.
#'
#' @param exprMatr matrix; A n x m expression matrix (n genes and m samples).
#' Rows must be named with gene symbols or IDs and columns with sample IDs.
#' @param regulators character; A vector of regulator genes (e.g. transcription
#' factors). If NULL, all genes in the expression matrix are used as regulators.
#' Regulators that are not listed the expression matrix are discarded.
#' @param nCores integer; Number of cores to run the network inference in
#' parallel.
#'
#' @return A tibble with regulatory links regulatoryGene-->targetGene sorted
#' by the likelihood score weight.
#'
#' @author Gregorio Alanis-Lobato \email{gregorio.alanis@crick.ac.uk}

genie3_net <- function(exprMatr, regulators = NULL, nCores = 1){
  if(is.null(regulators)){
    regulators <- rownames(exprMatr)
  }else{
    regulators <- regulators[regulators %in% rownames(exprMatr)]
  }
  weightMat <- GENIE3(as.matrix(exprMatr), regulators=regulators, nCores = nCores)
  linkList <- getLinkList(weightMat) %>%
    as_tibble() %>%
    arrange(desc(weight))
  return(linkList)
}

# Refine GENIE3 predictions with TF-gene binding data from ChIP-seq or ATAC-seq.
refined_genie3 <- function(exprMatr, regulators = NULL, tf_binding = NULL, nCores = 1){
  linkList <- genie3_net(exprMatr, regulators, nCores)
  if(!is.null(tf_binding)){
    linkList <- left_join(tf_binding, linkList,
                          by = c("regulatoryGene", "targetGene")) %>%
      filter(!is.na(weight)) %>%
      arrange(desc(weight))
  }
  return(linkList)
}
