options(repos = c(CRAN = "https://cran.r-project.org"))
if (!requireNamespace("parallel", quietly = TRUE))
  install.packages("parallel", quiet = TRUE)
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr", quiet = TRUE)

suppressPackageStartupMessages(library(parallel, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE))


#' Spearman's r regulatory network inference
#'
#' Given an n x m gene expression (n genes and m samples) and a list of
#' regulators assign likelihood scores to all possible regulator x gene
#' interactions using Spearman's correlation.
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
#' by correlation.
#'
#' @author Tom Bartlett \email{thomas.bartlett.10@ucl.ac.uk}
#' @author Gregorio Alanis-Lobato \email{gregorio.alanis@crick.ac.uk}

spearman <- function(exprMatr, regulators = NULL, nCores = 1){
  if(is.null(regulators)){
    regulators <- rownames(exprMatr)
  }else{
    regulators <- regulators[regulators %in% rownames(exprMatr)]
  }

  linkList <- as_tibble(expand.grid(regulators, rownames(exprMatr),
                                    stringsAsFactors = FALSE)) %>%
    filter(Var1 != Var2) %>%
    rename(regulatoryGene = Var1, targetGene = Var2)

  spearman_scores <- mclapply(1:nrow(linkList),
                              function(i){
                                cor(x = unlist(exprMatr[linkList$regulatoryGene[i],]),
                                    y = unlist(exprMatr[linkList$targetGene[i],]),
                                    method = "spearman")
                              },
                              mc.cores = nCores)
  linkList <- linkList %>%
    mutate(weight = unlist(spearman_scores)) %>%
    arrange(desc(weight))
  return(linkList)
}

# Refine Spearman's predictions with TF-gene binding data from ChIP-seq or ATAC-seq.
refined_spearman <- function(exprMatr, regulators = NULL, tf_binding = NULL,
                             nCores = 1){
  linkList <- spearman(exprMatr, regulators, nCores)
  if(!is.null(tf_binding)){
    linkList <- left_join(tf_binding, linkList,
                          by = c("regulatoryGene", "targetGene")) %>%
      filter(!is.na(weight)) %>%
      arrange(desc(weight))
  }
  return(linkList)
}
