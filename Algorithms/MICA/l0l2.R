options(repos = c(CRAN = "https://cran.r-project.org"))
if (!requireNamespace("L0Learn", quietly = TRUE))
  install.packages("L0Learn", quiet = TRUE)
if (!requireNamespace("parallel", quietly = TRUE))
  install.packages("parallel", quiet = TRUE)
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr", quiet = TRUE)
if (!requireNamespace("tidyr", quietly = TRUE))
  install.packages("tidyr", quiet = TRUE)

suppressPackageStartupMessages(library(L0Learn, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(parallel, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(dplyr,  warn.conflicts = FALSE, quietly = TRUE ))
suppressPackageStartupMessages(library(tidyr,  warn.conflicts = FALSE, quietly = TRUE ))


#' L0L2 sparse regression
#'
#' Given an n x m data matrix X (n predictors and m samples) and a m x 1 response vector y, such that
#' y = t(X)B + e, estimate B from the data (y, X).
#' @param folds integer; k for the internal k-fold cross-validation to tune
#' the hyperparameters lambda and gamma
#' @param scale boolean; Indicates whether the data should be scaled before
#' applying regression.
#'
#' @return A list with the coefficients that solve the regression problem.
#'
#' @author Tom Bartlett \email{thomas.bartlett.10@ucl.ac.uk}
#' @author Gregorio Alanis-Lobato \email{gregorio.alanis@crick.ac.uk}

solve_sparse_reg <- function(X, y, folds = 5, scale = FALSE, propor = 0.05){
  if(scale){
    y <- y / max(sd(y), 1e-6)
    X <- X / pmax(apply(X, 1, sd), 1e-6)
  }
  X <- t(X)

  # Use k-fold cross-validation to select the optimal values of the tuning
  # lambda and gamma
  cvfit <- L0Learn.cvfit(X, y, penalty = "L0L2", nFolds = folds,
                         algorithm = "CDPSI", maxSuppSize = propor*ncol(X))
  bestGamma_idx <- which.min(lapply(cvfit$cvMeans, min))
  if(length(bestGamma_idx)<1){
    print(bestGamma_idx)
  }
  bestLambda_idx <- which.min(cvfit$cvMeans[[bestGamma_idx]])
  bestLambda <- cvfit$fit$lambda[[bestGamma_idx]][bestLambda_idx]
  bestGamma <- cvfit$fit$gamma[bestGamma_idx]
  coeffs <- coef(cvfit, lambda = bestLambda, gamma = bestGamma)[, 1]
  coeffs <- coeffs[names(coeffs) != "Intercept"]

  names(coeffs) <- colnames(X)
  return(coeffs)
}

#' L0L2 regulatory network inference
#'
#' Given an n x m gene expression (n genes and m samples) and a list of
#' regulators assign likelihood scores to all possible regulator x gene
#' interactions using L0L2 sparse regression.
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
#' @author Tom Bartlett \email{thomas.bartlett.10@ucl.ac.uk}
#' @author Gregorio Alanis-Lobato \email{gregorio.alanis@crick.ac.uk}

l0l2 <- function(exprMatr, regulators = NULL, nCores = 1){
  if(is.null(regulators)){
    regulators <- rownames(exprMatr)
  }else{
    regulators <- regulators[regulators %in% rownames(exprMatr)]
  }
  targets <- rownames(exprMatr)

  reg_coeffs <- mclapply(seq_along(targets),
                        function(i){
                          solve_sparse_reg(X = as.matrix(exprMatr[regulators,]),
                                           y = unlist(exprMatr[targets[i],]),
                                           folds = 3, scale = TRUE,
                                           propor = 1)
                        },
                        mc.cores = nCores)

  linkList <- matrix(unlist(reg_coeffs),
                     nrow = length(regulators),
                     ncol = length(targets),
                     dimnames = list(regulators, targets)) %>%
    as_tibble(rownames = "regulatoryGene") %>%
    pivot_longer(cols = all_of(targets), names_to = "targetGene",
                 values_to = "weight") %>%
    filter(regulatoryGene != targetGene) %>%
    arrange(desc(weight))

  return(linkList)
}

# Refine L0L2 predictions with TF-gene binding data from ChIP-seq or ATAC-seq.
refined_l0l2 <- function(exprMatr, regulators = NULL, tf_binding = NULL,
                         nCores = 1){
  linkList <- l0l2(exprMatr, regulators, nCores)
  if(!is.null(tf_binding)){
    linkList <- left_join(tf_binding, linkList,
                          by = c("regulatoryGene", "targetGene")) %>%
      filter(!is.na(weight)) %>%
      arrange(desc(weight))
  }
  return(linkList)
}
