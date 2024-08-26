options(repos = c(CRAN = "https://cran.r-project.org"))
if (!requireNamespace("splines", quietly = TRUE))
  install.packages("splines", quiet = TRUE)
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr", quiet = TRUE)
if (!requireNamespace("parallel", quietly = TRUE))
  install.packages("parallel", quiet = TRUE)

suppressPackageStartupMessages(library(splines, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library(parallel,  warn.conflicts = FALSE, quietly = TRUE ))


mutual.information2.backend <- function(x1,x2=NULL,bins=10,spline.order=1,correct=FALSE,min.def=0)
{
  stopifnot(require(splines))
  if(spline.order > 1 & correct == TRUE){stop("The correction is only available for 'spline.order = 1'")}
  stopifnot(is.vector(x1))
  stopifnot(is.vector(x2))

  x1.def <- x1[complete.cases(x1,x2)]
  x2.def <- x2[complete.cases(x1,x2)]

  if( length(which(complete.cases(x1,x2))) / length(x1) < min.def )
  {
    (mi <- NA)
  } else
  {
    num.knots <- bins + spline.order
    knots     <- seq(0,(num.knots - 1),by=1)

    bspline.min <- knots[spline.order]
    bspline.max <- knots[length(knots) + 1 - spline.order]
    #print(paste("num.knots =",num.knots,"bsp.min =",bspline.min,"bsp.max =",bspline.max))

    # tranform data linearly into domain of B-spline functions
    x1.tr <- (x1.def - min(x1.def)) * (bspline.max - bspline.min) / (max(x1.def) - min(x1.def)) + bspline.min
    x2.tr <- (x2.def - min(x2.def)) * (bspline.max - bspline.min) / (max(x2.def) - min(x2.def)) + bspline.min

    if(all(is.na(x1.tr)) || all(is.na(x2.tr))){
      return(0)
    }else{

      x1.spline <- splineDesign(knots=knots,x=x1.tr,ord=spline.order,outer.ok = TRUE)
      x2.spline <- splineDesign(knots=knots,x=x2.tr,ord=spline.order,outer.ok = TRUE)

      p.x1    <- colSums(x1.spline) / length(x1.def)
      p.x2    <- colSums(x2.spline) / length(x2.def)
      p.x1.x2 <- as.vector( (t(x1.spline) %*% x2.spline) / length(x1.def) )

      h.x1    <- -sum(p.x1    * log2(p.x1)   ,na.rm=T)
      h.x2    <- -sum(p.x2    * log2(p.x2)   ,na.rm=T)
      h.x1.x2 <- -sum(p.x1.x2 * log2(p.x1.x2),na.rm=T)

      mi <- h.x1 + h.x2 - h.x1.x2

      if( correct == TRUE ){
        return(mi - (bins - 1) / (2 * length(x1.def)))
      }
      else{
        return(mi)
      }
    }
  }
}

mutual.information2 <- function(x1,x2=NULL,bins=10,spline.order=1,correct=FALSE,min.def=0)
{
  if(is.vector(x1))
  {
    stopifnot(is.vector(x2))
    mi <- mutual.information2.backend(x1,x2,bins=bins,spline.order=spline.order,correct=correct,min.def=min.def)
  } else
    if(is.data.frame(x1) || is.matrix(x1))
    {
      if (is.data.frame(x1)) { x1 <- as.matrix(x1) }
      if (is.null(x2)) { x2 <- x1 } else
      {
        if (is.data.frame(x2)) { x2 <- as.matrix(x2) }
      }

      ncx <- ncol(x1)
      ncy <- ncol(x2)

      mi <- matrix(0, nrow = ncx, ncol = ncy)
      rownames(mi) <- colnames(x1)
      colnames(mi) <- colnames(x2)

      for (i in 1:ncx)
      {
        for (j in 1:i)
        {
          mi[j,i] <- mi[i,j] <- mutual.information2.backend(x1[,i],x2[,j],bins=bins,spline.order=spline.order,correct=correct,min.def=min.def)
        }
      }
    }
  mi
}

#' MI regulatory network inference
#'
#' Given an n x m gene expression (n genes and m samples) and a list of
#' regulators assign likelihood scores to all possible regulator x gene
#' interactions using the spline estimation of mutual information (MI).
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

mutual_info <- function(exprMatr, regulators = NULL, nCores = 1){
  if(is.null(regulators)){
    regulators <- rownames(exprMatr)
  }else{
    regulators <- regulators[regulators %in% rownames(exprMatr)]
  }

  linkList <- as_tibble(expand.grid(regulators, rownames(exprMatr),
                                    stringsAsFactors = FALSE)) %>%
    filter(Var1 != Var2) %>%
    rename(regulatoryGene = Var1, targetGene = Var2)

  mi_scores <- mclapply(1:nrow(linkList),  function(i){
                          mutual.information2(x1 = unlist(exprMatr[linkList$regulatoryGene[i],]),
                                              x2 = unlist(exprMatr[linkList$targetGene[i],]),
                                              bins = floor(ncol(exprMatr)/5), spline.order = 3)
                          },
                        mc.cores = nCores)
  linkList <- linkList %>%
    mutate(weight = unlist(mi_scores)) %>%
    arrange(desc(weight))
  return(linkList)
}

# Refine MI predictions with TF-gene binding data from ChIP-seq or ATAC-seq.
refined_mi <- function(exprMatr, regulators = NULL, tf_binding = NULL,
                                nCores = 1){
  linkList <- mutual_info(exprMatr, regulators, nCores)
  if(!is.null(tf_binding)){
    linkList <- left_join(tf_binding, linkList,
                          by = c("regulatoryGene", "targetGene")) %>%
      filter(!is.na(weight)) %>%
      arrange(desc(weight))
  }
  return(linkList)
}
