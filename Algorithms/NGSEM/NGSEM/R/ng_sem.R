#' ng_sem
#'
#' a wrapper function to perform NG-SEM (CPU parallelization)
#'
#' @author Jiaying Zhao, \email{jyzhao@connect.hku.hk}
#'
#' @import parallel
#' @import foreach
#' @import MASS
#'
#' @export
#' @param dat single cell RNA-seq expression data matrix
#' @param nk number of Gaussian Components (default:3)
#' @param param.miter maximum iterations (default:5000)
#' @param param.error error bound (default:1e-6)
#' @param n.cores number of working CPU cores (default:24)
#' @param params Defalut is 500
#' @param g gene to estimate
#' @param g.seed seed for reproducibility
#' @param k Components
#' @return res: a list w.mat = w.mat,
#' \itemize{
#'  \item{w.mat} adjacency matrix for GRN
#'  \item{likelihood} a list of length for gene numbers containing likelihood recordings
#'  \item{iteration} a vector of length for gene numbers containing iteration recordings
#'  \item{Comp_Weight} a gene-by-component matrix containing component weights
#'  \item{Comp_Sigma2} a gene-by-component matrix containing component variance
#'  \item{time} CPU time in seconds
#'  }
#'
ng_sem <- function(dat,
                     nk=3,
                     param.miter=5000,
                     param.error=1e-6,
                     n.cores = 24,
                     seed_list = NULL
                     ){

  clusters <- parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(clusters)

  if((length(grep("matrix",class(dat)))==0)){
    dat <- as.matrix(dat)
  }

  X <- dat

  Data.Y <- X[,1:ncol(X)-1]
  n.G <- nrow(Data.Y)
  n.N <- ncol(Data.Y)

  if(is.null(seed_list)){
    seed_list <- rep(123456,n.G)
  }

  grn_W <- matrix(0,nrow=n.G,ncol=n.G,n.G)

  n.K <- nk

  params=list(n.K=n.K,n.N=n.N,Data.Y=Data.Y,param.error=param.error,param.miter=param.miter)

  start = Sys.time()

  res <- foreach(g=1:n.G, .inorder=TRUE)  %dopar%  {
    ng_sem_estimation(X=X,g=g,params=params,g.seed=seed_list[g])
  }
  end = Sys.time()
  time = as.numeric(difftime(end,start,units = "secs"))

  # for(i in seq_len(n.G)){
  #   grn_W[i,] <- res[[i]]$grn_w_g
  # }
  # w.mat <- grn_W

  Lik <- lapply(res, `[[`, 1)
  Iter <- unlist(lapply(res, `[[`, 2))
  Comp_Weight <- matrix(unlist(lapply(res, `[[`, 3)),nrow = n.G,byrow=T)
  Comp_Sigma2 <- matrix(unlist(lapply(res, `[[`, 4)),nrow = n.G,byrow=T)
  grn_W <- matrix(unlist(lapply(res, `[[`, 5)),nrow = n.G,byrow=T)

  w.mat <- grn_W
  rownames(w.mat) <- rownames(dat)
  colnames(w.mat) <- rownames(dat)

  res=list(
    w.mat = w.mat,
    likelihood = Lik,
    iteration = Iter,
    Comp_Weight = Comp_Weight,
    Comp_Sigma2 = Comp_Sigma2,
    time=time
  )

  parallel::stopCluster(clusters)

  return(res)
}




