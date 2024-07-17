#' ng_sem_estimation
#'
#' a wrapper function to perform NG-SEM on a specific gene
#'
#' @author Jiyaing Zhao, \email{jyzhao@connect.hku.hk}
#'
#' @export
#'
#' @param X scRNA-seq expression input
#' @param params a list of parameters n.K=n.K,n.N=n.N,Data.Y=Data.Y
#' \itemize{
#'  \item{n.K} number of Gaussian components
#'  \item{n.N} number of samples (T-1)
#'  \item{Data.Y} design matrix
#'  }
#' @param g gene to estimate
#' @param g.seed seed for reproducibility
#' @return res_for_gene
#' \itemize{
#'  \item{Lik_list} likelihood recordings
#'  \item{Iter_num} number of iterations
#'  \item{Comp_w} a vector for the weights of the Gaussian Mixture components
#'  \item{Comp_sig2} a vector for the variances of the Gaussian Mixture components
#'  \item{grn_w_g} grn weights for gene g
#'  }
#'

ng_sem_estimation <- function(X,g,params,g.seed){

  set.seed(g.seed)

  Data.Y <- params$Data.Y
  n.G <- params$n.G
  n.N <- params$n.N
  n.K <- params$n.K
  param.miter <- params$param.miter
  param.error <- params$param.error


  Data.y <- X[g,2:ncol(X)]
  Data.y <- as.numeric(Data.y)

  model.a <- as.numeric(MASS::ginv(t(Data.Y))%*%Data.y)
  model.mu <-  rep(0,n.K)
  model.Sigma2 <- runif(n.K)+1e-9

  prod.a.x <- rep(0,n.N)
  for(i in seq_len(n.N)){
    prod.a.x[i] <-sum(model.a*Data.Y[,i])
  }
  res_y_ax <- Data.y-prod.a.x

  Data <- list(Y=Data.Y,y=Data.y)

  # Random Initialization
  idx <- sample(1:n.N, n.K)
  m <- res_y_ax[idx]

  m1 <- matrix(unlist(lapply(m, function(x){
    x*res_y_ax
  })),nrow = as.integer(n.K),byrow = T)
  m2 <- sweep(m1,1,m*m/2,FUN = "-")
  idx.comp <- unlist(apply(m2,2,function(x){
    which(x==max(x))
  }))
  Gamma <- matrix(0,nrow=n.N,ncol=as.integer(n.K))
  for(i in seq_len(n.N)){
    Gamma[i,idx.comp[i]] <- 1
  }
  Nk <- apply(Gamma,2,sum)
  model.w <- Nk/sum(Nk)

  converged=F
  iter=1
  liklih <- c()
  model=list(mu=model.mu,Sigma2=model.Sigma2,w=model.w,a=model.a)

  liklih.pre <- -Inf

  while ((!converged) & iter < param.miter ) {

    iter=iter+1

    #--------E step
    e.step <- ng_sem_expectation(residual=res_y_ax,model = model)
    liklih.now <- e.step$liklih
    liklih <- c(liklih,e.step$liklih)

    #--------M step
    model <- ng_sem_maximization(e.step,Data,res_y_ax,model)

    prod.a.x <- rep(0,n.N)
    for(i in seq_len(n.N)){
      prod.a.x[i] <-sum(model$a*Data.Y[,i])
    }
    res_y_ax <- Data.y-prod.a.x
    converged <- ((liklih.now-liklih.pre) < as.numeric(param.error)*abs(liklih.now))
    liklih.pre <- liklih.now
  }

  # normalization to ensure the maximal weight is 1
  if(max(model$a)!=0){
    model$a=model$a/max(model$a)
  }

  res_for_gene <- list(Lik_list=liklih,Iter_num=iter,Comp_w=model$w,Comp_sig2=model$Sigma2,grn_w_g=as.numeric(model$a))
  res_for_gene
}
