#' Expectation step
#'
#' This function is used for the E step in the optimization
#'
#' @author Jiaying Zhao, \email{jyzhao@connect.hku.hk}
#'
#' @export
#'
#' @param residual residual for SEM y-model.a*x, a vecotr of length n.N(sample length)
#' @param model Gaussian Mixture model
#' \itemize{
#'  \item{mu} a vector for the mean of the Gaussian Mixture components
#'  \item{Sigma2} a vector for the variance of the Gaussian Mixture components
#'  \item{w} a vector for the weights of the Gaussian Mixture components
#'  }
#' @return res
#' \itemize{
#'  \item{Gamma} a matrix
#'  \item{liklih}
#'  }
ng_sem_expectation <- function(residual,model){
  mu <- model$mu
  Sigma2 <- model$Sigma2
  w <- model$w

  n.K <- length(mu)
  n.N <- length(residual)
  logRho <- matrix(0,nrow=n.N,ncol=n.K)

  for(i in seq_len(n.K)){
    logRho[,i] <- loggaussian(residual,mu[i],Sigma2[i])+log(w)[i]
  }

  logsumexp <- log(apply(exp(logRho),1,sum))
  if(any(is.infinite(logsumexp))){
    logsumexp <- log(apply(exp(logRho),1,sum)+1e-8)
  }

  liklih <- mean(logsumexp)
  logR = sweep(logRho,1,logsumexp,FUN = "-")
  Gamma = exp(logR)

  res=list(Gamma=Gamma,liklih=liklih)
  res
}
