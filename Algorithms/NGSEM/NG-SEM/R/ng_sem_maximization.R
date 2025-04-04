#' Maxinization step
#'
#' This function is used for the M step in the optimization
#'
#' @author Jiaying Zhao, \email{jyzhao@connect.hku.hk}
#' @export
#' @param e.step results of E-step
#' @param Data a list including Y, y
#' \itemize{
#'  \item{Y} designed matrix
#'  \item{y} dependent vector
#'  }
#' @param residual residual for SEM
#' @param model
#' \itemize{
#'  \item{mu} a vector for the mean of the Gaussian Mixture components
#'  \item{Sigma2} a vector for the variance of the Gaussian Mixture components
#'  \item{w} a vector for the weights of the Gaussian Mixture components
#'  }
#' @return model
#' \itemize{
#'  \item{mu} a vector for the mean of the Gaussian Mixture components
#'  \item{Sigma2} a vector for the variance of the Gaussian Mixture components
#'  \item{w} a vector for the weights of the Gaussian Mixture components
#'  \item{a} adjacent weights for GRN
#'  \item{n.K} number of components for the Gaussian Mixture
#'  }

ng_sem_maximization <- function(e.step,Data,residual,model){
  Data.Y <- Data$Y
  Data.y <- Data$y
  res_y_ax <- residual

  Gamma <- e.step$Gamma
  n.N <- nrow(Gamma)
  n.K <- ncol(Gamma)
  Nk <- apply(Gamma, 2, sum)
  w <- Nk/sum(Nk)

  mu <- rep(0,n.K)

  Sigma2 <- Gamma*matrix(unlist(lapply(mu,function(x){
    (res_y_ax-x)*(res_y_ax-x)
  })),ncol = n.K,byrow=F)

  Sigma2 <- sweep(Sigma2,2,1/Nk,FUN="*")
  Sigma2 <- colSums(Sigma2)+1e-6

  t <- sweep(Gamma,2,1/(2*Sigma2),FUN="*")
  wls.w <- rowSums(t)
  wls.w <- diag(wls.w)

  a <- MASS::ginv(Data.Y%*%wls.w%*%t(Data.Y)) %*% (Data.Y%*%wls.w%*%Data.y)
  model <- list(mu=mu,Sigma2=Sigma2,w=w,a=a,n.K=n.K)
  model
}



