#' Computing log-transformed 1-d Gaussian pdf
#'
#' This function is used to compute log-transformed Gaussian probability density function
#'
#' @author Jiaying Zhao \email{jyzhao@connect.hku.hk}
#' @export
#' @param x sample
#' @param mu mean of the Gaussian distribution
#' @param Sigma2 variance the Gaussian distribution
#' @return 1 numeric vector pdf of sample length
#' @examples
#' loggaussian(x=0,mu=0,Sigma2=1)
#'
#'
loggaussian <- function(x,mu,Sigma2){
  -(1/2)*log(2*pi*Sigma2)-(1/(2*Sigma2))*((x-mu)*(x-mu))
}
