####################
# PACKAGE required:
# pracma
####################

function(x,y){
  if(length(x)!=length(y)){
    stop("The input vectors should have the same size")
  }
  n <- length(x)
  xix <- order(x)
  x_ <- x[xix]
  y_ <- y[xix]
  NAidx <- is.na(x_)|is.na(y_)
  x_ <- x_[!NAidx]
  y_ <- y_[!NAidx]
  area <- trapz(x_,y_)
  return(area)
}