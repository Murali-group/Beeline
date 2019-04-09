function(x1,x2,alpha=0.05){
  
  if(alpha<0||alpha>1) stop("alpha should be between 0 and 1")
  
  sorted <- sort(c(x1,x2))
  binEdges <- c(c(min(sorted)-1,sorted),max(sorted)+1)
  binCounts1 <- c(hist(x1,breaks = binEdges, right = FALSE, plot = FALSE)$counts,0)
  binCounts2 <- c(hist(x2,breaks = binEdges, right = FALSE, plot = FALSE)$counts,0)
  sumCounts1 <- cumsum(binCounts1)/sum(binCounts1)
  sumCounts2 <- cumsum(binCounts2)/sum(binCounts2)
  sampleCDF1 <- sumCounts1[1:(length(sumCounts1)-1)]
  sampleCDF2 <- sumCounts2[1:(length(sumCounts2)-1)]
  N1 <- length(x1)
  N2 <- length(x2)
  N <- N1+N2
  
  #Compute test statistic of interest
  CMstatistic <- N1*N2/N^2*sum((sampleCDF1 - sampleCDF2)^2)
  
  #table of the limiting distribution
  z=c(0.00000,0.02480,0.02878,0.03177,0.03430,0.03656,0.03865,0.04061,0.04247,0.04427,0.04601,0.04772,0.04939,0.05103,0.05265,0.05426,0.05586,0.05746,0.05904,0.06063,0.06222,0.06381,0.06541,0.06702,0.06863,
    0.07025,0.07189,0.07354,0.07521,0.07690,0.07860,0.08032,0.08206,0.08383,0.08562,0.08744,0.08928,0.09115,0.09306,0.09499,0.09696,0.09896,0.10100,0.10308,0.10520,0.10736,0.10956,0.11182,0.11412,0.11647,
    0.11888,0.12134,0.12387,0.12646,0.12911,0.13183,0.13463,0.13751,0.14046,0.14350,0.14663,0.14986,0.15319,0.15663,0.16018,0.16385,0.16765,0.17159,0.17568,0.17992,0.18433,0.18892,0.19371,0.19870,0.20392,
    0.20939,0.21512,0.22114,0.22748,0.23417,0.24124,0.24874,0.25670,0.26520,0.27429,0.28406,0.29460,0.30603,0.31849,0.33217,0.34730,0.36421,0.38331,0.40520,0.43077,0.46136,0.49929,0.54885,0.61981,0.74346,
    1.16786)
  Pz=c(seq(0,0.99,0.01),0.999)
  
  #compute parameters of the statistic's distribution
  T_mean <- 1/6+1/6/N;
  T_var <- 1/45*(N+1)/N^2 * ( 4*N1*N2*N-3*(N1^2+N2^2)-2*N1*N2 ) / (4*N1*N2)
  #translate the T statistic into the limiting distribution
  CM_limiting_stat <- ( CMstatistic - T_mean ) / sqrt(45*T_var) + 1/6; 
  #interpolate
  if(CM_limiting_stat>z[length(z)]) pValue <- 1
  else if(CM_limiting_stat<z[1]) pValue <- 0
  else pValue <- approx(z,Pz,CM_limiting_stat)$y
  
  H <- 1*(alpha>(1-pValue))
  
  return(list(H=H,CM_limiting_stat=CM_limiting_stat,pValue=pValue))
}