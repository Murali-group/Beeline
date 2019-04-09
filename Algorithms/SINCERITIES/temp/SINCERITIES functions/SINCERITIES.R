###############################################################################
# SINCERITIES is a novel computational method for inferring
# gene regulatory network (GRN) from time-stamped cross-sectional
# single-cell expression data.
#
# result=SINCERITIES(DATA) 
# 
# DATA, a list containing the following information:
# - singleCELLdata: list of length n, where n is the number of
# capture time points. DATA$singleCELLdata[[k]] is a s_k by m matrix/dataframe
# containing observed expression levels of m genes in s_k single cells.
# - DATA.totDATA: S by m matrix, where S is the total number of single
# cells (i.e., S=s_1+s_2+...+s_n where n the number of capture time points) 
# and m is the number of genes.
# - DATA.time: vector of length n containing the cell capture time points or 
# time-stamps).
#
# result=SINCERITIES(DATA,distance) 
#
# distance: this parameter selects the distribution distance
# 1- for KS (Kolmogorov-Smirnov)  (* DEFAULT *)
# 2- for CM (Cramer-von Mises)
# 3- for AD (Anderson-Darling)
# 4- for Mean Expression Difference
#
# result=SINCERITIES(DATA,distance,method) 
#
# method: this parameter selects the regularization regression strategy
# 1- for RIDGE  (* DEFAULT *)
# 2- for ELASTIC-NET with automatic detection of optimal alpha parameter
# 3- for LASSO
# 4- for ELASTIC-NET with manual selection of alpha parameter
#
# result=SINCERITIES(DATA,distance,method,noDIAG)
#
# noDIAG: this parameter selects whether the auto-regulatory edge is
# inferred
# 0- GRN contains no auto-regulatory edge (* DEFAULT *)
# 1- GRN contain auto-regulatory edge
#
#
# result=SINCERITIES(DATA,distance,method,noDIAG,SIGN)
# SIGN: this parameter selects whether the sign / mode of the gene
# regulations is inferred
# 0- for unsigned GRN
# 1- for signed GRN (* DEFAULT *)
# SINCERITIES uses partial correlation analysis where a positive (negative) 
# correlation is taken as an indication of activation (repression). 
#
# OUTPUTS:
#
# result, a list containing the following information:
# -adj_matrix: m by m matrix containing the weights of regulatory edges.
# The larger adj_matrix(i,j) indicates higher confidence that the
# corresponding edge exists (i.e., gene i regulating gene j).
# -DISTANCE_matrix: n-1 by m matrix containing the (normalized)
# distribution distance (DD) computed during the network inference.
#
# Created by Nan Papili Gao and R version implemented by Ziyi Hua
# Institute for Chemical and Bioengineering 
# ETH Zurich
# E-mail:  nanp@ethz.ch
#
# Copyright. November 1, 2016.
#
# PACKAGES required:
# kSamples
# glmnet
# ppcor
###############################################################################

function(DATA,distance=1,method=1,noDIAG=0,SIGN=1){
  
  if(!distance%in%c(1,2,3,4)){
    stop("Choose distance metric with 1,2 and 3: 1 for Kolmogorov-Smirnov(KM), 2 for Cramer-von Mises(CM), 
         3 for Anderson-Darling(AD)")
  }
  
  if(!method%in%c(1,2,3,4)){
    stop("Choose regularization regression strategy with 1,2,3 and 4: 1 for RIGDE, 2 for ELASTIC NET with automatic
         detection of optimal alpha parameter, 3 for LASSO, 4 for ELASTIC NET with manual selection of alpha parameter")
  }
  
  
  if(!noDIAG%in%c(0,1)){
    stop("noDIAG should be either 0 or 1")
  }
  
  if(!SIGN%in%c(0,1)){
    stop("SIGN should be either 0 or 1")
  }
  
  #Initialization
  single_cell_data <- DATA$singleCELLdata
  time <- DATA$time
  numGENES <- DATA$numGENES
  num_time_points <- length(time)
  
  #Catch error
  if(num_time_points<5){
    stop('** DATA with a number of time points < 5. Please run SINCERITIES_CROSS_VALIDATION function **')
  }
  
  #Distribution Distance
  DISTANCE_matrix <- matrix(data=0,nrow=num_time_points-1,ncol = numGENES)
  #h <- matrix(data=0,nrow=num_time_points-1,ncol = numGENES)
  totalDATA <- single_cell_data[[1]]
  
  cmtest2 <- dget("SINCERITIES functions/cmtest2.R")
  
  for (ti in 1:(num_time_points-1)) {
    totalDATA <- rbind(totalDATA,single_cell_data[[ti+1]])
    data_ti <- t(single_cell_data[[ti]])
    data_ti_plus1 <- t(single_cell_data[[ti+1]])
    
    for (gi in 1:numGENES) {
      p1 <- data_ti[gi,]
      p2 <- data_ti_plus1[gi,]
      if(distance==1){
        test.stat <- ks.test(p1,p2)
        DISTANCE_matrix[ti,gi] <- test.stat$statistic 
      }else if(distance==2){
        DISTANCE_matrix[ti,gi] <- cmtest2(p1,p2)$CM_limiting_stat
      }else if(distance==3){
        test.stat <- ad.test(p1,p2)
        DISTANCE_matrix[ti,gi] <- test.stat$ad[2,1]
      }else if(distance==4){
        DISTANCE_matrix[ti,gi] <- abs(mean(p1)-mean(p2))
      }
    }
  }
  
  #normalization
  deltaT <- replicate(dim(DISTANCE_matrix)[2],time[2:length(time)]-time[1:(length(time)-1)])
  DISTANCE_matrix_normed <- DISTANCE_matrix/deltaT
  
  #Generate Y and X_matrix for glmnet
  if(method==1){
    alphas <- 0
  }else if(method==2){
    alphas <- seq(0,1,0.1)
  }else if(method==3){
    alphas <- 1
  }else{
    input <- readline(' *** Please input manually the alpha values (between 0 and 1) separated by comma: ')
    alphas <- as.numeric(unlist(strsplit(input,',')))
  }
  DISTANCE_matrix <- DISTANCE_matrix_normed
  X_matrix <- DISTANCE_matrix[1:(num_time_points-2),]
  
  #LOOCV settings
  nfold <- dim(X_matrix)[1]
  foldid <- 1:nfold
  keep <- TRUE
  pred_lambda_min <- matrix(0, nrow = numGENES, ncol = numGENES)
  
  lambda_res <- vector()
  alpha_res <- vector()
  
  for (gi in 1:numGENES) {
    
    lambda <-  vector()
    cvERROR <-  vector()
    beta <- matrix(data=0,nrow = dim(X_matrix)[2],ncol = length(alphas))
    
    for (test in 1:length(alphas)) {
      Y_vector <- DISTANCE_matrix[2:(num_time_points-1),gi]
      if(noDIAG==1){
        CV_results <- cv.glmnet(X_matrix,Y_vector,alpha=alphas[test],exclude=gi,nfolds = nfold, foldid = foldid,
                                keep = keep, lower.limits=0, upper.limits=Inf, grouped = FALSE)
      }else{
        CV_results <- cv.glmnet(X_matrix,Y_vector,alpha=alphas[test],nfolds = nfold, foldid = foldid,
                                keep = keep, lower.limits=0, upper.limits=Inf, grouped = FALSE)
      }
      lambda[test] <- CV_results$lambda.min
      cvERROR[test] <- CV_results$cvm[CV_results$lambda==CV_results$lambda.min]
      coef.CV_results <- coef.cv.glmnet(CV_results,s='lambda.min')
      beta[coef.CV_results@i[-1],test] = coef.CV_results@x[-1]
    }
    
    minIdx <- max(which(cvERROR==min(cvERROR)))
    lambda_res[gi] <- lambda[minIdx]
    alpha_res[gi] <- alphas[minIdx]
    pred_lambda_min[,gi] <- beta[,minIdx]
    
  }
  
  if(SIGN==1){
    parcorr_matrix <- pcor(DATA$totDATA,method = 'spearman')$estimate
    pred_lambda_min <- pred_lambda_min*sign(parcorr_matrix)
  }
  
  result <- list(DISTANCE_matrix=DISTANCE_matrix,adj_matrix=pred_lambda_min)
  return(result)
}