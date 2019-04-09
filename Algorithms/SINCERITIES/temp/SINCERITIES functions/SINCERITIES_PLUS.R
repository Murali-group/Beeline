###############################################################################
# SINCERITIES_PLUS is a novel computational method for inferring
# gene regulatory network (GRN) from time-stamped cross-sectional
# single-cell expression data.
#
# SINCERITIES_PLUS function is an extension of the default version of
# SINCERITIES, designed for single cell datasets with fewer than five time
# points.
# DEFAULT regularization regression strategy: RIDGE
# DEFAULT distribution distance: KS (Kolmogorov-Smirnov) 
#
# result=SINCERITIES_PLUS(DATA) 
# 
# DATA, a list containing the following information:
# - DATA$singleCELLdata: list of length n, where n is the number of
# capture time points. DATA$singleCELLdata[[k]] is a s_k by m matrix/dataframe
# containing observed expression levels of m genes in s_k single cells.
# - DATA$totDATA: S by m matrix, where S is the total number of single
# cells (i.e., S=s_1+s_2+...+s_n where n the number of capture time points) 
# and m is the number of genes.
# - DATA$time: vector of length n containing the cell capture time points or 
# time-stamps).
#
# result=SINCERITIES_PLUS(DATA,noDIAG)
#
# noDIAG: this parameter selects whether the auto-regulatory edge is
# inferred
# 0- GRN contains no auto-regulatory edge (* DEFAULT *)
# 1- GRN contain auto-regulatory edge
#
#
# result=SINCERITIES_PLUS(DATA,noDIAG,SIGN)
# SIGN: this parameter selects whether the sign / mode of the gene
# regulations is inferred
# 0- for unsigned GRN
# 1- for signed GRN (* DEFAULT *)
# SINCERITIES uses partial correlation analysis where a positive (negative) 
# correlation is taken as an indication of activation (repression). 
#
# result=SINCERITIES_PLUS(DATA,noDIAG,SIGN,CV_nfolds)
# CV_nfolds defines a partition of the data into CV_n_folds disjoint
# subsets for the cross validation
# CV_nfolds=5; (* DEFAULT *)
#
# OUTPUTS:
#
# result, a list containing the following information:
# -adj_matrix: m by m matrix containing the weights of regulatory edges.
# The larger adj_matrix[i,j] indicates higher confidence that the
# corresponding edge exists (i.e., gene i regulating gene j).
# -DISTANCE_matrix: n-1 by m matrix containing the (normalized)
# distribution distance (DD) computed during the network inference.
#
# Created by Nan Papili Gao and R version implemented by Ziyi Hua
# Institute for Chemical and Bioengineering 
# ETH Zurich
# E-mail:  nanp@ethz.ch
#
# Copyright. Apr 20, 2017
#
# PACKAGES required:
# cvTools
# glmnet
# ppcor
###############################################################################

function(DATA, noDIAG = 0, SIGN = 1, CV_nfolds = 5){
  
  #Initialization
  single_cell_data <- DATA$singleCELLdata
  time <- DATA$time
  numGENES <- DATA$numGENES
  num_time_points <- length(time)
  
  #Catch error
  if(num_time_points<3){
    stop('** The data must contain at least 3 time points **')
  }
  
  #K-fold CV
  library(cvTools)
  K <- CV_nfolds
  indices <- list()
  for(cv in 1:length(single_cell_data)){
    N <-  dim(single_cell_data[[cv]])[1]
    indices[[cv]] <- cvFolds(N,K = K)
  }
  
  error_CV <- array(0,dim = c(K,numGENES,100))
  
  for(cross in 1:K){
    data4test <- list()
    data4train <- list()
    for(cv in 1:length(single_cell_data)){
      test <- indices[[cv]]$subsets[indices[[cv]]$which==cross]
      train <- indices[[cv]]$subsets[indices[[cv]]$which!=cross]
      data4test[[cv]] <- single_cell_data[[cv]][test,]
      data4train[[cv]] <- single_cell_data[[cv]][train,]
    }
    
    #Distribution Distance
    DISTANCE_matrix_train <- matrix(data=0,nrow=num_time_points-1,ncol = numGENES)
    totalDATA <- data4train[[1]]
    
    for(ti in 1:(num_time_points-1)){
      totalDATA <- rbind(totalDATA,data4train[[ti+1]])
      data_ti <- t(data4train[[ti]])
      data_ti_plus1 <- t(data4train[[ti+1]])
      for(gi in 1:numGENES){
        p1 <- data_ti[gi,]
        p2 <- data_ti_plus1[gi,]
        test.stat <- ks.test(p1,p2)
        DISTANCE_matrix_train[ti,gi] <- test.stat$statistic 
      }
    }
    
    #Normalization
    deltaT <- replicate(dim(DISTANCE_matrix_train)[2],time[2:length(time)]-time[1:(length(time)-1)])
    DISTANCE_matrix_train_normed <- DISTANCE_matrix_train/deltaT
    X_matrix <- DISTANCE_matrix_train_normed[1:(num_time_points-2),]
    
    #Distribution Distance Test
    DISTANCE_matrix_test <- matrix(data=0,nrow=num_time_points-1,ncol = numGENES)
    totalDATA <- data4test[[1]]
    
    for(ti in 1:(num_time_points-1)){
      totalDATA <- rbind(totalDATA,data4test[[ti+1]])
      data_ti <- t(data4test[[ti]])
      data_ti_plus1 <- t(data4test[[ti+1]])
      for(gi in 1:numGENES){
        p1 <- data_ti[gi,]
        p2 <- data_ti_plus1[gi,]
        test.stat <- ks.test(p1,p2)
        DISTANCE_matrix_test[ti,gi] <- test.stat$statistic 
      }
    }
    
    #Normalization
    deltaT <- replicate(dim(DISTANCE_matrix_test)[2],time[2:length(time)]-time[1:(length(time)-1)])
    DISTANCE_matrix_test_normed <- DISTANCE_matrix_test/deltaT
    X_matrix_test <- DISTANCE_matrix_test_normed[1:(num_time_points-2),]
    
    #Generate Y and X_matrix for glmnet
    alphas <- 0 #Ridge Regression
    lambdas <- 10^seq(-2,2,length.out = 100)
    
    for(gi in 1:numGENES){
      Y_vector <- DISTANCE_matrix_train_normed[2:(num_time_points-1),gi]
      if(noDIAG==1){
        CV_results <- glmnet(X_matrix,Y_vector,alpha = alphas,exclude = gi,lambda = lambdas,
                             lower.limits = 0, upper.limits = Inf)
      }else{
        CV_results <- glmnet(X_matrix,Y_vector,alpha = alphas, lambda = lambdas,
                             lower.limits = 0, upper.limits = Inf)
      }
      Y_vector_test <- DISTANCE_matrix_test_normed[2:(num_time_points-1),gi]
      
      for(lambdacount in 1:length(CV_results$lambda)){
        beta_lambda <- as.matrix(CV_results$beta)[,lambdacount]
        error_CV[cross,gi,lambdacount] <- sum((Y_vector_test - X_matrix_test%*%beta_lambda)^2)
      }
    }
  }
  
  mean_error_CV <- apply(error_CV,c(2,3),mean)
  standard_error_mean <- apply(error_CV,c(2,3),sd)/sqrt(K)
  
  #Lambda min
  min_mean_error_CV <- apply(mean_error_CV,1,min)
  idx_lambda_min <- apply(mean_error_CV,1,which.min)
  lambda_min <- CV_results$lambda[idx_lambda_min]
  
  #Lambda 1SE
  idx_lambda_1SE <- vector(length = numGENES)
  for(gi in 1:numGENES){
    min_plu_1SE <- mean_error_CV[gi,idx_lambda_min[gi]]+standard_error_mean[gi,idx_lambda_min[gi]]
    idx_lambda_1SE[gi] <- which(mean_error_CV[gi,]<=min_plu_1SE)[1]
  }
  lambda_1SE <- CV_results$lambda[idx_lambda_1SE]
  
  SINCERITIES_final <- dget("SINCERITIES functions/SINCERITIES_final.R")
  pred_lambda <- SINCERITIES_final(single_cell_data,time,numGENES,num_time_points,lambdas,alphas,idx_lambda_min,noDIAG)
  
  if(SIGN==1){
    parcorr_matrix <- pcor(DATA$totDATA,method = 'spearman')$estimate
    pred_lambda <- pred_lambda*sign(parcorr_matrix)
  }
  
  result <- list(DISTANCE_matrix=DISTANCE_matrix_train_normed,adj_matrix=pred_lambda)
  return(result)
}