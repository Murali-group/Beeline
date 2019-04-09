function(single_cell_data,time,numGENES,num_time_points,lambdas,alphas,idx_lambda_min,noDIAG){
  #Distribution Distance
  DISTANCE_matrix <- matrix(data=0,nrow=num_time_points-1,ncol = numGENES)
  totalDATA <- single_cell_data[[1]]
  
  for(ti in 1:(num_time_points-1)){
    totalDATA <- rbind(totalDATA,single_cell_data[[ti+1]])
    data_ti <- t(single_cell_data[[ti]])
    data_ti_plus1 <- t(single_cell_data[[ti+1]])
    for(gi in 1:numGENES){
      p1 <- data_ti[gi,]
      p2 <- data_ti_plus1[gi,]
      test.stat <- ks.test(p1,p2)
      DISTANCE_matrix[ti,gi] <- test.stat$statistic 
    }
  }
  
  #Normalization
  deltaT <- replicate(dim(DISTANCE_matrix)[2],time[2:length(time)]-time[1:(length(time)-1)])
  DISTANCE_matrix_normed <- DISTANCE_matrix/deltaT
  X_matrix <- DISTANCE_matrix_normed[1:(num_time_points-2),]
  
  #glmnet
  pred_lambda_min <- matrix(0, nrow = numGENES, ncol = numGENES)
  for(gi in 1:numGENES){
    Y_vector <- DISTANCE_matrix_normed[2:(num_time_points-1),gi]
    if(noDIAG==1){
      CV_results <- glmnet(X_matrix,Y_vector,alpha = alphas,exclude = gi,lambda = lambdas,
                           lower.limits = 0, upper.limits = Inf)
    }else{
      CV_results <- glmnet(X_matrix,Y_vector,alpha = alphas, lambda = lambdas,
                           lower.limits = 0, upper.limits = Inf)
    }
    pred_lambda_min[,gi] <- as.matrix(CV_results$beta)[,idx_lambda_min[gi]]
  }
  
  return(pred_lambda_min)
}