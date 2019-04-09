####################
# PACKAGE required:
# pracma
####################

function(acc_ranked,adj,n_points){
  maxWEIGHT <- max(acc_ranked)
  acc_ranked <- acc_ranked/maxWEIGHT
  adj <- adj*3
  
  fp <- matrix(0,nrow = 1,ncol = n_points+1)
  tn <- matrix(0,nrow = 1,ncol = n_points+1)
  tp <- matrix(0,nrow = 1,ncol = n_points+1)
  fn <- matrix(0,nrow = 1,ncol = n_points+1)
  fpr <- matrix(0,nrow = 1,ncol = n_points+1)
  precision <- matrix(0,nrow = 1,ncol = n_points+1)
  recall <- matrix(0,nrow = 1,ncol = n_points+1)
  
  thr <- seq(0,1,by=1/n_points)
  sign_info <- sign(acc_ranked)
  sign_info[sign_info==0] <- 1
  for (i in 1:(n_points+1)) {
    adj_m <- (abs(acc_ranked)>=thr[i])*sign_info
    compare <- abs(adj_m+adj)
    fp[1,i] <- nnzero(compare<3&compare>0)
    tp[1,i] <- nnzero(compare==4)
    fn[1,i] <- nnzero(compare==3)
    tn[1,i] <- nnzero(compare==0)
    precision[1,i]=tp[1,i]/(tp[1,i]+fp[1,i]);
    recall[1,i]=tp[1,i]/(tp[1,i]+fn[1,i]);
    fpr[1,i]=fp[1,i]/(fp[1,i]+tn[1,i]);
  }
  precision <- rev(precision)
  recall <- rev(recall)
  precision_new <- c(1,precision)
  recall_new <- c(0,recall)
  fpr <- rev(fpr)
  
  auc <- dget("SINCERITIES functions/auc.R")
  AUROC <- auc(fpr,recall)
  AUPR <- auc(recall_new,precision_new)
  
  result <- list(AUROC=AUROC,AUPR=AUPR,fpr=fpr,recall=recall,recall_new=recall_new,precision_new=precision_new)
  return(result)
}