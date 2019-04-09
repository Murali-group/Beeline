#SIGN can be 0 or 1

function(connectivityMATRIX,genes,SIGN=0,fileNAME=NULL,saveFile=FALSE){
  
 
  if(is.null(fileNAME)){
    fileNAME <- 'GRNprediction.txt'
  }
  
  numGENES <- length(genes)
  interactions <- as.vector(connectivityMATRIX)
  
  edges <- vector(length = length(interactions))
  if(SIGN==1){
    edges[which(interactions<0)] <- "repression"
    edges[which(interactions>0)] <- "activation"
    edges[which(interactions==0)] <- "no regulation"
  }else{
    edges[which(interactions==0)] <- "no regulation"
    edges[which(interactions!=0)] <- "activation/repression"
  }
  
  interactions <- abs(interactions)
  
  targetGENES <- as.vector(replicate(numGENES,genes))
  sourceGENES <- as.vector(t(replicate(numGENES,genes)))
  
  df <- data.frame(sourceGENES,targetGENES,interactions,edges)
  colnames(df) <- c('SourceGENES','TargetGENES','Interaction','Edges')
  df <- df[order(-df$Interaction),]
  row.names(df) <- 1:dim(df)[1]
  if(saveFile) write.csv(df,file = fileNAME,row.names = FALSE,quote = FALSE)
  
  return(df)
}
