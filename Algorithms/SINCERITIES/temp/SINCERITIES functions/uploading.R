function(filename=NULL,header=TRUE){
  
  if(is.null(filename)){
    filename <- file.choose()
  }
  
  df <- read.csv(filename,header=header)
  NAIdx <- is.na(apply(df,1,sum))
  df <- df[!NAIdx,]
  
  totDATA <- df[,1:dim(df)[2]-1]
  timeline <- df[,dim(df)[2]]
  DATA.time <- sort(unique(timeline))
  DATA.num_time_points <- length(DATA.time)
  DATA.totDATA <- matrix(ncol = dim(df)[2]-1)
  DATA.timeline <- vector()
  
  for (k in 1:DATA.num_time_points) {
    I <- which(timeline==DATA.time[k])
    DATA.totDATA <- rbind(DATA.totDATA,as.matrix(totDATA[I,]))
    DATA.timeline <- c(DATA.timeline,timeline[I])
  }
  DATA.totDATA <- DATA.totDATA[-1,]
  DATA.totDATA[is.na(DATA.totDATA)] <- 0
  DATA.numGENES <- dim(DATA.totDATA)[2]
  
  if(!header){
    DATA.genes <- vector(length = DATA.numGENES)
    for (i in 1:DATA.numGENES) {
      DATA.genes[i] <- sprintf("Gene %d",i)
    }
  }else{
    DATA.genes <- colnames(df)[1:dim(df)[2]-1]
  }
  
  DATA.singleCELLdata <- by(DATA.totDATA,DATA.timeline,identity)
  DATA <- list(time=DATA.time,num_time_points=DATA.num_time_points,
               totDATA=DATA.totDATA,timeline=DATA.timeline,numGENES=DATA.numGENES,
               genes=DATA.genes,singleCELLdata=DATA.singleCELLdata)
  return(DATA)
}