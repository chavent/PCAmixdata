#' @export
splitmix<-function(data){
  data<-data.frame(data,check.names=T)
  data_desc <- sapply(data,class)
 
  X.quanti <- X.quali <- NULL
  col.quant <- c()
  col.qual <- c()

  for(i in 1:length(data_desc)){
    xdf <- as.data.frame(data_desc[i])
    if(xdf[1,1] == "factor" | xdf[1,1] == "character" | xdf[1,1] == "ordered"){
      col.qual <- c(col.qual, i)
    }
    else if(xdf[1,1] == "numeric" | xdf[1,1] =="integer"){
      col.quant <- c(col.quant, i)
    }
    else{
      stop("Undefined data type")
    }
  }
 
  if (length(col.quant)!=0) X.quanti<-data[,col.quant,drop=FALSE]
  if (length(col.qual)!=0) X.quali<-data[,col.qual,drop=FALSE]

  return(list(X.quanti=X.quanti,X.quali=X.quali,col.quant=col.quant,col.qual=col.qual))
}