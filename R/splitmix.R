#' @export
splitmix<-function(data){
  data<-data.frame(data,check.names=T)
  class.col <- unlist(lapply(data,class))
  col.quant <- which(class.col %in% c("numeric","integer"))
  col.qual <- which(class.col %in% c("factor","character"))
  if ("integer" %in% class.col) 
    warning("Columns of class integer are considered as quantitative")
  X.quanti <- X.quali <- NULL
  if (length(col.quant)!=0) X.quanti<-data[,col.quant,drop=FALSE]
  if (length(col.qual)!=0) X.quali<-data[,col.qual,drop=FALSE]
  
  return(list(X.quanti=X.quanti,X.quali=X.quali,col.quant=col.quant,col.qual=col.qual))
}