##' @export
##' @method predict PCAmix
predict.PCAmix <- function(object, X.quanti=NULL,X.quali=NULL,rename.level=FALSE,...)
{
  pca<-object
  if (!inherits(pca, "PCAmix")) 
    stop("use only with \"PCAmix\" objects")
#   data.new.split<-splitmix(data.new)
#   X.quanti<-data.new.split$X.quanti
#   X.quali<-data.new.split$X.quali
  
  if (pca$rename.level) 
    rename.level=TRUE
  if ((rename.level) & (!pca$rename.level))
    stop("perform PCAmix with argument rename.level=TRUE",call.=FALSE)
  
  rec <- recod(X.quanti,X.quali,rename.level=rename.level)
  Y <- rec$Y
  n <- rec$n
  beta <- pca$coef
  
  if (!is.null(X.quanti)) 
  {
    label <- rownames(X.quanti)
    n1 <- nrow(X.quanti)
    p1 <- ncol(X.quanti)
    if (p1 != pca$rec$p1) stop("The number of variables in X.quanti must be the same than in the learning set",call. = FALSE)
  }
  if (!is.null(X.quali))
  {
    label <- rownames(X.quali)
    n2 <- nrow(X.quali)
    p2 <- ncol(X.quali)
    if (p2 != pca$rec$p2) 
      stop("The number of variables in X.quali must be the same than in the learning set",call. = FALSE)
  }
  if (!is.null(X.quanti) && !is.null(X.quali))
  {
    if (n1 != n2) stop("The number of objects in X.quanti and X.quali must be the same",call. = FALSE)
    if (sum(rownames(X.quali)!=rownames(X.quanti))!=0) stop("The names of the objects in X.quanti and X.quali must be the same",call. = FALSE)
  }
  
  
  if (!setequal(colnames(pca$rec$X),colnames(rec$X)))
    stop("The colnames in the new sample are not appropiate",call. = FALSE)
  
  if (rec$p2 >0 )
  {
    #different levels in the new sample and in the learning sample
    if (!setequal(colnames(pca$rec$Y),colnames(Y)))
    {
      col.test <- is.element(colnames(Y),colnames(pca$rec$Y))
      Y <- Y[,col.test,drop=FALSE]
      col.beta <- c(TRUE,is.element(colnames(pca$rec$Y),colnames(Y)))
      beta <- lapply(beta,function(x) {x[col.beta,1,drop=FALSE]})
    }
  }
  
  coord <- matrix(,n,length(beta))
  for (g in 1: length(beta)) coord[,g] <-Y %*% beta[[g]][-1] +  beta[[g]][1]
  
  if (colnames(pca$sqload)[1]=="dim1.rot")  
    colnames(coord) <- paste("dim", 1:length(beta), sep = "",".rot")
  else
    colnames(coord) <- paste("dim", 1:length(beta), sep = "")
  rownames(coord) <- label
  return(coord)			
}
