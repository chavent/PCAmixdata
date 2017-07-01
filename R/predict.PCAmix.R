#' @export predict.MFAmix
#' @export
#' @name predict.PCAmix
#' @title  Prediction of new scores in PCAmix or PCArot
#' @description This function performs the scores of new observations 
#' on the principal components of PCAmix. If the components have been rotated,
#' this function performs the scores of the new observations on the 
#' rotated principal components. In other words, this function is projecting 
#' the new observations onto the principal components of PCAmix (or PCArot) 
#' obtained previoulsy on a separated dataset. Note that the new observations 
#' must be described with the same variables than those used in PCAmix (or PCArot).
#' @param object an object of class PCAmix obtained with the function 
#' \code{PCAmix} or \code{PCArot}.
#' @param X.quanti a numeric data matrix or an object that can be coerced to such 
#' a matrix (such as a numeric vector or a data frame with all numeric columns).
#' @param X.quali a categorical matrix of data, or an object that can be coerced 
#' to such a matrix (such as a character vector, a factor or a data frame with 
#' all factor columns).
#' @param rename.level boolean, if TRUE all the levels of the qualitative variables
#' are renamed as follows: "variable_name=level_name". This prevents to have 
#' identical names for the levels.
#' @param \ldots urther arguments passed to or from other methods. 
#' They are ignored in this function.
#' @return  Returns the matrix of the scores of the new observations on 
#' the principal components or on the rotated principal components of PCAmix.
#' @seealso \code{\link{PCAmix}},\code{\link{PCArot}}
#' @author Marie Chavent \email{marie.chavent@u-bordeaux.fr}, Amaury Labenne.
#' @references 
#' Chavent M., Kuentz-Simonet V., Labenne A., Saracco J., Multivariate analysis of mixed data: The PCAmixdata R package, arXiv:1411.4911 [stat.CO].
#' @examples 
#' # quantitative data
#' data(decathlon)
#' n <- nrow(decathlon)
#' sub <- sample(1:n,20)
#' pca<-PCAmix(decathlon[sub,1:10], graph=FALSE)
#' predict(pca,decathlon[-sub,1:10])
#' rot <- PCArot(pca,dim=4)
#' predict(rot,decathlon[-sub,1:10])
#' 
#' # quantitative and qualitative data
#' data(wine)
#' str(wine)
#' X.quanti <- splitmix(wine)$X.quanti
#' X.quali <- splitmix(wine)$X.quali
#' pca<-PCAmix(X.quanti[,1:27],X.quali,ndim=4,graph=FALSE)
#' n <- nrow(wine)
#' sub <- sample(1:n,10)
#' pca<-PCAmix(X.quanti[sub,1:27],X.quali[sub,],ndim=4)
#' pred <- predict(pca,X.quanti[-sub,1:27],X.quali[-sub,])
#' plot(pca,axes=c(1,2))
#' points(pred[,c(1,2)],col=2,pch=16)
#' text(pred[,c(1,2)], labels = rownames(X.quanti[-sub,1:27]), col=2,pos=3)


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
