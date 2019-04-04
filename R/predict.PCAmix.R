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
#' train <- sample(1:n,20)
#' pca <- PCAmix(decathlon[train,1:10], graph=FALSE)
#' predict(pca, decathlon[-train,1:10])
#' rot <- PCArot(pca,dim=4)
#' predict(rot,decathlon[-train,1:10])
#' 
#' # qualitative data
#' data(flower)
#' n <- nrow(flower)
#' train <- sample(1:n,10)
#' mca <- PCAmix(X.quali=flower[train,1:3], rename.level=TRUE, graph=FALSE)
#' predict(mca, X.quali=flower[-train,1:3])
#' 
#' # quantitative and qualitative data
#' data(wine)
#' X.quanti <- splitmix(wine)$X.quanti
#' X.quali <- splitmix(wine)$X.quali
#' n <- nrow(wine)
#' train <- sample(1:n, 10)
#' pca <-PCAmix(X.quanti[train,1:10], X.quali[train,], ndim=4)
#' pred <- predict(pca, X.quanti[-train,1:10], X.quali[-train,])
#' plot(pca,axes=c(1,2))
#' points(pred[,c(1,2)],col=2,pch=16)
#' text(pred[,c(1,2)], labels = rownames(X.quanti[-train,1:27]), col=2,pos=3)


predict.PCAmix <- function(object, X.quanti = NULL, X.quali = NULL,...)
{
  pca <- object
  if (!inherits(pca, "PCAmix")) 
    stop("use only with \"PCAmix\" objects")

  beta <- pca$coef
  p2 <- 0
  
  train.rec <- pca$rec
  train.rename.level <-  pca$rename.level
  
  if (!is.null(X.quanti)) 
  {
    label <- rownames(X.quanti)
    n1 <- nrow(X.quanti)
    if (is.null(train.rec$X.quanti))
      stop("No quantitative dataset for training PCAmix", call. = FALSE)       
    if (!setequal(colnames(train.rec$X.quanti), colnames(X.quanti)))
      stop("The names of the columns in X.quanti and in the learning dataset
           are different", call. = FALSE)
    Y1 <- as.matrix(X.quanti[,colnames(train.rec$X.quanti)])
    # imput missing values with mean values in the train dataset
    if (any(is.na(Y1))) 
    {
      for (v in 1:ncol(Y1))
      {
        ind <- which(is.na(Y1[,v])==TRUE)
        if(length(ind)>=1)
          Y1[ind,v] <- train.rec$g[v]
      }
    }
  } else Y1 <- NULL
  
  if (!is.null(X.quali))
  {
    label <- rownames(X.quali)
    n2 <- nrow(X.quali)
    p2 <- ncol(X.quali)
    if (is.null(train.rec$X.quali))
      stop("No qualitative dataset for training PCAmix", call. = FALSE)  
    if (!setequal(colnames(train.rec$X.quali), colnames(X.quali)))
      stop("The names of the columns in X.quali and in the learning dataset
           are different", call. = FALSE)
    GNA <- PCAmixdata::tab.disjonctif.NA(X.quali, 
                                         rename.level = train.rename.level)
    G <- as.matrix(replace(GNA,is.na(GNA),0))

    if (!setequal(colnames(train.rec$G),colnames(G)))
    {
      #levels in train not in test : levels are merged in test
      if (length(setdiff(colnames(train.rec$G), colnames(G))) > 0)
      {
        for (v in 1:p2)
          levels(X.quali[,v]) <- union(levels(X.quali[,v]),
                                              levels(train.rec$X.quali[,v]))
        GNA <- PCAmixdata::tab.disjonctif.NA(X.quali, 
                                             rename.level = train.rename.level)
        G <- as.matrix(replace(GNA,is.na(GNA),0))
      }
      
      #levels in test not in train : observations are removed
      if (length(setdiff(colnames(G), colnames(train.rec$G))) > 0)
      { 
        if (length(setdiff(colnames(G), colnames(train.rec$G))) == length(colnames(G)))
          stop("No level in common between the test and the training dataset",
               call. = FALSE)
        G2 <- G[,is.element(colnames(G),colnames(train.rec$G)), drop=FALSE]
        if (length(which(apply(G2,1,sum) < p2)) == nrow(G))
          stop("No observation in the test dataset can be predicted", 
               call. = FALSE)
        G <- G2
        if (any(apply(G2,1,sum) < p2))
        {
          warning("Predictions can not be preformed for some observations")
          G <- G2[-which(apply(G2,1,sum) < p2),, drop=FALSE]
          if (!is.null(Y1))
            Y1 <- Y1[-which(apply(G2,1,sum) < p2),,drop=FALSE]
        }
      }
    }
  G <- G[ , colnames(train.rec$G), drop = FALSE]
  } else G <- NULL
  
  if (!is.null(X.quanti) && !is.null(X.quali))
  {
    if (n1 != n2) 
      stop("The number of objects in X.quanti and X.quali must be 
                       the same", call. = FALSE)
    if (sum(rownames(X.quali)!=rownames(X.quanti))!=0) 
      stop("The names/order of the objects in X.quanti and X.quali must be the same", 
           call. = FALSE)
  }
  
  Y <- cbind(Y1,G)
  n <- nrow(Y)
  
  coord <- matrix(NA, n,length(beta))
  for (g in 1: length(beta)) 
    coord[,g] <- Y %*% beta[[g]][-1] +  beta[[g]][1]
  
  if (colnames(pca$sqload)[1]=="dim1.rot")  
    colnames(coord) <- paste("dim", 1:length(beta), sep = "",".rot")
  else
    colnames(coord) <- paste("dim", 1:length(beta), sep = "")
  rownames(coord) <- rownames(Y)
  return(coord)			
}
