#' @export
#' @name supvar
#' @title  Supplementary variables projection
#' @description \code{supvar} is a generic function for adding supplementary variables
#' in \code{PCAmix} or \code{MFAmix}. The function invokes invokes two methods which depend on
#' the class of the first argument.
#' @param obj an object of class \code{PCAmix} or \code{MFAmix}.
#' @param \dots further arguments passed to or from other methods.
#' @details This generic function has two methods \code{\link{supvar.PCAmix}} and
#' \code{\link{supvar.MFAmix}}

supvar <- function(obj,...) {
  UseMethod("supvar")
}

#' @export
#' @name supvar.PCAmix
#' @method supvar PCAmix
#' @title  Supplementary variables in PCAmix
#' @description Performs the coordinates of supplementary variables on the
#'   component of an object of class \code{PCAmix}.
#' @param obj an object of class \code{PCAmix}.
#' @param X.quanti.sup a numeric matrix of data.
#' @param X.quali.sup a categorical matrix of data.
#' @param rename.level boolean, if TRUE all the levels of the qualitative
#'   variables are renamed as follows: "variable_name=level_name". This
#'   prevents to have identical names of the levels.
#' @param \dots further arguments passed to or from other methods.
#' @seealso \code{\link{PCAmix}}
#' @examples
#' data(wine)
#' X.quanti <- splitmix(wine)$X.quanti[,1:5]
#' X.quali <- splitmix(wine)$X.quali[,1,drop=FALSE]
#' X.quanti.sup <-splitmix(wine)$X.quanti[,28:29]
#' X.quali.sup <-splitmix(wine)$X.quali[,2,drop=FALSE]
#' pca<-PCAmix(X.quanti,X.quali,ndim=4,graph=FALSE)
#' pcasup <- supvar(pca,X.quanti.sup,X.quali.sup)


supvar.PCAmix <- function(obj,X.quanti.sup=NULL,X.quali.sup=NULL,rename.level=FALSE,...)
{
  if (!inherits(obj, "PCAmix")) 
    stop("use only with \"PCAmix\" objects")
  pca <- obj
  rec.sup <- recod(X.quanti.sup, X.quali.sup,rename.level)
  W.sup <- rec.sup$W # standardized quantitative+centered dummy variables
  n <- rec.sup$n
  if (n!=pca$rec$n) stop("number of rows in X.quanti.sup or X.quali.sup is not correct",call.=FALSE)
  p1 <- rec.sup$p1
  p2 <- rec.sup$p - p1
  m <- ncol(W.sup) - p1
  U <- pca$scores.stand
  N <- rep(1/n, n)
  A <- t(W.sup)%*%diag(N)%*%U
  rownames(A) <- colnames(W.sup)
  colnames(A) <-  paste("dim", 1:ncol(A), sep = "")
  nor <- apply(W.sup,2,function(x){sqrt(sum(x^2)/n)})
  B <- sweep(A,1,STATS = nor,FUN = "/")^2 #cos2 !!!
  quanti.sup <- levels.sup <- NULL
  sqload.sup <- matrix(NA,rec.sup$p,ncol(A))
  rownames(sqload.sup) <- colnames(rec.sup$X)
  if (p1!=0)
  {
    quanti.sup <- list(coord=NULL,cos2=NULL)
    quanti.sup$coord <- A[1:p1,,drop=FALSE]
    quanti.sup$cos2 <- B[1:p1,,drop=FALSE]
    sqload.sup[1:p1,] <- A[1:p1,,drop=FALSE]^2
  }
  if (p2!=0)
  {
    levels.sup <- list(coord=NULL,cos2=NULL)
    ns <- apply(rec.sup$G, 2, sum)
    levels.sup$coord <- sweep(A[(p1 + 1):(p1 + m),],1,STATS=ns/n,FUN="/")
    levels.sup$cos2 <- B[(p1 + 1):(p1 + m),]
    C <-  sweep(A[(p1 + 1):(p1 + m),]^2,1,STATS=n/ns,FUN="*")
    for (j in 1:ncol(X.quali.sup))
    {
      #levelj <- levels(as.factor(X.quali.sup[,j]))
      levelj <- which(rec.sup$index==(j+p1)) -p1
      sqload.sup[(p1+j),] <- apply(C[levelj,,drop=FALSE],2,sum)
    }
  }
  if (!is.null(quanti.sup))
    pca$quanti.sup <- quanti.sup
  if (!is.null(levels.sup))
    pca$levels.sup <- levels.sup
  pca$sqload.sup <- sqload.sup
  pca$rec.sup <- rec.sup
  return(pca)
}