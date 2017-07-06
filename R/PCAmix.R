#' @importFrom stats cor sd
NULL
#' @importFrom utils combn
NULL
#' @export
#' @name PCAmix
#' @title Principal component analysis of mixed data 
#' @description Performs principal component analysis of  a set of 
#' individuals (observations) described by a mixture of qualitative and 
#' quantitative variables. PCAmix includes ordinary principal component 
#' analysis (PCA) and multiple correspondence analysis (MCA) as special cases.
#' @param X.quanti a numeric matrix of data, or an object that can be 
#' coerced to such a matrix (such as a numeric vector or a data frame with 
#' all numeric columns).
#' @param X.quali a categorical matrix of data, or an object that can be 
#' coerced to such a matrix (such as a character vector, a factor or a data 
#' frame with all factor columns).
#' @param ndim number of dimensions kept in the results (by default 5).
#' @param rename.level boolean, if TRUE all the levels of the qualitative 
#' variables are renamed as follows: "variable_name=level_name". 
#' This prevents to have identical names of the levels.
#' @param graph boolean, if TRUE the following graphics are displayed 
#' for the first two dimensions of PCAmix: component map of the individuals, 
#' plot of the squared loadings of all the variables (quantitative and 
#' qualitative), plot of the correlation circle (if  quantitative variables 
#' are available), component map of the levels (if qualitative variables 
#' are available).
#' @param weight.col.quanti vector  of weights for the quantitative variables.
#' @param weight.col.quali vector of the weights for the qualitative variables.
#' @details If X.quali is not specified (i.e. NULL), only quantitative 
#' variables are available and standard PCA is performed. If X.quanti 
#' is NULL, only qualitative variables are available and standard MCA 
#' is performed.
#' 
#' Missing values are replaced by means for quantitative variables
#' and by zeros in the indicator matrix for qualitative variables.
#' 
#' PCAmix performs squared loadings in (\code{sqload}). Squared loadings 
#' for a qualitative variable are correlation ratios between the variable 
#' and the principal components. For a quantitative variable, 
#' squared loadings are the squared correlations  between the variable 
#' and the principal components.
#' 
#' Note that when all the  p variables are qualitative, the factor 
#' coordinates (scores) of the n observations are equal to the factor 
#' coordinates (scores) of standard MCA times square root of p and the 
#' eigenvalues are then equal to the usual eigenvalues of MCA times p. 
#' When all the variables are quantitative, PCAmix gives exactly the same 
#' results as standard PCA.
#' @return  \item{eig}{a matrix containing the eigenvalues, the percentages of variance and the cumulative percentages of variance.}
#' \item{ind}{a list containing the results for the individuals (observations):
#'  \itemize{
#'     \item \code{$coord}: factor coordinates (scores) of the individuals,
#'     \item \code{$contrib}: absolute contributions of the individuals,
#'      \item \code{$contrib.pct}: relative contributions of the individuals,
#'      \item \code{$cos2}: squared cosinus of the individuals.
#'    }}
#' \item{quanti}{a list containing the results for the quantitative variables:
#'   \itemize{
#'     \item \code{$coord}: factor coordinates (scores) of the quantitative variables,
#'      \item \code{$contrib}: absolute contributions of the quantitative variables,
#'      \item \code{$contrib.pct}: relative contributions of the quantitative variables (in percentage),
#'     \item \code{$cos2}: squared cosinus of the quantitative variables.
#'   }}
#'\item{levels}{a list containing the results for the levels of the qualitative variables:
#'    \itemize{
#'      \item \code{$coord}: factor coordinates (scores) of the levels,
#'      \item \code{$contrib}: absolute contributions of the levels,
#'      \item \code{$contrib.pct}: relative contributions of the levels (in percentage),
#'     \item \code{$cos2}: squared cosinus of the levels.
#'   }}
#'\item{quali}{a list containing the results for the qualitative variables:
#'    \itemize{
#'      \item \code{$contrib}: absolute contributions of the qualitative variables (sum of absolute contributions of the levels of the qualitative variable),
#'      \item \code{$contrib.pct}: relative contributions (in percentage) of the qualitative variables (sum of relative contributions of the levels of the qualitative variable).
#'    }}
#'\item{sqload}{a matrix of dimension  (\code{p}, \code{ndim}) containing 
#' the squared loadings of the quantitative and qualitative variables.}
#'\item{coef}{the coefficients of the linear combinations used to 
#' construct the principal components of PCAmix, and to predict coordinates (scores) of new observations in the function \code{\link{predict.PCAmix}}.}
#'  \item{M}{the vector of the weights of the columns used in the Generalized Singular Value Decomposition.}
#' @seealso \code{\link{print.PCAmix}}, \code{\link{summary.PCAmix}}, \code{\link{predict.PCAmix}}, \code{\link{plot.PCAmix}}
#' @author Marie Chavent \email{marie.chavent@u-bordeaux.fr}, Amaury Labenne.
#' @references 
#' Chavent M., Kuentz-Simonet V., Labenne A., Saracco J., Multivariate analysis of mixed data: The PCAmixdata R package, arXiv:1411.4911 [stat.CO].
#' @examples 
#' #PCAMIX:
#' data(wine)
#' str(wine)
#' X.quanti <- splitmix(wine)$X.quanti
#' X.quali <- splitmix(wine)$X.quali
#' pca<-PCAmix(X.quanti[,1:27],X.quali,ndim=4)
#' pca<-PCAmix(X.quanti[,1:27],X.quali,ndim=4,graph=FALSE)
#' pca$eig
#' pca$ind$coord
#' 
#' #PCA:
#' data(decathlon)
#' quali<-decathlon[,13]
#' pca<-PCAmix(decathlon[,1:10])
#' pca<-PCAmix(decathlon[,1:10], graph=FALSE)
#' plot(pca,choice="ind",coloring.ind=quali,cex=0.8,
#'      posleg="topright",main="Scores")
#' plot(pca, choice="sqload",main="Squared correlations")
#' plot(pca, choice="cor",main="Correlation circle")
#' pca$quanti$coord
#' 
#' #MCA
#' data(flower)
#' mca <- PCAmix(X.quali=flower[,1:4],rename.level=TRUE)
#' mca <- PCAmix(X.quali=flower[,1:4],rename.level=TRUE,graph=FALSE)
#' plot(mca,choice="ind",main="Scores")
#' plot(mca,choice="sqload",main="Correlation ratios")
#' plot(mca,choice="levels",main="Levels")
#' mca$levels$coord
#' 
#' #Missing values
#' data(vnf)
#' PCAmix(X.quali=vnf,rename.level=TRUE)
#' vnf2<-na.omit(vnf)
#' PCAmix(X.quali=vnf2,rename.level=TRUE)

PCAmix<- function (X.quanti=NULL,X.quali=NULL,ndim=5,rename.level=FALSE,
                   weight.col.quanti=NULL,weight.col.quali=NULL,graph=TRUE)
{
  cl <- match.call()

  rec <- recod(X.quanti, X.quali,rename.level)
  n <- rec$n
  p <- rec$p
  p1 <- rec$p1
  p2 <- p - p1
  X <- rec$X
  W <- rec$W
  m <- ncol(W) - p1
 
  indexj <- rec$indexj
  
  #construction of the metrics
  N <- rep(1/n, n)
  M1 <- M2 <- NULL #standard metric for PCAmix
  P1 <- P2 <- NULL #supplementary metric for columns like the weights from MFAmix
  if (p1!=0) 
  {
      M1 <- rep(1,p1) 
      P1 <- rep(1,p1) 
      if (!(is.null(weight.col.quanti)))   
      {
        if (length(weight.col.quanti) != ncol(X.quanti))
          stop("the length of \"weight.col.quant\" is different from the number of columns in X.quanti",call. = FALSE)
        P1 <- weight.col.quanti
      }
  }
  
  if (p2!=0)
  {
    ns <- apply(rec$G, 2, sum)
    M2 <- n/ns
    P2 <- rep(1,m)
    if (!(is.null(weight.col.quanti)))   
    {
    if (length(weight.col.quali) != ncol(X.quali))
      stop("the length of \"weight.col.quali\" is different from the number of columns in X.quanti",call. = FALSE)
    P2 <- rep(weight.col.quali,rec$nbmoda) 
    }
  }

  M <- c(M1,M2)
  P <- c(P1,P2)
  M.col <- P*M
  names(M.col) <- colnames(W)
  
  #GSVD
  e <- svd.triplet(W, N, M.col)
  V.total.dim <- e$V
  U.total.dim <- e$U
  d.total.dim <- e$vs
  
  #explained inertia
  q <- qr(W)$rank
  eig <- matrix(0, q, 3)
  colnames(eig) <- c("Eigenvalue", "Proportion", "Cumulative")
  rownames(eig) <- paste("dim", 1:q, sep = " ")
  eig[, 1] <- e$vs[1:q]^2
  eig[, 2] <- 100 * eig[, 1]/sum(eig[, 1], na.rm = T)
  eig[1, 3] <- eig[1, 2]
  if (q > 1) 
  {
    for (j in 2:q) eig[j, 3] <- eig[j, 2] + eig[j - 1, 3]
  }
  
  #number of retained dimensions
  if (ndim <= 1) 
    stop("\"ndim\" must be an interger greater or equal to 2",call. = FALSE)
  ndim <- min(ndim, q)
  d <- e$vs[1:ndim]

  
  #scores
  U <- e$U[, 1:ndim,drop=FALSE]
  rownames(U) <- rownames(W)
  colnames(U) <- paste("dim", 1:ndim, sep = " ")

  F <- sweep(U,2,STATS=d,FUN="*")
  F.total.dim <- sweep(U.total.dim,2,STATS=d.total.dim,FUN="*")

  contrib.ind<-F^2/n 
  contrib.ind.pct <- sweep(contrib.ind, 2, STATS = d^2, FUN = "/")
  cos2.ind <- sweep(F^2, 1, STATS = apply(F.total.dim, 1, function(v) {return(sum(v^2))
  }), FUN = "/")
  
  result.ind <- list(coord = F, contrib = contrib.ind, 
                     contrib.pct = 100 * contrib.ind.pct,
                     cos2 = cos2.ind)
  
  #loadings and contributions
  A1 <- A2 <- NULL
  C1 <- C2 <- NULL
  contrib.quanti <- contrib.quali <- NULL
  
  V <- e$V[, 1:ndim,drop=FALSE]
  rownames(V) <- colnames(W)
  colnames(V) <- paste("dim", 1:ndim, sep = " ")
  
  if(p1 >0)
  {
    V1 <- V[1:p1, ,drop=FALSE]
    V1.total.dim <- V.total.dim[1:p1, ,drop=FALSE]
    A1 <- sweep(V1,2,STATS=d,FUN="*")
    A1.total.dim <- sweep(V1.total.dim,2,STATS=d.total.dim,FUN="*")
    #contrib.quanti <- sweep(A1^2, 1, STATS = M1, FUN = "*")
    contrib.quanti <- sweep(A1^2, 1, STATS = M1*P1, FUN = "*")
    contrib.quanti.pct <- sweep(contrib.quanti, 2, STATS = d^2, 
                                FUN = "/")
    cos2.quanti <- sweep(A1^2, 1, STATS = apply(A1.total.dim, 
                                                1, function(v) {
                                                  return(sum(v^2))
                                                }), FUN = "/")
  }
  if(p2 >0)
  {
    V2 <- V[(p1 + 1):(p1 + m), ,drop=FALSE]
    V2.total.dim <- V.total.dim[(p1 + 1):(p1+m), ,drop=FALSE]
    A2 <- sweep(V2,2,STATS=d,FUN="*")
    A2 <- sweep(A2,1,STATS=M2,FUN="*")
    A2.total.dim <- sweep(V2.total.dim,2,STATS=d.total.dim,FUN="*")
    A2.total.dim <- sweep(A2.total.dim,1,STATS=M2,FUN="*") 
    contrib.moda <- sweep(A2^2, 1, STATS = ns/n, FUN = "*")
    if (!is.null(weight.col.quali)) 
      contrib.moda <- sweep(contrib.moda, 1, STATS = P2, FUN = "*")
    contrib.moda.pct <- sweep(contrib.moda, 2, STATS = d^2, FUN = "/")
    cos2.moda <- sweep(A2^2, 1, STATS = apply(A2.total.dim, 1, function(v) {return(sum(v^2))}), FUN = "/")
    contrib.quali <-matrix(NA,p2,ndim)
    rownames(contrib.quali) <- colnames(X.quali)
    colnames(contrib.quali) <- paste("dim", 1:ndim, sep = " ")
    for (j in 1:p2)
    {
      contrib.quali[j,] <- apply(contrib.moda[which(indexj == (j+p1))-p1, ,drop=FALSE], 2, sum)
    }
    contrib.quali.pct<-sweep(contrib.quali, 2, STATS = d^2, FUN = "/")
  }
  
  sqload <- rbind(contrib.quanti, contrib.quali) #correlation ratio
  weight.col <- c(weight.col.quanti,weight.col.quali)
  if (!is.null(weight.col))
    sqload <- sweep(sqload,1,weight.col,"/")
  quali.eta2 <- NULL
  if (p2 > 0) quali.eta2 <- sqload[(p1+1):(p1+p2),,drop=FALSE]
  
  #organization of the results
  result.quanti <- result.levels <- result.quali <- NULL
  if (p1!=0) 
    result.quanti <- list(coord = A1, contrib= contrib.quanti, contrib.pct = 100 * contrib.quanti.pct, 
                          cos2 = cos2.quanti)
  if (p2!=0) 
  {
    result.levels <- list(coord = A2, contrib=contrib.moda, contrib.pct = 100 * contrib.moda.pct, 
                          cos2 = cos2.moda)
    result.quali<-list(contrib = contrib.quali, contrib.pct=contrib.quali.pct*100)
  }
  
 

  #coefficient of linear combinaisons defining PC
  
  coef <- list()
  for (i in 1:ndim)
  {
    beta <- V[,i]*M.col
    if (p1 > 0) beta[1:p1] <-  beta[1:p1]/rec$s[1:p1]
    beta0 <- -sum(beta*rec$g)
    coef[[i]] <- as.matrix(c(beta0,beta))
  }
  names(coef) <- paste("dim",1:ndim, sep = "")
  
  #matrix A for PCArot
  A2rot <- NULL
  if (p2 >0) A2rot <- sweep(A2,1,STATS=sqrt(ns/n) ,FUN="*")
  A <- rbind(A1,A2rot) 
  Z <- rec$Z
  res <- list(call = cl, 
              eig = eig, 
              ind = result.ind,
              quanti = result.quanti, 
              levels = result.levels,
              quali=result.quali,
              sqload = sqload, 
              coef = coef, Z = Z, 
              M = M.col,
              quanti.sup=NULL,
              levels.sup=NULL,
              sqload.sup=NULL,
              rec.sup=NULL,
              scores.stand = U, 
              scores = F, 
              V = V, 
              A = A, 
              categ.coord = A2, 
              quanti.cor = A1, 
              quali.eta2 = quali.eta2, 
              rec = rec, 
              ndim = ndim, 
              W = W, 
              rename.level=rename.level)
  class(res) <- "PCAmix"
  if (graph==TRUE) {
    plot.PCAmix(res)
    if (p1 != p) 
      plot.PCAmix(res, choice = "levels")
    if (p1 != 0) 
      plot.PCAmix(res, choice = "cor")
    plot.PCAmix(res, choice = "sqload")
  }
  return(res)
}