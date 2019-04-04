#' @export
#' @name predict.MFAmix
#' @title  Prediction of new scores in MFAmix
#' @description This function performs the scores of new observations 
#' on the principal components of MFAmix. In other words, this function 
#' is projecting the new observations onto the principal components 
#' of MFAmix obtained previoulsy on a separated dataset. 
#' Note that the new observations must be described with the 
#' same variables than those used in MFAmix. The groups of variables
#' must also be identical.
#' @param object an object of class MFAmix obtained with the function 
#' \code{MFAmix}.
#' @param data a data frame containing the description of the new observations 
#' on all the variables. This data frame  will be split into \code{G} groups according 
#' to the vector \code{groups}.
#' @param \ldots urther arguments passed to or from other methods. 
#' They are ignored in this function.
#' @return  Returns the matrix of the scores of the new observations on 
#' the principal components or on the rotated principal components of MFAmix.
#' @seealso \code{\link{MFAmix}}
#' @author Marie Chavent \email{marie.chavent@u-bordeaux.fr}, Amaury Labenne.
#' @references 
#' Chavent M., Kuentz-Simonet V., Labenne A., Saracco J.,
#' Multivariate analysis of mixed data: The PCAmixdata R package, 
#' arXiv:1411.4911 [stat.CO].
#' @examples 
#' data(gironde)
#' class.var<-c(rep(1,9),rep(2,5),rep(3,9),rep(4,4))
#' names<-c("employment","housing","services","environment")
#' dat<-cbind(gironde$employment,gironde$housing,
#'            gironde$services,gironde$environment)
#' n <- nrow(dat)
#' set.seed(10)
#' sub <- sample(1:n,520)
#' 
#' res<-MFAmix(data=dat[sub,],groups=class.var,
#'             name.groups=names, rename.level=TRUE, 
#'             ndim=3,graph=FALSE)
#' 
#' #Predict scores of new data
#' pred<-predict(res,data=dat[-sub,])
#' plot(res,choice="ind",cex=0.6,lim.cos2.plot=0.7)  
#' points(pred[1:5,c(1,2)],col=2,pch=16,cex=0.6)
#' text(pred[1:5,c(1,2)], labels = rownames(dat[-sub,])[1:5],
#'      col=2,pos=3,cex=0.6)


predict.MFAmix<-function (object, data,...) 
{
  mfa <- object
  if (!inherits(mfa, "MFAmix")) 
    stop("use only with \"MFAmix\" objects")
  
  X.quanti <- splitmix(data)$X.quanti
  X.quali <- splitmix(data)$X.quali
  
  beta <- mfa$coef
  p2 <- 0
  
  train.rec <- mfa$global.pca$rec
  train.rename.level <-  mfa$rename.level
  
  if (!is.null(X.quanti)) 
  {
    label <- rownames(X.quanti)
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
        for (v in 1:ncol(X.quali))
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
    G <- G[,colnames(train.rec$G), drop=FALSE]
  } else G <- NULL
  
  
  Y <- cbind(Y1,G)
  n <- nrow(Y)
  
  coord <- matrix(NA,n,length(beta))
  for (g in 1: length(beta)) 
    coord[,g] <- Y %*% beta[[g]][-1] +  beta[[g]][1]
  
  colnames(coord) <- paste("dim", 1:length(beta), sep = "")
  rownames(coord) <- rownames(Y)
  return(coord)		
}


