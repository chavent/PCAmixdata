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
#' @param rename.level boolean, if TRUE all the levels of the qualitative variables
#' are renamed as follows: "variable_name=level_name". This prevents to have 
#' identical names for the levels.
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


predict.MFAmix<-function (object, data, rename.level=FALSE,...) 
{
  mfa <- object
  if (!inherits(mfa, "MFAmix")) 
    stop("use only with \"MFAmix\" objects")
  
#   if (length(groups)!=ncol(data))
#     stop("\"groups\" must be a vector of size the number of variables in \"data\"",call.=FALSE)
#   
  if (mfa$rename.level) 
    rename.level=TRUE
  if ((rename.level) & (!mfa$rename.level))
    stop("perform MFAmix with argument rename.level=TRUE",call.=FALSE)
  
#  n <-nrow(data)
#   ngroup <-length(unique(groups))
#   if (length(name.groups)!=ngroup)
#     stop("Invalid length of \"name.groups\"",call.=FALSE)
#   
#   data.groups <- splitgroups(data=data,groups=groups,name.groups=name.groups)$data.groups 
#   listvar.group <- splitgroups(data=data,groups=groups,name.groups=name.groups)$listvar.groups
#   size.groups<-unlist(lapply(data.groups,ncol)) #number of variables in each group
#   
#   #reorder the data matrix by block 
#   data.ord <- data.groups[[1]]
#   for(g in 2:ngroup)
#   {
#     data.ord <- cbind(data.ord, data.groups[[g]])
#   }
#   #from now we use these data
#   groups<-rep(1:ngroup,size.groups)
#   data<-data.ord
#   
  X.quanti<-splitmix(data)$X.quanti
  X.quali<-splitmix(data)$X.quali
  
  rec <- recod(X.quanti, X.quali,rename.level=rename.level)
  Y <- rec$Y
  n <- rec$n
  beta <- mfa$global.pca$coef
  
#   ncol_tot<-length(beta[[1]])-1
#   #test if that all levels of the caregorical variable in the first dataset are observed in the new dataset.
#   if (ncol_tot!=ncol(Y)){
#     Ymodif<-matrix(0,nrow=n,ncol=ncol_tot)
#     colname.part<-colnames(Y)
#     colname.tot<-names(beta[[1]][-1,])
#     colnames(Ymodif)<-colname.tot
#     Ymodif[,colname.part]<-Y[,colname.part]
#     Y<-Ymodif  
#   }
  
  if (!is.null(X.quanti)) {
    label <- rownames(X.quanti)
    n1 <- nrow(X.quanti)
    p1 <- ncol(X.quanti)
    if (p1 != mfa$global.pca$rec$p1) 
      stop("The number of numerical variables in data must be the same than in the learning set",call.=FALSE)
  }
  if (!is.null(X.quali)) {
    label <- rownames(X.quali)
    n2 <- nrow(X.quali)
    p2 <- ncol(X.quali)
    if (p2 != mfa$global.pca$rec$p2) 
      stop("The number of categorical variables in data must be the same than in the learning set",call.=FALSE)
  }
  if (!is.null(X.quanti) && !is.null(X.quali)) {
    if (n1 != n2) 
      stop("The number of objects in X.quanti and X.quali must be the same",call.=FALSE)
    if (sum(rownames(X.quali) != rownames(X.quanti)) != 0) 
      stop("The names of the objects in X.quanti and X.quali must be the same",call.=FALSE)
  }
  
  if (!setequal(colnames(mfa$global.pca$rec$X),colnames(rec$X)))
    stop("The colnames in the new sample are not appropiate",call.=FALSE)
  
  if (rec$p2 >0 )
  {
    #different levels in the new sample and in the learning sample
    if (!setequal(colnames(mfa$global.pca$rec$Y),colnames(Y)))
    {
      col.test <- is.element(colnames(Y),colnames(mfa$global.pca$rec$Y))
      Y <- Y[,col.test,drop=FALSE]
      col.beta <- c(TRUE,is.element(colnames(mfa$global.pca$rec$Y),colnames(Y)))
      beta <- lapply(beta,function(x) {x[col.beta,1,drop=FALSE]})
    }
  }
  
  scores <- matrix(, n, length(beta))
  for (g in 1:length(beta)) scores[, g] <- Y %*% beta[[g]][-1] + 
    beta[[g]][1]
  colnames(scores) <- paste("dim", 1:length(beta), sep = "")
  rownames(scores) <- label
  return(scores)
}


