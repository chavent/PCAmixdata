##' @export
##' @method predict MFAmix
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


