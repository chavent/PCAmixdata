#' @export
#' @name supvar.MFAmix
#' @method supvar MFAmix
#' @title  Supplementary variables in MFAmix
#' @description Performs the coordinates of supplementary variables and groups on the component of an object of class \code{MFAmix}.
#' @param obj an object of class \code{MFAmix}.
#' @param data.sup a numeric matrix of data.
#' @param groups.sup a vector which gives the groups of the columns in \code{data.sup}.
#' @param name.groups.sup a vector which gives the names of the supplementary groups.
#' @param rename.level boolean, if TRUE all the levels of the qualitative variables
#' are renamed as follows: "variable_name=level_name". This prevents to have identical
#' names of the levels.
#' @param \dots further arguments passed to or from other methods.
#' @examples
#' data(wine)
#' X.quanti <- splitmix(wine)$X.quanti[,1:5]
#' X.quali <- splitmix(wine)$X.quali[,1,drop=FALSE]
#' X.quanti.sup <- splitmix(wine)$X.quanti[,28:29]
#' X.quali.sup <- splitmix(wine)$X.quali[,2,drop=FALSE]
#' data <- cbind(X.quanti,X.quali)
#' data.sup <- cbind(X.quanti.sup,X.quali.sup)
#' groups <-c(1,2,2,3,3,1)
#' name.groups <- c("G1","G2","G3")
#' groups.sup <- c(1,1,2)
#' name.groups.sup <- c("Gsup1","Gsup2")
#' mfa <- MFAmix(data,groups,name.groups,ndim=4,rename.level=TRUE,graph=FALSE)
#' mfa.sup <- supvar(mfa,data.sup,groups.sup,name.groups.sup,rename.level=TRUE)
#'
#'
supvar.MFAmix <- function(obj, data.sup, groups.sup, name.groups.sup,rename.level=FALSE,...)
{
  
  if (!inherits(obj, "MFAmix")) 
    stop("use only with \"MFAmix\" objects",call. = FALSE)
  mfa <- obj
  n <-nrow(data.sup)
  ngroup <-length(unique(groups.sup))
  if (length(name.groups.sup)!=ngroup)
    stop("Invalid length of \"name.groups.sup\"",call. = FALSE)
  
  data.groups <- splitgroups(data=data.sup,groups=groups.sup,name.groups=name.groups.sup)$data.groups 
  listvar.group <- splitgroups(data=data.sup,groups=groups.sup,name.groups=name.groups.sup)$listvar.groups
  size.groups<-unlist(lapply(data.groups,ncol)) #number of variables in each group
  
  #reorder the data matrix by block 
  data.ord <- data.groups[[1]]
  if (ngroup>1)
  {
    for(g in 2:ngroup)
    {
      data.ord <- cbind(data.ord, data.groups[[g]])
    }
  }
 
  
  #from now we use these data
  groups<-rep(1:ngroup,size.groups)
  data<-data.ord
  
  res.separate.pcamix <- list()
  data.groups.stand<-list()
  
  #PCAmix with each group
  ndim <- ncol(mfa$global.pca$scores)
  for(i in 1:ngroup){
    base.qt<-splitmix(data.groups[[i]])$X.quanti
    base.ql<-splitmix(data.groups[[i]])$X.quali
    res.separate.pcamix[[i]]<-PCAmix(X.quanti=base.qt, X.quali=base.ql, ndim=ndim, rename.level=rename.level, graph=FALSE)
    data.groups.stand[[i]]<-res.separate.pcamix[[i]]$Z
  }
  names(res.separate.pcamix) <- names(data.groups.stand) <- name.groups.sup
  
  #partial axes
  B <- mfa$global.pca$scores
  partial.axes.sup <- list()
  for (i in 1:ngroup)
  {
    A <-res.separate.pcamix[[i]]$scores
    partial.axes.sup[[i]] <- cor(A,B)
    rownames( partial.axes.sup[[i]]) <- paste(colnames(A),name.groups.sup[i],sep=".")
  }
  names(partial.axes.sup) <- name.groups.sup
 
  
  
  #projections of supplementary variables
  split <- splitmix(data.sup)
  X.quanti.sup <- split$X.quanti
  X.quali.sup <-  split$X.quali
  rec.sup <- recod(X.quanti.sup, X.quali.sup,rename.level)
  W.sup <- rec.sup$W # standardized quantitative+centered dummy variables
  n <- rec.sup$n
  if (n!=mfa$global.pca$rec$n) 
    stop("number of rows in data.sup is not correct",call. = FALSE)
  p1 <- rec.sup$p1
  p2 <- rec.sup$p - p1
  m <- ncol(W.sup) - p1
  U <- mfa$global.pca$scores.stand
  N <- rep(1/n, n)
  A <- t(W.sup)%*%diag(N)%*%U
  rownames(A) <- colnames(W.sup)
  colnames(A) <-  paste("dim", 1:ncol(A), sep = "")
  nor <- apply(W.sup,2,function(x){sqrt(sum(x^2)/n)})
  B <- sweep(A,1,STATS = nor,FUN = "/")^2 #cos2 !!!
  quanti.sup <- levels.sup <- NULL
  sqload.sup <- matrix(NA,rec.sup$p,ncol(A))
  rownames(sqload.sup) <- colnames(rec.sup$X)
  colnames(sqload.sup) <- colnames(A)
  if (p1!=0)
  {
    quanti.sup <- list(coord=NULL,cos2=NULL)
    quanti.sup$coord <- A[1:p1,]
    quanti.sup$cos2 <- B[1:p1,]
    sqload.sup[1:p1,] <- A[1:p1,]^2
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
  
  #indices
  index.col <- c(split$col.quant,split$col.qual) #in data ordered by block
  #res.global$sqload[index.col,]
  index.group <- groups.sup[index.col] #in data ordered by type
  index1 <- index2 <-  NULL
  if (rec.sup$p1>0) index1 <- index.group[1:rec.sup$p1] 
  if (rec.sup$p2>0) index2 <- rep(index.group[(rec.sup$p1+1):(rec.sup$p1+rec.sup$p2)],rec.sup$nbmoda)
  indexg <- c(index1,index2)
  
  # result for supplementary groups
  # the coordinate of a supp group is the sum of the squared loadings if its variables
  # divided by the first eigen value of its PCAmix. It is a weighted Lg coefficient
  
  coord.group.sup <- matrix(NA,ngroup,ndim)
  for (i in 1:ngroup)
    #coord.group.sup[i,] <- apply(sqload.sup[index.group==i,,drop=FALSE],2,sum)
    coord.group.sup[i,] <- apply(sqload.sup[index.group==i,,drop=FALSE],2,sum)/res.separate.pcamix[[i]]$eig[1,1]
  colnames(coord.group.sup) <- colnames(sqload.sup)
  rownames(coord.group.sup) <- name.groups.sup
  if (!is.null(levels.sup))
    mfa$levels.sup <- levels.sup
  if (!is.null(quanti.sup))
    mfa$quanti.sup <- quanti.sup
  mfa$sqload.sup <- sqload.sup
  mfa$partial.axes.sup <- partial.axes.sup
  mfa$listvar.group.sup <- listvar.group
  mfa$group.sup <- coord.group.sup
  mfa$index.groupsup <- index.group
  mfa$index.groupsup2 <- indexg
  mfa$rec.sup <- rec.sup
  
  return(mfa)
}