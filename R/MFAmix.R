#' @export
MFAmix<-function(data, groups, name.groups, ndim=5, rename.level=FALSE, 
                 graph=TRUE, axes=c(1,2))
{
  
  cl <- match.call()
  
  if (length(groups)!=ncol(data))
    stop("\"groups\" must be a vector of size the number of variables in \"data\"",call. = FALSE)
  
  n <-nrow(data)
  ngroup <-length(unique(groups))
  if (length(name.groups)!=ngroup)
    stop("Invalid length of \"name.groups\"",call. = FALSE)
  
  data.groups <- splitgroups(data=data,groups=groups,name.groups=name.groups)$data.groups 
  listvar.group <- splitgroups(data=data,groups=groups,name.groups=name.groups)$listvar.groups
  size.groups<-unlist(lapply(data.groups,ncol)) #number of variables in each group

  #reorder the data matrix by block 
  data.ord <- data.groups[[1]]
  if (ngroup > 1)
    for(g in 2:ngroup)
    {
      data.ord <- cbind(data.ord, data.groups[[g]])
    }

  
  #from now we use these data
  groups<-rep(1:ngroup,size.groups)
  data<-data.ord
  
  res.separate.pcamix <- list()
  data.groups.stand<-list()
  
  #PCAmix with each group
  for(i in 1:ngroup){
    base.qt<-splitmix(data.groups[[i]])$X.quanti
    base.ql<-splitmix(data.groups[[i]])$X.quali
    res.separate.pcamix[[i]]<-PCAmix(X.quanti=base.qt, X.quali=base.ql, ndim=ndim, rename.level=rename.level, graph=FALSE)
    data.groups.stand[[i]]<-res.separate.pcamix[[i]]$Z
  }
  names(res.separate.pcamix) <- names(data.groups.stand) <- name.groups
  
  eig.groups <- rep(NA,ngroup) #first eigenvalues of each PCAmix
  for(i in 1:ngroup) eig.groups[i]<-res.separate.pcamix[[i]]$eig[1,1]
  names(eig.groups) <- name.groups
  
  eig.separate <- matrix(NA,ndim,ngroup)
  colnames(eig.separate) <- name.groups
  rownames(eig.separate) <- paste("dim", 1:ndim, sep =" ")
  for (i in 1:ngroup)
  {
    v <- res.separate.pcamix[[i]]$eig[,1]
    eig.separate[1:min(ndim,length(v)),i] <-v[1:min(ndim,length(v))]
  }
   

  #global PCAmix
  split <- splitmix(data)
  data.quanti <- split$X.quanti
  data.quali <-  split$X.quali
  
  w.col <- rep(NA,ncol(data)) #columns are weighted with the inverse of the first eigenvalue of its group
  for (i in 1:ngroup)
  {
    sel <- which(groups==i)
    w.col[sel] <- 1/rep(eig.groups[i],length(sel))
  }
  
  w.col.quant <- w.col.qual <- NULL
  if (!is.null(data.quanti)) 
      w.col.quant <- w.col[split$col.quant]
  if (!is.null(data.quali)) 
      w.col.qual <- w.col[split$col.qual]
  
  res.global<-PCAmix(X.quanti=data.quanti, X.quali=data.quali,ndim=ndim,
                     rename.level=rename.level, graph=FALSE,
                     weight.col.quanti=w.col.quant, weight.col.quali=w.col.qual)
  
  #number of retained dimensions
  ndim<-ncol(res.global$V)
  
  #indices
  index.col <- c(split$col.quant,split$col.qual) #in data ordered by block
  #res.global$sqload[index.col,]
  index.group <- groups[index.col] #in data ordered by type
  rec <- recod(data.quanti,data.quali, rename.level=rename.level)
  index1 <- index2 <-  NULL
  if (rec$p1>0) index1 <- index.group[1:rec$p1] 
  if (rec$p2>0) index2 <- rep(index.group[(rec$p1+1):(rec$p1+rec$p2)],rec$nbmoda)
  indexg <- c(index1,index2)
  
  #partial individuals
 
  V <- res.global$V
  M <- res.global$M

  Fpart <- list()
  for (i in 1:ngroup) 
    {
    Wpart <- res.global$W
    Wpart[,-which(indexg==i)] <- 0
    Fpart[[i]] <- ngroup* sweep(Wpart,2,M,"*") %*% V
  }
  names(Fpart) <- name.groups
  
  
  list.inert.part <- lapply(Fpart,function(x) {apply(x^2,2,sum)})
  inert.part <- 0 #inertia of the n*ngroup partial individuals
  for (i in 1:ngroup)
    inert.part <- inert.part + list.inert.part[[i]]
  inert.b <- apply(res.global$scores^2,2,sum) #inertia of the gravity centers of partial individuals
  ratio.inert <- inert.b * ngroup/inert.part
  
  #partial axes
  B <- res.global$scores
  partial.axes <- list()
  for (i in 1:ngroup)
  {
    A <-res.separate.pcamix[[i]]$scores
    partial.axes[[i]] <- stats::cor(A,B)
    rownames( partial.axes[[i]]) <- paste(colnames(A),name.groups[i],sep=".")
  }
  names(partial.axes) <- name.groups
  
  #results for the groups
  contrib <- rbind(res.global$quanti$contrib, res.global$quali$contrib)
  contrib.group <- matrix(NA,ngroup,ndim)
  for (i in 1:ngroup)
    contrib.group[i,] <- apply(contrib[index.group==i,,drop=FALSE],2,sum)
  rownames(contrib.group) <- name.groups
  colnames(contrib.group) <- colnames(contrib)
  contrib.group.pct <- sweep(contrib.group,2,STATS=res.global$eig[1:ndim,1],FUN="/")*100
  Lg <- Lg.pond(data.groups.stand,1/eig.groups)
  RV <- RV.pond(data.groups.stand,1/eig.groups)
#   res.groups <- list(Lg = Lg, RV = RV, coord = contrib.group.pct/100,
#                      contrib=contrib.group, contrib.pct = contrib.group.pct)
  res.groups <- list(Lg = Lg, RV = RV, contrib=contrib.group, contrib.pct = contrib.group.pct)

  res<-list(call=cl,
            eig=res.global$eig,
            ind=res.global$ind,
            quanti=res.global$quanti,     
            levels=res.global$levels,
            quali=res.global$quali,
            sqload=res.global$sqload,
            coef=res.global$coef,
            eig.separate= eig.separate,
            separate.analyses=res.separate.pcamix,
            groups=res.groups,
            partial.axes=partial.axes,
            ind.partial=Fpart,
            listvar.group=listvar.group,
            global.pca=res.global,
            quanti.sup=NULL,
            levels.sup=NULL,
            sqload.sup=NULL,
            listvar.group.sup=NULL,
            partial.axes.sup=NULL,
            group.sup=NULL,
            rec.sup=NULL,
            inertia.ratio=ratio.inert,
            lst.groups=groups,
            index.group=index.group,
            index.group2=indexg,
            index.groupsup=NULL,
            index.groupsup2=NULL,
            rename.level=rename.level
  )
  
  
  class(res)<-c("MFAmix","list")
  #class(res)<-c("MFAmix")
  
  if (graph) {  
      plot.MFAmix(res,axes=axes,choice="axes",coloring.var="groups")
      plot.MFAmix(res,axes=axes,choice="groups",coloring.var="groups")
      plot.MFAmix(res,axes=axes,choice="ind",cex=0.8)
      plot.MFAmix(res,axes=axes,choice="sqload", coloring.var="groups")
     
     if (!is.null(res.global$quanti$coord)){
       plot.MFAmix(res,axes=axes,choice="cor", coloring.var="groups")
     }
     if (!is.null(res.global$levels$coord)){
       plot.MFAmix(res,axes=axes,choice="levels",coloring.var="groups")
     }
   }
 
  return(res)
}
