#' @export
#' @keywords internal
#' 
Lg<-function(liste.mat){
  n<-nrow(liste.mat[[1]])
  nbr.mat<-length(liste.mat)
  Z.liste<-list()
  for (i in 1:nbr.mat){
    Z.liste[[i]]<-recod(X.quanti=splitmix(liste.mat[[i]])$X.quanti,X.quali=splitmix(liste.mat[[i]])$X.quali,rename.level=TRUE)$Y 
    Z.liste[[i]]<-scale(Z.liste[[i]],scale=F)
  }
  names(Z.liste)<-names(liste.mat)
  
  Lg<-matrix(NA,ncol=nbr.mat,nrow=nbr.mat)
  rownames(Lg)<-colnames(Lg)<-names(Z.liste)
  for (i in 1: (nbr.mat)){
    for (j in 1:(nbr.mat)){
      Lg[i,j]<-sum(diag(t(Z.liste[[i]]) %*% Z.liste[[j]] %*% t(Z.liste[[j]]) %*% Z.liste[[i]]))*(1/n^2)   
    }
  }
  return(Lg)
}