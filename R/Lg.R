# @name Lg
# @title  Lg coefficient
# @description Computes a data fame with all the Lg coefficients between each matrix of a list.
# @param liste.mat a list of G matrices.
# @param ponde  a vector of size G with the ponderation associated to each matrix.
# @return a data frame with G rows and G colums with the Lg coefficients.
# @examples 
# V0<-c("a","b","a","a","b")
# V01<-c("c","d","e","c","e")
# V1<-c(5,4,2,3,6)
# V2<-c(8,15,4,6,5)
# V3<-c(4,12,5,8,7)
# V4<-c("vert","vert","jaune","rouge","jaune")
# V5<-c("grand","moyen","moyen","petit","grand")
# G1<-data.frame(V0,V01,V1)
# G2<-data.frame(V2,V3)
# G3<-data.frame(V4,V5)
# liste.mat<-list(G1,G2,G3)
# Lg(liste.mat)
# @keywords internal 
Lg<-function(liste.mat){
  n<-nrow(liste.mat[[1]])
  nbr.mat<-length(liste.mat)
  Z.liste<-liste.mat
 
  Lg<-matrix(NA,ncol=nbr.mat,nrow=nbr.mat)
  rownames(Lg)<-colnames(Lg)<-names(Z.liste)
  for (i in 1: (nbr.mat)){
    for (j in 1:(nbr.mat)){
      Lg[i,j]<-sum(diag(t(Z.liste[[i]]) %*% Z.liste[[j]] %*% t(Z.liste[[j]]) %*% Z.liste[[i]]))*(1/n^2)   
    }
  }
  return(Lg)
}