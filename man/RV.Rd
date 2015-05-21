\name{RV}
\alias{RV}
\title{Coefficient RV}
\usage{
  RV(liste.mat)
}
\arguments{
  \item{liste.mat}{a list of \code{k} matrices.}
}
\value{
  \item{RV}{a (\code{k}, \code{k}) matrix containing the RV coefficients.}
}
\description{
Computes RV coefficients between matrices.
}

\seealso{
\code{\link{RV.pond}}, \code{\link{Lg}}, \code{\link{Lg.pond}}, 
}

\examples{
V0<-c("a","b","a","a","b")
V01<-c("c","d","e","c","e")
V1<-c(5,4,2,3,6)
V2<-c(8,15,4,6,5)
V3<-c(4,12,5,8,7)
V4<-c("vert","vert","jaune","rouge","jaune")
V5<-c("grand","moyen","moyen","petit","grand")
G1<-data.frame(V0,V01,V1)
G2<-data.frame(V2,V3)
G3<-data.frame(V4,V5)
liste.mat<-list(G1,G2,G3)
RV(liste.mat)
}
\references{
  {Escofier B et Pages J (1998), Analyses factorielles simples et multiples, Dunod, 3e ed.
  
  Pages J, Analyse factorielle multiple appliquee aux variables qualitatives et aux donnees mixtes. Revue de statistique appliquee, tome 50, num 4 (2002), p. 5-37}
}