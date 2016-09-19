\name{PCArot}
\alias{PCArot}
\title{
Varimax rotation in PCAmix
}
\description{
Orthogonal rotation in PCAmix by maximization of the varimax function expressed in terms of PCAmix squared loadings (correlation ratios for qualitative variables and squared correlations for quantitative variables).
PCArot includes the ordinary varimax rotation in Principal Component Analysis (PCA) and a varimax-type rotation in Multiple Correspondence Analysis (MCA) as special cases.
}
\usage{
PCArot(obj, dim, itermax = 100, graph = TRUE)
}
\arguments{
  \item{obj}{
an object of class PCAmix.
}
  \item{dim}{
number of rotated Principal Components.
}
  \item{itermax}{
maximum number of iterations in the Kaiser's practical optimization algorithm based on successive pairwise planar rotations. 
}
  \item{graph}{
 boolean, if TRUE the following graphs are displayed for the first two dimensions after rotation: plot of the individuals (factor coordinates), plot of the variables (squared loadings)
 plot of the correlation circle (if quantitative variables are available), plot of the levels component map (if qualitative variables are available).
}
}

\value{

\item{eig}{variances of the ndim dimensions after rotation.}

\item{ind$coord}{a n by dim quantitative matrix which contains the coordinates (scores) of the n individuals on the dim rotated principal components.
}
\item{quanti$coord}{a p1  by dim quantitative matrix which contains the coordinates (loadings) of the p1  quantitative variables after rotation.
The coordinates of the quantitative variables after rotation are correlations with the rotated principal components.}

\item{levels$coord}{a m  by dim quantitative matrix which contains the coordinates  of the m levels on the dim rotated principal components.}

\item{quali$coord}{a p2 by dim quantitative matrix which contains the coordinates of the p2 qualitative variables on the dim rotated principal components. Coordinates of the qualitative variables after rotation are correlation ratio with the rotated principal components.}

\item{coef}{coefficients of the linear combinations used to construct the rotated principal components of PCAmix.}

\item{theta}{angle of rotation if dim is equal to 2. }

\item{T}{matrix of rotation.}
}
\references{
Chavent, M., Kuentz, V., Saracco, J. (2011), Orthogonal Rotation in PCAMIX. Advances in Classification and Data Analysis, Vol. 6, pp. 131-146.

Chavent M., Kuentz-Simonet V., Labenne A., Saracco J., Multivariate analysis of mixed data: The PCAmixdata R package, arXiv:1411.4911 [stat.CO].

Kiers, H.A.L., (1991), Simple structure in Component Analysis Techniques for mixtures of
qualitative and quantitative variables, Psychometrika, 56, 197-212.

}
\author{
Marie Chavent \email{marie.chavent@u-bordeaux.fr}, Vanessa Kuentz, Benoit Liquet, Jerome Saracco
}

\seealso{
\code{\link{plot.PCAmix}}, \code{\link{summary.PCAmix}}, \code{\link{PCAmix}}, \code{\link{predict.PCAmix}}
}

\details{
If X.quali is not specified (i.e. NULL) in the previous PCAmix step, only quantitative variables are available and  standard varimax rotation in PCA is performed. If X.quanti is NULL, only qualitative variables are available and  varimax-type rotation in MCA is performed. Note that p1 is the number of quantitative variables, p2 is the number of qualitative variables and m is the total number of levels of the p2 qualitative variables.
}
\examples{

#PCAMIX:
data(wine)
pca<-PCAmix(X.quanti=wine[,c(3:29)],X.quali=wine[,1:2],ndim=4,graph=FALSE)
pca

rot<-PCArot(pca,3)
rot
rot$eig #percentages of variances after rotation

plot(rot,choice="ind",coloring.ind=wine[,1],
	    posleg="bottomleft", main="Rotated scores")
plot(rot,choice="sqload",main="Squared loadings after rotation")
plot(rot,choice="levels",main="Levels after rotation")
plot(rot,choice="cor",main="Correlation circle after rotation")



#PCA:
data(decathlon)
quali<-decathlon[,13]
pca<-PCAmix(decathlon[,1:10], graph=FALSE)

rot<-PCArot(pca,3)
plot(rot,choice="ind",coloring.ind=quali,cex=0.8,
	    posleg="topright",main="Scores after rotation")
plot(rot, choice="sqload", main="Squared correlations after rotation")
plot(rot, choice="cor", main="Correlation circle after rotation")

#MCA
data(flower)
mca <- PCAmix(X.quali=flower[,1:4],rename.level=TRUE,graph=FALSE)

rot<-PCArot(mca,2)
plot(rot,choice="ind",main="Scores after rotation")
plot(rot, choice="sqload", main="Correlation ratios after rotation")
plot(rot, choice="levels", main="Levels after rotation")


}

\keyword{ multivariate }