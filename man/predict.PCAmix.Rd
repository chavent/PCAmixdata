\name{predict.PCAmix}
\alias{predict.PCAmix}
\title{Prediction of new scores in PCAmix or PCArot}

\description{
This function performs the scores of new observations on the principal components of PCAmix. If the components have been rotated, this function performs the scores of the
  new observations on the rotated principal components. In other words, this function is projecting the new observations onto the principal components of PCAmix (or PCArot) obtained previoulsy on a separated dataset. Note that the new observations must be described with the same variables than those used in PCAmix (or PCArot).
}

\usage{
  \method{predict}{PCAmix}(object, X.quanti = NULL, X.quali = NULL,...)
}
\arguments{
  \item{object}{an object of class PCAmix (output of the
  function PCAmix or PCArot).}

  \item{X.quanti}{a numerical data matrix with the description of the new
  observations on the quantitative variables.}

  \item{X.quali}{a categorical data matrix with the description of the new
  observations on the qualitative variables.}
  
  \item{...}{further arguments passed to or from other methods. They are ignored in this function.}
 } 
\value{
  Returns the matrix of the scores of the new observations on
  the principal components or on the rotated principal components of PCAmix.
}

\examples{
data(decathlon)
n <- nrow(decathlon)
sub <- sample(1:n,20)
pca<-PCAmix(decathlon[sub,1:10], graph=FALSE)
predict(pca,decathlon[-sub,1:10])
rot <- PCArot(pca,dim=4)
predict(rot,decathlon[-sub,1:10])
}

\seealso{
  \code{\link{PCAmix}}
}

\author{
  Marie Chavent \email{marie.chavent@math.u-bordeaux1.fr},
  Vanessa Kuentz, Benoit Liquet, Jerome Saracco
}

\references{
Chavent M., Kuentz-Simonet V., Labenne A., Saracco J., Multivariate analysis of mixed data: The PCAmixdata R package, arXiv:1411.4911 [stat.CO].
}

