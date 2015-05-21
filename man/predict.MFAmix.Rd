\name{predict.MFAmix}
\alias{predict.MFAmix}
\title{Prediction of new scores in MFAmix}

\description{
This function performs the scores of new observations 
on the principal components of MFAmix. In other words, this function is projecting the new observations onto the principal components of MFAmix obtained previoulsy on a separated dataset. Note that the new observations must be described with the same variables than those used in MFAmix. The groups of variables must also be identical.
}

\usage{
  \method{predict}{MFAmix}(object, data, groups, name.groups,...)
}
\arguments{
  \item{object}{an object of class MFAmix (output of the
  function MFAmix)}

 \item{data}{a data frame containing the description of the new observations on all the variables. This data frame  will be split into \code{G} groups according to the vector
  \code{groups}.}

  \item{groups}{a vector of size \code{p} whose values
  indicate at which group each variable belongs.}

  \item{name.groups}{a vector of size \code{G} which
  contains the names of the groups. Spaces and special characters are not allowed.}
  
  \item{...}{further arguments passed to or from other methods. They are ignored in this function.}
}
\value{
  Returns the matrix of the scores of the new observations on
  the principal components or on the rotated principal components of MFAmix.
}

\examples{
#MFAmix:
data(gironde)
class.var<-c(rep(1,9),rep(2,5),rep(3,9),rep(4,4))
nom.groupes<-c("employment","housing","services","environment")
data.mfamix<-cbind(gironde$employment,gironde$housing,
                gironde$services,gironde$environment)
res<-MFAmix(data=data.mfamix,groups=class.var,name.groups=nom.groupes, 
              rename.level=TRUE, ndim=3,graph=FALSE)

#Predict on the same data, we obtain the same scores
res.pred<-predict(res,data=data.mfamix,groups=class.var,
                name.groups=nom.groupes)
res$ind$coord/res.pred

}
\seealso{
  \code{\link{MFAmix}}
}

\references{
Chavent M., Kuentz-Simonet V., Labenne A., Saracco J., Multivariate analysis of mixed data: The PCAmixdata R package, arXiv:1411.4911 [stat.CO].
}

\author{
  Amaury Labenne \email{amaury.labenne@irstea.fr},
  Marie Chavent, Vanessa Kuentz, Benoit Liquet, Jerome Saracco
}

