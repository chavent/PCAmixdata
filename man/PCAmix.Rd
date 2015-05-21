\name{PCAmix}
\alias{PCAmix}
\title{Principal component analysis of mixed data}
\usage{
  PCAmix(X.quanti = NULL, X.quali = NULL, ndim = 5, rename.level=FALSE,
    weight.col = NULL, weight.row = NULL, graph = TRUE)
}

\description{
  Performs principal component analysis of  a set of individuals (observations) described by a mixture of qualitative and quantitative variables. PCAmix includes
  the ordinary principal component analysis (PCA) and multiple correspondence analysis (MCA) as special cases.
}

\arguments{
  \item{X.quanti}{a numerical data matrix of dimension (n,p1) with the description of n observations on p1 quantitative variables.}

  \item{X.quali}{a categorical data matrix of dimension (n,p2) with the description of n observations on p2 qualitative variables.}

  \item{ndim}{number of dimensions kept in the results (by
  default 5).}
  
 \item{rename.level}{boolean, if TRUE all the levels of the qualitative variables are renamed as follows: "variable_name=level_name". This prevents to have identical names of the levels.}

  \item{graph}{boolean, if TRUE the following graphics are
  displayed for the first two dimensions of PCAmix: component map of
  the individuals, plot of the squared loadings of all the variables (quantitative and qualitative), plot of the correlation circle (if
  quantitative variables are available), component map of the levels (if qualitative variables are available).}

  \item{weight.col}{the vector of the weights of the
  quantitative variables and the weights of the dummy variables of the levels. If \code{NULL}, the weights are all equal to 1.}

  \item{weight.row}{the vector of the weights of the
  individuals. If \code{NULL}, the weights are all equal to 1/\code{n}.}
}

\details{
If X.quali is not specified (i.e. NULL), only quantitative variables are available and standard PCA is performed. If X.quanti is NULL, only qualitative variables are available and standard MCA is performed.

Missing values are replaced by means for quantitative variables and by zeros in the indicator
matrix for qualitative variables.

PCAmix performs squared loadings in (\code{sqload}). Squared loadings for a qualitative variable are correlation ratios between the variable and the principal components. For a quantitative variable, squared loadings are the squared correlations  between the variable and the principal components.

Note that when all the  p variables are qualitative, the factor coordinates (scores) of the n observations are equal to the factor coordinates (scores) of standard MCA times square root of p and the eigenvalues are then equal to the usual eigenvalues of MCA times p. When all the variables are quantitative, PCAmix gives exactly the same results as standard PCA.
}


\value{
  \item{eig}{a matrix containing the eigenvalues, the percentages of variance and the cumulative percentages of variance.}

  \item{ind}{a list containing the results for the individuals (observations):
    \itemize{
  \item \code{$coord}: factor coordinates (scores) of the individuals,
  \item \code{$contrib}: absolute contributions of the individuals,
  \item \code{$contrib.pct}: relative contributions of the individuals,
  \item \code{$cos2}: squared cosinus of the individuals.
  }}

  \item{quanti}{a list containing the results for the quantitative variables:
  \itemize{
  \item \code{$coord}: factor coordinates (scores) of the quantitative variables,
  \item \code{$contrib}: absolute contributions of the quantitative variables,
  \item \code{$contrib.pct}: relative contributions of the quantitative variables (in percentage),
  \item \code{$cos2}: squared cosinus of the quantitative variables.
  }}

  \item{levels}{a list containing the results for the levels of the qualitative variables:
  \itemize{
   \item \code{$coord}: factor coordinates (scores) of the levels,
   \item \code{$contrib}: absolute contributions of the levels,
   \item \code{$contrib.pct}: relative contributions of the levels (in percentage),
   \item \code{$cos2}: squared cosinus of the levels.
  }}
  
  \item{quali}{a list containing the results for the qualitative variables:
  \itemize{
  \item \code{$contrib}: absolute contributions of the qualitative variables (sum of absolute contributions of the levels of the qualitative variable),
  \item \code{$contrib.pct}: relative contributions (in percentage) of the qualitative variables (sum of relative contributions of the levels of the qualitative variable).
  }}
 
  \item{sqload}{a matrix of dimension  (\code{p}, \code{ndim}) containing the squared loadings of the quantitative and qualitative variables.}

  \item{coef}{the coefficients of the linear combinations used to construct the principal components of PCAmix, and to predict coordinates (scores) of new observations in the function \code{\link{predict.PCAmix}}.}
 
  \item{M}{the vector of the weights of the columns used in the Generalized Singular Value Decomposition.}
  
}

\author{
  Marie Chavent \email{marie.chavent@u-bordeaux.fr},
  Vanessa Kuentz, Amaury Labenne, Benoit Liquet, Jerome Saracco
}

\references{
Chavent M., Kuentz-Simonet V., Labenne A., Saracco J., Multivariate analysis of mixed data: The PCAmixdata R package, arXiv:1411.4911 [stat.CO].

Chavent, M., Kuentz, V., Saracco, J. (2011), Orthogonal
  Rotation in PCAMIX. Advances in Classification and Data
  Analysis, Vol. 6, pp. 131-146.

Kiers, H.A.L., (1991), Simple structure in Component
  Analysis Techniques for mixtures of qualitative and
  quantitative variables, Psychometrika, 56, 197-212.

Le, S., Josse, J., and Husson, F. (2008). Factominer: an r package for multivariate analysis. Journal of statistical software, 25(1):1-18.

Pages, J. (2004). Analyse factorielle de donnees mixtes. 
Revue de statistique appliquee, 52(4):93-111.
}

\seealso{
\code{\link{print.PCAmix}}, \code{\link{summary.PCAmix}}, \code{\link{predict.PCAmix}}, \code{\link{plot.PCAmix}}
}

\examples{
#PCAMIX:
data(wine)
X.quanti <- wine[,c(3:29)]
X.quali <- wine[,c(1,2)]
pca<-PCAmix(X.quanti,X.quali,ndim=4)
pca<-PCAmix(X.quanti,X.quali,ndim=4,graph=FALSE)
pca$eig
pca$ind$coord

#PCA:
data(decathlon)
quali<-decathlon[,13]
pca<-PCAmix(decathlon[,1:10])
pca<-PCAmix(decathlon[,1:10], graph=FALSE)
plot(pca,choice="ind",quali=quali,cex=0.8,
    posleg="topright",main="Scores")
plot(pca, choice="sqload",main="Squared correlations")
plot(pca, choice="cor",main="Correlation circle")
pca$quanti$coord

#MCA
data(flower)
mca <- PCAmix(X.quali=flower[,1:4],rename.level=TRUE)
mca <- PCAmix(X.quali=flower[,1:4],rename.level=TRUE,graph=FALSE)
plot(mca,choice="ind",main="Scores")
plot(mca,choice="sqload",main="Correlation ratios")
plot(mca,choice="levels",main="Levels")
mca$levels$coord

#Missing values
data(vnf)
PCAmix(X.quali=vnf,rename.level=TRUE)
vnf2<-na.omit(vnf)
PCAmix(X.quali=vnf2,rename.level=TRUE)
}



\keyword{principal component analysis, mixed type data}

