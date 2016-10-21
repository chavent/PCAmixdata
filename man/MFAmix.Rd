\name{MFAmix}
\alias{MFAmix}
\title{Multiple factor analysis of mixed data}
\description{
Performs multiple factor analysis to analyze a set of individuals (observations) described by several groups of variables. Variables within a group can be a mixture of quantitative and qualitative variables. 
}
\usage{
  MFAmix(data, groups, name.groups, ndim=5, rename.level=FALSE, graph = TRUE,
    axes = c(1, 2))
}
\arguments{
  \item{data}{a data frame with \code{n} rows and \code{p}
  columns containing all the variables. This data frame
  will be split into \code{G} groups according to the vector
  \code{groups}.}

  \item{groups}{a vector which gives the groups of the columns in \code{data}.}

  \item{name.groups}{a vector of size \code{G} which gives the names of the groups.}

  \item{ndim}{number of dimensions kept in the results (by default 5).}
  
  \item{rename.level}{boolean, if TRUE all the levels of the qualitative variables are renamed as follows: "variable_name=level_name". This prevents to have identical names for the levels.}

  \item{graph}{boolean, if TRUE the following graphics are displayed for the first two dimensions of PCAmix: plot of the individuals coordinates, plot of the squared loadings of variables, plot of the partial axes, plot of the correlation circle (if quantitative variables are available), plot of the levels component map (if qualitative variables are available).}

  \item{axes}{a length 2 vector specifying the axes to
  plot.}
}


\value{

  \item{eig}{a matrix containing the eigenvalues, the percentages of variance and the cumulative percentages of variance.}

  \item{ind}{a list containing the results
  for the individuals (observations):
  \itemize{
  \item \code{$coord}: factor coordinates (scores) of the individuals,
  \item \code{$contrib}: absolute contributions of the individuals,
  \item \code{$contrib.pct}: relative contributions of the individuals,
  \item \code{$cos2}: squared cosinus of the individuals.
  }}

  \item{quanti}{a list containing the
  results for the quantitative variables:
  \itemize{
  \item \code{$coord}: factor coordinates (scores) of the quantitative variables,
  \item \code{$contrib}: absolute contributions of the quantitative variables,
  \item \code{$contrib.pct}: relative contributions of the quantitative variables (in percentage),
  \item \code{$cos2}: squared cosinus of the quantitative variables.
  }}

  \item{levels}{a list containing the
  results for the levels of the qualitative variables:
  \itemize{
  \item \code{$coord}: factor coordinates (scores) of the levels,
  \item \code{$contrib}: absolute contributions of the levels,
  \item \code{$contrib.pct}: relative contributions of the levels (in percentage),
  \item \code{$cos2}: squared cosinus of the levels.
  }}
  
  \item{quali}{a list containing the
  results for the qualitative variables:
  \itemize{
  \item \code{$contrib}: absolute contributions of the qualitative variables (sum of absolute contributions of the levels of the qualitative variable),
  \item \code{$contrib.pct}: relative contributions (in percentage) of the qualitative variables (sum of relative contributions of the levels of the qualitative variable).
  }}
 
  \item{sqload}{a matrix of dimension  (\code{p}, \code{ndim}) containing the squared loadings of the quantitative and qualitative variables.}
  
  \item{coef}{the coefficients of the linear combinations used to construct the principal components of MFAmix, and to predict coordinates (scores) of new observations in the function \code{\link{predict.MFAmix}}.}
  
  \item{eig.separate}{a matrix containing the \code{ndim} first eigenvalues of the separated analyses of each group.}

  \item{separate.analyses}{the results for the separated analyses of each group.}

  \item{groups}{a list containing the results for the groups:
  \itemize{
  \item \code{$Lg}: Lg coefficients between groups,
  \item \code{$RV}: RV coefficients between groups,
  \item \code{$contrib}: contributions of the groups (sum of variable contributions belonging to the group)
  \item \code{$contrib.pct}:  relative contributions of the groups times 100,
  }}

  \item{partial.axes}{a matrix containing the coordinates of the partial axes.}

   \item{ind.partial}{a list of \code{G} matrices containing the coordinates
  of the partial individuals.}
 
   \item{listvar.group}{list the variables in each group. It is usefull to check the adequacy between the vector \code{groups} and the vector \code{name.groups}.}
   
  \item{global.pca}{an object of class \code{PCAmix} containing the results of \code{MFAmix} considered as a unique \code{PCAmix}.}
}


\author{
  Amaury Labenne \email{amaury.labenne@irstea.fr},
  Marie Chavent, Vanessa Kuentz, Benoit Liquet, Jerome Saracco
}

\references{
Chavent M., Kuentz-Simonet V., Labenne A., Saracco J., Multivariate analysis of mixed data: The PCAmixdata R package, arXiv:1411.4911 [stat.CO].

Escofier, B. and Pages, J. (1994). Multiple factor analysis (afmult package). Computational statistics & data analysis, 18(1):121-140.
  
Le, S., Josse, J., and Husson, F. (2008). Factominer: an r package for multivariate analysis. Journal of statistical software, 25(1):1-18.
}

\seealso{
\code{\link{print.MFAmix}}, \code{\link{summary.MFAmix}}, \code{\link{predict.MFAmix}}, \code{\link{plot.MFAmix}}
}

\examples{
data(gironde)

class.var<-c(rep(1,9),rep(2,5),rep(3,9),rep(4,4))
names <- c("employment","housing","services","environment")

dat<-cbind(gironde$employment[1:20,],gironde$housing[1:20,],
      gironde$services[1:20,],gironde$environment[1:20,])
      
res<-MFAmix(data=dat,groups=class.var,
      name.groups=names, rename.level=TRUE, ndim=3,graph=FALSE)
      
summary(res)

}

\details{
Multiple Factor Analysis (MFA) developed by Escofier and Pages in 1983 is a method of factorial analysis to deal with multiple groups of variables collected on the same observations. The main idea of MFA is to normalize each group by dividing all the variables belonging to this group by the first eigenvalue coming from the Principal Component Analysis (PCA) of this group. Then, a usual PCA on all the weighted variables taken together is applied. Initially this method has been developed for groups only containing quantitative variables. Afterwards this method has been improved to deal simultaneously with groups of qualitative variables and groups of quantitative variables. The \code{MFAmix} method allows to perform MFA method for groups containing a mixture of quantitative and qualitative variables

One of the outputs available in the MFAmix method are the squared loadings (\code{sqload}). Squared loadings for a qualitative variable are correlation ratios between the variable and the principal components. For a quantitative variable, squared loadings are the squared correlation  between the variable and the principal components.

Some others outputs are specific to MFA:
\itemize{
  \item Coordinates of  groups are the sum of the absolute contributions of variables belonging to the groups,
  \item Partial individuals coordinates are factor coordinates of individuals according to a specific group. The partial coordinates can be achieved by projecting the data set of each group onto the principal component space of MFAmix,
  \item Partial axes of a group  are correlation between each principal components of the separated analyses of the group and the principal components of MFAmix.
}
}

\keyword{multivariate, mixed type data, data structured in groups}

