\name{plot.PCAmix}
\alias{plot.PCAmix}
\title{
Graphical outputs of PCAmix and PCArot
}
\description{
Displays the graphical outputs of PCAmix and PCArot. The individuals (observations), the quantitative variables and the levels of the qualitative variables are plotted as points using their factor coordinates (scores). All the variables (quantitative and qualitative) are plotted as points on the same graph using their squared loadings. 
}
\usage{
\method{plot}{PCAmix}(x,axes = c(1, 2), choice = "ind",label=TRUE,
                        coloring.ind=NULL,col.ind=NULL, coloring.var=NULL,
                        lim.cos2.plot=0,lim.contrib.plot=0, posleg="topleft",
                        xlim=NULL,ylim=NULL, cex=1,leg=TRUE,main=NULL,cex.leg=1, ...)
}
\arguments{
  \item{x}{
  an object of class PCAmix obtained with the function \code{PCAmix} or \code{PCArot}.
}
  \item{axes}{
 a length 2 vector specifying the components to plot.
}
  \item{choice}{
 the graph to plot: 
  \itemize{
  \item "ind" for the individuals component map,
  \item "cor" for the correlation circle if quantitative variables are available in the data,
  \item "levels" for the levels components map (if qualitative variables are available in the data),
    \item "sqload" for the plot of the squared loadings of all the variables.
  }
}

\item{label}{
 boolean, if FALSE the labels of the points are not plotted.}
  
\item{coloring.ind}{a qualitative variable such as a character vector or a factor of size n (the number of individuals). The individuals are colored according to the levels of this variable. If NULL, the individuals are not colored.}

\item{col.ind}{a vector of colors, of size the number of levels of \code{coloring.ind}. If NULL, colors are chosen automatically.}

\item{coloring.var}{If "type", the variables in the plot of the squared loadings are colored according to their type (quantitative or qualitative). If NULL, variables are not colored.
  }
  
\item{lim.cos2.plot}{a value between 0 and 1. Points with squared cosinus below this value are not plotted.}
  
\item{lim.contrib.plot}{a value between 0 and 100. Points with relative contributions (in percentage) below this value  are not plotted.}


  
  \item{posleg}{
position of the legend.
}
  \item{xlim}{
a numeric vectors of length 2, giving the x coordinates range.
}
  \item{ylim}{
a numeric vectors of length 2, giving the y coordinates range.
}

  \item{main}{a string corresponding to the title of the
  graph to draw.}

  \item{cex}{cf. function \code{par} in the \bold{graphics}
package}
  

 \item{leg}{
if TRUE a legend is displayed for choice="sqload" and coloring.var="type".
}

 \item{cex.leg}{
a numerical value giving the amount by which the legend should be magnified. Default is 0.8.
}

  \item{\dots}{
arguments to be passed to methods, such as graphical parameters.
}
}

\author{
  Amaury Labenne \email{amaury.labenne@irstea.fr},
  Marie Chavent, Vanessa Kuentz, Benoit Liquet, Jerome Saracco
}

\seealso{
\code{\link{summary.PCAmix}},\code{\link{PCAmix}},\code{\link{PCArot}}, 
}

\details{The observations can be colored according to the levels of a qualitative variable. The observations, the quantitative variables and the levels can be selected according to their squared cosine (lim.cos2.plot) or their relative contribution (lim.contrib.plot) to the component map. Only points with squared cosine or relative contribution greater than a given threshold are plotted. Note that the relative contribution of a point to the component map (a plan) is the sum of the absolute contributions to each dimension, divided by the sum of the corresponding eigenvalues.}

\references{
Chavent M., Kuentz-Simonet V., Labenne A., Saracco J., Multivariate analysis of mixed data: The PCAmixdata R package, arXiv:1411.4911 [stat.CO].

Le, S., Josse, J., and Husson, F. (2008). Factominer: an r package for multivariate analysis. Journal of statistical software, 25(1):1-18.
}
\examples{
data(gironde)

X.quanti<-splitmix(gironde$housing[1:20,])$X.quanti
X.quali<-splitmix(gironde$housing[1:20,])$X.quali

res<-PCAmix(X.quanti,X.quali, rename.level=TRUE, ndim=3,graph=FALSE)

#----quantitative variables on the correlation circle
plot(res,choice="cor",cex=0.6)

#----individuals component map
plot(res,choice="ind",cex=0.6)

#----individuals colored with the qualitative variable "houses"
plot(res,choice="ind",cex=0.6,coloring.ind=X.quali$houses) 

#----individuals selected according to their squared cosine
plot(res,choice="ind",cex=0.6,lim.contrib.plot=0.8)

#----all the variables plotted with the squared loadings
plot(res,choice="sqload",cex=0.6)

#----variables colored according to their type (quanti or quali)
plot(res,choice="sqload",cex=0.6,coloring.var="type") 

#----levels component map
plot(res,choice="levels",cex=0.6)

}

