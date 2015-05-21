\name{plot.MFAmix}
\alias{plot.MFAmix}
\title{Graphical outputs of MFAmix}

\description{
Displays the graphical outputs of MFAmix. Individuals (observations), quantitative variables and levels of the qualitative variables are plotted as points using their factor coordinates (scores) in MFAmix. All the variables (quantitative and qualitative) are plotted on the same graph as points using their squared loadings. The groups of variables are plotted using their contributions to the component coordinates. Partial axes  and partial individuals of separated analyses can also be plotted.
}

\usage{
\method{plot}{MFAmix}(x, axes = c(1, 2), choice = "axes",
label=TRUE, coloring.var = NULL, coloring.ind=NULL,
col.ind=NULL, col.groups=NULL,  partial = NULL, chrono=FALSE,
lim.cos2.plot = 0, lim.contrib.plot=0, xlim = NULL,  ylim = NULL,
cex = 1, main = NULL, new.plot = FALSE,leg=TRUE,
posleg="topleft",cex.leg=0.8, ...)
}

\arguments{
  \item{x}{an object of class MFAmix.}

  \item{axes}{a length 2 vector specifying the components to plot.}

  \item{choice}{the graph to plot: 
  \itemize{
  \item "ind" for the individuals, 
  \item "cor" for the correlation circle of the quantitative variables, 
  \item "levels" for the levels of of the qualitative variables,
  \item "sqload" for the plot of the squared loadings of all the variables, 
  \item "groups" for the plot of the contributions of the groups, 
  \item "axes" for the correlation circle of the partial axes.  
  }}

  \item{label}{boolean,  if FALSE the labels of the points are not plotted.}
  
   \item{coloring.ind}{a qualitative variable such as a character vector or a factor of size n (the number of individuals). The individuals are colored according to the levels of this variable. If NULL, the individuals are not colored.}
  
  \item{col.ind}{a vector of colors, of size the number of levels of \code{coloring.ind}. If NULL, colors are chosen automatically.}
  
  \item{coloring.var}{a value to choose among: 
  \itemize{
  \item "type": the variables in the plot of the squared loadings are colored according to their type (quantitative or qualitative), 
  \item "groups": the variables are colored according to their group.
  \item NULL: variables are not colored.
  }} 
  \item{col.groups}{a vector of colors, of size the number of groups. If NULL, colors are chosen automatically.}


  \item{partial}{a vector with the row names of the individuals, for which the partial points should be drawn.
  (by default, partial = NULL and no partial points are
  drawn). Partial points are colored according to \code{col.groups}}
  
  \item{chrono}{boolean, if TRUE, the partial points are linked (useful when groups correspond to different times).}


\item{lim.cos2.plot}{a value between 0 and 1. Points with squared cosinus below this value are not plotted.}
  
\item{lim.contrib.plot}{a value between 0 and 100. Points with relative contributions (in percentage) below this value  are not plotted.}

  \item{posleg}{position of the legend.}
  
  \item{xlim}{a numeric vectors of length 2, giving the x
  coordinates range.}

  \item{ylim}{a numeric vectors of length 2, giving the y
  coordinates range.}
  
 \item{main}{a string corresponding to the title of the
  graph to draw.}

  \item{cex}{cf. function \code{par} in the \bold{graphics}
  package.}
 
  \item{new.plot}{boolean, if TRUE, a new graphical device
  is created.}
  
  \item{leg}{boolean, if TRUE, a legend is displayed.}
  
  \item{cex.leg}{a numerical value giving the amount by which
  the legend should be magnified. Default is 0.8.}
  
  \item{...}{arguments to be passed to methods, such as graphical parameters.}
}

\details{The observations can be colored according to the levels of a qualitative variable. The observations, the quantitative variables and the levels can be selected according to their squared cosine (lim.cos2.plot) or their relative contribution (lim.contrib.plot) to the component map. Only points with squared cosine or relative contribution greater than a given threshold are plotted. Note that the relative contribution of a point to the component map (a plan) is the sum of the absolute contributions to each dimension, divided by the sum of the corresponding eigenvalues.}

\references{
Chavent M., Kuentz-Simonet V., Labenne A., Saracco J., Multivariate analysis of mixed data: The PCAmixdata R package, arXiv:1411.4911 [stat.CO].

Escofier, B. and Pages, J. (1994). Multiple factor analysis (afmult package). Computational statistics & data analysis, 18(1):121-140.
  
Le, S., Josse, J., and Husson, F. (2008). Factominer: an r package for multivariate analysis. Journal of statistical software, 25(1):1-18.
}

\examples{
data(gironde)

class.var<-c(rep(1,9),rep(2,5),rep(3,9),rep(4,4))
names <- c("employment","housing","services","environment")
dat<-cbind(gironde$employment[1:20,],gironde$housing[1:20,],
      gironde$services[1:20,],gironde$environment[1:20,])
res<-MFAmix(data=dat,groups=class.var,
      name.groups=names, rename.level=TRUE, ndim=3,graph=FALSE)

#---- quantitative variables
plot(res,choice="cor",cex=0.6)
plot(res,choice="cor",cex=0.6,coloring.var="groups")
plot(res,choice="cor",cex=0.6,coloring.var="groups",
    col.groups=c("red","yellow","pink","brown"),leg=TRUE)

#----partial axes
plot(res,choice="axes",cex=0.6)
plot(res,choice="axes",cex=0.6,coloring.var="groups")
plot(res,choice="axes",cex=0.6,coloring.var="groups",
    col.groups=c("red","yellow","pink","brown"),leg=TRUE)


#----groups
plot(res,choice="groups",cex=0.6)   #no colors for groups
plot(res,choice="groups",cex=0.6,coloring.var="groups") 
plot(res,choice="groups",cex=0.6,coloring.var="groups",
    col.groups=c("red","yellow","pink","blue")) 

#----squared loadings
plot(res,choice="sqload",cex=0.8)    #no colors for groups
plot(res,choice="sqload",cex=0.8,coloring.var="groups",
    posleg="topright") 
plot(res,choice="sqload",cex=0.6,coloring.var="groups",
    col.groups=c("red","yellow","pink","blue"),ylim=c(0,1)) 

plot(res,choice="sqload",cex=0.8,coloring.var="type",
    cex.leg=0.9,posleg="topright")  

#----individuals 
plot(res,choice="ind",cex=0.6) 

#----individuals with squared cosine greater than 0.5
plot(res,choice="ind",cex=0.6,lim.cos2.plot=0.5)  

#----individuals colored with a qualitative variable
nbchem <- gironde$services$chemist[1:20]
plot(res,choice="ind",cex=0.6,coloring.ind=nbchem,
    posleg="topright")   
plot(res,choice="ind",coloring.ind=nbchem,
  col.ind=c("pink","brown","darkblue"),label=FALSE,posle="topright")     

#----partial individuals colored by groups
plot(res,choice="ind",partial="all",cex=0.6)
plot(res,choice="ind",partial=c("AUBIAC","ARCACHON"),
    cex=0.6,posleg="bottomright")

#----levels of qualitative variables
plot(res,choice="levels",cex=0.8)
plot(res,choice="levels",cex=0.8,coloring.var="groups")

#levels with squared cosine greater than 0.6
plot(res,choice="levels",cex=0.8, lim.cos2.plot=0.6)
}

\author{
  Amaury Labenne \email{amaury.labenne@irstea.fr},
  Marie Chavent, Vanessa Kuentz, Benoit Liquet, Jerome Saracco
}