\name{summary.MFAmix}
\alias{summary.MFAmix}
\title{
  Summary of a 'MFAmix' object
}
\description{
  This is a method for the function summary for objects of the class \code{MFAmix}.
}
\usage{
  \method{summary}{MFAmix}(object, ...)
}
\arguments{
  \item{object}{
    an object of class MFAmix obtained with the function \code{MFAmix}.
  }
  \item{\dots}{
    further arguments passed to or from other methods.
  }
}

\value{
  Returns the total number of observations, the number of quantitative variables, the number of qualitative variables with the total number of levels. And all those values are also given by groups.
}


\seealso{
  \code{\link{plot.MFAmix}},\code{\link{MFAmix}} 
}


