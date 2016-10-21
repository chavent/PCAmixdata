\name{splitmix}
\alias{splitmix}
\title{splitmix}
\usage{
  splitmix(data)
}
\description{
  Splits a mixed data matrix in two data sets: one with the quantitative variables and one with the qualitative variables. Here, the columns of class "integer are considered quantitative. If you want this column to be considered as qualitative, it must be of class character of factor.
}

\arguments{
  \item{data}{a data matrix or a data.frame with a mixture of quantitative and qualitative variables.}
}
\value{
  \item{X.quanti}{a data matrix containing only the quantitative variables.}

  \item{X.quali}{A data.frame containing only the qualitative variables.}
}

\examples{
data(decathlon)
data.split <- splitmix(decathlon)
data.split$X.quanti
data.split$X.quali
}

