\name{splitmix}
\alias{splitmix}
\title{splitmix}
\usage{
  splitmix(base)
}
\description{
  Splits a mixed data matrix in two data sets: one with the quantitative variables and one with the qualitative variables.
}

\arguments{
  \item{base}{a mixed data matrix (a data.frame).}
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

