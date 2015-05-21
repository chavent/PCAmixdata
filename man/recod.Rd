\name{recod}
\alias{recod}
\title{Recoding of the data matrices}
\usage{
  recod(X.quanti, X.quali,rename.level=FALSE)
}
\arguments{
  \item{X.quanti}{a numerical data matrix.}

  \item{X.quali}{a categorical data matrix.}
  
  \item{rename.level}{boolean, if TRUE all the levels of the qualitative variables are renamed as follows: "variable_name=level_name".}
}
\value{
  \item{X}{X.quanti and X.quali concatenated in a single
  matrix.}

  \item{Y}{X.quanti with missing values replaced with mean
  values concatenated with the indicator
  matrix of X.quali with missing values replaced by zeros.}

  \item{Z}{X.quanti standardized (centered and reduced by
  standard deviations) concatenated with the indicator
  matrix of X.quali centered and reduced with the square
  roots of the relative frequencies of the categories.}

  \item{W}{X.quanti standardized (centered and reduced by
  standard deviations) concatenated with the indicator
  matrix of X.quali centered.}

  \item{n}{the number of observations.}

  \item{p}{the total number of variables}

  \item{p1}{the number of quantitative variables}

  \item{p2}{the number of qualitative variables}

  \item{g}{the means of the columns of Y}

  \item{s}{the standard deviations of the columns of
  Y}

  \item{G}{The indicator matix of X.quali with missing
  values replaced by 0}

  \item{Gcod}{The indicator matix G reduced with the square
  roots of the relative frequencies of the categories}
}
\description{
  Recoding of the quantitative and of the qualitative data matrix.
}

