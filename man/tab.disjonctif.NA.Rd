\name{tab.disjonctif.NA}
\alias{tab.disjonctif.NA}
\title{Built an indicator matrix}
\usage{
  tab.disjonctif.NA(tab,rename.level=FALSE)
}
\description{
  This function built the indicator matrix of a qualitative
  data matrix. Missing observations are indicated as NA.
}

\arguments{
  \item{tab}{a categorical data matrix.}
   \item{rename.level}{boolean, if TRUE all the levels of the qualitative variables are renamed as follows: "variable_name=level_name".}

}
\value{
 Returns the indicator matrix with NA for missing observations.
}

\details{
This function uses the code of the function "tab.disjonctif" implemented in the package  \bold{FactoMineR} but is different. Here, a NA value appears when a category has not been observed in a row. In the function "tab.disjonctif" of the package  \bold{FactoMineR}, a new column is created in that case. In the output of the function recodqual, the NA are replaced with 0.
}

\examples{
data(vnf)
X <- vnf[1:10,9:12]
tab.disjonctif.NA(X)
recodqual(X)
}

