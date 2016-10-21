\name{splitgroups}
\alias{splitgroups}
\title{splitgroups}
\usage{
  splitgroups(data, groups, name.groups)
}

\description{
If the p variables of a data matrix of dimension (n,p) are separated into G groups, this functions splits this data matrix into \code{G} datasets according the groups membership. 
}

\arguments{
  \item{data}{the a data matrix into \code{G} datasets with \code{n} rows
  and \code{p} columns.}

  \item{groups}{a vector of size \code{p} whose values
  indicate at which group belongs each variable.}

  \item{name.groups}{a vector of size \code{G} which
  contains names for each group we want to create. The names are given in the increasing order of the numbers of the groups.}
}
\value{
  \item{data.groups}{a list of G data matrix: one matrix for each group.}

  \item{listvar.groups}{The list of the variables in each group.}
}


\examples{
data(decathlon)
split <- splitgroups(decathlon,groups=c(rep(1,10),2,2,3),
          name.groups=c("Epreuve","Classement","Competition"))
split$data.groups$Epreuve
}

