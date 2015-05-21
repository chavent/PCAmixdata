\name{gironde}

\alias{gironde}

\docType{data}

\title{gironde}

\description{A list of 4 datasets caracterizing conditions of life of 542 cities in Gironde. The four datasets correspond to four thematics relative to conditions of life. Each dataset contains a different number of variables (quantitative and/or qualitative). The first three datasets come from the 2009 population census realized in Gironde by INSEE (Institut National de la Statistique et des Etudes Economiques). The fourth come from an IGN (Institut National de l'Information Geographique et forestiere) database.}

\value{
\item{gironde$employment}{This data frame contains the description of 542 cities by 9 quantitative variables. These variables are related to employment conditions like, for instance, the average income (income), the percentage of farmers (farmer).}

\item{gironde$housing}{This data frame contains the description of 542 cities by 5 variables (2 qualitative variables and 3 quantitative variables). These variables are related to housing conditions like, for instance, the population density (density), the percentage of counsil housing within the cities (council).}

\item{gironde$services}{This data frame contains the description of 542 cities by 9 qualitative variables. These variables are related to the number of services within the cities, like, for instance, the number of bakeries (baker) or the number of post office (postoffice).}

\item{gironde$environment}{This data frame contains the description of 542 cities by 4 quantitative variables. These variables are related to the natural environment of the cities, like, for instance the percentage of agricultural land (agricul) or the percentage of buildings (building).}
}

\format{
A list of 4 data frames.
}
  
  \usage{data(gironde)}

  \source{
  www.INSEE.fr
  
  www.ign.fr
  
  http://siddt.grenoble.cemagref.fr/
  
Multivariate analysis of mixed data: The PCAmixdata R package, M. Chavent, V.  Kuentz-Simonet, A. Labenne, J. Saracco, arXiv:1411.4911 [stat.CO]
  }
  
  
  \keyword{datasets}
  