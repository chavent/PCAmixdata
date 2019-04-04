#' @export
#' @name tab.disjonctif.NA
#' @title  Built an indicator matrix
#' @description This function built the indicator matrix of a qualitative data matrix. 
#' Missing observations are indicated as NAs.
#' @param tab a categorical data matrix..
#' @param rename.level boolean, if TRUE all the levels of the qualitative variables 
#' are renamed as follows: variable_name=level_name.
#' @details This function uses the code of the function tab.disjonctif implemented 
#' in the package FactoMineR but is different. Here, a NA value appears when 
#' a category has not been observed in a row. In the function tab.disjonctif of the package FactoMineR, 
#' a new column is created in that case.
#' @return Returns the indicator matrix with NA for missing observations.
#' @examples 
#' data(vnf)
#' X <- vnf[1:10,9:12]
#' tab.disjonctif.NA(X)
#' 
#' 
tab.disjonctif.NA <- function (tab,rename.level=FALSE) {
  tab <- as.data.frame(tab)
  modalite.disjonctif <- function(i) {
    moda <- tab[, i]
    nom <- names(tab)[i]
    n <- length(moda)
    moda <- as.factor(moda)
    x <- matrix(0, n, length(levels(moda)))
    ind<-(1:n) + n * (unclass(moda) - 1)
    indNA<-which(is.na(ind))
    
    x[(1:n) + n * (unclass(moda) - 1)] <- 1
    x[indNA,]<-NA 
    if (rename.level==TRUE){
      dimnames(x) <- list(row.names(tab), paste(nom, levels(moda), 
                                                sep = "="))
    }    else{
      dimnames(x) <- list(row.names(tab), levels(moda))
    }
    
    
    return(x)
  }
  if (ncol(tab) == 1) 
    res <- modalite.disjonctif(1)
  else {
    res <- lapply(1:ncol(tab), modalite.disjonctif)
    res <- as.matrix(data.frame(res, check.names = FALSE))
  }
  return(res)
}

