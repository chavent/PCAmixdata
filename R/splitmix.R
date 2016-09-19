splitmix<-function(base){
  type<-NULL
  base<-data.frame(base,check.names=T)
  testint <- FALSE
  for (v in 1:ncol(base)) {
    if (is.numeric(base[, v])) 
      type = c(type, "QT")
    if (is.integer(base[, v])) 
    {
     #warning("Columns of class integer are considered as quantitative")
      testint <- TRUE
      
    }
    if (is.character(base[, v])) 
    {
      type = c(type, "QL")
      #warning("Columns of class character are considered as qualititative")
    } 
    if (is.factor(base[, v])) 
    {
      type = c(type, "QL")
    } 
  }
  if (testint) warning("Columns of class integer are considered as quantitative")
  ind.QT<-which(type=="QT")
  ind.QL<-which(type=="QL")
  ind.tot<-factor(type)
  X.quanti<-data.frame(base[,ind.QT])
  names(X.quanti)<-names(base[ind.QT])
  X.quali<-data.frame(base[,ind.QL])
  names(X.quali)<-names(base[ind.QL])
  if (ncol(X.quanti)==0) {X.quanti<-NULL}
  if (ncol(X.quali)==0) {X.quali<-NULL}
  
  if (nb.level(ind.tot)==2) {
    typ.group<-"MIX"
  } else if (levels(ind.tot)=="QL") {
    typ.group<-"QL"
  } else 
    typ.group<-"QT"
  
  return(list(X.quanti=X.quanti,X.quali=X.quali,typ.group=typ.group))
}