#' @export
#' @keywords internal
Lg.pond<-function(liste.mat,ponde){
  Lg.pond<-Lg(liste.mat)
  Lg.pond<-sweep(Lg.pond,1,STATS=ponde,FUN="*")
  Lg.pond<-sweep(Lg.pond,2,STATS=ponde,FUN="*")
  return(Lg.pond)
}