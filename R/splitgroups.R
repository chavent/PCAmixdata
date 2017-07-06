#' @export
splitgroups<-function(data,groups,name.groups){
  num.group <- sort(unique(groups))
  if (length(name.groups) != length(num.group))
    stop("groups is not in adequation with name.groups",call.=FALSE)
  data.group <- list()
  for (i in 1:length(num.group))
    data.group[[i]] <- data[,which(groups == num.group[i]),drop=FALSE]
  names(data.group) <- name.groups
  listvar.groups <- lapply(data.group,colnames)
  return(list(data.groups=data.group,listvar.groups=listvar.groups))
}
