#' @export
# recodqual <- function(X,rename.level=FALSE)
# {
#   nmoda <- apply(X, 2, 
#                  function(x){length(levels(as.factor(x)))})
#   if (sum(nmoda==1)!=0)
#   {
#     listmoda <- toString(paste(colnames(X)[nmoda==1]))
#     stop(sprintf("The column(s) %s in X.quali contain(s) only one category", 
#                  listmoda), call.=FALSE)
#   }
#   GNA <- tab.disjonctif.NA(X,rename.level)
#   G <- replace(GNA,is.na(GNA),0)
#   return(G)	
# }


recodqual <-
function(X,rename.level=FALSE)
	{
		#X <- as.matrix(X)
		GNA <- tab.disjonctif.NA(X,rename.level)
		G <- replace(GNA,is.na(GNA),0)
		n <- nrow(GNA)
		if (n > 1)
		{
		  ns <- apply(G,2,sum)
		  nmiss <- apply((apply(GNA,2,is.na)),2,sum)
		  if(sum((n-nmiss)==ns)!=0)
		  {
		    levelsunique <- colnames(G)[(n-nmiss)==ns ]
		    stop("There are columns in X.quali where all
		          the categories are identical", call.=FALSE)

		  }
		}
		return(G)
	}

