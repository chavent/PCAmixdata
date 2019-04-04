#' @export
#' @title  Graphical outputs of MFAmix
#' @description Displays the graphical outputs of MFAmix. 
#' Individuals (observations), quantitative variables and 
#' levels of the qualitative variables are plotted as points 
#' using their factor coordinates (scores) in MFAmix. 
#' All the variables (quantitative and qualitative) are plotted 
#' on the same graph as points using their squared loadings. 
#' The groups of variables are plotted using their contributions 
#' to the component coordinates. Partial axes  and partial 
#' individuals of separated analyses can also be plotted.
#' @param x an object of class MFAmix obtained with the function \code{MFAmix}.
#' @param axes a length 2 vector specifying the components to plot.
#' @param choice the graph to plot: 
#' \itemize{
#' \item "ind" for the individuals, 
#' \item "cor" for the correlation circle of the quantitative variables, 
#' \item "levels" for the levels of of the qualitative variables,
#' \item "sqload" for the plot of the squared loadings of all the variables, 
#' \item "groups" for the plot of the contributions of the groups, 
#' \item "axes" for the correlation circle of the partial axes.  
#' }
#' @param label boolean, if FALSE the labels of the points are not plotted.
#' @param coloring.ind a qualitative variable such as a character vector 
#' or a factor of size n (the number of individuals). The individuals 
#' are colored according to the levels of this variable. If NULL, the 
#' individuals are not colored.
#' @param nb.partial.axes f choice="axes", the maximum number of partial axes 
#' related to each group to plot on the correlation circle. By default equal to 3.
#' @param col.ind a vector of colors, of size the number of levels of 
#' \code{coloring.ind}. If NULL, colors are chosen automatically.
#' @param coloring.var a value to choose among: 
#' \itemize{
#' \item "type": the variables in the plot of the squared loadings are colored 
#' according to their type (quantitative or qualitative), 
#' \item "groups": the variables are colored according to their group.
#' \item NULL: variables are not colored.
#' }
#' @param col.groups a vector of colors, of size the number of groups. 
#' If NULL, colors are chosen automatically.
#' @param partial a vector of class character with the row names of the individuals,
#' for which the partial individuals should be drawn. 
#' By default partial = NULL and no partial points are drawn. 
#' Partial points are colored according to \code{col.groups}
#' @param lim.cos2.plot a value between 0 and 1. Points with squared 
#' cosinus below this value are not plotted.
#' @param lim.contrib.plot a value between 0 and 100. Points with relative contributions
#' (in percentage) below this value  are not plotted.
#' @param posleg position of the legend.
#' @param xlim a numeric vectors of length 2, giving the x coordinates range. If NULL (by default) 
#' the range is defined automatically (recommended). 
#' @param ylim a numeric vectors of length 2, giving the y coordinates range. 
#' If NULL (by default) the range is defined automatically (recommended).
#' @param main a string corresponding to the title of the graph to draw.
#' @param cex cf. function \code{par} in the \bold{graphics} package
#' @param leg boolean, if TRUE, a legend is displayed..
#' @param cex.leg a numerical value giving the amount by which the legend should 
#' be magnified. Default is 0.8.
#' @param posleg.sup position of the legend for the supplementary groups.
#' @param col.groups.sup a vector of colors, of size the number of supplementary groups. 
#' If NULL, colors are chosen automatically.
#' @param nb.paxes.sup if choice="axes", the maximum number of partial axes 
#' of supplementary groups ploted on the correlation circle. By default equal to 3.
#' @param \ldots arguments to be passed to methods, such as graphical parameters.
#' @seealso \code{\link{summary.PCAmix}},\code{\link{PCAmix}},\code{\link{PCArot}}
#' @author , \email{marie.chavent@u-bordeaux.fr}, Amaury Labenne
#' @details The observations can be colored according to the levels of a qualitative
#' variable. The observations, the quantitative variables and the levels can be 
#' selected according to their squared cosine (lim.cos2.plot) or their relative 
#' contribution (lim.contrib.plot) to the component map. Only points with squared 
#' cosine or relative contribution greater than a given threshold are plotted. 
#' Note that the relative contribution of a point to the component map (a plan) 
#' is the sum of the absolute contributions to each dimension, divided by the 
#' sum of the corresponding eigenvalues.
#' @references 
#' Chavent M., Kuentz-Simonet V., Labenne A., Saracco J., Multivariate analysis of mixed data: The PCAmixdata R package, arXiv:1411.4911 [stat.CO].
#' @examples
#' data(gironde)
#' class.var<-c(rep(1,9),rep(2,5),rep(3,9),rep(4,4))
#' names <- c("employment","housing","services","environment")
#' dat <- cbind(gironde$employment[1:20,],gironde$housing[1:20,],
#'            gironde$services[1:20,],gironde$environment[1:20,])
#' res <- MFAmix(data=dat,groups=class.var,
#'             name.groups=names, rename.level=TRUE, ndim=3,graph=FALSE)
#' 
#' #---- quantitative variables
#' plot(res,choice="cor",cex=0.6)
#' plot(res,choice="cor",cex=0.6,coloring.var="groups")
#' plot(res,choice="cor",cex=0.6,coloring.var="groups",
#'      col.groups=c("red","yellow","pink","brown"),leg=TRUE)
#' 
#' #----partial axes
#' plot(res,choice="axes",cex=0.6)
#' plot(res,choice="axes",cex=0.6,coloring.var="groups")
#' plot(res,choice="axes",cex=0.6,coloring.var="groups",
#'      col.groups=c("red","yellow","pink","brown"),leg=TRUE)
#' 
#' #----groups
#' plot(res,choice="groups",cex=0.6)   #no colors for groups
#' plot(res,choice="groups",cex=0.6,coloring.var="groups") 
#' plot(res,choice="groups",cex=0.6,coloring.var="groups",
#'      col.groups=c("red","yellow","pink","blue")) 

#' #----squared loadings
#' plot(res,choice="sqload",cex=0.8)    #no colors for groups
#' plot(res,choice="sqload",cex=0.8,coloring.var="groups",
#'      posleg="topright") 
#' plot(res,choice="sqload",cex=0.6,coloring.var="groups",
#'      col.groups=c("red","yellow","pink","blue"),ylim=c(0,1)) 

#' plot(res,choice="sqload",cex=0.8,coloring.var="type",
#'      cex.leg=0.9,posleg="topright")  
#' 
#' #----individuals 
#' plot(res,choice="ind",cex=0.6) 
#' 
#' #----individuals with squared cosine greater than 0.5
#' plot(res,choice="ind",cex=0.6,lim.cos2.plot=0.5)  
#' 
#' #----individuals colored with a qualitative variable
#' nbchem <- gironde$services$chemist[1:20]
#' plot(res,choice="ind",cex=0.6,coloring.ind=nbchem,
#'      posleg="topright")   
#' plot(res,choice="ind",coloring.ind=nbchem,
#'      col.ind=c("pink","brown","darkblue"),label=FALSE,posleg="topright")     
#' 
#' #----partial individuals colored by groups
#' plot(res,choice="ind",partial=c("AUBIAC","ARCACHON"),
#'     cex=0.6,posleg="bottomright")
#' 
#' #----levels of qualitative variables
#' plot(res,choice="levels",cex=0.8)
#' plot(res,choice="levels",cex=0.8,coloring.var="groups")
#' 
#' #levels with squared cosine greater than 0.6
#' plot(res,choice="levels",cex=0.8, lim.cos2.plot=0.6)
#' 
#' #supplementary groups
#' data(wine)
#' X.quanti <- splitmix(wine)$X.quanti[,1:5]
#' X.quali <- splitmix(wine)$X.quali[,1,drop=FALSE]
#' X.quanti.sup <- splitmix(wine)$X.quanti[,28:29]
#' X.quali.sup <- splitmix(wine)$X.quali[,2,drop=FALSE]
#' data <- cbind(X.quanti,X.quali)
#' data.sup <- cbind(X.quanti.sup,X.quali.sup)
#' 
#' groups <-c(1,2,2,3,3,1)
#' name.groups <- c("G1","G2","G3")
#' groups.sup <- c(1,1,2)
#' name.groups.sup <- c("Gsup1","Gsup2")
#' mfa <- MFAmix(data,groups,name.groups,ndim=4,rename.level=TRUE,graph=FALSE)
#' mfa.sup <- supvar(mfa,data.sup,groups.sup,name.groups.sup,rename.level=TRUE)
#' plot(mfa.sup,choice="sqload",coloring.var="groups")
#' plot(mfa.sup,choice="axes",coloring.var="groups")
#' plot(mfa.sup,choice="groups",coloring.var="groups")
#' plot(mfa.sup,choice="levels",coloring.var="groups")
#' plot(mfa.sup,choice="levels")
#' plot(mfa.sup,choice="cor",coloring.var = "groups")
#' 
plot.MFAmix <- function(x, axes = c(1, 2), choice = "ind", label=TRUE, coloring.var = "not", coloring.ind=NULL, nb.partial.axes=3,
                       col.ind=NULL, col.groups=NULL, partial = NULL, lim.cos2.plot = 0, lim.contrib.plot=0, xlim = NULL,  ylim = NULL,
                       cex = 1, main = NULL, leg=TRUE,posleg="topleft",cex.leg=0.8, 
                       col.groups.sup=NULL,posleg.sup="topright",nb.paxes.sup=3,...) 
{
  
  cl<-match.call()
  if (!inherits(x, "MFAmix")) 
    stop("use only with \"MFAmix\" objects")
  
  n <- nrow(x$ind$coord)
  if (is.null(x$sqload.sup)) sup <- FALSE else sup <- TRUE
  
  if (!(choice %in% c("ind", "sqload", "levels", "cor", "axes", "groups"))) 
    stop("\"choice\" must be either \"ind\",\"sqload\",\"cor\", \"levels\",\"axes\" or \"groups\"",call.=FALSE)
  if ((choice=="levels") & is.null(x$levels) & is.null(x$levels.sup))
    stop("\"choice=levels\" is not possible with pure quantitative data",call. = FALSE)
  if ((choice=="cor") & is.null(x$quanti) & is.null(x$quanti.sup))
    stop("\"choice=cor\" is not possible with pure qualitative data",call. = FALSE)
  
  if (!is.logical(leg))
    stop("argument \"leg\" must be TRUE or FALSE",call. = FALSE)
  
  if (lim.cos2.plot != 0 & lim.contrib.plot!=0)
    stop("use either \"lim.cos2.plot\" OR \"lim.contrib.plot\"",call. = FALSE)
  
  if (!is.null(partial))
    if (!is.character(partial) | length(which(rownames(x$ind$coord) %in% partial))==0)
      stop("invalid values in \"partial\"",call. = FALSE)
  
  if (!is.null(coloring.ind))
  {
    if (choice!="ind")
      warning("use \"coloring.ind\" only if choice=\"ind\"")
    if (!is.null(partial)) 
      warning("use \"coloring.ind\" only if partial=NULL")
  }
  
  if (!is.null(coloring.ind))
  {
    if (!is.factor(coloring.ind) | length(coloring.ind)!=n)
    {
      warning("\"coloring.ind\" must be either NULL or a qualitative variable of class factor. Its length must be equal to the number of individuals")
      coloring.ind=NULL
    }
  } 
  
  if (coloring.var!="groups" &  coloring.var!="type" & coloring.var!="not")
    stop("'coloring.var' must be one of the following: 'not', 'type', 'groups'",call. = FALSE)
  
  if (coloring.var=="type")
    if (choice=="ind" | choice=="cor" | choice=="levels"| choice=="axes")
    {
      warning("coloring.var=\"type\" is not used if choice=\"ind\", \"cor\",\"levels\" or \"axes\"",call. = FALSE)
      coloring.var <- "not"
    }

  if (coloring.var=="groups")
    if (choice=="ind") {
      warning("\"coloring.var\" is not used if choice=\"ind\"")
      coloring.var <- "not"
    }
  
 
  eig.axes <- x$eig[axes,1]
  
  dim1 <- axes[1]
  dim2 <- axes[2]
  lab.x <- paste("Dim ", dim1, " (", signif(x$eig[axes[1],2], 4), " %)", sep = "")
  lab.y <- paste("Dim ",dim2, " (", signif(x$eig[axes[2],2], 4), " %)", sep = "")
  
  #group <- x$lst.groups
  ngroup <- nrow(x$groups$contrib)
  name.groups <- names(x$partial.axes)
  name.groups.sup <- names(x$partial.axes.sup)
  if (sup) ngroupsup <- nrow(x$group.sup)
  
  if (is.null(col.groups)) 
    col.groups <- 2:(ngroup+1) else
      if (length(col.groups)!=ngroup)
      {
        warning("invalid length of \"col.groups\"")
        col.groups <- 2:(ngroup+1)
      }
  if (sup)
    if (is.null(col.groups.sup)) 
      col.groups.sup <- (ngroup+2):(ngroup+ngroupsup+1) else
        if (length(col.groups.sup)!=ngroupsup)
        {
          warning("invalid length of \"col.groups.sup\"")
          col.groups.sup <- (ngroup+2):(ngroup+ngroupsup+1)
        }
 
  p1 <- x$global.pca$rec$p1
  p <- x$global.pca$rec$p
  p2<-x$global.pca$rec$p2
  m <- ncol(x$global.pca$rec$W)-p1
 
  if (sup)
  {
    p.sup <- x$rec.sup$p
    p1.sup <- x$rec.sup$p1
    p2.sup <- x$rec.sup$p2
  }
 
  
  # plot of the partial axes on a correlation circle
  
  if (choice == "axes") 
  {
    if (is.null(main)) main <- "Partial axes"
    if (is.null(xlim)) xlim <- c(-1.1, 1.1)
    if (is.null(ylim)) ylim <- c(-1.1, 1.1)
    
    graphics::plot(0, 0, xlab = lab.x, ylab = lab.y, xlim = xlim, ylim = ylim, col = "white", asp = 1, 
         cex = cex, main = main,...)
    x.cercle <- seq(-1, 1, by = 0.01)
    y.cercle <- sqrt(1 - x.cercle^2)
    graphics::lines(x.cercle, y = y.cercle)
    graphics::lines(x.cercle, y = -y.cercle)
    graphics::abline(v = 0, lty = 2, cex = cex)
    graphics::abline(h = 0, lty = 2, cex = cex)
    
    coord.paxes <- NULL
    col.paxes <- NULL
    for (i in 1:ngroup)
    {
      nmax <- min(nrow(x$partial.axes[[i]]),nb.partial.axes)
      coord.paxes <- rbind(coord.paxes, x$partial.axes[[i]][1:nmax,c(dim1,dim2)])
      col.paxes <- c(col.paxes,rep(col.groups[i],nmax))
    }
    if (coloring.var != "groups") col.paxes <- rep("black",nrow(coord.paxes))
    
    for (v in 1:nrow(coord.paxes)) 
    {
      graphics::arrows(0, 0, coord.paxes[v, 1], coord.paxes[v, 2], length = 0.1, angle = 15, code = 2, col = col.paxes[v], cex = cex)
      if (abs(coord.paxes[v, 1]) > abs(coord.paxes[v, 2])) {
        if (coord.paxes[v, 1] >= 0) 
          pos <- 4
        else pos <- 2
      }
      else {
        if (coord.paxes[v, 2] >= 0) 
          pos <- 3
        else pos <- 1
      }
      graphics::text(coord.paxes[v, 1], y = coord.paxes[v, 2], labels = rownames(coord.paxes)[v], 
           pos = pos, col = col.paxes[v], cex = cex,...)
    }
       if ((coloring.var == "groups") & (leg==TRUE)) 
      graphics::legend((posleg), legend = name.groups, text.col = col.groups, cex = cex.leg)
    if (sup)
    {
      coord.paxes <- NULL
      col.paxes <- NULL
      for (i in 1:ngroupsup)
      {
        nmax <- min(nrow(x$partial.axes.sup[[i]]),nb.paxes.sup)
        coord.paxes <- rbind(coord.paxes, x$partial.axes.sup[[i]][1:nmax,c(dim1,dim2)])
        col.paxes <- c(col.paxes,rep(col.groups.sup[i],nmax))
      }
      if (coloring.var != "groups") col.paxes <- rep("black",nrow(coord.paxes))
      
      for (v in 1:nrow(coord.paxes)) 
      {
        graphics::arrows(0, 0, coord.paxes[v, 1], coord.paxes[v, 2], length = 0.1, angle = 15, lty=5, code = 2, col = col.paxes[v], cex = cex)
        if (abs(coord.paxes[v, 1]) > abs(coord.paxes[v, 2])) {
          if (coord.paxes[v, 1] >= 0) 
            pos <- 4
          else pos <- 2
        }
        else {
          if (coord.paxes[v, 2] >= 0) 
            pos <- 3
          else pos <- 1
        }
        graphics::text(coord.paxes[v, 1], y = coord.paxes[v, 2], labels = rownames(coord.paxes)[v], 
             pos = pos, col = col.paxes[v], cex = cex,...)
      }
      if ((coloring.var == "groups") & (leg==TRUE)) 
        graphics::legend((posleg.sup), legend = name.groups.sup, text.col = col.groups.sup, cex = cex.leg)
    }
    
  }   
  
  #plot of the groups according to their contribution
  if (choice == "groups") 
  {
    if (is.null(main)) main <- "Groups contributions"
    coord.groups <- x$groups$contrib[, axes, drop = FALSE]
    
    xmax <- max(coord.groups[,1],x$groups.sup[,dim1], xlim)
    xlim <- c(0, xmax * 1.2)
    
    ymax <- max(coord.groups[,2],x$groups.sup[,dim1],ylim)
    ylim <- c(0, ymax * 1.2)
    

    if (coloring.var != "groups") 
    {
      col.groups = rep("darkred", nrow(coord.groups))
      if (sup) col.groups.sup <- rep("blue", ngroupsup)
    }
     
    
    graphics::plot(coord.groups, xlab = lab.x, ylab = lab.y, xlim = xlim, ylim = ylim, pch = 17, col = col.groups, 
         cex = cex, main = main, cex.main = cex * 1.2, asp = 1,...)
    graphics::abline(v = 0, lty = 2,cex=cex)
    graphics::abline(h = 0, lty = 2,cex=cex)
    
    if (label) 
      graphics::text(coord.groups[, 1], y = coord.groups[, 2], labels = rownames(coord.groups), 
           pos = 3, col = col.groups, cex = cex)
    if (sup)
    {
      graphics::points(x$group.sup[,axes,drop=FALSE], xlab = lab.x, ylab = lab.y, xlim = xlim, ylim = ylim, pch = 2, col = col.groups.sup, 
           cex = cex, main= main, cex.main = cex * 1.2, asp = 1,...)
      if (label) 
        graphics::text(x=x$group.sup[,dim1], y=x$group.sup[,dim2],  labels = rownames(x$group.sup), 
             pos = 3, col = col.groups.sup, cex = cex,...)
    }
  }
  
  # plot of the variables according to their squared loadings
  if (choice=="sqload") 
  {
    if (is.null(main)) main <- "Squared loadings"
    xmax <- max(x$sqload[, dim1],x$sqload.sup[, dim1],xlim)
    xlim <- c(-0.1, xmax * 1.2)
    
    ymax <- max(x$sqload[, dim2],x$sqload.sup[, dim2],ylim)
    ylim <- c(-0.1, ymax * 1.2)
    
    graphics::plot(0, 0, type="n",xlab = lab.x, ylab = lab.y, xlim = xlim, ylim = ylim,cex=cex,main=main,...)
    graphics::abline(v = 0, lty = 2,cex=cex)
    graphics::abline(h = 0, lty = 2,cex=cex)
    
    col.var <- rep(1,p) #color of a point-variable
    if (coloring.var == "groups") 
    {
      for (i in 1:ngroup)
        col.var[which(x$index.group==i)] <- col.groups[i]
    }
    if (coloring.var=="type")
    {
      if (p1 >0) col.var[1:p1] <- "blue"
      if (p2 >0) col.var[(p1+1):p] <- "red"
    }
    
    for (j in 1:p)
    {
      graphics::arrows(0,0,x$sqload[j,dim1],x$sqload[j,dim2], 
             length = 0.1, angle = 15, code = 2,cex=cex,col=col.var[j])
      if (label)
      {
        if (x$sqload[j,dim1] > x$sqload[j,dim1]) 
        { 
          pos <- 4
        } 
        else  pos <- 3
        graphics::text(x$sqload[j,dim1],x$sqload[j,dim2],
             labels = rownames(x$sqload)[j], pos = pos,cex=cex,col=col.var[j],...)  
      }
    }
    
    if ((coloring.var == "groups") & (leg==TRUE)) 
      graphics::legend((posleg), legend = name.groups, text.col = col.groups, cex = cex.leg)
    if (coloring.var=="type" & (leg==TRUE)) 
      graphics::legend(posleg, legend = c("numerical","categorical"), text.col = c("blue","red"), cex=cex.leg)
    
    if (sup)
    {
      col.var.sup <- rep(4,)
      #color of a supp point-variable
      if (coloring.var == "groups") 
        for (i in 1:ngroupsup)
          col.var.sup[which(x$index.groupsup==i)] <- col.groups.sup[i]
      if (coloring.var=="type")
      {
        if (p1.sup >0) col.var.sup[1:p1.sup] <- "blue"
        if (p2.sup >0) col.var.sup[(p1.sup+1):p.sup] <- "red"
      }
      for (j in 1:nrow(x$sqload.sup)) 
      {
        graphics::arrows(0, 0, x$sqload.sup[j, dim1], x$sqload.sup[j, dim2], 
               length = 0.1, angle = 15, code = 2, lty=5, col=col.var.sup[j],cex = cex,...)
        if (label) 
        {
          if (x$sqload.sup[j, dim1] > x$sqload.sup[j, dim2]) 
          {
            pos <- 4
          } else pos <- 3
          graphics::text(x$sqload.sup[j, dim1], x$sqload.sup[j, dim2], labels = rownames(x$sqload.sup)[j], 
               pos = pos, cex = cex,col=col.var.sup[j],...)
        }
      }
      if ((coloring.var == "groups") & (leg==TRUE)) 
        graphics::legend((posleg.sup), legend = name.groups.sup, text.col = col.groups.sup, cex = cex.leg)
    }
    
     }    
  
  #plot of the quantitative variables on a correlation circle
  if (choice == "cor") 
  {
    if (is.null(main))  main <- "Correlation circle"
    if (is.null(xlim)) xlim = c(-1.1, 1.1)
    if (is.null(ylim)) ylim = c(-1.1, 1.1)
    
    graphics::plot(0, 0, main = main, xlab = lab.x, ylab = lab.y, 
         xlim = xlim, ylim = ylim, col = "white", 
         asp = 1, cex = cex,...)
    x.cercle <- seq(-1, 1, by = 0.01)
    y.cercle <- sqrt(1 - x.cercle^2)
    graphics::lines(x.cercle, y = y.cercle)
    graphics::lines(x.cercle, y = -y.cercle)
    graphics::abline(v = 0, lty = 2, cex = cex)
    graphics::abline(h = 0, lty = 2, cex = cex)
    
    #plot active quantitative variables
    if (!is.null(x$quanti))
    {
      if (lim.cos2.plot == 0 & lim.contrib.plot==0)
      {
        lim.plot<-0
        base.lim<-x$quanti$cos2[,axes]
      }
      
      if (lim.cos2.plot != 0)
      {
        lim.plot<-lim.cos2.plot
        base.lim<-x$quanti$cos2[,axes]
      }
      
      if(lim.contrib.plot!=0)
      {
        lim.plot<-lim.contrib.plot
        base.lim<-x$quanti$contrib[,axes]
        base.lim<-100*(base.lim/sum(eig.axes))     
      }
    }
    
    coord.var <- x$quanti$coord[, axes, drop = FALSE]
    
    col.var <- rep(1,p) #color of a the point-variable (all quantitative)
    if (coloring.var == "groups") 
    {
      for (i in 1:ngroup)
        col.var[which(x$index.group==i)] <- col.groups[i]
    }
    col.var <- col.var[1:p1]
    
    test.empty.plot<-c()  
    
    for (v in 1:nrow(coord.var)) 
    {
      if (sum(base.lim[v, ] , na.rm = TRUE) >= lim.plot && !is.na(sum(base.lim[v, ], na.rm = TRUE))) 
      {
        test.empty.plot<-c(test.empty.plot,1)
        graphics::arrows(0, 0, coord.var[v, 1], coord.var[v,2], length = 0.1, angle = 15, code = 2,cex = cex,col=col.var[v])
        
        if (label) 
        {
          if (abs(coord.var[v, 1]) > abs(coord.var[v, 2])) 
          {
            if (coord.var[v, 1] >= 0) 
              pos <- 4
            else pos <- 2
          }
          else {
            if (coord.var[v, 2] >= 0) 
              pos <- 3
            else pos <- 1
          }
          graphics::text(coord.var[v, 1], y = coord.var[v, 2], 
               labels = rownames(coord.var)[v], pos = pos, cex = cex,col=col.var[v])
        }
      }
    }
    
    if(is.null(test.empty.plot))
      warning("\"lim.cos.plot\" (or \"lim.contrib.plot\") is too large. No variable can be plotted")
    if ((coloring.var == "groups") & (leg==TRUE)) 
      graphics::legend(posleg, legend = name.groups[unique(x$index.group[1:p1])], text.col = unique(col.var), cex = cex.leg)
    
    if (!is.null(x$quanti.sup))
    {
      coord.var.sup <- x$quanti.sup$coord[, axes, drop = FALSE]
      
      col.var.sup <- rep(1,p1.sup) #color of a the point-variable (all quantitative)
      if (coloring.var == "groups") 
      {
        for (i in 1:ngroupsup)
          col.var.sup[which(x$index.groupsup==i)] <- col.groups.sup[i]
      }
      col.var.sup <- col.var.sup[1:p1.sup]
      for (v in 1:nrow(coord.var.sup)) 
      {
        graphics::arrows(0, 0, coord.var.sup[v, 1], coord.var.sup[v,2], length = 0.1, 
               angle = 15, code = 2,cex = cex,col= col.var.sup[v],lty=5)
        
        if (label) 
        {
          if (abs(coord.var.sup[v, 1]) > abs(coord.var.sup[v,2])) 
          {
            if (coord.var.sup[v, 1] >= 0) 
              pos <- 4
            else pos <- 2
          }
          else {
            if (coord.var.sup[v, 2] >= 0) 
              pos <- 3
            else pos <- 1
          }
          graphics::text(coord.var.sup[v, 1], y = coord.var.sup[v, 2], 
               labels = rownames(coord.var.sup)[v], pos = pos, cex = cex,col=col.var.sup[v])
        }
      }
      if ((coloring.var == "groups") & (leg==TRUE)) 
        graphics::legend(posleg.sup, legend = name.groups.sup[unique(x$index.groupsup[1:p1.sup])], 
                         text.col = unique(col.var.sup), cex = cex.leg)
    }
    
  }
  
  # plot of the individuals
  if (choice == "ind") 
  {
    if (is.null(main)) 
      main <- "Individuals component map"
    
    coord.ind <- x$ind$coord
  
    if (lim.cos2.plot == 0 & lim.contrib.plot==0)
    {
      lim.plot<-0
      select.ind <- 1:nrow(coord.ind)
    }
    
    if (lim.cos2.plot != 0 & lim.contrib.plot==0)
    {
      lim.plot <- lim.cos2.plot
      base.lim <- x$ind$cos2[,axes]
      select.ind <- which(apply(base.lim[,],1,sum)>=lim.plot)    
    }
    
    if (lim.cos2.plot == 0 & lim.contrib.plot!=0)
    {
      lim.plot <- lim.contrib.plot
      base.lim <- x$ind$contrib[,axes]
      base.lim <- 100*(base.lim/sum(eig.axes))
      select.ind <- which(apply(base.lim[,],1,sum)>=lim.plot)
    }
    
    
    if (is.null(partial))
    {
      if (length(select.ind)==0)
        stop("\"lim.cos.plot\" (or \"lim.contrib.plot\") is too large. No individuals can be plotted",call. = FALSE)
      
      coord.ind <- coord.ind[select.ind, , drop=FALSE]
  
      xmin <- min(xlim,coord.ind[, dim1])
      xmax <- max(xlim,coord.ind[, dim1])
      xlim <- c(xmin, xmax) * 1.2
      
      ymin <- min(ylim,coord.ind[, dim2])
      ymax <- max(ylim,coord.ind[, dim2])
      ylim <- c(ymin, ymax) * 1.2
      
      if (is.null(col.ind) | is.null(coloring.ind))
      {
        col.plot.ind <- rep("black",nrow(coord.ind))
      }
      
      if (is.factor(coloring.ind))
      { 
        quali<-coloring.ind
        if (!is.null(col.ind))
        { 
          levels(quali) <- col.ind
          col.plot.ind <- quali
        }
        if (is.null(col.ind))
          col.plot.ind <- as.numeric(quali)
      }
      col.plot.ind.total<-col.plot.ind
      col.plot.ind <- col.plot.ind[select.ind]
      
      graphics::plot(coord.ind[, axes,drop=FALSE], xlim = xlim, ylim = ylim, xlab = lab.x, 
           ylab = lab.y, pch = 20, col = as.character(col.plot.ind), 
           cex = cex, main=main,...)
      graphics::abline(h = 0, lty = 2, cex = cex)
      graphics::abline(v = 0, lty = 2, cex = cex)
      
      if (leg==TRUE & is.factor(coloring.ind))
        graphics::legend(posleg, legend =paste(cl["coloring.ind"],levels(coloring.ind),sep="="), text.col = levels(as.factor(col.plot.ind.total)), 
               cex =cex.leg)
      
      if (label) 
        graphics::text(coord.ind[, axes], labels = rownames(coord.ind), 
             pos = 3, col = as.character(col.plot.ind), cex = cex,...)
      
    }
    
    if (!is.null(partial))
    {
      select.partial <- which(rownames(coord.ind[select.ind,,drop=FALSE]) %in% partial)
      if (length(select.partial)==0)
        stop("\"lim.cos.plot\" (or \"lim.contrib.plot\") is too large. No partial individuals can be plotted",call. = FALSE)
      
      coord.ind.part <- coord.ind[select.partial, , drop=FALSE]
      xmin <- min(xlim,coord.ind[, dim1])
      xmax <- max(xlim,coord.ind[, dim1])
      ymin <- min(ylim,coord.ind[, dim2])
      ymax <- max(ylim,coord.ind[, dim2])
      
      for (i in 1:ngroup)
      { 
        t <- x$ind.partial[[i]][select.partial,axes,drop=FALSE]
        xmin <- min(xmin, t[,1])
        xmax <- max(xmax, t[,1])
        ymin <- min(ymin,t[,2])
        ymax <- max(ymax, t[,2])
      }
      
      xlim <- c(xmin, xmax) * 1.2
      ylim <- c(ymin, ymax) * 1.2
      
      col.plot.ind <- rep("black",nrow(coord.ind))
      
      graphics::plot(as.matrix(coord.ind[,axes]), xlim = xlim, ylim = ylim, xlab = lab.x, 
           ylab = lab.y, pch = 20, cex = cex,main=main,...)
      graphics::abline(h = 0, lty = 2, cex = cex)
      graphics::abline(v = 0, lty = 2, cex = cex)
      
      if (label) 
        graphics::text(coord.ind.part[, axes,drop=FALSE], labels = rownames(coord.ind.part), 
                       pos = 3, col = as.character(col.plot.ind), cex = cex)
      
      for (i in 1:ngroup)
      {
        t <- x$ind.partial[[i]][select.partial,axes,drop=FALSE]
        graphics::points(t,col=col.groups[i],pch=20,...)
        for (j in 1:length(select.partial))
        {
          m <- list(x=c(coord.ind.part[j,dim1],t[j,1]), y=c(coord.ind.part[j,dim2],t[j,2]))
          graphics::lines(m,col=col.groups[i]) 
        }
      }
      
      if(leg==TRUE)
        graphics::legend(posleg, legend = name.groups, text.col = col.groups, cex = cex.leg)
    }
  }
  
  #plot of the levels
  if (choice == "levels") 
    {
    if (is.null(main))  main <- "Levels component map"
    xmin <- min(xlim,x$levels$coord[, dim1],x$levels.sup$coord[, dim1])
    xmax <- max(xlim,x$levels$coord[, dim1],x$levels.sup$coord[, dim1])
    xlim <- c(xmin, xmax) * 1.2
    
    ymin <- min(ylim,x$levels$coord[, dim2],x$levels.sup$coord[, dim2])
    ymax <- max(ylim,x$levels$coord[, dim2],x$levels.sup$coord[, dim2])
    ylim <- c(ymin, ymax) * 1.2
    
    graphics::plot(0,0, xlim = xlim, ylim = ylim,
         xlab = lab.x, ylab = lab.y, type="n", cex = cex,main=main, ...)
    graphics::abline(h = 0, lty = 2, cex = cex)
    graphics::abline(v = 0, lty = 2, cex = cex)
    
    #plot levels of active variables
    if (!is.null(x$levels))
    {
      if (lim.cos2.plot == 0 & lim.contrib.plot==0)
      {
        lim.plot<-0
        base.lim<-x$levels$cos2[,axes]
      }
      
      if (lim.cos2.plot != 0)
      {
        lim.plot<-lim.cos2.plot
        base.lim<-x$levels$cos2[,axes]
      }
      
      if (lim.contrib.plot!=0)
      {
        lim.plot<-lim.contrib.plot
        base.lim<-x$levels$contrib[,axes]
        base.lim<-100*(base.lim/sum(eig.axes))    
      }
      
      coord.lev <- x$levels$coord[, axes, drop = FALSE]
      
      col.lev <- rep(1,p1+m) #color of a the point-level 
      if (coloring.var == "groups") 
      {
        for (i in 1:ngroup)
          col.lev[which(x$index.group2==i)] <- col.groups[i]
      }
      col.lev <- col.lev[(p1+1):(p1+m)]
      
      test.empty.plot<-c()
      for (v in 1:nrow(coord.lev)) 
      {
        if (sum(base.lim[v, ], na.rm = TRUE) >= lim.plot && !is.na(sum(base.lim[v, ], na.rm = TRUE))) {
          test.empty.plot<-c(test.empty.plot,1)
          graphics::points(coord.lev[v, 1], coord.lev[v,2], col=col.lev[v], pch=20,cex = cex,...)
          
          if (label) 
          {
            if (abs(coord.lev[v, 1]) > abs(coord.lev[v,2])) 
            {
              if (coord.lev[v, 1] >= 0) 
                pos <- 4
              else pos <- 2
            }
            else {
              if (coord.lev[v, 2] >= 0) 
                pos <- 3
              else pos <- 1
            }
            graphics::text(coord.lev[v, 1], y = coord.lev[v, 2], col=col.lev[v],
                 labels = rownames(coord.lev)[v], pos = pos, cex = cex)
          }
        }
      }
      if (is.null(test.empty.plot))
        warning("\"lim.cos.plot\" (or \"lim.contrib.plot\") is too large. No level can be plotted")
      if ((coloring.var == "groups") & (leg==TRUE)) 
        graphics::legend(posleg, legend = name.groups[unique(x$index.group2[(p1+1):(p1+m)])], 
               text.col = unique(col.lev), cex = cex.leg)
    }
    #plot levels of supplementary variables
    if (!is.null(x$levels.sup))
    {
      coord.lev.sup <- x$levels.sup$coord[, axes, drop = FALSE]
      m.sup <- nrow(coord.lev.sup)
      col.lev.sup <- rep(4,p1.sup+m.sup) #color of a the point-level 
      if (coloring.var == "groups") 
      {
        for (i in 1:ngroupsup)
          col.lev.sup[which(x$index.groupsup2==i)] <- col.groups.sup[i]
      }
      col.lev.sup <- col.lev.sup[(p1.sup+1):(p1.sup+m.sup)]
      for (v in 1:nrow(coord.lev.sup)) 
      {
        
        graphics::points(coord.lev.sup[v, 1], coord.lev.sup[v,2], pch=1,cex = cex,
               col=col.lev.sup[v],...)
        
        if (label) 
        {
          if (abs(coord.lev.sup[v, 1]) > abs(coord.lev.sup[v,2])) 
          {
            if (coord.lev.sup[v, 1] >= 0) 
              pos <- 4
            else pos <- 2
          }
          else {
            if (coord.lev.sup[v, 2] >= 0) 
              pos <- 3
            else pos <- 1
          }
          graphics::text(coord.lev.sup[v, 1], y = coord.lev.sup[v, 2], 
               labels = rownames(coord.lev.sup)[v], pos = pos, 
               col=col.lev.sup[v], cex = cex)
        }
      }
      if ((coloring.var == "groups") & (leg==TRUE)) 
        graphics::legend(posleg.sup, legend = name.groups[unique(x$index.groupsup2[(p1+1):(p1+m)])], 
               text.col = unique(col.lev.sup), cex = cex.leg)
    }
  }
}
