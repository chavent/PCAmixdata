##' @export
##' @method plot MFAmix

plot.MFAmix<-function (x, axes = c(1, 2), choice = "axes", label=TRUE, coloring.var = "not", coloring.ind=NULL, nb.partial.axes=3,
                       col.ind=NULL, col.groups=NULL, partial = NULL, lim.cos2.plot = 0, lim.contrib.plot=0, xlim = NULL,  ylim = NULL,
                       cex = 1, main = NULL, new.plot = FALSE,leg=TRUE,posleg="topleft",cex.leg=0.8, 
                       col.groups.sup=NULL,posleg.sup="topright",nb.paxes.sup=3,...) 
{
  
  cl<-match.call()
  if (!inherits(x, "MFAmix")) 
    stop("use only with \"MFAmix\" objects")
  
  res.mfa <- x
  n <- nrow(res.mfa$ind$coord)
  if (is.null(res.mfa$sqload.sup)) sup <- FALSE else sup <- TRUE
  
  if (!(choice %in% c("ind", "sqload", "levels", "cor", "axes", "groups"))) 
    stop("\"choice\" must be either \"ind\",\"sqload\",\"cor\", \"levels\",\"axes\" or \"groups\"",call.=FALSE)
  if ((choice=="levels") & is.null(res.mfa$levels) & is.null(res.mfa$levels.sup))
    stop("\"choice=levels\" is not possible with pure quantitative data",call. = FALSE)
  if ((choice=="cor") & is.null(res.mfa$quanti) & is.null(res.mfa$quanti.sup))
    stop("\"choice=cor\" is not possible with pure qualitative data",call. = FALSE)
  
  if (!is.logical(leg))
    stop("argument \"leg\" must be TRUE or FALSE",call. = FALSE)
  
  if (lim.cos2.plot != 0 & lim.contrib.plot!=0)
    stop("use either \"lim.cos2.plot\" OR \"lim.contrib.plot\"",call. = FALSE)
  
  if (!is.null(partial))
    if (!is.character(partial) | length(which(rownames(res.mfa$ind$coord) %in% partial))==0)
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
  
 
  eig.axes <- res.mfa$eig[axes,1]
  
  dim1 <- axes[1]
  dim2 <- axes[2]
  lab.x <- paste("Dim ", dim1, " (", signif(res.mfa$eig[axes[1],2], 4), " %)", sep = "")
  lab.y <- paste("Dim ",dim2, " (", signif(res.mfa$eig[axes[2],2], 4), " %)", sep = "")
  
  #group <- res.mfa$lst.groups
  ngroup <- nrow(res.mfa$groups$contrib)
  name.groups <- names(res.mfa$partial.axes)
  name.groups.sup <- names(res.mfa$partial.axes.sup)
  if (sup) ngroupsup <- nrow(res.mfa$group.sup)
  
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
 
  p1 <- res.mfa$global.pca$rec$p1
  p <- res.mfa$global.pca$rec$p
  p2<-res.mfa$global.pca$rec$p2
  m <- ncol(res.mfa$global.pca$rec$W)-p1
 
  if (sup)
  {
    p.sup <- res.mfa$rec.sup$p
    p1.sup <- res.mfa$rec.sup$p1
    p2.sup <- res.mfa$rec.sup$p2
  }
 
  
  # plot of the partial axes on a correlation circle
  
  if (choice == "axes") 
  {
    if (new.plot) 
      dev.new()
    if (is.null(main)) main <- "Partial axes"
    if (is.null(xlim)) xlim <- c(-1.1, 1.1)
    if (is.null(ylim)) ylim <- c(-1.1, 1.1)
    
    plot(0, 0, xlab = lab.x, ylab = lab.y, xlim = xlim, ylim = ylim, col = "white", asp = 1, 
         cex = cex, main = main,...)
    x.cercle <- seq(-1, 1, by = 0.01)
    y.cercle <- sqrt(1 - x.cercle^2)
    lines(x.cercle, y = y.cercle)
    lines(x.cercle, y = -y.cercle)
    abline(v = 0, lty = 2, cex = cex)
    abline(h = 0, lty = 2, cex = cex)
    
    coord.paxes <- NULL
    col.paxes <- NULL
    for (i in 1:ngroup)
    {
      nmax <- min(nrow(res.mfa$partial.axes[[i]]),nb.partial.axes)
      coord.paxes <- rbind(coord.paxes, res.mfa$partial.axes[[i]][1:nmax,c(dim1,dim2)])
      col.paxes <- c(col.paxes,rep(col.groups[i],nmax))
    }
    if (coloring.var != "groups") col.paxes <- rep("black",nrow(coord.paxes))
    
    for (v in 1:nrow(coord.paxes)) 
    {
      arrows(0, 0, coord.paxes[v, 1], coord.paxes[v, 2], length = 0.1, angle = 15, code = 2, col = col.paxes[v], cex = cex)
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
      text(coord.paxes[v, 1], y = coord.paxes[v, 2], labels = rownames(coord.paxes)[v], 
           pos = pos, col = col.paxes[v], cex = cex)
    }
    
    if (sup)
    {
      coord.paxes <- NULL
      col.paxes <- NULL
      for (i in 1:ngroupsup)
      {
        nmax <- min(nrow(res.mfa$partial.axes.sup[[i]]),nb.paxes.sup)
        coord.paxes <- rbind(coord.paxes, res.mfa$partial.axes.sup[[i]][1:nmax,c(dim1,dim2)])
        col.paxes <- c(col.paxes,rep(col.groups.sup[i],nmax))
      }
      if (coloring.var != "groups") col.paxes <- rep("black",nrow(coord.paxes))
      
      for (v in 1:nrow(coord.paxes)) 
      {
        arrows(0, 0, coord.paxes[v, 1], coord.paxes[v, 2], length = 0.1, angle = 15, lty=5, code = 2, col = col.paxes[v], cex = cex)
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
        text(coord.paxes[v, 1], y = coord.paxes[v, 2], labels = rownames(coord.paxes)[v], 
             pos = pos, col = col.paxes[v], cex = cex)
      }
    }
    
  }   
  
  #plot of the groups according to their contribution
  if (choice == "groups") 
  {
    if (new.plot) dev.new()
    if (is.null(main)) main <- "Groups contributions"
    coord.groups <- res.mfa$groups$contrib[, axes, drop = FALSE]
    
    xmax <- max(coord.groups[,1],res.mfa$groups.sup[,dim1], xlim)
    xlim <- c(0, xmax * 1.2)
    
    ymax <- max(coord.groups[,2],res.mfa$groups.sup[,dim1],ylim)
    ylim <- c(0, ymax * 1.2)
    

    if (coloring.var != "groups") 
    {
      col.groups = rep("darkred", nrow(coord.groups))
      if (sup) col.groups.sup <- rep("blue", ngroupsup)
    }
     
    
    plot(coord.groups, xlab = lab.x, ylab = lab.y, xlim = xlim, ylim = ylim, pch = 17, col = col.groups, 
         cex = cex, main = main, cex.main = cex * 1.2, asp = 1,...)
    abline(v = 0, lty = 2,cex=cex)
    abline(h = 0, lty = 2,cex=cex)
    
    if (label) 
      text(coord.groups[, 1], y = coord.groups[, 2], labels = rownames(coord.groups), 
           pos = 3, col = col.groups, cex = cex)
    if (sup)
    {
      points(res.mfa$group.sup[,axes,drop=FALSE], xlab = lab.x, ylab = lab.y, xlim = xlim, ylim = ylim, pch = 2, col = col.groups.sup, 
           cex = cex, main= main, cex.main = cex * 1.2, asp = 1,...)
      if (label) 
        text(x=res.mfa$group.sup[,dim1], y=res.mfa$group.sup[,dim2],  labels = rownames(res.mfa$group.sup), 
             pos = 3, col = col.groups.sup, cex = cex)
    }
  }
  
  # plot of the variables according to their squared loadings
  if (choice=="sqload") 
  {
    if (is.null(main)) main <- "Squared loadings"
    if (new.plot) dev.new()
    xmax <- max(res.mfa$sqload[, dim1],res.mfa$sqload.sup[, dim1],xlim)
    xlim <- c(-0.1, xmax * 1.2)
    
    ymax <- max(res.mfa$sqload[, dim2],res.mfa$sqload.sup[, dim2],ylim)
    ylim <- c(-0.1, ymax * 1.2)
    
    plot(0, 0, type="n",xlab = lab.x, ylab = lab.y, xlim = xlim, ylim = ylim,cex=cex,main=main,...)
    abline(v = 0, lty = 2,cex=cex)
    abline(h = 0, lty = 2,cex=cex)
    
    col.var <- rep(1,p) #color of a point-variable
    if (coloring.var == "groups") 
    {
      for (i in 1:ngroup)
        col.var[which(res.mfa$index.group==i)] <- col.groups[i]
    }
    if (coloring.var=="type")
    {
      if (p1 >0) col.var[1:p1] <- "blue"
      if (p2 >0) col.var[(p1+1):p] <- "red"
    }
    
    for (j in 1:p)
    {
      arrows(0,0,res.mfa$sqload[j,dim1],res.mfa$sqload[j,dim2], 
             length = 0.1, angle = 15, code = 2,cex=cex,col=col.var[j])
      if (label)
      {
        if (res.mfa$sqload[j,dim1] > res.mfa$sqload[j,dim1]) 
        { 
          pos <- 4
        } 
        else  pos <- 3
        text(res.mfa$sqload[j,dim1],res.mfa$sqload[j,dim2],
             labels = rownames(res.mfa$sqload)[j], pos = pos,cex=cex,col=col.var[j])  
      }
    }
    
    if ((coloring.var == "groups") & (leg==TRUE)) 
      legend((posleg), legend = name.groups, text.col = col.groups, cex = cex.leg)
    if (coloring.var=="type" & (leg==TRUE)) 
      legend(posleg, legend = c("numerical","categorical"), text.col = c("blue","red"), cex=cex.leg)
    
    if (sup)
    {
      col.var.sup <- rep(4,)
      #color of a supp point-variable
      if (coloring.var == "groups") 
        for (i in 1:ngroupsup)
          col.var.sup[which(res.mfa$index.groupsup==i)] <- col.groups.sup[i]
      if (coloring.var=="type")
      {
        if (p1.sup >0) col.var.sup[1:p1.sup] <- "blue"
        if (p2.sup >0) col.var.sup[(p1.sup+1):p.sup] <- "red"
      }
      for (j in 1:nrow(res.mfa$sqload.sup)) 
      {
        arrows(0, 0, res.mfa$sqload.sup[j, dim1], res.mfa$sqload.sup[j, dim2], 
               length = 0.1, angle = 15, code = 2, lty=5, col=col.var.sup[j],cex = cex,...)
        if (label) 
        {
          if (res.mfa$sqload.sup[j, dim1] > res.mfa$sqload.sup[j, dim2]) 
          {
            pos <- 4
          } else pos <- 3
          text(res.mfa$sqload.sup[j, dim1], res.mfa$sqload.sup[j, dim2], labels = rownames(res.mfa$sqload.sup)[j], 
               pos = pos, cex = cex,col=col.var.sup[j],...)
        }
      }
      if ((coloring.var == "groups") & (leg==TRUE)) 
        legend((posleg.sup), legend = name.groups.sup, text.col = col.groups.sup, cex = cex.leg)
    }
    
     }    
  
  #plot of the quantitative variables on a correlation circle
  if (choice == "cor") 
  {
    if (is.null(main))  main <- "Correlation circle"
    if (new.plot) dev.new()
    
    if (is.null(xlim)) xlim = c(-1.1, 1.1)
    if (is.null(ylim)) ylim = c(-1.1, 1.1)
    
    plot(0, 0, main = main, xlab = lab.x, ylab = lab.y, 
         xlim = xlim, ylim = ylim, col = "white", 
         asp = 1, cex = cex,...)
    x.cercle <- seq(-1, 1, by = 0.01)
    y.cercle <- sqrt(1 - x.cercle^2)
    lines(x.cercle, y = y.cercle)
    lines(x.cercle, y = -y.cercle)
    abline(v = 0, lty = 2, cex = cex)
    abline(h = 0, lty = 2, cex = cex)
    
    #plot active quantitative variables
    if (!is.null(res.mfa$quanti))
    {
      if (lim.cos2.plot == 0 & lim.contrib.plot==0)
      {
        lim.plot<-0
        base.lim<-res.mfa$quanti$cos2[,axes]
      }
      
      if (lim.cos2.plot != 0)
      {
        lim.plot<-lim.cos2.plot
        base.lim<-res.mfa$quanti$cos2[,axes]
      }
      
      if(lim.contrib.plot!=0)
      {
        lim.plot<-lim.contrib.plot
        base.lim<-res.mfa$quanti$contrib[,axes]
        base.lim<-100*(base.lim/sum(eig.axes))     
      }
    }
    
    coord.var <- res.mfa$quanti$coord[, axes, drop = FALSE]
    
    col.var <- rep(1,p) #color of a the point-variable (all quantitative)
    if (coloring.var == "groups") 
    {
      for (i in 1:ngroup)
        col.var[which(res.mfa$index.group==i)] <- col.groups[i]
    }
    col.var <- col.var[1:p1]
    
    test.empty.plot<-c()  
    
    for (v in 1:nrow(coord.var)) 
    {
      if (sum(base.lim[v, ] , na.rm = TRUE) >= lim.plot && !is.na(sum(base.lim[v, ], na.rm = TRUE))) 
      {
        test.empty.plot<-c(test.empty.plot,1)
        arrows(0, 0, coord.var[v, 1], coord.var[v,2], length = 0.1, angle = 15, code = 2,cex = cex,col=col.var[v])
        
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
          text(coord.var[v, 1], y = coord.var[v, 2], 
               labels = rownames(coord.var)[v], pos = pos, cex = cex,col=col.var[v])
        }
      }
    }
    
    if(is.null(test.empty.plot))
      warning("\"lim.cos.plot\" (or \"lim.contrib.plot\") is too large. No variable can be plotted")
    if ((coloring.var == "groups") & (leg==TRUE)) 
      legend(posleg, legend = name.groups[unique(res.mfa$index.group[1:p1])], text.col = unique(col.var), cex = cex.leg)
    
    if (!is.null(res.mfa$quanti.sup))
    {
      coord.var.sup <- res.mfa$quanti.sup$coord[, axes, drop = FALSE]
      
      col.var.sup <- rep(1,p1.sup) #color of a the point-variable (all quantitative)
      if (coloring.var == "groups") 
      {
        for (i in 1:ngroupsup)
          col.var.sup[which(res.mfa$index.groupsup==i)] <- col.groups.sup[i]
      }
      col.var.sup <- col.var.sup[1:p1.sup]
      for (v in 1:nrow(coord.var.sup)) 
      {
        arrows(0, 0, coord.var.sup[v, 1], coord.var.sup[v,2], length = 0.1, 
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
          text(coord.var.sup[v, 1], y = coord.var.sup[v, 2], 
               labels = rownames(coord.var.sup)[v], pos = pos, cex = cex,col.var.sup[v])
        }
      }
      if ((coloring.var == "groups") & (leg==TRUE)) 
        legend(posleg.sup, legend = name.groups.sup[unique(res.mfa$index.groupsup[1:p1.sup])], text.col = unique(col.var.sup), cex = cex.leg)
    }
    
  }
  
  # plot of the individuals
  if (choice == "ind") 
  {
    if (is.null(main)) main <- "Individuals component map"
    if (new.plot) dev.new()
    
    coord.ind <- res.mfa$ind$coord
    
    if (lim.cos2.plot == 0 & lim.contrib.plot==0)
    {
      lim.plot<-0
      select.ind<-1:nrow(coord.ind)
    }
    
    if (lim.cos2.plot != 0 & lim.contrib.plot==0)
    {
      lim.plot <- lim.cos2.plot
      base.lim <- res.mfa$ind$cos2[,axes]
      select.ind <- which(apply(base.lim[,],1,sum)>=lim.plot)    
    }
    
    if (lim.cos2.plot == 0 & lim.contrib.plot!=0)
    {
      lim.plot <- lim.contrib.plot
      base.lim <- res.mfa$ind$contrib[,axes]
      base.lim <- 100*(base.lim/sum(eig.axes))
      select.ind <- which(apply(base.lim[,],1,sum)>=lim.plot)
    }
    
    
    if (is.null(partial))
    {
      if (length(select.ind)==0)
        stop("\"lim.cos.plot\" (or \"lim.contrib.plot\") is too large. No individuals can be plotted",call. = FALSE)
      
      coord.ind <- coord.ind[select.ind, , drop=FALSE]
  
      xmin <- min(xlim, coord.ind[, dim1])
      xmax <- max(xlim, coord.ind[, dim1])
      xlim <- c(xmin, xmax) * 1.2
      
      ymin <- min(ylim,coord.ind[, dim2])
      ymax <- max(ylim, coord.ind[, dim2])
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
      
      plot(coord.ind[, axes], xlim = xlim, ylim = ylim, xlab = lab.x, 
           ylab = lab.y, pch = 20, col = as.character(col.plot.ind), 
           cex = cex, main=main, ...)
      abline(h = 0, lty = 2, cex = cex)
      abline(v = 0, lty = 2, cex = cex)
      
      if (leg==TRUE & is.factor(coloring.ind))
        legend(posleg, legend =paste(cl["coloring.ind"],levels(coloring.ind),sep="="), text.col = levels(as.factor(col.plot.ind.total)), 
               cex =cex.leg)
      
      if (label) 
        text(coord.ind[, axes], labels = rownames(coord.ind), 
             pos = 3, col = as.character(col.plot.ind), cex = cex,...)
      
    }
    
    if (!is.null(partial))
    {
      select.partial <- which(rownames(coord.ind[select.ind,,drop=FALSE]) %in% partial)
      if (length(select.partial)==0)
        stop("\"lim.cos.plot\" (or \"lim.contrib.plot\") is too large. No partial individuals can be plotted",call. = FALSE)
      
      coord.ind <- coord.ind[select.partial, , drop=FALSE]
      xmin <- min(xlim, coord.ind[, dim1])
      xmax <- max(xlim, coord.ind[, dim1])
      ymin <- min(ylim,coord.ind[, dim2])
      ymax <- max(ylim, coord.ind[, dim2])
      
      for (i in 1:ngroup)
      { 
        t <- res.mfa$ind.partial[[i]][select.partial,axes,drop=FALSE]
        xmin <- min(xmin, t[,1])
        xmax <- max(xmax, t[,1])
        ymin <- min(ymin,t[,2])
        ymax <- max(ymax, t[,2])
      }
      
      xlim <- c(xmin, xmax) * 1.2
      ylim <- c(ymin, ymax) * 1.2
      
      plot(coord.ind[,axes], xlim = xlim, ylim = ylim, xlab = lab.x, 
           ylab = lab.y, pch = 20, cex = cex,main=main,...)
      abline(h = 0, lty = 2, cex = cex)
      abline(v = 0, lty = 2, cex = cex)
      
      for (i in 1:ngroup){
        t <-res.mfa$ind.partial[[i]][select.partial,axes,drop=FALSE]
        points(t,col=col.groups[i],pch=20)
        for (j in 1:length(select.partial)){
          m<-list(x=c(coord.ind[j,dim1],t[j,1]), y=c(coord.ind[j,dim2],t[j,2]))
          lines(m,col=col.groups[i]) 
        }
      }
      
      if(leg==TRUE)
        legend(posleg, legend = name.groups, text.col = col.groups, cex = cex.leg)
      
      if (label) 
        text(coord.ind[, axes], labels = rownames(coord.ind), 
             pos = 3, cex = cex,...)
    }
  }
  
  #plot of the levels
  if (choice == "levels") 
    {
    if (is.null(main))  main <- "Levels component map"
    if (new.plot) dev.new()
    
    xmin <- min(xlim,res.mfa$levels$coord[, dim1],res.mfa$levels.sup$coord[, dim1])
    xmax <- max(xlim,res.mfa$levels$coord[, dim1],res.mfa$levels.sup$coord[, dim1])
    xlim <- c(xmin, xmax) * 1.2
    
    ymin <- min(ylim,res.mfa$levels$coord[, dim2],res.mfa$levels.sup$coord[, dim2])
    ymax <- max(ylim,res.mfa$levels$coord[, dim2],res.mfa$levels.sup$coord[, dim2])
    ylim <- c(ymin, ymax) * 1.2
    
    plot(0,0, xlim = xlim, ylim = ylim,
         xlab = lab.x, ylab = lab.y, type="n", cex = cex,main=main, ...)
    abline(h = 0, lty = 2, cex = cex)
    abline(v = 0, lty = 2, cex = cex)
    
    #plot levels of active variables
    if (!is.null(res.mfa$levels))
    {
      if (lim.cos2.plot == 0 & lim.contrib.plot==0)
      {
        lim.plot<-0
        base.lim<-res.mfa$levels$cos2[,axes]
      }
      
      if (lim.cos2.plot != 0)
      {
        lim.plot<-lim.cos2.plot
        base.lim<-res.mfa$levels$cos2[,axes]
      }
      
      if (lim.contrib.plot!=0)
      {
        lim.plot<-lim.contrib.plot
        base.lim<-res.mfa$levels$contrib[,axes]
        base.lim<-100*(base.lim/sum(eig.axes))    
      }
      
      coord.lev <- res.mfa$levels$coord[, axes, drop = FALSE]
      
      col.lev <- rep(1,p1+m) #color of a the point-level 
      if (coloring.var == "groups") 
      {
        for (i in 1:ngroup)
          col.lev[which(res.mfa$index.group2==i)] <- col.groups[i]
      }
      col.lev <- col.lev[(p1+1):(p1+m)]
      
      test.empty.plot<-c()
      for (v in 1:nrow(coord.lev)) 
      {
        if (sum(base.lim[v, ], na.rm = TRUE) >= lim.plot && !is.na(sum(base.lim[v, ], na.rm = TRUE))) {
          test.empty.plot<-c(test.empty.plot,1)
          points(coord.lev[v, 1], coord.lev[v,2], col=col.lev[v], pch=20,cex = cex,...)
          
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
            text(coord.lev[v, 1], y = coord.lev[v, 2], col=col.lev[v],
                 labels = rownames(coord.lev)[v], pos = pos, cex = cex)
          }
        }
      }
      if (is.null(test.empty.plot))
        warning("\"lim.cos.plot\" (or \"lim.contrib.plot\") is too large. No level can be plotted")
      if ((coloring.var == "groups") & (leg==TRUE)) 
        legend(posleg, legend = name.groups[unique(res.mfa$index.group2[(p1+1):(p1+m)])], 
               text.col = unique(col.lev), cex = cex.leg)
    }
    #plot levels of supplementary variables
    if (!is.null(res.mfa$levels.sup))
    {
      coord.lev.sup <- res.mfa$levels.sup$coord[, axes, drop = FALSE]
      m.sup <- nrow(coord.lev.sup)
      col.lev.sup <- rep(4,p1.sup+m.sup) #color of a the point-level 
      if (coloring.var == "groups") 
      {
        for (i in 1:ngroupsup)
          col.lev.sup[which(res.mfa$index.groupsup2==i)] <- col.groups.sup[i]
      }
      col.lev.sup <- col.lev.sup[(p1.sup+1):(p1.sup+m.sup)]
      for (v in 1:nrow(coord.lev.sup)) 
      {
        
        points(coord.lev.sup[v, 1], coord.lev.sup[v,2], pch=1,cex = cex,
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
          text(coord.lev.sup[v, 1], y = coord.lev.sup[v, 2], 
               labels = rownames(coord.lev.sup)[v], pos = pos, 
               col=col.lev.sup[v], cex = cex)
        }
      }
      if ((coloring.var == "groups") & (leg==TRUE)) 
        legend(posleg.sup, legend = name.groups[unique(res.mfa$index.groupsup2[(p1+1):(p1+m)])], 
               text.col = unique(col.lev.sup), cex = cex.leg)
    }
  }
}
