###############################################################################
# FONCTIONS POUR L'ACP                                                        #
###############################################################################

# ACP
# Paramètre :	Table R des données
# Sortie :	un objet de type ACP
ACP<-function(d,transpose=T,scale=F,center=T) {
	if(transpose) d<-t(d);
	resacp <-prcomp(x = d,retx = T,center = center,scale = scale);
	resacp$n.obs<-dim(d)[1];
	resacp$percentVar<- resacp$sdev^2 / sum( resacp$sdev^2 )
	return(resacp);
}

	
# Fonctions présentes dans R pour interpréter une ACP :
# soient "acp" un objet de type acp et T la table des données
#
# coefficients de corrélation entre les variables cor(T)
# Inertie des axes, axe i : acpn$sdev^2 , acpn$sdev[i]^2
# Nouvelles coordonnées des individus : acpn$scores
# De l'individu i sur CP j : acpn$scores[i,j]
# Graphique des inerties (valeurs propres) : plot(acp)
# Plan principal : biplot(acp)
# Plan i x j : biplot(acp,c(c1=i, c2=j))


#Tableau des inerties et des pourcentages cumulés
# Paramètre : résultat ACP
# Sortie : tableau des inerties, pourcentages et pourcentages cumulés
VP <- function(resacp) {
	N <- length(resacp$sdev);
	tab <- matrix(nrow=N,ncol=3);
	s <- sum(resacp$sdev^2);
	s1 <- 0;
	for (i in 1:N) {
		tab[i,1] <- resacp$sdev[i]^2;
		tab[i,2] <- resacp$sdev[i]^2/s;
		s1 <- s1+resacp$sdev[i]^2;
		tab[i,3] <- 100*s1/s;
	};
	return(tab)}


# Corrélations entre les axes et les variables initiales
# Paramètres :	table R des données
#		résultat ACP (produit par princomp)
# Sortie :	la matrice des corrélations
AXEVAR <- function(resacp) {return(resacp$rotation)}


# Corrélations entre les k premiers axes et les variables initiales
# Paramètres :	table R des données
#		résultat ACP (produit par princomp)
#		nombre axes
# Sortie :	la matrice des corrélations
AXEVARk <- function(d,resacp,k) {
	return(resacp$rotation[,1:k])
}

# Contribution de la ligne i à l'inertie de l'axe j
# Paramètres :	résultat ACP
#		numéro ligne
#		numéro axe
# Sortie :	pourcentage de la contribution
CTRij <- function(resacp,i,j) {
	x <- resacp$rotation[i,j]^2/(resacp$n.obs * resacp$sdev[j]^2);
        x <- 100*x;
	return(x)}


# Tableau des contribution des lignes aux axes
# Paramètres :	résultat ACP
#		nombre axes
# Sortie :	tableau des pourcentages des contributions
CTR <- function(resacp, nbax) {
      	matrice <- matrix(nrow=resacp$n.obs,ncol=nbax);
        row.names(matrice) <- row.names(resacp$x);
        for (j in 1:nbax) 
		for (i in 1:resacp$n.obs) matrice[i,j] <- CTRij(resacp,i,j);
     
       return(matrice)}

# Fonction utilitaire
SOMME2 <- function(resacp) {
	N <- resacp$n.obs ; M <- ncol(resacp$x);
	s2 <- vector (mode = "numeric", N);
	for (i in 1:N)
		for (j in 1:M) s2[i] <- s2[i] + resacp$x[i,j]^2;
	return(s2)
}

# Cosinus ** 2 des angles de projection
# Paramètres :	résultat ACP
#		nombre axes
# Sortie :	tableau des cos2 des angles de projection
COS2TETA <- function(resacp, nbax) {
	N <- resacp$n.obs ; 
	c2teta <- matrix(nrow=N,ncol=nbax);
	row.names(c2teta) <- row.names(resacp$x);
	s2 <- SOMME2(resacp);
	for (i in 1:N)
		for (j in 1:nbax) c2teta[i,j] <- resacp$x[i,j]^2 / s2[i];
	return(c2teta)
}

############################################################
#PLOTS													   #
############################################################
# Raccourci pour faire afficher un plan de projection
# Paramètres :	résultat ACP
#		premier axe choisi
#		deuxième axe choisi	
PLAN <- function(resacp,i,j) {biplot(resacp,c(c1=i,c2=j))}


# visualisation 3d ggplot

acp3d<-function(pca, comp=1:3, group=rep(1,pca$n.obs), plotVars = FALSE, pointSize=2, plotText=FALSE){
	if(!require("rgl")) stop("You must install rgl");
	if(length(comp)!=3) stop("You must give a vector of 3 integer for comp parameter")
	if(!plotVars){
		x<-pca$x
	}else{
		x<-pca$rotation
	}
	if(is.null(levels(group))){ colors="black"}
	else{
		hashCol<-rainbow(nlevels(group))
		names(hashCol)<-levels(group)
		colors<-hashCol[group]
	}

	
	percentVar <- pca$percentVar
	plot3d(x[,comp[1]],x[,comp[2]],x[,comp[3]],
		xlab=paste0("PC",comp[1],": ",round(percentVar[comp[1]] * 100),"% variance"), 
		ylab=paste0("PC",comp[2],": ",round(percentVar[comp[2]] * 100),"% variance"), 
		zlab=paste0("PC",comp[3],": ",round(percentVar[comp[3]] * 100),"% variance"),
	col=colors,size=pointSize,type=ifelse(plotText,"n","p"))
	
	legend3d("topright", legend = names(hashCol), pch = 16, col = hashCol, cex=1, inset=c(0.02))
	
	if(plotText) text3d(x[,comp[1]],x[,comp[2]],x[,comp[3]],texts=rownames(x),cex=pointSize,col=colors)
	if(plotVars) spheres3d(x=0,y=0,z=0, radius = 1,alpha=0.5,color="white")
	spheres3d(x=0,y=0,z=0, radius = 0.005,alpha=1,color="red")
}

# visualisation 2d ggplot
acp2d<-function(pca, comp=1:2,group=NULL, plotVars = FALSE, pointSize=2, plotText=FALSE,fixedCoord=FALSE,main=NULL,ellipse=FALSE,color=NULL){
	if(!require("ggplot2")) stop("You must install ggplot2");
	if(length(comp)!=2) stop("You must give a vector of 2 integer for comp parameter");

	percentVar <- pca$percentVar
	functPlot<-ifelse(plotText,geom_text,geom_point)
	coor=ifelse(plotVars,"rotation","x")
	
	if(is.null(group)){
		d <- data.frame(PC1=pca[[coor]][,comp[1]], PC2=pca[[coor]][,comp[2]]);
		graph<-ggplot(data=d, mapping = aes(x=PC1, y=PC2, label = rownames(d)))
	}else{
		d <- data.frame(PC1=pca[[coor]][,comp[1]], PC2=pca[[coor]][,comp[2]], group=group);
		graph<-ggplot(data=d, mapping = aes(x=PC1, y=PC2,colour=group, label = rownames(d)))
	}

	graph<-graph+functPlot(size=pointSize)+
		xlab(paste0("PC",comp[1],": ",round(percentVar[comp[1]] * 100),"% variance")) +
		ylab(paste0("PC",comp[2],": ",round(percentVar[comp[2]] * 100),"% variance")) 
	if(fixedCoord)graph <- graph + coord_fixed(ratio=percentVar[comp[2]]/percentVar[comp[1]])
	if(ellipse) graph<-graph+stat_ellipse()
	if(!is.null(main)) graph <- graph + ggtitle(main)
	if(!is.null(color)) graph <- graph + scale_color_manual(values=color)
	
	print(graph)
}

###############################################################################
# Projections d'une partition (obtenue ici avec hclust) sur un plan factoriel #
# (ici : projection sur le plan principal)                                    #
###############################################################################

# Paramètres : résultat de l'acp (renvoyé par princomp),
#              résultat de la cah (renvoyé par hclust),
#              nombre de classes de la cah

# Sortie : les éléments des classes projetés sur le plan principal
#          une couleur par classes

CAHsurACP<-function(acp,cah,k){
    n<-acp$obs
    classe<-cutree(cah,k)
    couleur<-vector("numeric",n)
    # pour avoir la liste des couleurs, on peut raper sous R : colors()
    liste_coul<-c("blue","red", "green", "grey", "orange", "turquoise", "yellow")
    for(i in 1:n){
        couleur[i]<-liste_coul[classe[i]]
    }
    plot (acp$scores[,1],acp$scores[,2], col=couleur)
}

##########################
#Autre type de projection#
##########################

merge0dist<-function(disMat){
	mat<-as.matrix(disMat)
	merged<-list()
	found<-TRUE
	while(found==TRUE){
	  found<-FALSE
	  for(i in 2:nrow(mat)){
		for(j in 1:(i-1)){
		  if(mat[i,j]==0){
			newNames<-rownames(mat)
			newNames<-newNames[-i]
			newMat<-mat[-i,-i]
			colnames(newMat)<-rownames(newMat)<-newNames
			merged[[rownames(mat)[j]]]<-c(merged[[rownames(mat)[j]]],rownames(mat)[i])
			mat<-newMat
			found<-TRUE
			break
		  }
		}
		if(found) break
	  }
	}
	return(list(distMat=as.dist(mat),merged=merged))
}

NMDS<-function(data,transpose=TRUE,scale=FALSE,center=FALSE,metric=dist,ndim=2,maxit=100){
	merged<-FALSE
	require(MASS)
	if(transpose) data <- t(data)
	d <- metric(data)  # euclidean distances between the rows
	if(min(d,na.rm=TRUE)==0){
		merged<-TRUE
		md<-merge0dist(d)
		d<-md$distMat
		mergedSample<-md$merged	
	}
	fit <- isoMDS(d, k=ndim, maxit=maxit) # k is the number of dim
	fit$coord<-fit$points
	fit$points<-NULL
	if(merged){
		for(sple in names(mergedSample)){
			values<-matrix(rep(fit$coord[sple,],length(mergedSample[[sple]])),nrow=length(mergedSample[[sple]]),byrow = TRUE)
			rownames(values)<-mergedSample[[sple]]
			fit$coord<-rbind(fit$coord,values)
		}
	}
	return(fit)
}

proj2d<-function(obj,coord=NULL, axis=1:2,group=NULL, pointSize=2, plotText=FALSE,main=NULL,ellipse=FALSE){
	if(is.null(coord)) coord<-obj$coord
	if(!require("ggplot2")) stop("You must install ggplot2");
	if(length(axis)!=2) stop("You must give a vector of 2 integer for axis parameter");
	functPlot<-ifelse(plotText,geom_text,geom_point)
	if(is.null(group)){
		d <- data.frame(Axis1=coord[,axis[1]], Axis2=coord[,axis[2]]);
		rownames(d)<-rownames(coord)
		graph<-ggplot(data=d, mapping = aes(x=Axis1, y=Axis2, label = rownames(d)))
	}else{
		d <- data.frame(Axis1=coord[,axis[1]], Axis2=coord[,axis[2]], group=group);
		rownames(d)<-rownames(coord)
		graph<-ggplot(data=d, mapping = aes(x=Axis1, y=Axis2,colour=group, label = rownames(d)))
	}
	graph<-graph+functPlot(size=pointSize)+
		xlab(paste0("Axis",axis[1])) +
		ylab(paste0("Axis",axis[2])) 
	if(!is.null(main)) graph <- graph + ggtitle(main)
	if(ellipse) graph<-graph+stat_ellipse()
	print(graph)
}

addIndivACP<-function(acp,indivs,transpose=TRUE,combineMat=TRUE){
  indivs<-as.matrix(indivs)
  if(transpose)indivs<-t(indivs)
  indivs<-indivs[,rownames(acp$rotation)]
  newTab<-matrix(nrow = nrow(indivs),ncol=ncol(acp$rotation),data = 0,dimnames = list(rownames(indivs),colnames(acp$rotation)))
  for(i in 1:nrow(indivs)) newTab[i,]<-colSums(indivs[i,]*acp$rotation)
  if(combineMat){
    acp$x<-rbind(acp$x,newTab)
    return(acp)
  }else{
    return(newTab)
  }
}

#' Add grids to a scatterplot3d
#' 
#' @description The goal of this function is to add grids on an existing
#'  plot created using the package scatterplot3d
#' @param x,y,z numeric vectors specifying the x, y, z coordinates of points.
#'  x can be a matrix or a data frame containing 3 columns corresponding to
#'  the x, y and z coordinates. In this case the arguments y and z are optional
#' @param grid specifies the facet(s) of the plot on which grids should be drawn.
#'  Possible values are the combination of "xy", "xz" or "yz".
#'  Example: grid = c("xy", "yz"). The default value is TRUE to add grids only on xy facet.
#' @param col.grid,lty.grid color and line type to be used for grids
#' @param lab a numerical vector of the form c(x, y, len).
#'  The values of x and y give the (approximate) number of tickmarks on the x and y axes.
#' @param lab.z the same as lab, but for z axis
#' @param scale.y of y axis related to x- and z axis
#' @param angle angle between x and y axis
#' @param "xlim, ylim, zlim" the x, y and z limits (min, max) of the plot.
#' 
#' @note
#' Users who want to extend an existing scatterplot3d graphic with the
#'  function addgrids3d, should consider to set the arguments scale.y, angle, ...,
#'  to the value used in scatterplot3d.
#' 
#' @author Alboukadel Kassambara \email{alboukadel.kassambara@@gmail.com}
#' @references http://www.sthda.com
#' 
#' @example
#' library(scatterplot3d)
#' data(iris)
#' scatterplot3d(iris[, 1:3], pch = 16, grid=T, box=F)
#' addgrids3d(iris[, 1:3], grid = c("xy", "xz", "yz"))
addgrids3d <- function(x, y=NULL, z=NULL, grid = TRUE,
                    col.grid = "grey", lty.grid = par("lty"),
                    lab = par("lab"), lab.z = mean(lab[1:2]),
                    scale.y = 1, angle = 40,
                    xlim=NULL, ylim=NULL, zlim=NULL){
  
  
  if(inherits(x, c("matrix", "data.frame"))){
    x <- as.data.frame(x)
    y <- unlist(x[,2])
    z <- unlist(x[,3])
    x <- unlist(x[,1])
  }
  
  p.lab <- par("lab")
  
  angle <- (angle%%360)/90
  yz.f <- scale.y * abs(if (angle < 1) angle else if (angle >3) angle - 4 else 2 - angle)
  yx.f <- scale.y * (if (angle < 2) 1 - angle else angle - 3)
  
  
  # x axis range
  x.range <- range(x[is.finite(x)], xlim)
  x.prty <- pretty(x.range, n = lab[1], min.n = max(1, min(0.5 *lab[1], p.lab[1])))
  x.scal <- round(diff(x.prty[1:2]), digits = 12)
  x <- x/x.scal
  x.range <- range(x.prty)/x.scal
  x.max <- ceiling(x.range[2])
  x.min <- floor(x.range[1])
  if (!is.null(xlim)) {
    x.max <- max(x.max, ceiling(xlim[2]/x.scal))
    x.min <- min(x.min, floor(xlim[1]/x.scal))
  }
  x.range <- range(x.min, x.max)
  
  # y axis range
  y.range <- range(y[is.finite(y)], ylim)
  y.prty <- pretty(y.range, n = lab[2], min.n = max(1, min(0.5 *lab[2], p.lab[2])))
  y.scal <- round(diff(y.prty[1:2]), digits = 12)
  y.add <- min(y.prty)
  y <- (y - y.add)/y.scal
  y.max <- (max(y.prty) - y.add)/y.scal
  if (!is.null(ylim))
    y.max <- max(y.max, ceiling((ylim[2] - y.add)/y.scal))
  
  # Z axis range
  z.range <- range(z[is.finite(z)], zlim)
  z.prty <- pretty(z.range, n = lab.z, min.n = max(1, min(0.5 *lab.z, p.lab[2])))
  z.scal <- round(diff(z.prty[1:2]), digits = 12)
  z <- z/z.scal
  z.range <- range(z.prty)/z.scal
  z.max <- ceiling(z.range[2])
  z.min <- floor(z.range[1])
  if (!is.null(zlim)) {
    z.max <- max(z.max, ceiling(zlim[2]/z.scal))
    z.min <- min(z.min, floor(zlim[1]/z.scal))
  }
  z.range <- range(z.min, z.max)
  
  # Add grid
  if ("xy" %in% grid || grid == TRUE) {
    i <- x.min:x.max
    segments(i, z.min, i + (yx.f * y.max), yz.f * y.max + 
               z.min, col = col.grid, lty = lty.grid)
    i <- 0:y.max
    segments(x.min + (i * yx.f), i * yz.f + z.min, x.max + 
               (i * yx.f), i * yz.f + z.min, col = col.grid, lty = lty.grid)
  }
   
  if ("xz" %in% grid) {
    i <- x.min:x.max
    segments(i + (yx.f * y.max), yz.f * y.max + z.min, 
             i + (yx.f * y.max), yz.f * y.max + z.max, 
             col = col.grid, lty = lty.grid)
    temp <- yx.f * y.max
    temp1 <- yz.f * y.max
    i <- z.min:z.max
    segments(x.min + temp,temp1 + i, 
             x.max + temp,temp1 + i , col = col.grid, lty = lty.grid)
    
  }
  
  if ("yz" %in% grid) {
    i <- 0:y.max
    segments(x.min + (i * yx.f), i * yz.f + z.min,  
             x.min + (i * yx.f) ,i * yz.f + z.max,  
             col = col.grid, lty = lty.grid)
    temp <- yx.f * y.max
    temp1 <- yz.f * y.max
    i <- z.min:z.max
    segments(x.min + temp,temp1 + i, 
             x.min, i , col = col.grid, lty = lty.grid)
    }
  
}