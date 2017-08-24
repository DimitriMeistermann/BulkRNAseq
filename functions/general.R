#!/usr/bin/R

#<Alias
rn<-rownames;
cn<-colnames;
len<-length;
inter<-intersect;
#>

lire<-function(x, character=FALSE){
	if(character){
		d<-read.table(file = x,sep = "\t",header=T,row.names = 1,colClasses = "character")
	}else{
		d<-read.table(file = x,sep = "\t",header=T,row.names = 1)
	}
	return(d)
}

ecrire<-function(x,file="default.tsv",headRow="Name"){
	options(warn=-1) #Supress unecessary warning about append 
	write.table(x = "NAME\t",file = file,sep = "\t",eol="",quote=F,row.names=F,col.names=F)
	write.table(x=x,file=file,sep="\t", row.names = T, col.names = T, quote = FALSE,append=T)
	options(warn=0)
}

gmean<-function(x, keepZero=F){ #moyenne géometrique
	if(!keepZero){
		x<-x[which(x!=0)]
	}
	return( prod(x)^(1/length(x)) )
}

plotText<-function(cols,label=colnames(cols),cex=.4,xlab="",ylab="",main=""){
	plot(cols,type="n",xlab=xlab,ylab=ylab,main=main)
	text(cols,labels=label,cex=cex)
}

plotMatrixText<-function(mat,main="",cex=1,rowtitle="",rownames=TRUE,colnames=TRUE,insertBlank=c()){
	#Convertie une matrice de caractère en un plot sur la sortie graphique
	matTxt<-matrix("",nrow(mat),ncol(mat))
	for(j in 1:ncol(mat)) matTxt[,j]<-as.character(mat[,j])
	for(j in insertBlank) matTxt<-cbind(matTxt[,1:j],rep("",nrow(matTxt)),matTxt[,(j+1):ncol(matTxt)])
	if(colnames){
		cnMat<-colnames(mat)
		for(j in insertBlank) cnMat<-c(cnMat[1:j],"",cnMat[(j+1):length(cnMat)])
		matTxt<-rbind(cnMat,matTxt)
	}
	if(rownames) matTxt<-cbind(c(rowtitle,rownames(mat)),matTxt)
	#coord calculation
	coordMatrixX<- matrix(0,nrow(matTxt),ncol(matTxt))
	coordMatrixY<-coordMatrixX
	for(i in nrow(coordMatrixX):1){
	  for(j in ncol(coordMatrixY):1){
		coordMatrixX[i,j]<-nrow(coordMatrixX)-i
		coordMatrixY[i,j]<-j
	  }
	}
	plot(1,xlim=c(1,ncol(matTxt)+2),ylim=c(0,nrow(matTxt)),type="n", axes = F, main=main,ylab="",xlab="")
	text(coordMatrixY,coordMatrixX,labels = matTxt,cex=cex,adj = c(0,0))
}

genColWithGrep<-function(nameVector,patternVector){
	#give a color vector from a character vector with a pattern vector
	cols=rainbow(length(patternVector))
	colVector<-rep("#FFFFFFFF",length(nameVector))
	i<-1
	for(pattern in patternVector){
		colVector[grep(pattern,nameVector)]<-cols[i]
		i<-i+1
	}
	return(colVector)
}

genColWithFactors<-function(factorVector){
	levels<-levels(as.factor(factorVector))
	cols=rainbow(length(levels))
	names(cols)<-levels
	return(cols)
	
}

genCOlWithTab<-function(x){
  col<-list()
  for(l in names(x)) col[[l]]<-rgb(x[[l]]$r/255,x[[l]]$g/255,x[[l]]$b/255,names = rownames(x[[l]]))
  return(col)
}


cv<-function(x){ #Coefficient of variation
	return(sd(x)/mean(x));
}

cv2<-function(x){ #Coefficient of variation of Lanner
	return(sd(x)^2/mean(x)^2);
}

se<-function(x){ #Standard mean error
	return(sd(x)/sqrt(length(x)));
}

uncenter<-function(x){
	#transform vector to have no negative value
	return(x+abs(min(x)));
}

takefirst<-function(x){
	uniqDat<-unique(x)
	caseUniq<-c()
	for(i in uniqDat) caseUniq<-c(caseUniq,which(i==x)[1])
	return(x[caseUniq])
}


LevenshteinDist<-function(vectA,vectB){
	lenA<-length(vectA);
	ncols<-lenA+1;
	lenB<-length(vectB);
	nrows<-lenB+1;
	CompMat<-data.frame(matrix(0,nrows,ncols));
	print(dim(CompMat))
	CompMat[,1]<-0:lenB;
	CompMat[1,]<-0:lenA;
	print("mat initialized")
	for(i in 2:nrows){
		for(j in 2:ncols){
			cost<-if(vectA[j-1]==vectB[i-1]) 0 else 1;
			CompMat[i,j]<-min(CompMat[i-1,j]+1, CompMat[i,j-1]+1, CompMat[i-1,j-1]+cost)
		}
		print("i=")
		print(i)
	}
	return(CompMat[nrows,ncols])
}

#Aggregation de ligne/colonne selon une variable qualitative
aggregCols<-function(dataframe,vector,fun=mean){ 
  vector<-as.factor(as.character(vector))
  #pour forcer le mise à jour de l'attribut levels
  lvl<-levels(vector)
  res<-matrix(0,nrow(dataframe),length(lvl))
  rownames(res)<-rownames(dataframe)
  colnames(res)<-lvl
  for(i in lvl){
	if(length(which(vector==i))>1){
		res[,i]<-apply(dataframe[,which(vector==i)],1,fun)
	}
	else{
		res[,i]<-dataframe[,which(vector==i)]
	}
  }
  return(res)
}

aggregRows<-function(dataframe,vector,fun=mean){
  vector<-as.factor(as.character(vector))
  lvl<-levels(vector)
  res<-matrix(0,length(lvl),ncol(dataframe))
  colnames(res)<-colnames(dataframe)
  rownames(res)<-lvl
  for(i in lvl){
	if(length(which(vector==i))>1){
		res[i,]<-apply(dataframe[which(vector==i),],2,fun)
	}else if(length(which(vector==i))==1){
		res[i,]<-unlist(dataframe[which(vector==i),])
	}
  }
  return(res)
}

autoGparFontSizeMatrix<-function(n){ #Calcule automatiquement la taille de police selon le nombre de colonnes ou lignes (empirique)
	n=max(n,50)
	n=min(n,1000)
	return(gpar(fontsize=1/n*600))
}

ConvertKey<-function(keyList, tabKey,colOldKey=1,colNewKey=2){
	#@param keyList: vector of id to be converted
	#@param tabKey : dataframe of x columns with at least two contain old and new keys correspondance
	#@param colOldKey : index of old keys column in tabKey 
	#@param colNewKey : index of new keys column in tabKey 
	hashCorr<-tabKey[,colNewKey]
	names(hashCorr)<-tabKey[,colOldKey]
	returned<-hashCorr[keyList]
	names(returned)<-NULL
	return(as.character(returned))
}

ConvertKeyMatrix<-function(tab,tabKey,colOldKey=1,colNewKey=2,first=TRUE,dim=1,fun){
	if(dim==2) tab<-t(tab)
	keyList<-rownames(tab)
	newkey<-ConvertKey(keyList, tabKey,colOldKey,colNewKey)
	names(newkey)<-as.character(1:length(newkey))
	if(first){
	  newkey<-newkey[which(!is.na(newkey))]
	  newkey<-takefirst(newkey)
	  tab<-tab[as.numeric(names(newkey)),]
	  rownames(tab)<-newkey
	}else{
	  tab<-tab[which(!is.na(newkey)),]
	  newkey<-newkey[which(!is.na(newkey))]
	  tab<-aggregRows(tab,newkey,fun=fun)
	}
	if(dim==2) tab<-t(tab)
	return(tab)
}

supprNAnames<-function(x ,side=1){
	#Vector case
	if(is.null(dim(x))){
		return(x[which(!is.na(names(x)))]);
	}
	#Col case
	if(side==2){
		return(x[,which(!is.na(colnames(x)))]);
	}
	#Row case
	else{
		return(x[which(!is.na(rownames(x))),]);
	}
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

mobileMean<-function(datVector, step=2){
#Lissage de signal à partir de moyenne mobile
  pond<-datVector
  for(i in 1:length(datVector)){
    pond[i]=mean(datVector[max(1,(i-step)):min((i+step),length(datVector))])
  }
  return(pond)
}

OrderBestUpDown<-function(x){ #-3,4,-2,2,6 --> 6,4,2,-3,-2
#x : named vector with positive and negative values  
	up<-sort(x[which(x>0)],decreasing =TRUE)
	down<-sort(x[which(x<=0)])
	return(c(up,down))
}


rowScale<-function(data,center=TRUE,scaled=FALSE){
	data<-t(data)
	data<-t(scale(data,center=center,scale=scaled))
	return(data)
}

corrDist<-function(x){
	return(as.dist(1-cor(x)))
}

chead<-function(x, n=5){
	print(x[1:n,1:n])
}

inflaslope<-function(x, ordering=TRUE, decreasing = FALSE){
#Calculer point d'ébouli des valeur
	if(ordering) x<-sort(x, decreasing=decreasing)
	l1<-length(x)
	t1<-cumsum(x)
	t2<-abs(t1[1:(l1-1)]-t1[2:l1])[1:(l1-1)]
	t3<-which(t2>sd(x))[1]
	return(t3)
}