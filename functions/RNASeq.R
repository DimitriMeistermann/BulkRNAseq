#!/usr/bin/R

#need general.R

filterMostExprimedGenesBySample <-function(data, numberOfGenes=nrow(data),minFreqOfGene=1,maxFreqOfGene=ncol(data),threshold=min(data)){
	#Fonction de filtrage des counts, on prend les x gènes les plus exprimés pour chache échantillon puis on effectue une jointure complète
	#puis on filtre selon la fréquence du gène
	#data: dataframe des counts
	#numberOfGenes: cutoff sur les rank d'expression des gènes, plus il est élevé moins le filtre est stringeant
	#minFreqOfGene : nombre de fois minimum où le gène revient dans la liste des x gènes les plus exprimés (où le gène > threshold)
    #maxFreqOfGene : nombre de fois maximale où le gène est exprimé
	#threshold: l'expression du gène doit être plus grande que ce paramètre pour que la fréquence soit comptée comme 1
	mostExprimedGenes<-list()
	for(i in 1:ncol(data)){
		col<- data[,i]
		names(col)<-rownames(data)
		col<-col[which(col>threshold)];
		mostExprimedGenes[[i]]<-names(sort(col,decreasing = T))[1:min(numberOfGenes,length(col))]
	}
	rm(col)

	freqTable<-summary(as.factor(unlist(mostExprimedGenes)),maxsum=nrow(data))
	mostExprimedgenesVector<-names(freqTable[which(freqTable>=minFreqOfGene & freqTable<=maxFreqOfGene)])
	return(data[mostExprimedgenesVector,])
}

vectSup0<-function(x){
	a<-x[which(x>0)]
	return(length(a))
}

nbGeneSup0<-function(data){
	nS<-ncol(data)
	i<-apply(data,1,vectSup0)
	return(length(i[which(i>0)]))
}

filterMostExprimedGenesBySampleThres <-function(data, numberOfGenes=3000,maxGenes=nrow(data)/2){
	#Fonction de filtrage des counts, on prend les x gènes les plus exprimés pour chache échantillon puis on effectue une jointure complète
	#puis on filtre selon la fréquence du gène
	#data: dataframe des counts
	#numberOfGenes: cutoff des gènes, plus il est élevé moins le filtre est stringeant
	#maxGenes : nombre maximal de gène que l'on veut en sortie
	mostExprimedGenes<-data.frame(matrix(ncol = ncol(data), nrow=numberOfGenes)) #création d'un dataframe vide
	for(i in 1:ncol(data)){
	col<- data[,i]
	names(col)<-rownames(data)
	mostExprimedGenes[,i]<-names(sort(col,decreasing = T,method="shell"))[1:numberOfGenes]
	}
	rm(col)
	colnames(mostExprimedGenes)<-colnames(data)
	freqTable<-sort(summary(as.factor(unlist(mostExprimedGenes)),maxsum=nrow(data)),decreasing=TRUE,method="shell")
	maxGenes<-min(maxGenes,length(freqTable))
	mostExprimedgenesVector<-names(freqTable[1:maxGenes])
	return(data[mostExprimedgenesVector,])
}

interceptDistThres<-function(data,threshold=nrow(data)/3){
	print("construction of database");
	nrows<-nrow(data);
	ncols<-ncol(data);
	listOfgenesBySample=matrix("",threshold,ncols);
	for(i in 1:ncols){
		col<-data[,i];
		names(col)<-1:nrows;
		#tri alp ?
		listOfgenesBySample[,i]<-names(sort(col,decreasing = T,method="shell"))[1:threshold];
	}
	print("construction of dist matrix");
	distMat<-matrix(0,ncols,ncols)
	rownames(distMat)<-colnames(data);
	colnames(distMat)<-rownames(distMat);
	c<-1;
	for(i in 2:ncols){
		for(j in 1:c){
			distMat[i,j]<-length(setdiff(listOfgenesBySample[,i],listOfgenesBySample[,j]));
		}
		c<-c+1;
	}
	distMat<-as.dist(distMat);
	return(distMat);
}

interceptDistanceLevenstein<-function(data){
	print("construction of database");
	nrows<-nrow(data);
	ncols<-ncol(data);
	listOfgenesBySample=matrix("",nrows,ncols);
	for(i in 1:ncols){
		col<-data[,i];
		names(col)<-1:nrows;
		listOfgenesBySample[,i]<-names(sort(col,decreasing = T,method="shell"));
	}
	print("construction of dist matrix");
	distMat<-matrix(0,ncols,ncols)
	rownames(distMat)<-colnames(data);
	colnames(distMat)<-rownames(distMat);
	c<-1;
	for(i in 2:ncols){
		for(j in 1:c){
			distMat[i,j]<-LevenshteinDist(listOfgenesBySample[,i],listOfgenesBySample[,j]);
		}
		c<-c+1;
	}
	distMat<-as.dist(distMat);
	return(distMat);
}

compareVectReplace<-function(data,verbose=FALSE){
	#return a dist matrix
	#@param data : dataframe or matrx (rows = genes, cols = samples)
	# tri alph
	print("construction of database");
	nrows<-nrow(data);
	ncols<-ncol(data);
	listOfgenesBySample=matrix(0,nrows,ncols);
	for(i in 1:ncols){
		col<-data[,i];
		names(col)<-1:nrows;
		listOfgenesBySample[,i]<-as.integer(names(sort(col,decreasing = T,method="shell")));
	}
	print("construction of dist matrix");
	distMat<-matrix(0,ncols,ncols)
	rownames(distMat)<-colnames(data);
	colnames(distMat)<-rownames(distMat);
	tot<-length(as.dist(distMat))
	c<-1;
	out<-0;
	comptT<-0
	for(i in 2:ncols){
		for(j in 1:c){
			distMat[i,j]<-.Call("orderDist",as.integer(listOfgenesBySample[,i]),as.integer(listOfgenesBySample[,j]))
			comptT<-comptT+1;
		}
		if(verbose) print(paste0(comptT,"/",tot));
		c<-c+1;
	}
	distMat<-as.dist(distMat);
	return(distMat);
}

genesRankVSExpr<-function(data, numberOfGenes=nrow(data)){
	#Stats par gènes
	nsamples<-ncol(data)
	mostExprimedGenes<-data.frame(matrix(ncol = nsamples, nrow=numberOfGenes));
	for(i in 1:nsamples){
		col<- data[,i];
		names(col)<-rownames(data);
		mostExprimedGenes[,i]<-names(sort(col,decreasing = T))[1:numberOfGenes];
	}
	colnames(mostExprimedGenes)<-colnames(data);
	ranks<-data.frame(matrix(ncol = nsamples, nrow=numberOfGenes));
	rownames(ranks)<-mostExprimedGenes[,1]
	colnames(ranks)<-colnames(data);
	for(i in 1:nsamples){
		col<-1:numberOfGenes;
		names(col)<-mostExprimedGenes[,i];
		ranks[,i]<-col[rownames(ranks)]
	}
	meanRank<-apply(ranks,1,mean)
	sdRank<-apply(ranks,1,sd)
	minRank<-apply(ranks,1,min)
	maxRank<-apply(ranks,1,max)
	data<-data[names(meanRank),]
	sdExpr<-apply(data,1,mean)
	minExpr<-apply(data,1,sd)
	maxExpr<-apply(data,1,max)
	meanExpr<-apply(data,1,min)
	res<-data.frame(meanRank = meanRank,SDRank=sdRank, minRank=minRank, maxRank=maxRank, 
		meanExpr = meanExpr,sdExpr=sdExpr,maxExpr=maxExpr,minExpr=minExpr) ;
	res<-res[order(res$meanRank),]
	return(res);
}

examineRNAseqSamples<-function(x, uncenter= FALSE){
	#Donne différentes stat par échantillon 
	if(uncenter){
		x<-x-min(x)
		zero<-0
	}else{
		zero<-min(x)
	}
	mean<-apply(x,2,mean)
	sd<-apply(x,2,sd)
	count<-colSums(x)
	CV<-apply(x,2,cv)
	noGenEx<-rep(0,ncol(x))
	for(i in 1:ncol(x)) noGenEx[i]<-length(which(x[,i]>zero))

	return(data.frame(mean=mean, sd=sd,CV=CV,TotalGenEx=noGenEx,TotalCount=count))
}

retrieveSexHumanEmbryoKmeans<-function(d,patternLen=3){
	#return a matrix of count of "male" and "female" predicted cells in each embryo (Kmeans method)
	maleGene<-c("DDX3Y","EIF1AY","TTTY15","RPS4Y1")
	k<-kmeans(t(d[maleGene,]),2)
	mORf<-rowSums(k$centers)
	if(mORf[1]<mORf[2]){
		mf<-c("F","M")
	}else{
		mf<-c("M","F")
	}
	count<-as.factor(mf[k$cluster])
	names(count)<-substr(colnames(d),1,patternLen)
	embryos<-as.factor(unique(substr(colnames(d),1,patternLen)))
	res<-list()
	res$count<-data.frame(matrix(ncol=2,nrow=length(embryos)))
	colnames(res$count)<-c("Male","Female")
	rownames(res$count)<-embryos
	
	for(embryo in embryos){
		res$count[embryo,1]<-length(which(count[which(names(count)==embryo)]=="M"))
		res$count[embryo,2]<-length(which(count[which(names(count)==embryo)]=="F"))
	}
	res$freq<-res$count/rowSums(res$count)
	return(res)
}

retrieveSexHumanEmbryoACP<-function(d,patternLen=3){
	#return a freq matrix of "male" and "female" predicted cells in each embryo (ACP method)
	maleGene<-c("DDX3Y","EIF1AY","TTTY15","RPS4Y1")
	z<-colSums(d[maleGene,])
	#voir ici pour renvoyer tout mâle ou tout femelle
	#normer acp sans réduire ? et prendre seuil
	
	a<-ACP(t(d[maleGene,]))
	M<-as.factor(substr(names(a$x[which(a$x[,1]>0),1]),1,patternLen)) #sélection du nom de l'embryon seul, pas des cellules (substr)
	F<-as.factor(substr(names(a$x[which(a$x[,1]<0),1]),1,patternLen))
	embryos<-as.factor(unique(substr(colnames(d),1,patternLen)))
	res<-list()
	res$count<-data.frame(matrix(ncol=2,nrow=length(embryos)))
	colnames(res$count)<-c("Male","Female")
	rownames(res$count)<-embryos
	
	for(embryo in embryos){
		res$count[embryo,1]<-length(which(M==embryo))
		res$count[embryo,2]<-length(which(F==embryo))
	}
	res$freq<-res$count/rowSums(res$count)
	return(res)
}

#' Easy and  quick view of expression
#' @param expr dataframe or matrix. Expression.
#' @param conditions dataframe of factors.  Sample table with sample in row and annotations in column.
#' @param legendName character. Custom legend name.
#' @param errorBar character. What represent the bars around each point. "se" : standard error mean, "ci": confidance interval (distribution must follow a normal law) "na":none
#' @param ciRate numeric. Confidance interval rate if errorBar = ci, 0.95 = CI 95%
#' @param geom character. GGplot function name for representation, Ex : geom="point" if geom_point is wanted.
#' @param addLine logical. Add extra line to the plot ?
#' @param xaxis character. 'gene' for representing genes and 'annot' for annotation on x axis (remaining parameter will be represented with color)
#' @param negValue logical. error bars sub zero ?
#' @param scale character. 'identity' or 'log' scaled ?
#' @param breaks numeric. Position of breaks on y axis.
#' @param xRotationLab numeric. Rotation of x axis labels
#' @param hjust numeric. Horizontal justification
#' @param main character. Main title
#' @param xLabelSize numeric. Size of x-axis labels
#' @param colorType character. How ggplot color graph, 'fill', 'contour' or 'both'
#' @param returnTab logical. Return processed data ready to plot by GGplot ?
#' @param ... Additional parameter passed to ggplot2
#' @return Nothing
#' @examples
#' genes<-c("geneA","geneB")
#' samples<-paste0("sample",as.character(1:10))
#' #Generation of random expression
#' expr<-rbind(abs(c(rnorm(5,6,4),rnorm(5,8,2))),abs(c(rnorm(5,100,15),rnorm(5,50,20))))
#' rownames(expr)<-genes; colnames(expr)<-samples
#'
#' #Generation of sample table
#' sampleAnnot<-data.frame(group=as.factor(c(rep("g1",5),rep("g2",5))),row.names=samples)
#'
#' plotExpr(expr=expr,conditions=sampleAnnot,scale="log")
plotExpr<-function(expr,conditions=NULL,legendName="gene",errorBar="se", ciRate=0.95,  geom="point", addLine=FALSE,  xaxis="gene", negValue=FALSE, 
	scale="identity",breaks = waiver(),xRotationLab=0, hjust=0.5,main=NULL,xLabelSize=10, colorType="contour",returnTab=FALSE,...){
	
	require(grid)
	require(ggplot2)
		
	if(!errorBar%in%c("se","ci","na")) stop("errorBar must be 'ci' or 'se'")
	if(!xaxis%in%c("gene","annot")) stop("xaxis must be 'gene' or 'annot'")
	if(!scale%in%c("identity","log")) stop("scale must be 'identity' or 'log'")
	if(!colorType%in%c("contour","fill","both")) stop("scale must be 'identity' or 'log'")
	
	colorContour<-FALSE
	colorFill<-FALSE
	if(colorType=="contour" | colorType=="both") colorContour=TRUE
	if(colorType=="fill" | colorType=="both") colorFill=TRUE
	
	if(is.vector(expr))  expr<-t(data.frame(expr=expr,row.names = names(expr)))
	if(is.null(conditions)) conditions<-data.frame(samples=as.factor(colnames(expr)))
	if(is.vector(conditions))  conditions<-data.frame(cond=conditions,row.names = colnames(expr))
	numPlots = ncol(conditions)
	cols<-floor(sqrt(numPlots))
	layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
	ncol = cols, nrow = ceiling(numPlots/cols))
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

	ciFactor<-1
	if(errorBar=="ci") ciFactor<-qnorm(ciRate+(1-ciRate)/2)
	
	for(condIndex in 1:ncol(conditions)){
		conditionName<-colnames(conditions)[condIndex]
		cond<-unlist(conditions[,condIndex])
		names(cond)<-rownames(conditions)
		if(!is.factor(cond)) stop("You must give a factor dataframe/matrix") #Qualitatif
			tabGraph<-data.frame(matrix(nrow = length(levels(cond))*nrow(expr),ncol=5))
			colnames(tabGraph)<-c(conditionName,"means","errBarMin","errBarMax","gene")
			tabGraph[,conditionName]<-rep(levels(cond),nrow(expr))
			tabGraph$gene<-rep(rownames(expr),each=length(levels(cond)))

		for(lvl in levels(cond)){
			for(exprIndex in 1:nrow(expr)){
				values<-unlist(expr[exprIndex,which(colnames(expr)%in%names(cond[which(cond==lvl)]))])
				nameExpr<-rownames(expr)[exprIndex]
				tabGraph[which(tabGraph[,conditionName]==lvl & tabGraph$gene == nameExpr),"means"]<-mean(values)
				tabGraph[which(tabGraph[,conditionName]==lvl & tabGraph$gene == nameExpr),"errBarMin"]<-  se(values)*ciFactor
				tabGraph[which(tabGraph[,conditionName]==lvl & tabGraph$gene == nameExpr),"errBarMax"]<-  se(values)*ciFactor
			}
		}
		tabGraph[,conditionName]<-factor(tabGraph[,conditionName],levels=levels(cond))
		tabGraph$gene<-factor(tabGraph$gene,levels=rownames(expr))

		if(!negValue)tabGraph$errBarMin[which(tabGraph$means-tabGraph$errBarMin<0)]<-tabGraph$means[which(tabGraph$means-tabGraph$errBarMin<0)]
		
		if(xaxis=="annot"){
			g<-"gene"
			x<-conditionName
		}else{
			g<-conditionName
			x<-"gene"
		}
		paramAES<-list(y="means",x=x,group=g)
			
		if(colorContour) paramAES$colour= g
		if(colorFill) paramAES$fill= g
		
		graph<-ggplot(data=tabGraph,mapping=do.call("aes_string",paramAES))
		
		if(geom=="bar"){
			graph<-graph+geom_bar(stat = "identity",width=.4,position = "dodge",...)
		}else{
			graph<-graph+do.call(paste0("geom_",geom),list(...))
		}
		if(addLine) graph<-graph+geom_line()
		
		if(errorBar != "na"){
			if(geom=="bar"){
				graph<-graph+geom_errorbar(width=.1, aes(ymin=means-errBarMin, ymax=means+errBarMax),position = position_dodge(width=.4))
			}else{
				graph<-graph+geom_errorbar(width=.1, aes(ymin=means-errBarMin, ymax=means+errBarMax))
			}
		}

		if(scale=="log"){
			graph<-graph+scale_y_log10(breaks=breaks)
		}else{
			graph<-graph+scale_y_continuous(breaks=breaks)
		}

		if(! is.null(main)) graph<-graph+ggtitle(main)
		graph<-graph+theme(axis.text.x = element_text(angle = xRotationLab, hjust = hjust,vjust=.3,size=xLabelSize))
		matchidx <- as.data.frame(which(layout == condIndex, arr.ind = TRUE))
		print(graph, vp = viewport(layout.pos.row = matchidx$row,layout.pos.col = matchidx$col))
	}
	if(returnTab) return(tabGraph)
}

#Utile pour les packages qui demandes des fonctions comme argument
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}
selectGene<- function(x){
  return(x==1)
}

###

UMI2UPM<-function(data){ #Normalisation UPM
	data.UPM <- sweep(data, 2, colSums(data),`/`)
	data.UPM <-data.UPM * 1000000
	return(data.UPM)
}

plotDistrib<-function(data,type="boxplot",conditions=NULL,main=NULL){
  require(grid)
  require(ggplot2)
  #data : matrix of expression data
  #type: 'violin' or 'boxplot' ?
  #conditions : vector of factors to color plot
  if(!type%in%c("violin","boxplot")) stop("type must be 'violin', 'hist' or 'boxplot'")
  
  vectorDat<-as.vector(as.matrix(data))
  tabGraph<-data.frame(val=vectorDat,sample=rep(cn(data),each=nrow(data)))
  if(!is.null(conditions)){
    tabGraph[,conditionName]=rep(conditions,each=nrow(data))
  }
  
  tabGraph$sample<-factor(tabGraph$sample,levels=cn(data))
  
  if(is.null(conditions)){
    graph<-ggplot(data = tabGraph,mapping = aes(sample,val))
  }else{
    graph<-ggplot(data = tabGraph,mapping = aes(sample,val,color=cond))
  }
  
    graph<-graph+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.3))
  if(type=="violin") graph<-graph+geom_violin()
  if(type=="boxplot") graph<-graph+geom_boxplot()
  if(! is.null(main)) graph<-graph+ggtitle(main)
  
  print(graph)
}

