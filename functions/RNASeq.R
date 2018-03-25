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

retrieveSexHumanEmbryoKmeans<-function(d,group){
	group<-droplevels(as.factor(group))
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
	embryos<-as.factor(levels(group))
	names(count)<-group
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
#' @param axis.names character. A two element vector containing axis names
#' @param colorScale list or vector. If condition has one column, a character specifying color, else a list of character. 
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
plotExpr<-function(expr,conditions=NULL,legendName=NULL,errorBar="se", ciRate=0.95,  geom="point", addLine=FALSE,  xaxis="gene", negValue=FALSE, 
	scale="identity",breaks = waiver(),xRotationLab=0, hjust=0.5,main=NULL,xLabelSize=10, colorType="contour",returnTab=FALSE,colorScale=NULL,axis.names=NULL,...){
	
	require(grid)
	require(ggplot2)
		
	if(!errorBar%in%c("se","ci","na")) stop("errorBar must be 'na', 'ci' or 'se'")
	if(!xaxis%in%c("gene","annot")) stop("xaxis must be 'gene' or 'annot'")
	if(!scale%in%c("identity","log")) stop("scale must be 'identity' or 'log'")
	if(!colorType%in%c("contour","fill","both")) stop("scale must be 'identity' or 'log'")
	
	colorContour<-FALSE
	colorFill<-FALSE
	nullConditions<-is.null(conditions)
	nullColorScale<-is.null(colorScale)
	nullLegendName<-is.null(legendName)
	colorScaleList<-NULL
	if(colorType=="contour" | colorType=="both") colorContour=TRUE
	if(colorType=="fill" | colorType=="both") colorFill=TRUE
	
	if(is.vector(expr))  expr<-t(data.frame(expr=expr,row.names = names(expr)))
	if(nullConditions){
		conditions<-as.factor(colnames(expr))
		errorBar<-"na"
	}
	if(!(is.data.frame(conditions)|is.matrix(conditions)))  conditions<-data.frame(cond=conditions,row.names = colnames(expr))
	numPlots = ncol(conditions)
	if(numPlots>1 & is.list(colorScale)) colorScaleList<-colorScale
	cols<-floor(sqrt(numPlots))
	layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
	ncol = cols, nrow = ceiling(numPlots/cols))
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

	ciFactor<-1
	if(errorBar=="ci") ciFactor<-qnorm(ciRate+(1-ciRate)/2)
	
	for(condIndex in 1:ncol(conditions)){
		conditionName<-colnames(conditions)[condIndex]
		if(!is.null(colorScaleList)) colorScale<-colorScaleList[[conditionName]]
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
		
		if(xaxis=="gene"){
			g<-"gene"
			x<-conditionName
		}else{
			g<-conditionName
			x<-"gene"
		}
		paramAES<-list(x=x,group=g)
		if(!nullConditions){
			if(colorContour) paramAES$colour= g
			if(colorFill) paramAES$fill= g
		}
		addParam<-list(...)
		if(geom=="bar"){
			addParam$stat = "identity";addParam$width=.4;addParam$position = "dodge"
		}
		paramAES$y="means"
		graph<-ggplot(data=tabGraph,mapping=do.call("aes_string",paramAES))+do.call(paste0("geom_",geom),addParam)
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
		
		if(nullLegendName) legendName=g
		if(nullColorScale) colorScale<-ggplotColours(nlevels(tabGraph[,g]))
		if(colorFill) graph<-graph+scale_fill_manual(name=legendName,values=colorScale)
		if(colorContour) graph<-graph+scale_color_manual(name=legendName,values=colorScale)
		
		if(! is.null(main)) graph<-graph+ggtitle(main)
		if(is.null(axis.names)) axis.names<-c("","Expression")
		graph<-graph+xlab(axis.names[1]) + ylab(axis.names[2])
		graph<-graph+theme(axis.text.x = element_text(angle = xRotationLab, hjust = hjust,vjust=.3,size=xLabelSize),
			panel.background = element_rect(fill = NA,colour="black"),
			panel.grid.major = element_line(colour = "grey50"),
			panel.grid.minor = element_line(colour = "grey50"))
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

TPMfullLength<-function(data, gene.length){
	gene.length.kb <- gene.length[rn(data)]/1000
	data<-sweep(data, 1, gene.length.kb,`/`)
	return(UMI2UPM(data))
}

plotDistrib<-function(data,type="boxplot",conditions=NULL,main=NULL,conditionName="Batch"){
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
    graph<-ggplot(data = tabGraph,mapping = aes_string("sample","val",color=conditionName))
  }
  
    graph<-graph+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.3))
  if(type=="violin") graph<-graph+geom_violin()
  if(type=="boxplot") graph<-graph+geom_boxplot()
  if(! is.null(main)) graph<-graph+ggtitle(main)
  
  print(graph)
}

#Return gene id correspondance, GO species code and KEGG species code
getSpeciesData<-function(sample.species="Human",genes,updateSpeciesPackage=FALSE){
	require(gage)
	data(bods)
	species<-list()
	species.data<-data.frame(bods)
	species.index<-which(species.data$species==sample.species)
	if(len(species.index)!=1) stop("Wrong species name, type \"species.data$species\" to see available species")
	species$package<-as.character(species.data[species.index,"package"])
	species$kegg<-as.character(species.data[species.index,"kegg.code"])
	species$go<-strsplit(as.character(species$package),split = ".",fixed = TRUE)[[1]][2]
	if(updateSpeciesPackage) biocLite(species$package)
	library(species$package,character.only = TRUE)
	species$GeneIdTable<-select(get(species$package),genes, "ENTREZID","SYMBOL")
	return(species)
}

#Take a gene list, add transcription factor data and long Gene name
detailOnGenes<-function(x,tfDat,speciesDat){
	require(AnnotationDbi)
	if(is.data.frame(x) | is.matrix(x)){
		geneSym<-rn(x)
		res<-data.frame(x)
	}else{
		if(is.null(names(x))){
			res<-data.frame(row.names=x)
			geneSym<-x
		}else{
			geneSym<-names(x)
			res<-data.frame(val=x,row.names=names(x))
		}
	}
	geneNames<-select(get(speciesDat$package),geneSym,"GENENAME","SYMBOL")
	res$LongName<-ConvertKey(geneSym,tabKey = geneNames,colOldKey = "SYMBOL",colNewKey = "GENENAME")
	genesTF<-intersect(rn(tfDat),geneSym)
	res$TFdegree<-0
	res[genesTF,"TFdegree"]<-tfDat[genesTF,"tf_degree"]
	return(res)
}


#Optimal cutree on hclust
best.cutree <- function(hc, min=3, max=20, loss=FALSE, graph=FALSE, ...){
  if (class(hc)!="hclust") hc <- as.hclust(hc)
  max <- min(max, length(hc$height))
  inert.gain <- rev(hc$height)
  intra <- rev(cumsum(rev(inert.gain)))
  relative.loss = intra[min:(max)]/intra[(min - 1):(max - 1)]
  best = which.min(relative.loss)
  names(relative.loss) <- min:max
  if (graph) {
	temp <- relative.loss
	temp[best] <- NA
	best2 <- which.min(temp)
	pch <- rep(1, max-min+1)
	pch[best] <- 16
	pch[best2] <- 21
	plot(min:max, relative.loss, pch=pch, bg="grey75", ...)
  } else {
	if (loss)
	  relative.loss
	else
	  best + min - 1
  }
}

Enrich<-function(x, corrIdGenes=NULL,database=c("reactom","goBP"),way="upreg",maxSize=500,nperm=1000,customAnnot=NULL,returnLeadingEdge=FALSE,
  keggDisease=FALSE,species="Human",returnTerm=FALSE,...){
	require(fgsea)
	require(AnnotationDbi)
	select<-AnnotationDbi::select
	validDBs<-c("kegg","reactom","goBP","goCC","goMF","custom")
	if(!(combineLogical(database%in%validDBs))) stop(paste0("Error, valid values for database are: ",paste0(validDBs,collapse=", ")))
	if(!(way%in%c("upreg","downreg","both"))) stop("Error, valid values for way are: upreg, downreg or both") 
	if(is.null(customAnnot) & "custom"%in%database) stop("You must give a value a list in customAnnot if database=custom")
	
	if(is.data.frame(x) | is.matrix(x)){
		tempx<-x
		x<-tempx[,1]
		names(x)<-rownames(tempx)
	}
	options(warn=-1)
	if(is.null(corrIdGenes)) corrIdGenes<-getSpeciesData(species,names(x))$GeneIdTable
	geneSym<-names(x)
	geneEntrez<-ConvertKey(names(x),tabKey = corrIdGenes,colOldKey = "SYMBOL",colNewKey = "ENTREZID")
	new_x<-x[!is.na(geneEntrez)]
	names(new_x)<-geneEntrez[!is.na(geneEntrez)]
	geneEntrez<-names(new_x)
	
	if(way=="downreg") new_x <- -new_x
	if(way=="both")    new_x <- abs(new_x)
	
	db_terms<-list()
	if(is.list(customAnnot)) db_terms$custom<-customAnnot
	
	if(!(length(database)<=1 & database[1]=="custom")){
		if("reactom"%in%database){ 
			require("reactome.db")
			db_terms$reactom<- reactomePathways(geneEntrez)
		}
		if("kegg"%in%database){
			require("gage")
			kg.species <- kegg.gsets(tolower(species), id.type="entrez")
			db_terms$kegg<- if(keggDisease) kg.species$kg.sets else kg.species$kg.sets[kg.species$sigmet.idx]
		}
		if("go"%in%substr(database,1,2)){
			require("gage")
			go.species <- go.gsets(tolower(species), id.type="entrez")
			if("goBP"%in%database) db_terms$goBP<-go.species$go.sets[go.species$go.subs$BP]
			if("goMF"%in%database) db_terms$goMF<-go.species$go.sets[go.species$go.subs$MF]
			if("goCC"%in%database) db_terms$goCC<-go.species$go.sets[go.species$go.subs$CC]
		}
	}
	options(warn=0)
	if(length(db_terms)==0) stop("Error, no term in any database was found")
	if(returnTerm) return(db_terms)
	res<-list()
	for(db in names(db_terms)){
		db_terms[[db]]<-sapply(db_terms[[db]],function(term){
		  term[term %in% geneEntrez]
		})
		res[[db]]<-fgsea(db_terms[[db]], new_x,nperm=nperm,maxSize=maxSize,...)
		res[[db]]<-res[[db]][order(res[[db]]$NES),]
		res[[db]]$leadingEdge<- if(returnLeadingEdge) sapply(res[[db]]$leadingEdge,ConvertKey,tabKey=corrIdGenes,colOldKey = "ENTREZID",colNewKey = "SYMBOL") else NULL
		res[[db]]$database<-db
	}
	res<-do.call("rbind", res)
	return(res)
}

EnrichFisher<-function(x, corrIdGenes=NULL,database=c("reactom","goBP"),maxSize=500,minSize=2,returnGene=FALSE, keggDisease=FALSE,species="Human"){
	require(AnnotationDbi)
	select<-AnnotationDbi::select
	validDBs<-c("kegg","reactom","goBP","goCC","goMF","custom")
	if(!(combineLogical(database%in%validDBs))) stop(paste0("Error, valid values for database are: ",paste0(validDBs,collapse=", ")))
	if(is.null(customAnnot) & "custom"%in%database) stop("You must give a value a list in customAnnot if database=custom")

	if(is.data.frame(x) | is.matrix(x)){
		tempx<-x
		x<-tempx[,1]
		names(x)<-rownames(tempx)
	}

	if(is.null(corrIdGenes)) corrIdGenes<-getSpeciesData(species,names(x))$GeneIdTable
	geneSym<-names(x)
	geneEntrez<-ConvertKey(names(x),tabKey = corrIdGenes,colOldKey = "SYMBOL",colNewKey = "ENTREZID")
	new_x<-x[!is.na(geneEntrez)]
	names(new_x)<-geneEntrez[!is.na(geneEntrez)]
	geneEntrez<-names(new_x)

	db_terms<-list()
	if(is.list(customAnnot)) db_terms$custom<-customAnnot

	if(!(length(database)<=1 & database[1]=="custom")){
		if("reactom"%in%database){ 
			require("reactome.db")
			db_terms$reactom<- reactomePathways(geneEntrez)
		}
		if("kegg"%in%database){
			require("gage")
			kg.species <- kegg.gsets(tolower(species), id.type="entrez")
			db_terms$kegg<- if(keggDisease) kg.species$kg.sets else kg.species$kg.sets[kg.species$sigmet.idx]
			db_terms$kegg
		}
		if("go"%in%substr(database,1,2)){
			require("gage")
			go.species <- go.gsets(tolower(species), id.type="entrez")
			if("goBP"%in%database) db_terms$goBP<-go.species$go.sets[go.species$go.subs$BP]
			if("goMF"%in%database) db_terms$goMF<-go.species$go.sets[go.species$go.subs$MF]
			if("goCC"%in%database) db_terms$goCC<-go.species$go.sets[go.species$go.subs$CC]
		}
	}

	nInterest<-length(which(new_x))
	nNotInterest<-length(which(!new_x))

	results<-list()
	for(db in names(db_terms)){
		db_terms[[db]]<-sapply(db_terms[[db]],function(term){
			term[term %in% geneEntrez]
		})

		len_term<-sapply(db_terms[[db]],length)
		db_terms[[db]]<-db_terms[[db]][len_term>=minSize & len_term<=maxSize]

		nGeneByterm<-sapply(db_terms[[db]],length)
		nGeneOfInterestByterm<-sapply( db_terms[[db]],function(term){
			return(length(which(new_x[term])))
		})
		results[[db]]<-data.frame(row.names = names(db_terms[[db]]))
		results[[db]]$pval<-phyper(q = nGeneOfInterestByterm-0.5, m = nInterest,n = nNotInterest, k = nGeneByterm, lower.tail=FALSE)
		results[[db]]$padj<-p.adjust(results[[db]]$pval,method = "BH")
		results[[db]]$nGeneOfInterest<-nGeneOfInterestByterm
		results[[db]]$nGene<-nGeneByterm
		if(returnGene){
			results[[db]]$Genes<- sapply(db_terms[[db]],function(term){
				term[term%in%geneEntrez[new_x]]
			});
			results[[db]]$Genes<- sapply(results[[db]]$Genes,ConvertKey,tabKey=corrIdGenes,colOldKey = "ENTREZID",colNewKey = "SYMBOL")
		}
	}
	return(results)
}


###ViewKEGG####
viewKEGG<-function(x,pathway,corrIdGenes=NULL,species="Human"){	
	if(is.data.frame(x) | is.matrix(x)){
		tempx<-x
		x<-tempx[,1]
		names(x)<-rownames(tempx)
	}
	if(is.null(corrIdGenes)) corrIdGenes<-getSpeciesData(species,names(x))$GeneIdTable
	entrezId<-corrIdGenes[corrIdGenes$SYMBOL%in%names(x),"ENTREZID"];
	notNA<-which(!is.na(entrezId))
	dat<-x[notNA];
	names(dat)<-entrezId[notNA]
	pathview(gene.data = dat, pathway.id = pathway, species = "hsa",kegg.native=TRUE,low="blue",mid="white",high="red",na.col="black")
}

####Convert 2 Entrez
Sym2Entrez<-function(x, corrIdGenes=NULL,species="Human"){
	require(AnnotationDbi)
	select<-AnnotationDbi::select
	validDBs<-c("kegg","reactom","goBP","goCC","goMF","custom")
	options(warn=-1)
	if(is.null(corrIdGenes)) corrIdGenes<-getSpeciesData(species,x)$GeneIdTable
	return(ConvertKey(x,tabKey = corrIdGenes,colOldKey = "SYMBOL",colNewKey = "ENTREZID"))
}
	

###Fonction pour retrouver les branch "Heatmap branch" de monocle
retrieveBranch<-function(cds,branch_point){
  require(monocle)
  require(igraph)
  pr_graph_cell_proj_mst <- minSpanningTree(cds)
  
  root_cell <- cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell
  root_state <- pData(cds)[root_cell, ]$State
  pr_graph_root <- subset(pData(cds), State == root_state)
  
  closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  root_cell_point_in_Y <- closest_vertex[row.names(pr_graph_root), ]
  
  root_cell <- names(which(igraph::degree(pr_graph_cell_proj_mst, v = root_cell_point_in_Y, 
                                  mode = "all") == 1, useNames = T))[1]
  paths_to_root <- list()
  
  pr_graph_cell_proj_mst <- minSpanningTree(cds)
  
  mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
  branch_cell <- mst_branch_nodes[branch_point] #Récupérer nom échantillon branch point
  mst_no_branch_point <- pr_graph_cell_proj_mst - V(pr_graph_cell_proj_mst)[branch_cell]
  path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst, 
                                     branch_cell, root_cell)
  path_to_ancestor <- names(unlist(path_to_ancestor$vpath))
  for (backbone_nei in V(pr_graph_cell_proj_mst)[suppressWarnings(nei(branch_cell))]$name) {
    descendents <- bfs(mst_no_branch_point, V(mst_no_branch_point)[backbone_nei], 
                       unreachable = FALSE)
    descendents <- descendents$order[!is.na(descendents$order)]
    descendents <- V(mst_no_branch_point)[descendents]$name
    if (root_cell %in% descendents == FALSE) {
      path_to_root <- unique(c(path_to_ancestor, branch_cell, 
                               descendents))
      closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
      path_to_root <- row.names(closest_vertex)[which(V(pr_graph_cell_proj_mst)[closest_vertex]$name %in% path_to_root)]
      closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
      path_to_root <- intersect(path_to_root, colnames(cds))
      paths_to_root[[backbone_nei]] <- path_to_root
    }
  }
  return(paths_to_root)
}

unsupervisedClustering<-function(x,transpose=TRUE,method.dist="pearson",method.hclust="ward.D2",bootstrap=FALSE,nboot=10){
	if(transpose) x<-t(x)
	if(bootstrap){
		require(pvclust)
		resClust<-pvclust(t(x),nboot=nboot,method.hclust = method.hclust,parallel = TRUE,method.dist = method.dist)$hclust
	}else{
		if(method.dist=="pearson"){
			resDist<-corrDist(x)
		}else{
			resDist<-dist(x, method = method.dist)
		}
		resClust<-stats::hclust(resDist,method = method.hclust)
	}
	return(resClust)
}