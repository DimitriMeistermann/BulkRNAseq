#!/usr/bin/env Rscript

setwd("D:/PHDwork/Scripts/BulkRNAseq/") #Put your Working Directory

library(BiocParallel)
library(foreach)
library(data.table)
library(gage)
library(AnnotationDbi)
library(stats4)
library(BiocGenerics)
library(parallel)
library(DESeq2)
library(sva)
library(ggplot2)
library(AnnotationDbi)
library(fgsea)
library(ggrepel)
library(rjson)
library(circlize)
library(ComplexHeatmap)
library(ggbeeswarm)
library(uwot)
library(dbscan)
library(MASS)
library(simplifyEnrichment)
library(WGCNA)
library(qualpalr)
library(reshape2)
library(RJSONIO)
library(stringr)
library(doParallel)
library(scales)
library(pathview)

#multithread init
BPPARAM = bpparam()
register(BPPARAM = BPPARAM)
registerDoParallel(cores = BPPARAM$workers)
wd=getwd()


######################################
#### Configuration of parameters #####
######################################

# /!\ Please verify carefully this section /!\ #

#expression file path
expressionData<-"data/rawCounts.tsv" 

#sample annotation table path
sampleTable<-"data/sampleAnnot.tsv" 

#Name of column in sample annotation to use for batch correction, can be set to NULL if there is no batch correction to perform.
batchColumn<-"run" 

#[Optional] Path of file that is containing 
comparisonToDoFile<-"data/ComparisonToDo.tsv"

#name of conditions column in sample table for DE genes (design of the experiment). Several values can be provided
#example with test data : c("culture_media",line")
condColumns<-NULL

#[Optional] other columns used for plots and analysis but not included in experimental design, can be set to NULL
otherInterestingColumn<-c("passage")

#[Optional] seed of the Random Nuomber Genrator (RNG). If it is fixed by an integer, results of the script will be reproducible (See ?Random).
# Can be set to NULL. If so, some results will be slightly random.
# It can be useful to turn it to NULL for being sure that your results are biased by an exceptional event in the analysis.
randomSeed<-666

#[Optional] json file with color scales, can be set to NULL.
# It can be incomplete: color scales that are not provided in this file will be computed automatically.
# For factors, order of the levels will be kept from this file
colorScaleFile="data/colorScales.json" 

#data(bods); print(bods) #to see available species
sample.species<-"Human"

#Benjamini & Hochberg pvalue threshold to consider a gene as differentially expressed
padjThreshold<-0.05

#adjusted p-val threshold for enrichment test
enrichThreshold<-0.05 

#Absolute Log2(Fold-Change) threshold (if LFCthreshold=1, gene is differentially expressed if expressed 2 time more or less between folds)
LFCthreshold<-1 

#update org.species database ? It is recommended to change this parameter to TRUE one time by semester
updateSpeciesPackage<-FALSE

#Minimum number of UMI in a sample to be kept in the analysis
minimumTotalCount<-200000 

#Minimum number of expressed genes in a sample to be kept in the analysis
minimumExpressedGenes<-5000 

#Minimum mean count of a gene to be kept in the analysis
minimumMeanCounts <- 0.1

topGeneShown = 50

#bootstrap unsupervised clustering ? bootstrap increase a lot computing time and may generate errors if N is too small
bootstrap<-FALSE 

#if bootsrap, number of bootstrap replications
nboot=30 

#Number of max pathway heatmap by comparison
pathwayInFig<-50


######################
#### Loading Data ####
######################
print("Loading data...")
source("https://raw.githubusercontent.com/DimitriMeistermann/veneR/main/loadFun.R") #Importing home made functions

### prepocess ###
if(!is.null(randomSeed)) set.seed(randomSeed)
if(is.null(condColumns) & is.null(comparisonToDoFile)) stop("You must provide a value for condColumn or ComparisonToDoFile")
#Creating dirs
dir.create("results",showWarnings=F)
dir.create("figs",showWarnings=F)
dir.create("rsave",showWarnings=F)
dir.create("resPerComparison",showWarnings=F)
#Loading files
rawCounts<-fastRead(expressionData,as.matrix = TRUE)
sampleAnnot<-fastRead(sampleTable,stringsAsFactors = TRUE)
rownames(sampleAnnot)<-make.names(rownames(sampleAnnot),unique = TRUE)
batch<-TRUE; if(is.null(batchColumn)) batch<-FALSE
for(column in c(condColumns,batchColumn)) sampleAnnot[,column]<-as.factor(as.character(sampleAnnot[,column]))# to be sure that condColums are factors

#Check if there is problems in sample names
if(sum(!rn(sampleAnnot)%in%cn(rawCounts))>0){warning("Error, these samples don't exist in expression data, please check sample names"); print(rn(sampleAnnot)[!rn(sampleAnnot)%in%cn(rawCounts)])}
if(sum(!cn(rawCounts)%in%rn(sampleAnnot))>0){warning("these samples don't exist in sample table, they will be removed from expression data"); print(cn(rawCounts)[!cn(rawCounts)%in%rn(sampleAnnot)])}

print("Processing quality control...")

### QC of raw counts table ###
sampleQCstats<-examineRNAseqSamples(rawCounts);fastWrite(sampleQCstats,"results/QualityControlOfSamples.tsv")
pdf("figs/QualityControlOfSamples.pdf",width = 9.5,height = 9);print(
	ggplot(sampleQCstats,mapping = aes(x=TotalCount,y=TotalGenEx,label=rn(sampleQCstats)))+
		geom_vline(xintercept = minimumTotalCount,color="red",lwd=1)+
		geom_hline(yintercept = minimumExpressedGenes,color="red",lwd=1)+
		geom_point()+
		xlab("Number of counts (Log10 scale)")+
		ylab("Number of expressed genes")+
		scale_x_log10()+
		geom_text_repel()+
		ggtitle("QC of samples, choosen thresholds in red")
);dev.off()

geneQCstats<-data.frame(mean=rowMeans(rawCounts),cv2=apply(rawCounts,1,cv2))
fastWrite(geneQCstats,"results/QualityControlOfGenes.tsv")

pdf("figs/QualityControlOfGenes.pdf",width = 9.5,height = 9);print(
	ggplot(geneQCstats[geneQCstats$mean>0,],mapping = aes(x=mean,y=cv2))+
		geom_vline(xintercept = minimumMeanCounts,color="red",lwd=1)+
		geom_point()+
		scale_x_log10()+
		xlab("Mean (log10 scale)")+
		ylab("Squared coefficient of variation")+
		ggtitle("QC of genes with at least one count, minimum mean counts in red")
);dev.off()

sampleAnnot <- cbind(sampleAnnot,sampleQCstats)

sample2keep<-rn(sampleQCstats)[sampleQCstats$TotalGenEx>minimumExpressedGenes & sampleQCstats$TotalCount>minimumTotalCount]
sample2rm<-rn(sampleAnnot)[!rn(sampleAnnot)%in%sample2keep]; if(len(sample2rm)>0){warning("these samples don't pass the quality control, they will be removed from analysis"); print(sample2rm)}
sampleAnnot<-sampleAnnot[sample2keep,]
rawCounts<-rawCounts[geneQCstats$mean>minimumMeanCounts,rn(sampleAnnot)]

print("Creating various stuff for analyses...")

### process of comparisons matrix (compMatrix) ###
if(!is.null(comparisonToDoFile)){
	comparisonToDoTable<-fastRead(comparisonToDoFile,row.names = NULL)
	compMatrix<-apply(comparisonToDoTable,1,function(row){
		if(!row["condColumn"]%in%cn(sampleAnnot)) stop(paste0(row["condColumn"]," does not exist in columns of sample annotation"))
		for(value in c(row["downLevel"],row["upLevel"])){
			if(!value%in%levels(sampleAnnot[,row["condColumn"]])) stop(paste0(value," is not a level of ",row["condColumn"]))
		}
		return(c(row["condColumn"],row["downLevel"],row["upLevel"]))
	})
	condColumns<-unique(comparisonToDoTable$condColumn)
	rm(comparisonToDoTable)
}else{
	compMatrix<-do.call("cbind",lapply(condColumns,function(condColumn){
		mat<-combn(levels(sampleAnnot[,condColumn]),2)
		return(rbind(rep(condColumn,ncol(mat)),mat))
	}))
	rownames(compMatrix)<-c("condColumn","downLevel","upLevel")
}
colnames(compMatrix)<-make.names(apply(compMatrix,2,function(x) paste(x["condColumn"],x["downLevel"],x["upLevel"],sep = "_")))
comparisonToDo<-colnames(compMatrix)
for(comp in comparisonToDo) dir.create(paste0("resPerComparison/",comp),showWarnings=F)

### color scales generation ###
annot2Plot <- c(condColumns, batchColumn,otherInterestingColumn)
colorScalesToGen <- annot2Plot
colorScales<-list()

if(!is.null(colorScaleFile)){
	colorScalesFromFile <- lapply(rjson::fromJSON(file = colorScaleFile),function(lvl) unlist(lvl))
	for(colorScaleName in names(colorScalesFromFile)){
		if(!colorScaleName%in%colorScalesToGen) stop("Condition '",colorScaleName,"' does not match with existing condition names in ",colorScaleFile)
		colorScale<-colorScalesFromFile[[colorScaleName]]
		annotVect<-sampleAnnot[,colorScaleName]
		if(!is.null(names(colorScale))){ #factors
			if(is.numeric(annotVect)){
				warning(colorScaleName," is numeric but encoded as factors (color vector has names) in ",colorScaleFile,". It will be converted to factors.")
				sampleAnnot[,colorScaleName]<-as.factor(as.character(sampleAnnot[,colorScaleName]))
				annotVect<-sampleAnnot[,colorScaleName]
			}else if(!is.factor(annotVect)){
				stop(colorScaleName," is not factors or numeric, please check the sample annotation table.")
			}
			if(sum(!levels(annotVect) %in% names(colorScale))>0) stop("Levels of ",colorScaleName," are existing in sample annotation table but not in ",colorScaleFile)
			sampleAnnot[,colorScaleName]<-factor(annotVect,levels = names(colorScale)) #reordering levels
			colorScales[[colorScaleName]]<-colorScale
		}else{ #numeric
			if(!is.numeric(annotVect)) stop(colorScaleName," is not numeric but encoded as numeric (color vector has no names) in ",colorScaleFile)
			colorScales[[colorScaleName]]<-colorScale
		}
	}
	colorScalesToGen<-setdiff(colorScalesToGen,names(colorScalesFromFile))
	rm(colorScalesFromFile)
}

continuousPalettes<-list(
	c("#440154","#6BA75B","#FDE725"),
	c("#2EB538","#1D1D1B","#DC0900"),
	c("#FFFFC6","#FF821B","#950961")
)

i<-1;for(colorScaleName in colorScalesToGen){
	annotVect<-sampleAnnot[,colorScaleName]
	if(is.numeric(annotVect)){
		colorScales[[colorScaleName]]<-continuousPalettes[[i]]
		i<-i+1
	}else{
		sampleAnnot[,colorScaleName]<-as.factor(as.character(sampleAnnot[,colorScaleName]))
		annotVect<-sampleAnnot[,colorScaleName]
		colorScales[[colorScaleName]]<-mostDistantColor2(nlevels(annotVect))
		names(colorScales[[colorScaleName]])<-levels(annotVect)
	}
}

#Species preparation
species.data<-getSpeciesData2(sample.species)


###########################################
#### DESeq, generation of counts table ####
###########################################
print("Computing DESeq2 model...")

formulaChar<-paste0("~",paste(condColumns,collapse = "+"))
if(batch){
	sampleAnnot[,batchColumn]<-as.factor(as.character(sampleAnnot[,batchColumn]))
	formulaChar<-paste0(formulaChar,"+",batchColumn)
}


dds <- DESeqDataSetFromMatrix(countData = rawCounts,
                                    colData = sampleAnnot,
                                    design = formula(formulaChar)) 

dds <- DESeq(dds,parallel=TRUE)

print("Correcting batch effect & saving counts tables...")

normCounts<-counts(dds,normalize=TRUE)
logCounts<-assay(vst(dds))

logCounts<-logCounts-min(logCounts) #so 0 is 0

fastWrite(normCounts,"results/normCounts.tsv")
fastWrite(logCounts,"results/logCounts.tsv")

if(batch){
  uncorrectedCounts<-logCounts
	logCounts<-ComBat(uncorrectedCounts,batch = sampleAnnot[,batchColumn],	
										mod = model.matrix(data=sampleAnnot,formula(paste0("~",paste(condColumns,collapse="+")))),BPPARAM = bpparam())
  
	logCounts<-scaleRangePerRow(logCounts,uncorrectedCounts)
  #each gene expression of the corrected values was subtracted by the minimum of the gene 
  #expression before the batch correction. 
  #This step does not change the relative expression of genes; 
  #however, it permits an easier interpretation of the expression values as minimums cannot be less than zero.
	
  normCounts<-2^logCounts - 1
  fastWrite(normCounts,"results/normCorrectedCounts.tsv")
  fastWrite(logCounts,"results/logCorrectedCounts.tsv")
}

fastWrite(CPM(rawCounts),"results/countsPerMillion.tsv")

print("Checkpoint #1")
save.image("rsave/Step1.RData")


#############################
####Unsupervised analysis####
#############################
#Over-dispersion plot
print("Computing overdispersion plot...")

geneAnnot<-getMostVariableGenes4(normCounts,minCount = 0,plot = FALSE)

genesOrderedByDispersion<-rn(geneAnnot)[order(geneAnnot$residuals,decreasing = TRUE)]

pdf("figs/overDispersionPlot.pdf",width = 16,height = 9)
ggplot(geneAnnot,aes(x=mu,y=cv2,label=rownames(dispTable),fill=residuals))+
	geom_point(stroke=1/8,colour = "black",shape=21)+geom_line(aes(y=fitted),color="red",size=1.5)+
	scale_x_log10()+scale_y_log10()+xlab("mean")+ylab("Squared coefficient of variation")+
	computeColorScaleFun(c("#440154FF","white","#FF7D1B"),geneAnnot$residuals,useProb = F,
											 midColorIs0 = TRUE,returnGGscale = TRUE,geomAes = "fill")+
	geom_text_repel(data = geneAnnot[genesOrderedByDispersion[1:topGeneShown],],inherit.aes = FALSE,
								mapping = aes(x=mu,y=cv2,label=genesOrderedByDispersion[1:topGeneShown]),size=3,fontface="bold.italic")
dev.off()

### Principal Component Analysis (PCA) ###
print("Computing PCA...")
pca<-PCA(logCounts)
if(batch){ #batch effect control
	condColumn<-condColumns[1] #ploT only 1st condition column for QC
	pcaUncorrected <- PCA(uncorrectedCounts)
	g1<-pca2d(pcaUncorrected,comp=c(1,2),group = sampleAnnot[condColumn],colorScales =  colorScales[[condColumn]],
						main="PCA on normalized/transformed counts",returnGraph = TRUE)
	g2<-pca2d(pcaUncorrected,comp=c(1,2),group = sampleAnnot[batchColumn],colorScales = colorScales[[batchColumn]],
						main="PCA on normalized/transformed counts",returnGraph = TRUE)
	g3<-pca2d(pca,comp=c(1,2),group = sampleAnnot[condColumn],colorScales = colorScales[[condColumn]],
						main="PCA on normalized/transformed/corrected counts",returnGraph = TRUE)
	g4<-pca2d(pca,comp=c(1,2),group = sampleAnnot[batchColumn],colorScales = colorScales[[batchColumn]],
						main="PCA on normalized/transformed/corrected counts",returnGraph = TRUE)
	pdf("figs/BatchEffectCorrectionControl.pdf",width = 10,height = 9)
	multiplot(g1,g2,g3,g4,cols = 2)
	dev.off()
	rm(pcaUncorrected)
}


visualizedComponent<-4
fixedCoord=FALSE #is X vs Y scale proportional to explained variance ?
plotLabelRepel=TRUE #show sample names ?

pdf(file = "figs/PrincipalComponentAnalysis.pdf",width=10,height=10)
barplotPercentVar(pca)
for(i in 1:(visualizedComponent-1)){
  for(j in (i+1):visualizedComponent){
  	for(annot in annot2Plot){
  		pca2d(pca,group = sampleAnnot[annot],pointSize = 4,comp = c(i,j),main="Principal Component Analysis",
  					fixedCoord = fixedCoord,colorScales = colorScales[[annot]],plotLabelRepel=plotLabelRepel)
  	}
  	pca2d(pca,pointSize = 1.2,comp = c(i,j),plotVars = TRUE, outlierLabel = TRUE,
  				fixedCoord = fixedCoord,colorScales = condColors,main="Correlation circle")
  	#change outlierLabelThres (between 0 and 1 to display names of more or less genes)
  }
}

dev.off()

fastWrite(pca$rotation,file="results/contribGenesPCA.tsv")

###Principal Component Regression (PCR)###
pcReg<-PCR(pca,sampleAnnot[annot2Plot],nComponent = 10)

pdf("figs/PrincipalComponentRegression.pdf",width = 8,height = 5);print(
ggBorderedFactors(ggplot(pcReg,aes(x=PC,y=Rsquared,fill=Annotation))+
	geom_beeswarm(pch=21,size=4,cex = 3)+
	xlab("Principal component")+ylab("R²")+
	scale_fill_manual(values = mostDistantColor2(length(annot2Plot)))+
	theme(
		panel.grid.major.y = element_line(colour = "grey75"),
		panel.grid.minor.y = element_line(colour = "grey75"),
		panel.background = element_rect(fill = NA,colour="black")
	)
)
);dev.off()

### Correlation heatmap ###

print("Computing correlation heatmap of samples...")
sampleCorrelations<-cor(logCounts)

pdf(file = "figs/HeatmapCorPearson.pdf",width=10,height=9)
heatmap.DM3(sampleCorrelations,preSet = "correlation",center = FALSE,midColorIs0 = FALSE,
						sampleAnnot = sampleAnnot[annot2Plot],colorAnnot = colorScales)

if(batch) heatmap.DM3(cor(uncorrectedCounts),preSet = "correlation",center = FALSE,midColorIs0 = FALSE,
											sampleAnnot = sampleAnnot[annot2Plot],colorAnnot = colorScales,
											name = "Pearson\ncorrelation\n(no batch\neffect correction)")
dev.off()

### UMAP & Clustering ###
print("Computing UMAP and unsupervised gene/sample clusters...")
selectedGenes<-rn(geneAnnot)[geneAnnot$residuals>0]

umap<-make.umap2(logCounts[selectedGenes,],n_epochs=2000)
umapHighDim<-make.umap2(logCounts[selectedGenes,],
												n_components = 10,n_epochs=2000,n_neighbors = min(20,nrow(sampleAnnot)))

sampleClustering=dbscan::hdbscan(umapHighDim,minPts = 2)
sampleAnnot$sampleClusters <- as.factor(paste0("k",sampleClustering$cluster))

if(is.null(colorScales$sampleClusters)) colorScales$sampleClusters <- genColorsForAnnots(sampleAnnot["sampleClusters"])$sampleClusters

pdf("figs/UMAP.pdf",width = 8,height = 7)
proj2d(umap,group=sampleAnnot["sampleClusters"],colorScale = colorScales$sampleClusters,plotFactorsCentroids = TRUE,fixedCoord = T)
for(annot  in annot2Plot) proj2d(umap,group=sampleAnnot[annot],colorScale = colorScales[[annot]],fixedCoord = T)
dev.off()

geneClustering=hierarchicalClustering(logCounts[selectedGenes,],transpose = FALSE,method.dist = "bicor")
geneCluster<-cutree(geneClustering,k = best.cutree2(geneClustering,min = 4))

geneAnnot$Module<-"M0";geneAnnot[names(geneCluster),"Module"]<-paste0("M",alphabetOrderByAdding0(geneCluster))
geneAnnot$Module<-as.factor(geneAnnot$Module)
trueGeneModules<-levels(geneAnnot$Module)[levels(geneAnnot$Module)!="M0"]


moduleActivationScore<-data.frame(sapply(trueGeneModules,function(m){
	eigengenes(logCounts,genes = rn(geneAnnot)[geneAnnot$Module==m ] )
}));names(moduleActivationScore)<-trueGeneModules

moduleMembership <- suppressWarnings(bicor(t(logCounts), moduleActivationScore));
colnames(moduleMembership) <- trueGeneModules

geneAnnot$Membership<-NA
for(gene in rn(geneAnnot)[geneAnnot$Module!="M0"]) {
	geneAnnot[gene,"Membership"]<-moduleMembership[gene,as.character(geneAnnot[gene,"Module"])]
}

pdf("figs/moduleActivationScore.pdf",width = 10,height = 10)
heatmap.DM3(t(moduleActivationScore),scale = F,preSet = "default",
						midColorIs0 = TRUE,sampleAnnot = sampleAnnot[annot2Plot],colorAnnot = colorScales[annot2Plot],
						name="gene\nmodule\nactivation",column_split = sampleAnnot$sampleClusters)
dev.off()

AUCs<-getMarkers3(logCounts,sampleAnnot$sampleClusters,BPPARAM=BPPARAM)

bestMarkersPerCluster<-VectorListToFactor(lapply(data.frame(AUCs),function(auc){
	rn(AUCs)[order(auc,decreasing = TRUE)][1:topGeneShown]
}))

pdf("figs/markersOfUnsupervisedClusters.pdf",width = 10,height = 10)
heatmap.DM3(logCounts[names(bestMarkersPerCluster),],sampleAnnot = sampleAnnot[annot2Plot],colorAnnot = colorScales[annot2Plot],
						column_split = sampleAnnot$sampleClusters,row_split=bestMarkersPerCluster,clustering_distance_columns = "pearson")
dev.off()

fastWrite(AUCs,"results/markerAUCsPerSampleCluster.tsv")
fastWrite(moduleActivationScore,"results/moduleActivationScore.tsv")
fastWrite(moduleMembership,"results/moduleMembership.tsv")
fastWrite(geneAnnot,"results/geneAnnotations.tsv")
fastWrite(sampleAnnot,"results/sampleAnnotations.tsv")

###Super Heatmap###
print("Computing the super heatmap...")
htList<-list()

genes<-rn(geneAnnot)[geneAnnot$Module=="M0"]
ht<-heatmap.DM3(logCounts[genes,],showGrid = F,
													 use_raster = TRUE,raster_quality = 5,returnHeatmap = TRUE,
													 column_split = sampleAnnot$sampleClusters,
													 height = unit(log2(len(genes))/5, "cm"),
													 cluster_row_slices = FALSE,sampleAnnot = sampleAnnot[annot2Plot],colorAnnot = colorScales[annot2Plot],
													 row_title_gp=gpar(fontsize=9),column_title_gp=gpar(fontsize=9),
													 cluster_column_slices = FALSE,border = TRUE,name="expression",row_title_rot = 0,
													 show_row_names = F,show_column_names = F,row_title = paste0("M0\n",length(genes)," genes"))

for(mod in trueGeneModules){
	genes<-rn(geneAnnot)[geneAnnot$Module==mod]
	htList[[mod]]<-heatmap.DM3(logCounts[genes,],showGrid = F,
									 use_raster = TRUE,raster_quality = 5,returnHeatmap = TRUE,
									 column_split = sampleAnnot$sampleClusters,
									 height = unit(log2(len(genes))/5, "cm"),cluster_row_slices = FALSE,
									 row_title_gp=gpar(fontsize=9),column_title_gp=gpar(fontsize=9),
									 cluster_column_slices = FALSE,border = TRUE,name="expression",row_title_rot = 0,
									 show_row_names = F,show_column_names = F,row_title = paste0(mod,"\n",length(genes)," genes")
	)
}


for(i in 1:length(htList)){
	suppressWarnings(ht<- ht %v% htList[[i]])
}

pdf("figs/superHeatmap.pdf",height = length(trueGeneModules)+1.5,width = 10)
draw(ht,clustering_distance_columns = "pearson",clustering_method_columns="ward.D2")
dev.off()

rm(htList,ht)

print("Checkpoint #2")
save.image("rsave/Step2.RData")



######################################################
#### Differentially expressed (DE) genes analysis ####
######################################################


print("Computing differentially expressed genes...")
DEresults<-foreach(comp=colnames(compMatrix)) %dopar%{ #change 'dopar' by 'do' to run in single thread mode
	library(DESeq2)
	library(ggplot2)
	condColumn<-compMatrix[1,comp]
	upLevel<-compMatrix[2,comp]
	downLevel<-compMatrix[3,comp]
	print(paste0("Computing differentially expressed genes for comparison [",condColumn,": ",downLevel," vs ",upLevel,"]"))
	
	
	DEresult<-data.frame(results(dds, contrast=c(condColumn,downLevel,upLevel), 
																				independentFiltering = TRUE,alpha = padjThreshold))
	DEresult<-cbind(data.frame(gene=rownames(DEresult),stringsAsFactors = FALSE),DEresult)
	DEresult$isDE<-"NONE"
	DEresult[!is.na(DEresult$padj) & DEresult$padj<padjThreshold & DEresult$log2FoldChange > LFCthreshold,"isDE"]<-"UPREG"
	DEresult[!is.na(DEresult$padj) & DEresult$padj<padjThreshold & DEresult$log2FoldChange < -LFCthreshold,"isDE"]<-"DOWNREG"
	DEresult$isDE<-as.factor(DEresult$isDE)
	fastWrite(DEresult,paste0("resPerComparison/",comp,"/DESeqResults.tsv"),row.names = FALSE)
	
	### P val histogram ###
	pdf(paste0("resPerComparison/",comp,"/pvalHistogram.pdf"),width = 10,height = 8)
	print(ggplot(data.frame(DEresult),mapping=aes(pvalue))+
		geom_histogram(binwidth = 0.05)+
		geom_hline(yintercept = median(table(cut(DEresult$pvalue,breaks = seq(0,1,0.05),include.lowest = FALSE))))
	)
	dev.off()
	
	### Volcano plots ###
	pdf(paste0("resPerComparison/",comp,"/volcanoPlot.pdf"),width = 9,height = 8)
	volcanoPlot.DESeq2(DEresult = DEresult,formula = formulaChar,condColumn = condColumn,downLevel = downLevel,upLevel = upLevel,
										 padjThreshold = padjThreshold,LFCthreshold = LFCthreshold,topGene = topGeneShown)
	dev.off()
		
	### MA plot ###
	gene2Plot<-order(DEresult$padj)
	gene2Plot<-gene2Plot[!is.na(DEresult[gene2Plot,"isDE"])][1:topGeneShown]
	
	pdf(paste0("resPerComparison/",comp,"/MA-plot.pdf"),width = 10,height = 8)
	print(ggplot(data=DEresult,mapping = aes(x=baseMean,y=log2FoldChange,colour=isDE))+
		geom_point(size=1)+
		scale_x_log10()+
		scale_color_manual(values = c("#3AAA35","grey75","#E40429"),guide=FALSE)+
		theme(
			panel.background = element_rect(fill = NA,colour="black"),
			panel.grid.major = element_line(colour = NA),
			panel.grid.minor = element_line(colour = NA)
		)+
		geom_text_repel(data = DEresult[gene2Plot,],mapping = aes(x=baseMean,y=log2FoldChange,label=gene),
										size=3,inherit.aes = FALSE,fontface="bold.italic")+
		ggtitle(paste0("MA-plot for comparison ",compMatrix[1,comp],": ",compMatrix[2,comp]," vs ",compMatrix[3,comp]))+
		xlab("Mean expression")+ylab("log2(Fold-Change)"))
	dev.off()
	
	pdf(paste0("resPerComparison/",comp,"/focusOnTop10genes.pdf"),width = 12,height = 8)
	plotExpression(normCounts[gene2Plot[1:10],],sampleAnnot[condColumn],colorScale = colorScales[[condColumn]],printGraph = TRUE)
	dev.off()
	
	return(DEresult)
};names(DEresults)<-colnames(compMatrix)

DEgenesPerComparison<-lapply(DEresults,function(DEresult){ DEresult$gene[DEresult$isDE!="NONE"] })

### DE genes heatmaps###

print("Computing heatmaps and upset plot of differentially expressed genes...")

#Select comparisons with number of DE > 0
significantComparisons<-c(); for(comp in comparisonToDo) if(sum(DEresults[[comp]]$isDE!="NONE")>1) significantComparisons<-c(significantComparisons,comp)

for (comp in significantComparisons){
	DEgenes<-DEgenesPerComparison[[comp]]
	pdf(paste0("resPerComparison/",comp,"/heatmapDEgenes.pdf"),width = 10,height = 10)
	heatmap.DM3(logCounts[DEgenes,],
							column_split= ifelse(sampleAnnot[,compMatrix[1,comp]] %in% c(compMatrix[2,comp],compMatrix[3,comp]),comp,"other samples"),
							sampleAnnot = sampleAnnot[annot2Plot],colorAnnot = colorScales[annot2Plot],row_split=DEresults[[comp]][DEgenes,"isDE"])
	dev.off()
}


####UpSet plot####
pdf("figs/upsetPlot.pdf",width = 10,height = 6)
customUpsetPlot(DEgenesPerComparison,universe = rownames(logCounts))
dev.off()

print("Checkpoint #3")
save.image("rsave/Step3.RData")

###############################
#### Functional Enrichment ####
###############################

print("Downloading gene set databases...")

# geneSetsDatabase=list(msigdb=gmtPathways("data/msigdb.v7.2.symbols.gmt")) #alternative for human
geneSetsDatabase<-getDBterms2(rn(logCounts),speciesData = species.data,returnGenesSymbol = TRUE,database=c("kegg","goBP"))

### fGSEA ####
enrichResultsGSEAPerComp<-foreach(comp=significantComparisons) %do%{
	condColumn<-compMatrix[1,comp]
	upLevel<-compMatrix[2,comp]
	downLevel<-compMatrix[3,comp]
	print(paste0("Computing GSEA enrichment for comparison [",condColumn,": ",downLevel," vs ",upLevel,"]"))
	
	DEres<-DEresults[[comp]]
  scoreEnrich<-calConsensusRanking(rownames(DEres),DEres$pvalue,DEres$log2FoldChange)
  
  enrichResultsGSEA<-Enrich.gsea(scoreEnrich,speciesData = species.data,db_terms = geneSetsDatabase,returnGenes = TRUE)
  

  enrichResultsGSEA$significant<-FALSE
  enrichResultsGSEA[enrichResultsGSEA$padj<enrichThreshold,"significant"]<-TRUE
  
  enrichResultsGSEA<-data.frame(enrichResultsGSEA[,1:(ncol(enrichResultsGSEA)-2)],enrichResultsGSEA[,"significant"],enrichResultsGSEA[,"genes"]) #reorder columns
  
  exportEnrich.2(enrichResultsGSEA,paste0("resPerComparison/",comp,"/GSEA_enrichment.tsv"))
  
  pdf(paste0("resPerComparison/",comp,"/GSEA_pvalHistogram.pdf"),width = 6 ,height = 5)
  print(ggplot(enrichResultsGSEA,mapping=aes(pval))+
  	geom_histogram(binwidth = 0.05)+
  	geom_hline(yintercept = median(table(cut(enrichResultsGSEA$pval,breaks = seq(0,1,0.05),include.lowest = FALSE)))))
  dev.off()
  
  
  pdf(paste0("resPerComparison/",comp,"/GSEA_volcano.pdf"),width = 11 ,height = 10)
  print(ggplot(enrichResultsGSEA,aes(x=NES,y=-log10(padj),color=significant))+
  	geom_point()+theme_bw()+scale_color_manual(values=c("grey75","black"))+
  	geom_text_repel(data = enrichResultsGSEA[order(enrichResultsGSEA$pval)[1:15],],aes(x=NES,y=-log10(padj),label=pathway),
  									inherit.aes = FALSE,color="grey50"))
  dev.off()
  
  
  if(sum(enrichResultsGSEA$significant)>0){
  	pathwayToPlot<-enrichResultsGSEA[enrichResultsGSEA$significant,]
  	pathwayToPlot<-pathwayToPlot[order(pathwayToPlot$pval),]
  	pathwayToPlot<-pathwayToPlot[1:min(nrow(pathwayToPlot),pathwayInFig),]
  	
  	pdf(paste0("resPerComparison/",comp,"/GSEA_heatmapPerPathway.pdf"),width = 10,height = 15,useDingbats = FALSE)
  	for(i in 1:nrow(pathwayToPlot)){
  		term<-pathwayToPlot[i,"pathway"]
  		genesOfPathway<-pathwayToPlot[i,]$genes[[1]]
  		genesNotNA<-genesOfPathway[genesOfPathway%in%rownames(logCounts)]
  		db<-pathwayToPlot[i,"database"]
  		if(!length(genesNotNA)>2) next
  		heatmap.DM3(logCounts[genesNotNA,],
				column_split = ifelse(sampleAnnot[,compMatrix[1,comp]] %in% c(compMatrix[2,comp],compMatrix[3,comp]),comp,"other samples"),
				sampleAnnot = sampleAnnot[annot2Plot],colorAnnot = colorScales[annot2Plot],
				row_title=paste0(term,"\nFrom ",db,", ",length(genesNotNA),"/",length(genesOfPathway)," genes detected")
  		)
  	}
  	dev.off()
  	
  	if("kegg"%in%enrichResultsGSEA$database){
  		dirOfPathview<-paste0("resPerComparison/",comp,"/pathview")
  		dir.create(dirOfPathview,showWarnings = FALSE)
  		keggPathway<-unlist(enrichResultsGSEA[enrichResultsGSEA$database=="kegg" & enrichResultsGSEA$significant,"pathway"])
  		logFC<-DEres$log2FoldChange;names(logFC)<-rownames(DEres)
  		setwd(dirOfPathview);for(pathway in keggPathway) viewKEGG(logFC,pathway,speciesData=species.data);setwd(wd)
  	}
  }
  return(enrichResultsGSEA)
}


### Gene Set Differential Activation (experimental) ####
print("Preparing GSDA enrichment...")

geneSetsEigens<-computeActivationScore(logCounts,db_terms = geneSetsDatabase)

for(db in names(geneSetsEigens)){
	fastWrite(geneSetsEigens[[db]]$eigen,paste0("results/geneSetActivations_",db,".tsv"))
	write.vectorList(geneSetsEigens[[db]]$contribution,paste0("results/geneSetsContributions_",db,".tsv"))
}

for(comp in significantComparisons){
	condColumn<-compMatrix[1,comp]
	upLevel<-compMatrix[2,comp]
	downLevel<-compMatrix[3,comp]
	print(paste0("Computing GSDA enrichment for comparison [",condColumn,": ",downLevel," vs ",upLevel,"]"))
	
	enrichResultsGSDA<-GSDA(geneSetEigens = geneSetsEigens,colData = sampleAnnot,contrast = compMatrix[,comp],
													db_terms = geneSetsDatabase)
	pdf(paste0("resPerComparison/",comp,"/GSDA_pvalHistogram.pdf"),width = 6 ,height = 5)
	print(ggplot(enrichResultsGSDA,mapping=aes(pval))+
		geom_histogram(binwidth = 0.05)+
		geom_hline(yintercept = median(table(cut(enrichResultsGSDA$pval,breaks = seq(0,1,0.05),include.lowest = FALSE))))
	)
	dev.off()
	
	enrichResultsGSDA$significant<-FALSE
	enrichResultsGSDA[enrichResultsGSDA$padj<enrichThreshold,"significant"]<-TRUE
	
	fastWrite(enrichResultsGSDA,paste0("resPerComparison/",comp,"/GSDA_enrichment.tsv"),row.names = FALSE)
	
	pdf(paste0("resPerComparison/",comp,"/GSDA_Volcano.pdf"),width = 11 ,height = 10)
	print(ggplot(enrichResultsGSDA,aes(x=log2FoldChange,y=-log10(padj),color=significant))+
		geom_point()+theme_bw()+scale_color_manual(values=c("grey75","black"))+
		geom_text_repel(data = enrichResultsGSDA[order(enrichResultsGSDA$pval)[1:15],],aes(x=log2FoldChange,y=-log10(padj),label=term),
										inherit.aes = FALSE,color="grey50")
	)
	dev.off()
	
	if(sum(enrichResultsGSDA$significant)>1){
		selectedTerms<-enrichResultsGSDA[enrichResultsGSDA$significant,]
		selectedTerms<-selectedTerms[order(selectedTerms$pval),]
		selectedTerms<-selectedTerms[1:min(nrow(selectedTerms),pathwayInFig),] #use top gene shown 
		
		selectedContrib<-apply(selectedTerms,1,function(term){
			geneSetsEigens[[ term["database"] ]]$contribution[[ term["term"] ]]
		});names(selectedContrib)<-selectedTerms$term
		selectedActivations<-t(apply(selectedTerms,1,function(term){
			geneSetsEigens[[ term["database"] ]]$eigen[term["term"],]
		}));rownames(selectedActivations)<-selectedTerms$term
		
		
		pdf(paste0("resPerComparison/",comp,"/GSDA_heatmap.pdf"),width = 25,height = 12)
		heatmap.DM3(selectedActivations,midColorIs0 = T,center=F,
								column_split= ifelse(sampleAnnot[,compMatrix[1,comp]] %in% c(compMatrix[2,comp],compMatrix[3,comp]),comp,"other samples"),
								column_title_gp=gpar(fontsize=8),
								name = "Activation score",preSet = "default",
								right_annotation=rowAnnotation("gene contribution" = GSDA.HeatmapAnnot(contributions = selectedContrib,width = unit(4,"inches"))),
								row_names_side ="left",row_dend_side ="left",sampleAnnot = sampleAnnot[annot2Plot],colorAnnot = colorScales[annot2Plot],
								row_names_max_width = unit(8, "inches"),autoFontSizeRow=FALSE,row_names_gp=gpar(fontsize=1/nrow(selectedActivations)*400))
		dev.off()
	}
}

## Over representation analysis of gene modules ##
print("Computing enrichment of gene modules...")

dir.create("results/EnrichmentPerGeneModule",showWarnings = F)
null<-foreach(module=trueGeneModules) %dopar%{
	isInModule<-geneAnnot$Module==module;names(isInModule)<-rownames(geneAnnot)
	exportEnrich.2(
		Enrich.simple(isInModule,speciesData = species.data,db_terms = geneSetsDatabase,returnGenes = TRUE),
		paste0("results/EnrichmentPerGeneModule/",module,".tsv")
	)
}

print("Checkpoint #4")
save.image("rsave/Step4.RData")

###############################
#### Generation of web app ####
###############################
print("Building web application...")

webAppFolder<-"webApp/"

dir.create(webAppFolder,showWarnings = F)
dir.create(paste0(webAppFolder,"data/geneCor"),showWarnings = F,recursive = T)
dir.create(paste0(webAppFolder,"data/gene"),showWarnings = F,recursive = T)


exprForUI<-normCounts

rownames(exprForUI)<-str_replace_all(rownames(exprForUI),"/","_")

print("Web app: writing gene expressions...")

null<-foreach(gene=rn(exprForUI)) %dopar% {
	newDat<-data.frame(expression=exprForUI[gene,])
	fastWrite(newDat,fileName =paste0(webAppFolder,"data/gene/",gene,".tsv"),row.names = F)
}

write(toJSON(rn(exprForUI)),paste0(webAppFolder,"data/geneList.json"))

print("Web app: writing gene per gene correlations...")

#export top 50 correlation
exprLogForUI<-logCounts[rownames(normCounts),]
rownames(exprLogForUI)<-rownames(exprForUI)

null<-foreach(gene = rn(exprLogForUI)) %dopar%{
	corr<-t(cor(exprLogForUI[gene,],t(exprLogForUI)))
	res<-corr[order(corr[,1],decreasing = T)[c(2:51,(nrow(corr)-50):nrow(corr))],,drop=F]
	res<-data.frame(Gene=rn(res),pearsonCor=res[,1])
	fastWrite(res,
						file=paste0(webAppFolder,"data/geneCor/",gene,".tsv"), row.names = FALSE)
}

rm(exprForUI)
rm(exprLogForUI)


print("Web app: writing some files...")

colnames(umap)<-c("x","y")

fastWrite(sampleAnnot[,names(colorScales)],paste0(webAppFolder,"data/sampleAnnot.tsv"))
fastWrite(umap,paste0(webAppFolder,"data/coorProjection.tsv"))
fastWrite(sampleCorrelations,paste0(webAppFolder,"data/corSamples.tsv"))


#auto
htmlOptionMain=""
htmlOptionViolin=""


for(col in names(colorScales)){
	htmlOptionMain=paste0(htmlOptionMain,strrep("\t",8),"<option value=\"",col,"\">",col,"</option>\n")
	if(is.factor(sampleAnnot[,col])){
		htmlOptionViolin=paste0(htmlOptionViolin,strrep("\t",10),"<option value=\"",col,"\">",col,"</option>\n")
	}
}

colorScalesUI<-colorScales[sapply(colorScales,function(x) !is.null(names(x)))]

listExportJSON<-list(domain=lapply(colorScalesUI,function(x){names(x)}),
										 range=lapply(colorScalesUI,function(x){names(x)<-NULL; x}))

write(toJSON(listExportJSON),paste0(webAppFolder,"data/colorScales.json"))

mainAnnot<-"sampleClusters"
backGroundAnnot<-"sampleClusters"
selectableAnnotation<-names(colorScalesUI[1])

htmlSkeleton<-readChar(paste0(webAppFolder,"prepare/htmlSkeleton.txt"), file.info(paste0(webAppFolder,"prepare/htmlSkeleton.txt"))$size)
htmlOut<-str_replace(htmlSkeleton,"<\\?optionMain\\?>",htmlOptionMain)
htmlOut<-str_replace(htmlOut,"<\\?optionViolin\\?>",htmlOptionViolin)
write(htmlOut,paste0(webAppFolder,"index.html"))

runUIskeleton<-readChar(paste0(webAppFolder,"prepare/runUIskeleton.txt"), file.info(paste0(webAppFolder,"prepare/runUIskeleton.txt"))$size)
runUIOut<-str_replace(runUIskeleton,"<\\?backGroundAnnot\\?>",backGroundAnnot)
runUIOut<-str_replace(runUIOut,"<\\?mainAnnot\\?>",mainAnnot)
runUIOut<-str_replace(runUIOut,"<\\?selectableAnnotation\\?>",selectableAnnotation)
write(runUIOut,paste0(webAppFolder,"js/runUI.js"))


####Gene sets
fastWrite(moduleActivationScore,paste0(webAppFolder,"data/GenesClustEigen.tsv"))
fastWrite(geneAnnot[,c("Module","Membership")],paste0(webAppFolder,"data/GenesClust.tsv"))
write(toJSON(cn(moduleActivationScore)),paste0(webAppFolder,"data/GenesClustList.json"))

####

print("Analyses are over !")
