#!/usr/bin/R

library("Biobase")
library("BiocGenerics")
library("BiocInstaller")
library(shiny)
library(flashClust)
library(DESeq2)
library(gage)
library(ggplot2)
library(fgsea)
library(Rcpp)
library(reactome.db)
library(pvclust)
library(ComplexHeatmap)
library(rgl)
library(pathview)
library(GO.db)
library(fdrtool)
library(genefilter)
library(plyr)
library(limma)
library(RColorBrewer)
library(circlize)
library(fgsea)
library(AnnotationDbi)
library(reactome.db)
library(ggrepel)

#### Config param  #####
#basic parameters (please verify carefully)
setwd("~/PHDwork/Scripts/DGE-Seq_Analyzer/") #Put your Working Directory
RequiredfunctionsDir<-"~/PHDwork/functions/" #functions used by script path
expressionData<-"data/exprDat.tsv" #expression file path
sampleTable<-"data/sampleAnnot.tsv" #sample sheet path
condCol<-"Milieu" #name of condition column in sample table for DE genes
batchCol<-"Run" #Name of batch column in sample table for batch correction, 'NULL' if no batch

#advanced parameters
sample.species<-"Human" #data(bods); print(bods) #to see available species
AdjPVaLthreshold<-0.05 #Benjamini & Hochberg pvalue threshold to consider a gene as differentially expressed
LocalFDRthreshold<-0.05 #local adjusted fdr threshold to consider a gene as differentially expressed (package fdrtool)
Enrichthreshold<-0.05 #adjusted p-val threshold for enrichment test
logFCthreshold<-1 #Absolute Log2(Fold-Change) threshold (if logFCthreshold=1, gene is differentially expressed if expressed 2 time more or less between folds)
nTopGo<-100 #n TOp Go Term/Ontology
updateSpeciesPackage<-FALSE #update org.species database ? It is recommended to change this parameter to TRUE one time by semester
minimumTotalCount<-200000 #Mininim number of UMI in a sample to keep the sample in analysis
minimumGeneExpressed<-5000 #Mininim number of expressed genes in a sample to keep the sample in analysis
conditionSeparator<-"_" #separator for comparison name generation, ex: for condition A and B with underscore, comp A vs B = A_B
bootstrap<-FALSE #bootstrap unsupervised clustering ? bootstrap increase a lot computing time and may generate errors if N is too small
nboot=30 #if bootsrap, number of bootstrap replications
allSampleInHt<-FALSE #for heatmaps based on differentialy expression, include all samples or only samples involved in comparisons ?
GenesInFig<-100 #here tape a number of genes that will correspond to top DE genes or provide a liste of genes
PathwayInFig<-10 #Number of max pathway heatmap by comparison

#####LOADING DATA#####
#Importation des fonctions maison
source(paste0(RequiredfunctionsDir,"/general.R"))
source(paste0(RequiredfunctionsDir,"/projection.R"))
source(paste0(RequiredfunctionsDir,"/RNASeq.R"))
#Create dirs
dir.create("results",showWarnings=F)
dir.create("figs",showWarnings=F)
dir.create("rsave",showWarnings=F)
#Load files
exprDat<-lire(expressionData)
sampleAnnot<-lire(sampleTable)
#Species preparation
species.data<-getSpeciesData(sample.species,rn(exprDat))

#####DESeq#####
#prepocess
sampleExprInAnnot<-rn(sampleAnnot)[!rn(sampleAnnot)%in%cn(exprDat)]; if(len(sampleExprInAnnot)>0){warning("Error, these samples don't exist in expression data, please check sample names"); print(sampleExprInAnnot)}
sampleAnnotInExpr<-cn(exprDat)[!cn(exprDat)%in%rn(sampleAnnot)]; if(len(sampleAnnotInExpr)>0){warning("these samples don't exist in sample table, they will be removed from expression data"); print(sampleAnnotInExpr)}
absSamples<-examineRNAseqSamples(exprDat);ecrire(absSamples,"results/SamplesAbstract.tsv")
sample2keep<-rn(absSamples)[absSamples$TotalGenEx>minimumGeneExpressed & absSamples$TotalCount>minimumTotalCount]
sample2rm<-rn(sampleAnnot)[!rn(sampleAnnot)%in%sample2keep]; if(len(sample2rm)>0){warning("these samples don't pass the quality control, they will be removed from analysis"); print(sample2rm)}
sampleAnnot<-sampleAnnot[sample2keep,]
exprDat<-exprDat[,rn(sampleAnnot)]

batch<-TRUE; if(is.null(batchCol)) batch<-FALSE
exprDat<-filterMostExprimedGenesBySample(exprDat) # Gene filtering
sampleAnnot[,condCol]<-as.factor(as.character(sampleAnnot[,condCol]))

formulaChar<-paste0("~",condCol)
if(batch){
	sampleAnnot[,batchCol]<-as.factor(as.character(sampleAnnot[,batchCol]))
	formulaChar<-paste0(formulaChar,"+",batchCol)
}


dds <- DESeqDataSetFromMatrix(countData = exprDat,
                                    colData = sampleAnnot,
                                    design = formula(formulaChar)) 

dds <- DESeq(dds)

exprDatN<-counts(dds,normalize=TRUE)

if(ncol(exprDat)<=20){
  exprDatT<-assay(rlog(dds))
}else{
  exprDatT<-assay(vst(dds))
}

exprDatTm<-exprDatT-min(exprDatT)

ecrire(exprDatN,"results/exprNormalized.tsv")
ecrire(exprDatT,"results/exprTransformed.tsv")

if(batch){
  exprDatNa<-removeBatchEffect(exprDatN,batch = sampleAnnot[,batchCol],design = model.matrix(data=sampleAnnot,formula(paste0("~",condCol))))
  exprDatTa<-removeBatchEffect(exprDatT,batch = sampleAnnot[,batchCol],design = model.matrix(data=sampleAnnot,formula(paste0("~",condCol))))
  exprW<-exprDatTa-min(exprDatTa)
  ecrire(exprDatNa,"results/exprNormalizedAdjusted.tsv")
  ecrire(exprDatTa,"results/exprTransformedAdjusted.tsv")
}else{
  exprW<-exprDatTm
}

samples<-colnames(exprW)
sampleHt<-samples
genes<-rownames(exprW)

save.image("rsave/Step1.R.RData")

##### Quality control#####
pdf(file = "figs/DistribCountPerSample.pdf",width=10,height=10)
if(batch){
  plotDistrib(exprDatT,main="Barplots of normalized/transformed data",conditions = sampleAnnot[,batchCol])
  plotDistrib(exprDatTa,main="Barplots of normalized/transformed/adjusted data",conditions = sampleAnnot[,batchCol])
}else{
  plotDistrib(exprDatT,main="Barplots of normalized/transformed data")
}
acpT<-ACP(exprDatT)
acp2d(acpT,comp=c(1,2),plotText = TRUE,group = sampleAnnot[,condCol],main="PCA on normalized/transformed data")
if(batch){
  acp2d(acpT,comp=c(1,2),plotText = TRUE,group = sampleAnnot[,batchCol],main="PCA on normalized/transformed data")
  acpTa<-ACP(exprDatTa)
  acp2d(acpTa,comp=c(1,2),plotText = TRUE,group = sampleAnnot[,condCol],main="PCA on normalized/transformed/adjusted data")
  acp2d(acpTa,comp=c(1,2),plotText = TRUE,group = sampleAnnot[,batchCol],main="PCA on normalized/transformed/adjusted data")
}
dev.off()


exprDat.UPM <- UMI2UPM(exprDat)
ecrire(exprDat.UPM,"results/exprDatUPM.tsv")

rowScaledExpr<-rowScale(exprW,center=TRUE,scale=TRUE)
quantile.expr<-quantile(unlist(rowScaledExpr),seq(0,1,.01))
colHA<-colorRamp2(c(quantile.expr[2],0,quantile.expr[100]),c("blue","white","red"))
#####PCA#####
compo<-4

acp<-ACP(exprW) 
pdf(file = "figs/PCA.pdf",width=10,height=10)
barplot(acp$percentVar,names.arg = round(acp$percentVar*100,2),main = "Contribution of each componant in PCA")

for(i in 1:(compo-1)){
  for(j in (i+1):compo){
    acp2d(acp,group = sampleAnnot[,condCol],plotText = TRUE,pointSize = 4,comp = c(i,j),main="PCA",fixedCoord = F)
    acp2d(acp,group = sampleAnnot[,condCol],pointSize = 2,comp = c(i,j),main="PCA",fixedCoord = F)
    acp2d(acp,pointSize = 2,comp = c(i,j),plotVars = TRUE, plotText = TRUE,fixedCoord = F)
  }
}

dev.off()

ecrire(acp$rotation,file="results/contribGenesPCA.tsv")


#####Clustering & Heatmap#####
corSample<-cor(exprW,method = "pearson")
clustSamples<-unsupervisedClustering(corSample,bootstrap = bootstrap,nboot=nboot,method.dist = "euclidean")
if(batch){
  corSampleUnadjust<-cor(exprDatT,method = "pearson")
  clustSamplesUnadjust<-unsupervisedClustering(corSampleUnadjust,nboot=nboot,bootstrap = bootstrap,method.dist = "euclidean")
}

#cah
sampleGroup<-sampleAnnot[,condCol]
names(sampleGroup)<-rownames(sampleAnnot)

if(batch){
  annot<-sampleAnnot[,c(condCol,batchCol)]
}else{
  annot<-sampleAnnot[condCol]
}

colTopAnnot<-vector("list", ncol(annot))
names(colTopAnnot)<-cn(annot)
colFun<-c(ggplotColours,rainbow)
i<-1
for(col in cn(annot)){
  colTopAnnot[[col]]<-colFun[[i]](nlevels(annot[,col]))
  names(colTopAnnot[[col]])<-levels(annot[,col])
  i<-i+1
  if(i==4) i<-1
}
ha<-HeatmapAnnotation(annot,col = colTopAnnot)
Ht<-Heatmap(matrix = corSample, cluster_rows = clustSamples, cluster_columns = clustSamples, top_annotation = ha, 
            name="Pearson\ncorrelation",row_names_gp = autoGparFontSizeMatrix(ncol(corSample)),
            column_names_gp = autoGparFontSizeMatrix(ncol(corSample)),
            row_dend_reorder = FALSE, column_dend_reorder = FALSE,col=heat.colors(100))

if(batch){
  Htna<-Heatmap(matrix = corSampleUnadjust, cluster_rows = clustSamplesUnadjust, cluster_columns = clustSamplesUnadjust, top_annotation = ha, 
              name="Pearson\ncorrelation\nnon adjusted",row_names_gp = autoGparFontSizeMatrix(ncol(corSampleUnadjust)),
              column_names_gp = autoGparFontSizeMatrix(ncol(corSampleUnadjust)),
              row_dend_reorder = FALSE, column_dend_reorder = FALSE,col=heat.colors(100))
}
Htc<-Heatmap(absSamples[clustSamples$labels,"TotalCount"], name = "Total\nexpression", 
             col=colorRamp2(c(0,max(absSamples$TotalCount)),c("black","green")), show_row_names = FALSE, width = unit(5, "mm"))

pdf(file = "figs/HeatmapCorPearson.pdf",width=10,height=9)
print(Ht+Htc)
if(batch) print(Htna+Htc)
dev.off()

save.image("rsave/Step2.R.RData")


#####Differentially expressed (DE) genes analysis#####
#Generation of each comparison
conds<-levels(sampleAnnot[,condCol])
compMatrix<-combn(conds,2)
comps<-apply(compMatrix,2,function(x)paste0(x[1],conditionSeparator,x[2]))
colnames(compMatrix)<-comps

res<-list()
sampleByComp<-list()
for(comp in comps){
  cond1<-compMatrix[1,comp]
  cond2<-compMatrix[2,comp]
  sampleByComp[[comp]]<-samples[sampleAnnot[,condCol]%in% c(cond1,cond2)]
  #DESeq results for each comparison
  res[[comp]]<-results(dds, contrast=c(condCol,cond1,cond2), independentFiltering = TRUE,alpha = min(LocalFDRthreshold,AdjPVaLthreshold))
  res[[comp]]$meanInComp<-rowMeans(exprDatN[,sampleAnnot[,condCol]%in%compMatrix[,comp]])
}


#gathering each results in a dataframe
DE<-list()
FDR.res<-list()
for(comp in comps){
  DE[[comp]]<-data.frame(res[[comp]])
  rownames(DE[[comp]])<-genes
  FDR.res[[comp]]<-fdrtool(res[[comp]]$pvalue[which(!is.na(res[[comp]]$pvalue))], statistic= "normal", plot = F)
  for(item in c("pval","qval","lfdr")){
    stat<-FDR.res[[comp]][[item]]
    FDR.res[[comp]][[item]]<-rep(NA,length(genes))
    FDR.res[[comp]][[item]][which(!is.na(res[[comp]]$pvalue))]<-stat
  }
  DE[[comp]]$localFDR<-FDR.res[[comp]]$qval
}

#sample/group correspondence 
sampleGroupVector<-colData(dds)[[condCol]]
names(sampleGroupVector)<-colnames(dds)

DE.sel<-list()
for(comp in comps){
  DE.sel[[comp]]<-list()
  DE.sel[[comp]]$up<-DE[[comp]][which(DE[[comp]]$padj<AdjPVaLthreshold & DE[[comp]]$log2FoldChange > logFCthreshold &  DE[[comp]]$localFDR<LocalFDRthreshold),]
  DE.sel[[comp]]$down<-DE[[comp]][which(DE[[comp]]$padj<AdjPVaLthreshold & DE[[comp]]$log2FoldChange < -logFCthreshold &  DE[[comp]]$localFDR<LocalFDRthreshold),]
  DE.sel[[comp]]$isDE<-rbind(DE.sel[[comp]]$up,DE.sel[[comp]]$down)
  DE.sel[[comp]]$notDE<-DE[[comp]][setdiff(rn(DE[[comp]]),rn(DE.sel[[comp]]$isDE)),]
  DE[[comp]]$DE<-"NONE"
  DE[[comp]][rn(DE.sel[[comp]]$up),"DE"]<-"UP"
  DE[[comp]][rn(DE.sel[[comp]]$down),"DE"]<-"DOWN"
  DE[[comp]]$DE<-factor(DE[[comp]]$DE,levels = c("DOWN","NONE","UP"))
}

for(comp in comps) ecrire(DE[[comp]],paste0("results/DESeqRes_",comp,".tsv"))

###P val histogram###
for(comp in comps){
	histList<-list()
	for(ptype in c("pvalue","padj","localFDR")){
		histList[[ptype]]<-ggplot(data.frame(DE[[comp]]),mapping=aes_string(ptype))+
			geom_histogram(binwidth = 0.05)+
			geom_hline(yintercept = median(table(cut(DE[[comp]][,ptype],breaks = seq(0,1,0.05),include.lowest = FALSE))))
	}
	pdf(paste0("figs/pvalHistogram_",comp,".pdf"),width = 8,height = 10)
	multiplot(plotlist = histList,cols = 1)
	dev.off()
}


#####Volcano plots#####

for(comp in comps){
  cond1<-compMatrix[1,comp]
  cond2<-compMatrix[2,comp]
  
  pdf(paste0("figs/VolcanoPlot_",comp,".pdf"),width=10,height=10)
  plot(DE[[comp]]$log2FoldChange,-log10(DE[[comp]]$padj),type="n",ylab="-log10(DEseq padj)",xlab="log2(Fold-Change)",col="black",
       main=paste0("Volcano plot of comparison ",cond1," (pos FC) VS ",cond2," (neg FC)\nBenjamini & Hochberg method (",nrow(DE.sel[[comp]]$isDE),"/",nrow(DE[[comp]])," DE genes)"))
  abline(h=-log10(AdjPVaLthreshold))
  abline(v = -logFCthreshold)
  abline(v = logFCthreshold)
  if(nrow(DE.sel[[comp]]$up)>0) text(DE.sel[[comp]]$up$log2FoldChange,-log10(DE.sel[[comp]]$up$padj),labels = rownames(DE.sel[[comp]]$up),cex=.2,col="red")
  if(nrow(DE.sel[[comp]]$down)>0) text(DE.sel[[comp]]$down$log2FoldChange,-log10(DE.sel[[comp]]$down$padj),labels = rownames(DE.sel[[comp]]$down),cex=.2,col="green")
  points(DE.sel[[comp]]$notDE$log2FoldChange,-log10(DE.sel[[comp]]$notDE$padj),cex=.2,col="black",pch=16)
  
  plot(DE[[comp]]$log2FoldChange,-log10(DE[[comp]]$localFDR),type="n",ylab="-log10(qval)",xlab="log2(Fold-Change)",col="black",
       main=paste0("Volcano plot of comparison ",cond1," (pos FC) VS ",cond2," (neg FC)\nLocal adjusted fdr method (",nrow(DE.sel[[comp]]$isDE),"/",nrow(DE[[comp]])," DE genes)"))
  abline(h=-log10(LocalFDRthreshold))
  abline(v = -logFCthreshold)
  abline(v = logFCthreshold)
  if(nrow(DE.sel[[comp]]$up)>0) text(DE.sel[[comp]]$up$log2FoldChange,-log10(DE.sel[[comp]]$up$localFDR),labels = rownames(DE.sel[[comp]]$up),cex=.2,col="red")
  if(nrow(DE.sel[[comp]]$down)>0) text(DE.sel[[comp]]$down$log2FoldChange,-log10(DE.sel[[comp]]$down$localFDR),labels = rownames(DE.sel[[comp]]$down),cex=.2,col="green")
  points(DE.sel[[comp]]$notDE$log2FoldChange,-log10(DE.sel[[comp]]$notDE$localFDR),cex=.2,col="black",pch=16)
  
  dev.off()
}

####MA plot####
graphList<-list()
for(comp in comps){
  data<-data.frame(DE[[comp]][res[[comp]]$meanInComp>0.5,])
  data$NAME<-rn(data)
  data<-data[order(abs(data$pval)),]
  dataDE<-data[1:GenesInFig,]
  graphList[[comp]]<-ggplot(data=data,mapping = aes(x=meanInComp,y=log2FoldChange,colour=DE))+
    geom_point()+
    scale_x_log10()+
    scale_color_manual(values = c("green","grey75","red"),guide=FALSE)+
    geom_text_repel(data = dataDE,mapping = aes(x=meanInComp,fontface="bold.italic",
                                                y=log2FoldChange,label=NAME),inherit.aes = FALSE,max.iter = 8000)+
    theme(
      panel.background = element_rect(fill = NA,colour="black"),
      panel.grid.major = element_line(colour = NA),
      panel.grid.minor = element_line(colour = NA)
    )+
    ggtitle(paste0("MA-plot for comparison ",compMatrix[1,comp]," vs ",compMatrix[2,comp]))+
    xlab("Mean expression")+ylab("log2(Fold-Change)")
}

pdf("figs/MA-plot.pdf",width = 10,height = 8)
for(comp in comps) print(graphList[[comp]])
dev.off()

#### DE genes clustering#####
#Select comparisons with number of DE > 0
compsDE<-c(); for(comp in comps) if(nrow(DE.sel[[comp]]$isDE)>1) compsDE<-c(compsDE,comp)
#DE genes storage in dataframes/vectors
DEgenes.names<-list()

for (comp in compsDE) DEgenes.names[[comp]]<-rn(DE.sel[[comp]]$isDE)

allDEgenes.names<-unique(unlist(DEgenes.names))

exprDE<-list()
exprDE.scaled<-list()
haByComp<-list()
for (comp in compsDE){
  if(!allSampleInHt){
    sampleHt<-sampleByComp[[comp]]
    haByComp[[comp]]
    if(batch){ 
      tempAnnot<-droplevels(sampleAnnot[sampleHt,c(condCol,batchCol)]) }else{
      tempAnnot<-droplevels(sampleAnnot[sampleHt,condCol,drop=FALSE])}
    
    colTopAnnotTemp<-vector("list", ncol(tempAnnot))
    names(colTopAnnotTemp)<-cn(tempAnnot)
    i<-1
    for(col in cn(tempAnnot)){
      colTopAnnotTemp[[col]]<-colFun[[i]](nlevels(tempAnnot[,col]))
      names(colTopAnnotTemp[[col]])<-levels(tempAnnot[,col])
      i<-i+1
      if(i==4) i<-1
    }
    haByComp[[comp]]<-HeatmapAnnotation(tempAnnot,col = colTopAnnotTemp)
  }else{
    haByComp[[comp]]<-ha
  }
  exprDE[[comp]]<-exprW[DEgenes.names[[comp]],sampleHt,drop=FALSE]
  exprDE.scaled[[comp]]<-rowScale(exprDE[[comp]],center = TRUE,scaled = TRUE)
}

#! slow if bootsrap
hclustGeneDE<-list()
hclustSampleDE<-list()
for (comp in compsDE){
  bootTemp<-bootstrap; if(nrow(exprDE.scaled[[comp]])>10) bootTemp<-FALSE
  hclustGeneDE[[comp]]<-unsupervisedClustering(exprDE.scaled[[comp]],transpose = F,nboot=nboot,bootstrap = bootTemp)
  hclustSampleDE[[comp]]<-unsupervisedClustering(exprDE.scaled[[comp]],transpose = T,nboot=nboot,bootstrap = bootTemp,method.dist = "euclidean")
}

pdf("figs/clustDEgene.pdf",width = 10,height=10)
for(comp in compsDE){
  print(
    Heatmap(exprDE.scaled[[comp]],top_annotation = haByComp[[comp]], row_names_gp = autoGparFontSizeMatrix(nrow(exprDE.scaled[[comp]])),
      cluster_rows = hclustGeneDE[[comp]],col = colHA,
      cluster_columns = hclustSampleDE[[comp]],name=comp,column_names_gp = autoGparFontSizeMatrix(ncol(exprDE.scaled[[comp]])) )
  )
}
dev.off();

####DE secondary analysis####
#DE genes digital matrix, 1=gene is DE, 0=gene not DE
if(len(conds)>2){
  DE.df<-matrix(0,nrow(exprDat),length(comps))
  colnames(DE.df)<-comps
  rownames(DE.df)<-genes
  DE.df<-data.frame(DE.df)
  for(comp in compsDE) DE.df[DEgenes.names[[comp]],comp]<-1
  ecrire(DE.df,"results/DEgenesByComparison.tsv")
  
  #Gene specificity
  specificGenesCond<-matrix(0,nrow = len(genes),ncol=nlevels(sampleAnnot[,condCol]),dimnames = list(genes,conds))
  addFreq<-1/(len(conds)-1)
  for(i in 1:len(conds)){
      cond<-conds[i]
      InComp<-which(compMatrix[1,]==cond | compMatrix[2,]==cond)
      notInComp<-which(!(compMatrix[1,]==cond | compMatrix[1,]==cond))
      specificGenesCond[,i]<-rowSums(DE.df[,InComp,drop=FALSE])-rowSums(DE.df[,notInComp,drop=FALSE])
  }
  ecrire(specificGenesCond,"results/geneSpecificityByConditions.tsv")
  
  abstract.DE<-matrix(0,nrow = len(conds),ncol = len(conds),dimnames = list(conds,conds))
  for(comp in comps){
    cond1<-compMatrix[1,comp]
    cond2<-compMatrix[2,comp]
    abstract.DE[cond1,cond2]<-abstract.DE[cond2,cond1]<-nrow(DE.sel[[comp]]$isDE)
  }
  
  pdf("figs/summaryDE.pdf")
  print(Heatmap(abstract.DE, rect_gp = gpar(col= "white"),
      cell_fun =function(j, i, x, y, w, h, col) {grid.text(as.character(abstract.DE[i,j]),x,y) },
      name="Number of\nDE gene",column_title = "Summary of differential\nexpression genes analysis",
      row_dend_reorder = F,column_dend_reorder = F))
  
  dev.off()
}
save.image("rsave/Step3.R.RData")

####Functionnal Enrichment####

logFC<-list()
scoreEnrich<-list()
enrichRes<-list()
compWtPathway<-c()
### Several GO term ###
enrichDBs<-getDBterms(rn(exprW), speciesData=species.data,database=c("kegg","reactom","goBP","goCC","goMF"))

for(comp in compsDE){
  logFC[[comp]] = res[[comp]]$log2FoldChange;names(logFC[[comp]])<-rn(res[[comp]])
  scoreEnrich[[comp]]<-calConsensusRanking(rn(res[[comp]]),res[[comp]]$pvalue,logFC[[comp]])
  enrichRes[[comp]]<-Enrich(scoreEnrich[[comp]],corrIdGenes = species.data$GeneIdTable,
    returnLeadingEdge = TRUE,db_terms = enrichDBs)
  exportEnrich(enrichRes[[comp]],paste0("results/Enrich_",comp,".tsv"))
  if(length(which(enrichRes[[comp]]$padj<Enrichthreshold))>0) compWtPathway<-c(compWtPathway,comp)
  enrichRes[[comp]]<-enrichRes[[comp]][enrichRes[[comp]]$padj<Enrichthreshold,]
}


dir.create("figs/enrich",showWarnings = FALSE)
for(comp in compWtPathway){
  pdf(paste0("figs/enrich/EnrichPathway",comp,".pdf"),width = 10,height = 15,useDingbats = FALSE)
  for(i in 1:min(PathwayInFig,nrow(enrichRes[[comp]]))){
    selGene<-intersect(rn(exprW),enrichRes[[comp]]$leadingEdge[i][[1]])
    if(! length(selGene)>2) next
    data<-na.omit(data.frame(rowScale(exprW[selGene,sampleByComp[[comp]]],center = TRUE,scaled = TRUE)))
    pathway<-enrichRes[[comp]][i,enrichRes[[comp]][i,"pathway"]]
    db<-enrichRes[[comp]][i,enrichRes[[comp]][i,"database"]]
    print(Heatmap(data,top_annotation = haByComp[[comp]],clustering_distance_rows = "pearson",
      row_title = paste0(pathway," (",db,")"),
      row_names_gp = autoGparFontSizeMatrix(nrow(data)),
      column_names_gp = autoGparFontSizeMatrix(ncol(data)),
      clustering_distance_columns = "euclidean",clustering_method_columns = "ward.D2",clustering_method_rows = "ward.D2",
      name="z-scored\nexpression",col = colHA,rect_gp = gpar(col= "black")))
  }
  dev.off()
}




for(comp in compWtPathway){
  dir.create(paste0("figs/enrich/pathview_",comp),showWarnings = FALSE)
  setwd(paste0("figs/enrich/pathview_",comp))
  keggPathway<-unlist(enrichRes[[comp]][enrichRes[[comp]]$database=="kegg","pathway"])
  for(pathway in keggPathway){
    viewKEGG(logFC[[comp]],pathway,corrIdGenes = species.data$GeneIdTable)
  }
  setwd("../../..")
}


save.image("rsave/Step4.R.RData")

