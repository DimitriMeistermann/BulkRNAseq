#!/usr/bin/R

library("Biobase")
library("BiocGenerics")
library("BiocInstaller")
library(shiny)
library(flashClust)
library(DESeq2)
library(gage)
library(ggplot2)
library(topGO)
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

#### Config param  #####
setwd("~/PHDwork/07-DGE-Seq/RUN7/DEgenesAhmed/") #Change Working Directory
RequiredfunctionsDir<-"functions/" 
expressionData<-"data/exprDat_Ahmed.tsv"
sampleTable<-"data/sampleAnnot_Ahmed.tsv"
condCol<-"Group.compact" #name of condition column in sample table for DE genes
batchCol<-NULL #Name of batch column in sample table for batch correction, 'NULL' if no batch
p.valDEthreshold<-0.05
logFCthreshold<-1
sample.species<-"Rat" #data(bods); print(bods) #to see available species
nTopGo<-100 #n TOp Go Term/Ontology
updateSpeciesPackage<-FALSE #update org.species database ?

#####LOADING DATA#####
#Importation des fonctions maison
source(paste0(RequiredfunctionsDir,"/general.R"))
source(paste0(RequiredfunctionsDir,"/projection.R"))
source(paste0(RequiredfunctionsDir,"/RNASeq.R"))


exprDat<-lire(expressionData)
sampleAnnot<-lire(sampleTable)
sampleAnnot[,condCol]<-as.factor(sampleAnnot[,condCol])

#Species preparation
data(bods)
species.data<-data.frame(bods)
species.index<-which(species.data$species==sample.species)
if(len(species.index)!=1) stop("Wrong species name, type \"species.data$species\" to see available species")
species.package<-as.character(species.data[species.index,"package"])
species.KEGG<-as.character(species.data[species.index,"kegg.code"])
if(updateSpeciesPackage) biocLite(species.package)
library(species.package,character.only = TRUE)
corrIdGenes<-AnnotationDbi::select(get(species.package),rn(exprDat), "ENTREZID","SYMBOL")


#####Create dirs#####
dir.create("results",showWarnings=F)
dir.create("figs",showWarnings=F)
dir.create("rsave",showWarnings=F)

#####DESeq#####
batch<-TRUE
if(is.null(batchCol)) batch<-FALSE

exprDat<-filterMostExprimedGenesBySample(exprDat) # Gene filtering

formulaChar<-paste0("~",condCol)
if(batch){
	sampleAnnot[,batchCol]<-as.factor(sampleAnnot[,batchCol])
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

save.image("rsave/Step1.R.RData")

##### Quality control#####
absSamples<-examineRNAseqSamples(exprDat)

pdf(file = "figs/QC.pdf",width=10,height=10)
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
clustSamples<-pvclust(exprW,nboot=30,method.hclust="ward.D2",parallel = TRUE)
corSample<-cor(exprW,method = "pearson")
if(batch){
  clustSamplesUnadjust<-pvclust(exprDatT,nboot=30,method.hclust="ward.D2",parallel = TRUE)
  corSampleUnadjust<-cor(exprDatT,method = "pearson")
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
colFun<-c(rainbow,topo.colors,terrain.colors)
i<-1
for(col in cn(annot)){
  colTopAnnot[[col]]<-colFun[[i]](nlevels(annot[,col]))
  names(colTopAnnot[[col]])<-levels(annot[,col])
  i<-i+1
  if(i==4) i<-1
}
ha<-HeatmapAnnotation(annot,col = colTopAnnot)
Ht<-Heatmap(matrix = corSample, cluster_rows = clustSamples$hclust, cluster_columns = clustSamples$hclust, top_annotation = ha, 
            name="Pearson\ncorrelation",row_names_gp = autoGparFontSizeMatrix(ncol(corSample)),
            column_names_gp = autoGparFontSizeMatrix(ncol(corSample)),
            row_dend_reorder = FALSE, column_dend_reorder = FALSE,col=heat.colors(100))

if(batch){
  Htna<-Heatmap(matrix = corSampleUnadjust, cluster_rows = clustSamplesUnadjust$hclust, cluster_columns = clustSamplesUnadjust$hclust, top_annotation = ha, 
              name="Pearson\ncorrelation\nnon adjusted",row_names_gp = autoGparFontSizeMatrix(ncol(corSampleUnadjust)),
              column_names_gp = autoGparFontSizeMatrix(ncol(corSampleUnadjust)),
              row_dend_reorder = FALSE, column_dend_reorder = FALSE,col=heat.colors(100))
}
Htc<-Heatmap(absSamples[clustSamples$hclust$labels,"TotalCount"], name = "Total\nexpression", 
             col=colorRamp2(c(0,max(absSamples$TotalCount)),c("black","green")), show_row_names = FALSE, width = unit(5, "mm"))

pdf(file = "figs/HeatmapCorPearson.pdf",width=10,height=9)
print(Ht+Htc)
if(batch) print(Htna+Htc)
dev.off()

save.image("rsave/Step2.R.RData")


#####DE analysis#####
#On défini les condition de départ d'interêt
conds<-levels(sampleAnnot[,condCol])

comps<-c()
for (i in 1:length(conds))
  for(j in setdiff(1:length(conds),1:i))
    comps<-c(comps,paste0(conds[i],"_",conds[j]))

res<-list()
for(comp in comps){
  cond1<-unlist(strsplit(comp,"_"))[1]
  cond2<-unlist(strsplit(comp,"_"))[2]
  #résultats DESeq2 pour chaque comparaison
  res[[comp]]<-results(dds, contrast=c(condCol,cond1,cond2), independentFiltering = TRUE)
}

samples<-colnames(exprW)
genes<-rownames(exprW)


#On récupère les statistiques de résultat pour chaque comparaison dans des dataframe
DE<-list()
FDR.res<-list()


for(comp in comps){
  DE[[comp]]<-data.frame(res[[comp]])
  rownames(DE[[comp]])<-genes
  FDR.res[[comp]]<-fdrtool(res[[comp]]$stat[which(!is.na(res[[comp]]$stat))], statistic= "normal", plot = F)
  for(item in c("pval","qval","lfdr")){
    stat<-FDR.res[[comp]][[item]]
    FDR.res[[comp]][[item]]<-rep(NA,length(genes))
    FDR.res[[comp]][[item]][which(!is.na(res[[comp]]$stat))]<-stat
  }
  DE[[comp]]$fdr.qval<-FDR.res[[comp]]$qval
}

#Tableau correspondance échantillon/groupe
sampleGroupVector<-colData(dds)[[condCol]]
names(sampleGroupVector)<-colnames(dds)

for(comp in comps) ecrire(DE[[comp]],paste0("results/DESeqRes_",comp,".tsv"))


#####Volcano plots#####
DE.sel<-list()
for(comp in comps){
  DE.sel[[comp]]<-list()
  DE.sel[[comp]]$up<-DE[[comp]][which(DE[[comp]]$padj<p.valDEthreshold & DE[[comp]]$log2FoldChange > logFCthreshold),]
  DE.sel[[comp]]$down<-DE[[comp]][which(DE[[comp]]$padj<p.valDEthreshold & DE[[comp]]$log2FoldChange < -logFCthreshold),]
  DE.sel[[comp]]$isDE<-rbind(DE.sel[[comp]]$up,DE.sel[[comp]]$down)
  DE.sel[[comp]]$notDE<-DE[[comp]][setdiff(rn(DE[[comp]]),rn(DE.sel[[comp]]$isDE)),]
}


for(comp in comps){
  cond1<-unlist(strsplit(comp,"_"))[1]
  cond2<-unlist(strsplit(comp,"_"))[2]
  
  pdf(paste0("figs/VolcanoPlot_",comp,".pdf"),width=10,height=10)
  plot(DE[[comp]]$log2FoldChange,-log10(DE[[comp]]$padj),type="n",ylab="-log10(DEseq padj)",xlab="log2(Fold-Change)",col="black",
       main=paste0("Volcano plot of comparison ",cond1," (pos FC) VS ",cond2," (neg FC)\nDESeq padj method (",nrow(DE.sel.deseq[[comp]]$isDE),"/",nrow(DE[[comp]])," DE genes)"))
  abline(h=-log10(p.valDEthreshold))
  abline(v = -logFCthreshold)
  abline(v = logFCthreshold)
  if(nrow(DE.sel.deseq[[comp]]$up)>0) text(DE.sel.deseq[[comp]]$up$log2FoldChange,-log10(DE.sel.deseq[[comp]]$up$padj),labels = rownames(DE.sel.deseq[[comp]]$up),cex=.2,col="red")
  if(nrow(DE.sel.deseq[[comp]]$down)>0) text(DE.sel.deseq[[comp]]$down$log2FoldChange,-log10(DE.sel.deseq[[comp]]$down$padj),labels = rownames(DE.sel.deseq[[comp]]$down),cex=.2,col="green")
  points(DE.sel.deseq[[comp]]$notDE$log2FoldChange,-log10(DE.sel.deseq[[comp]]$notDE$padj),cex=.2,col="black",pch=16)
  dev.off()
}


####Clustering DE#####
#Select comparisons with number of DE > 0
compsDE<-c(); for(comp in comps) if(nrow(DE.sel[[comp]]$isDE)>1) compsDE<-c(compsDE,comp)


#Stockage des gènes DE et de leur expression dans des vecteurs/dataframe
DEgenes.names<-list()

for (comp in compsDE) DEgenes.names[[comp]]<-rn(DE.sel[[comp]]$isDE)

allDEgenes.names<-unique(unlist(DEgenes.names))

exprDE<-list()
exprDE.scaled<-list()
for (comp in compsDE){
  exprDE[[comp]]<-exprW[DEgenes.names[[comp]],,drop=FALSE]
  exprDE.scaled[[comp]]<-exprDE[[comp]]-rowMeans(exprDE[[comp]])
}

#! long
hclustGeneDE<-list()
hclustSampleDE<-list()
for (comp in compsDE){
    hclustGeneDE[[comp]]<-pvclust(t(exprDE.scaled[[comp]]),nboot=10,method.hclust = "ward.D2",parallel = TRUE)
    hclustSampleDE[[comp]]<-pvclust(exprDE.scaled[[comp]],nboot=30,method.hclust = "ward.D2",parallel = TRUE)
}

pdf("figs/clustDEgene.pdf",width = 10,height=10)
for(comp in compsDE){
  print(
    Heatmap(exprDE.scaled[[comp]],top_annotation = ha, row_names_gp = autoGparFontSizeMatrix(nrow(exprDE.scaled[[comp]])),cluster_rows = hclustGeneDE[[comp]]$hclust,
                 cluster_columns = hclustSampleDE[[comp]]$hclust,name=comp,column_names_gp = autoGparFontSizeMatrix(ncol(exprDE.scaled[[comp]])) )
  )
}
dev.off();

#Pour chaque comparaison les gène DE sont séparées en 2
hclustDE.cuted<-list()
for(comp in compsDE)
  hclustDE.cuted[[comp]]<-cutree(hclustGeneDE[[comp]]$hclust, k = 2)

#Création de la matrice DE, contient tous les gènes, avec 1 DE positif, 0 DE négatif pour la condition
DE.df<-matrix(0,nrow(exprDat),length(comps))
colnames(DE.df)<-comps
rownames(DE.df)<-genes
DE.df<-data.frame(DE.df)

for(comp in compsDE){
  DE.df[DEgenes.names[[comp]],comp]<-1
}



#Etude logFC, s'assurer d'à quoi correspond un logFC négatif
#On écrit la moyenne du Fold chanage par cluster de condition pour savoir qui est activé pour quelle groupe
meanFC<-matrix(0,length(compsDE)*2,3)
rownames(meanFC)<-paste0(rep(compsDE,each=2),".",rep(1:2,length(compsDE)))
colnames(meanFC)<-c("MeanLogFC","MeanExprCondLeft","MeanExprCondRight")

for(comp in compsDE){
  meanFC[paste0(comp,".1"),1]<- mean(DE[[comp]][names(hclustDE.cuted[[comp]][which(hclustDE.cuted[[comp]]==1)]),"log2FoldChange"])
  meanFC[paste0(comp,".2"),1]<- mean(DE[[comp]][names(hclustDE.cuted[[comp]][which(hclustDE.cuted[[comp]]==2)]),"log2FoldChange"])
  condLeft<-unlist(strsplit(comp,"_"))[1]
  condRight<-unlist(strsplit(comp,"_"))[2]
  meanFC[paste0(comp,".1"),2]<- mean(exprW[names(hclustDE.cuted[[comp]][which(hclustDE.cuted[[comp]]==1)]),rownames(sampleAnnot)[which(sampleAnnot[,condCol]==condLeft)]])
  meanFC[paste0(comp,".2"),2]<- mean(exprW[names(hclustDE.cuted[[comp]][which(hclustDE.cuted[[comp]]==2)]),rownames(sampleAnnot)[which(sampleAnnot[,condCol]==condLeft)]])
  meanFC[paste0(comp,".1"),3]<- mean(exprW[names(hclustDE.cuted[[comp]][which(hclustDE.cuted[[comp]]==1)]),rownames(sampleAnnot)[which(sampleAnnot[,condCol]==condRight)]])
  meanFC[paste0(comp,".2"),3]<- mean(exprW[names(hclustDE.cuted[[comp]][which(hclustDE.cuted[[comp]]==2)]),rownames(sampleAnnot)[which(sampleAnnot[,condCol]==condRight)]])
}


ecrire(meanFC,file = "results/clusterLogFC.tsv")
abstract.DE<-matrix(0,nrow = len(conds),ncol = len(conds),dimnames = list(conds,conds))
for(comp in comps){
  cond1<-unlist(strsplit(comp,"_"))[1]
  cond2<-unlist(strsplit(comp,"_"))[2]
  abstract.DE[cond1,cond2]<-abstract.DE[cond2,cond1]<-nrow(DE.sel[[comp]]$isDE)
}

pdf("figs/summaryDE.pdf")
print(Heatmap(abstract.DE, rect_gp = gpar(col= "white"),
    cell_fun =function(j, i, x, y, w, h, col) {grid.text(as.character(abstract.DE[i,j]),x,y) },
    name="Number of\nDE gene",column_title = "Summary of differentially\nexpressed genes analysis",
    row_dend_reorder = F,column_dend_reorder = F))

dev.off()
save.image("rsave/Step3.R.RData")

####TopGO####
alg<-list() #alg: gene vector list with ensembl id as names
for(comp in compsDE){
  alg[[comp]]<-DE.df[,comp]
  names(alg[[comp]])<-ConvertKey(rownames(DE.df),tabKey = corrIdGenes,colOldKey = "SYMBOL",colNewKey = "ENTREZID")#convertir gene name en ensembl id
  alg[[comp]]<-supprNAnames(alg[[comp]]) #supprimer les genes ne correspondant pas à un ensembl id
}


onts = c( "MF", "BP", "CC" ) #liste des Go ontology
tab = as.list(onts)
names(tab) = onts
topGOResults<-list()
GOdata<-list()
resultTopGO.elim   <-list()
resultTopGO.classic<-list()

for(comp in compsDE){ #for each comparison
  print(comp)
  GOdata[[comp]]<-list()
  resultTopGO.elim[[comp]]<-list()
  resultTopGO.classic[[comp]]<-list()
  pdf(paste0("figs/TopGO_",comp,".pdf"),onefile = TRUE)
  
  for(ont in onts){ #For each ontology (molecular function,bio process, cell componant)
    
    ## prepare data
    GOdata[[comp]][[ont]] <- new("topGOdata", ontology = ont, allGenes = alg[[comp]], geneSel =selectGene,
                                  annot = annFUN.org, mapping=species.package, ID="entrez")
    
    ## run tests
    resultTopGO.elim[[comp]][[ont]]  <- runTest(GOdata[[comp]][[ont]], algorithm = "elim", statistic = "Fisher" )
    resultTopGO.classic[[comp]][[ont]]<-runTest(GOdata[[comp]][[ont]], algorithm = "classic", statistic = "Fisher" )
    
    
    pValue.classic <- score(resultTopGO.classic[[comp]][[ont]])
    pValue.elim <-    score(resultTopGO.elim[[comp]][[ont]])[names(pValue.classic)]
    gstat <- termStat(GOdata[[comp]][[ont]], names(pValue.classic))
    gSize <- gstat$Annotated / max(gstat$Annotated) * 4
    gCol <- colMap(gstat$Significant)
    plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
         pch = 19, cex = gSize, col = gCol, main = paste0("ontology:",ont," comp:",comp))
    
    showSigOfNodes(GOdata[[comp]][[ont]], score(resultTopGO.elim[[comp]][[ont]]), firstSigNodes = 5, useInfo = 'all')
    
    
    ## look at results
    tab[[ont]] <-GenTable( GOdata[[comp]][[ont]], Fisher.elim = resultTopGO.elim[[comp]][[ont]], 
                           Fisher.classic = resultTopGO.classic[[comp]][[ont]],
                           orderBy = "Fisher.classic" , topNodes = nTopGo)
    tab[[ont]]$ont<-rep(ont,nrow(tab[[ont]]))
    print(paste0("ontology ",ont," done"))
  }
  
  topGOResults[[comp]] <- rbind.fill(tab)
  dev.off()
  ecrire(topGOResults[[comp]],paste0("results/","TopGO_",comp,".tsv"))
}

####HEATMAP TOPGO####

#Création du tableau de conversion ensembl id vers gene Name
hashNameToEntrez<-genes
names(hashNameToEntrez)<-ConvertKey(genes,tabKey = corrIdGenes,colOldKey = "SYMBOL",colNewKey = "ENTREZID")
hashNameToEntrez<-supprNAnames(hashNameToEntrez)

### un GO terme ###
#PARAMETRE
compsDE
comp<-compsDE[6]
topGOResults[[comp]][,c("Term","Rank in Fisher.classic","ont")] #visualisation pour choisir le rang d'interêt
index<-1 #index d'interêt

#BEGIN (à executer une fois les paramètres défini jusqu'à END
allRes<-topGOResults[[comp]]
goID <- allRes[index, "GO.ID"] # on recupere ici le GO term avec l'index voulu
term<-  allRes[index, "Term"] #même chose pour la définition du Go
goID.genes.ensemblID <- unlist(genesInTerm(GOdata[[comp]][[allRes[index,"ont"]]], whichGO=goID)) # on recupere l'ensemble des genes associes a ce GO terme 

#preparation de la matrice d'expression
ExpressionDataDE = exprW[alg[[comp]]==1,]
ExpressionData.selecGO=ExpressionDataDE[intersect(hashNameToEntrez[goID.genes.ensemblID],rownames(ExpressionDataDE)),] # matrice d'expression restreinte aux genes associes au terme
data = ExpressionData.selecGO-min(ExpressionData.selecGO)
data = data-rowMeans(data)


#Heatmap
Heatmap(data, column_title = paste0("Disregulated gene from ",goID,":\n",term,"\nin comparison ",comp),
        name = "expression", top_annotation = ha, show_heatmap_legend=T,
        row_names_gp = autoGparFontSizeMatrix(nrow(data)),column_names_gp = autoGparFontSizeMatrix(ncol(data)),
        clustering_distance_rows = "pearson",
        clustering_method_rows = "ward.D2",
        clustering_distance_columns = "pearson",
        clustering_method_columns = "ward.D2"
)
# END

### plusieurs terme GO ###

index=c(1:10,nTopGo:(nTopGo+10),(2*nTopGo):(2*nTopGo+10))
pdf(file = "figs/HeatmapGOterms.pdf",width = 10,height=15)
for(comp in compsDE){
  cond1<-unlist(strsplit(comp,"_"))[1]
  cond2<-unlist(strsplit(comp,"_"))[2]

  allRes<-topGOResults[[comp]]
  goIDs <- allRes[index, "GO.ID"];  names(goIDs)<-index
  goTerms <- allRes[index, "Term"]; names(goTerms)<-index
  
  goID.genes<-list()
  for(i in index){
    genes.selGo<-hashNameToEntrez[unlist(genesInTerm(GOdata[[comp]][[allRes[i,"ont"]]], whichGO=goIDs[as.character(i)]))]
    names(genes)<-c()
    goID.genes[[as.character(i)]]<-genes.selGo
  }
  goID<-unlist(goID.genes); names(goID)<-c()
  
  goAnnotations <-c(); 
  for(i in index){
    goAnnotations<-c(goAnnotations,rep(goTerms[as.character(i)], length(goID.genes[[as.character(i)]])));
  }
  names(goAnnotations) = goID
  
  #prepraration de la matrice d'expression
  ExpressionDataDE = exprW[alg[[comp]]==1,]
  ExpressionData.selecGO=ExpressionDataDE[intersect(goID,rownames(ExpressionDataDE)),] # matrice d'expression restreinte aux genes associes au terme
  data = ExpressionData.selecGO-min(ExpressionData.selecGO)
  data = data-rowMeans(data)
  
  sampleClustGO<-pvclust(data,method.hclust = "ward.D2",method.dist = "euclidean",parallel = TRUE,nboot = 10)
  
  print(Heatmap(data, column_title = paste0("Disregulated gene from several GO terms\nin comparison ",cond1," VS ",cond2),
          name = "Expr",  top_annotation = ha, show_heatmap_legend=TRUE,
          column_names_gp = autoGparFontSizeMatrix(ncol(data)),
          cluster_columns = sampleClustGO$hclust,
          split = goAnnotations[rownames(data)], row_title_gp = gpar(fontsize = 6), 
          row_title_rot =  0, row_names_gp = autoGparFontSizeMatrix(nrow(data)) 
  ))
}
dev.off()

####GAGE####
#Database preparation
data(kegg.gs)
data(go.gs)


go.species <- go.gsets(tolower(sample.species), id.type="entrez")
kg.species <- kegg.gsets(tolower(sample.species), id.type="entrez")
go.gs.bp <- go.species$go.sets[go.species$go.subs$BP]
kegg.gs <- kg.species$kg.sets[kg.species$sigmet.idx]

groupIndex<-list()
for(cond in conds){
  groupIndex[[cond]]<-which(sampleGroup==cond)
}

gene.entrezID = ConvertKey(rownames(exprW),corrIdGenes,"SYMBOL","ENTREZID")
ExpressionData.entrez<-exprW
rownames(ExpressionData.entrez)<-gene.entrezID
ExpressionData.entrez<-ExpressionData.entrez[which(!is.na(rownames(ExpressionData.entrez))),]

dataDE.entrez<-list()
for(comp in comps){
  dataDE.entrez[[comp]]<-DEgenes.names[[comp]]
  dataDE.entrez[[comp]] <- ConvertKey(DEgenes.names[[comp]],corrIdGenes,"SYMBOL","ENTREZID")
  dataDE.entrez[[comp]]<-dataDE.entrez[[comp]][which(!is.na(dataDE.entrez[[comp]]))]
}

#ici compare = "unpaired" car ind non paire entre les 2 conditions (aller voir compare = "unparied" et same.dir = TRUE)
RESgage.kegg<-list()
RESgage.go.bp<-list()
diffExpr<-list()

for(comp in comps){
  ref<- groupIndex[[unlist(strsplit(comp,"_"))[1]]]
  samp<-groupIndex[[unlist(strsplit(comp,"_"))[2]]]
  RESgage.kegg[[comp]] <-gage(exprs=ExpressionData.entrez,gsets = kegg.gs,ref = ref, samp = samp, compare = "unpaired",same.dir = FALSE) #on prend les pathway les plus perturbés
}
# possibilité de faire fonctionner avec sameDirection=FALSE  !!

select<-AnnotationDbi::select #Conflit entre WGCNA et pathview (au cas où le package WGCNA est chargée)

sel<-list()
logFC<-list()
for(comp in comps){
  logFC[[comp]] = res[[comp]]$log2FoldChange #! log2fc sert de valeur par pour les pathview mais toujours dans le même sens
}

for(comp in comps){
  names(logFC[[comp]])<-ConvertKey(genes,corrIdGenes,"SYMBOL","ENTREZID")
  logFC[[comp]]<-logFC[[comp]][which(!is.na(names(logFC[[comp]])))]
  #Greater (up regulated)
  sel[[comp]] <- rownames(RESgage.kegg[[comp]]$greater)[which(RESgage.kegg[[comp]]$greater[,"q.val"] < 0.01)]
  dir.create(paste0("pathviews/pathview_",comp),recursive = TRUE,showWarnings = FALSE)
  setwd(paste0("pathviews/pathview_",comp))
  for(pathway in sel[[comp]]){
    pathview(gene.data = logFC[[comp]], pathway.id = pathway, species = species.KEGG,kegg.native=FALSE)
    pathview(gene.data = logFC[[comp]], pathway.id = pathway, species = species.KEGG,kegg.native=TRUE )
  }
  ecrire(RESgage.kegg[[comp]]$greater[sel[[comp]],],"statPathway.tsv")
  setwd("../..")
}


####Heatmap pathway#### 

### un Kegg terme ###
##PARAMETRE
#comp<-comps[1]
#RESgage.kegg[[comp]]$greater[sel[[comp]],c("stat.mean","q.val")] #visualisation pour choisir le rang d'interêt
#kegg.name<-"hsa04660 T cell receptor signaling pathway" #nom du pathway d'interêt
#
##BEGIN (à executer une fois les paramètres défini jusqu'à END
#allRes<-topGOResults[[comp]]
#kegg.genes<- ConvertKey(unique(unlist(kegg.gs[kegg.name])),corrIdGenes,colOldKey = "ENTREZID",colNewKey = "SYMBOL")
#
##preparation de la matrice d'expression
#ExpressionData.selecKEGG=exprW[intersect(kegg.genes,rownames(exprW)),] # matrice d'expression restreinte aux genes associes au terme
#data = ExpressionData.selecKEGG-min(ExpressionData.selecKEGG)
#data = data-rowMeans(data)
#
#
##Heatmap
#Heatmap(data, column_title = paste0("Disregulated gene from\n",kegg.name,"\nin comparison ",comp),
#        name = "expression", top_annotation = ha, show_heatmap_legend=T,
#        row_names_gp = autoGparFontSizeMatrix(nrow(data)),column_names_gp = autoGparFontSizeMatrix(ncol(data))
#)
## END


####annotation stats####
pathAbstract<-list()
statFC<-list()
geneCompound<-list()
logFC.sym<-list()

for(comp in comps){
  pathAbstract[[comp]]<-data.frame(RESgage.kegg[[comp]]$greater)
  pathAbstract[[comp]]<-pathAbstract[[comp]][which(pathAbstract[[comp]][,"q.val"]<0.01),c("stat.mean","q.val")]
  statFC[[comp]]<-matrix(0,nrow(pathAbstract[[comp]]),7)
  colnames(statFC[[comp]])<-c("n_up","mean_up","n_down","mean_down","n","mean","ratioUpDown")
  pathAbstract[[comp]]<-cbind(pathAbstract[[comp]],statFC[[comp]])
  
  geneCompound[[comp]]<-kegg.gs[rownames(pathAbstract[[comp]])]
  
  logFC.sym[[comp]]<-res[[comp]][,"log2FoldChange"]
  names(logFC.sym[[comp]])<-rownames(res[[comp]])
  
  for(pathway in names(geneCompound[[comp]])){
    geneCompound[[comp]][[pathway]]<-ConvertKey(geneCompound[[comp]][[pathway]],tabKey = corrIdGenes,colOldKey = "ENTREZID",colNewKey = "SYMBOL")
    geneCompound[[comp]][[pathway]]<-geneCompound[[comp]][[pathway]][which(!is.na(geneCompound[[comp]][[pathway]]))]
    geneCompound[[comp]][[pathway]]<-intersect(geneCompound[[comp]][[pathway]],genes)
    path.logFC<-logFC.sym[[comp]][geneCompound[[comp]][[pathway]]]
    
    pathAbstract[[comp]][pathway,"n_up"]<-length(path.logFC[which(path.logFC>0)])
    pathAbstract[[comp]][pathway,"mean_up"]<-mean(path.logFC[which(path.logFC>0)])
    pathAbstract[[comp]][pathway,"n_down"]<-length(path.logFC[which(path.logFC<0)])
    pathAbstract[[comp]][pathway,"mean_down"]<-mean(path.logFC[which(path.logFC<0)])
    pathAbstract[[comp]][pathway,"n"]<-length(path.logFC)
    pathAbstract[[comp]][pathway,"mean"]<-mean(path.logFC,na.rm=TRUE)
  }
  pathAbstract[[comp]][,"ratioUpDown"]<-pathAbstract[[comp]][,"n_up"]/pathAbstract[[comp]][,"n_down"]
  
}

for(comp in comps){
  ecrire(pathAbstract[[comp]],paste0("results/statpath_",comp,".tsv"))
}


save.image("rsave/Step4.R.RData")

