
webAppFolder<-"webApp/"

dir.create(webAppFolder,showWarnings = F)
dir.create(paste0(webAppFolder,"data/geneCor"),showWarnings = F,recursive = T)
dir.create(paste0(webAppFolder,"data/gene"),showWarnings = F,recursive = T)


exprForUI<-normCounts

rownames(exprForUI)<-str_replace_all(rownames(exprForUI),"/","_")

print("Web app: writing gene expressions...")

null<-foreach(gene=rn(exprForUI)) %dopar% {
    library(oob)
    newDat<-data.frame(expression=exprForUI[gene,])
    fastWrite(newDat,fileName =paste0(webAppFolder,"data/gene/",gene,".tsv"),row.names = F)
}

write(toJSON(rn(exprForUI)),paste0(webAppFolder,"data/geneList.json"))

print("Web app: writing gene per gene correlations...")

#export top 50 correlation
exprLogForUI<-logCounts[rownames(normCounts),]
rownames(exprLogForUI)<-rownames(exprForUI)

null<-foreach(gene = rn(exprLogForUI)) %dopar%{
    library(oob)
    corr<-t(cor(exprLogForUI[gene,],t(exprLogForUI)))
    res<-corr[order(corr[,1],decreasing = T)[c(2:51,(nrow(corr)-50):nrow(corr))],,drop=F]
    res<-data.frame(Gene=rn(res),pearsonCor=res[,1])
    fastWrite(res,
              file=paste0(webAppFolder,"data/geneCor/",gene,".tsv"), row.names = FALSE)
}

rm(exprForUI)
rm(exprLogForUI)


print("Web app: writing some files...")

colnames(trimap)<-c("x","y")

fastWrite(sampleAnnot[,names(colorScales)],paste0(webAppFolder,"data/sampleAnnot.tsv"))
fastWrite(trimap,paste0(webAppFolder,"data/coorProjection.tsv"))
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