viewKEGG <-
    function(x,
             pathway,
             corrIdGenes = NULL,
             species = "Human",
             speciesData = NULL,
             directory = getwd(),
             ...) {
        require(pathview)
        
        blacklist <- c("hsa04215 Apoptosis - multiple species")
        if (pathway %in% blacklist) {
            warning(
                pathway,
                " is blacklisted as it contains issues in vizualisation, it will not be rendered."
            )
            return(NULL)
        }
        
        if (is.data.frame(x) | is.matrix(x)) {
            tempx <- x
            x <- tempx[, 1]
            names(x) <- rownames(tempx)
        }
        if (is.null(speciesData))
            speciesData <- getSpeciesData2(species)
        if (is.null(corrIdGenes))
            corrIdGenes <- speciesData$GeneIdTable
        entrezId <-
            ConvertKey(
                keyList = names(x),
                tabKey = corrIdGenes,
                colOldKey = "SYMBOL",
                colNewKey = "ENTREZID"
            )
        
        notNA <- which(!is.na(entrezId))
        if (length(notNA) > 0) {
            dat <- x[notNA]
            
            names(dat) <- entrezId[notNA]
            dat <- dat[takefirst(names(dat), returnIndex = T)]
            
            pathway.id <- strsplitNth(pathway, " ")
            
            pathview(
                gene.data = dat,
                pathway.id = pathway.id,
                species = speciesData$kegg,
                kegg.native = TRUE,
                low = "#4B9AD5",
                mid = "white",
                high = "#FAB517",
                na.col = "grey75",
                kegg.dir = directory,
                ...
            )
        } else{
            warning("no entrez id were found")
        }
    }


plotPCsPairs <-
    function(pca,
             sampleAnnot,
             annot2Plot,
             colorScales,
             visualizedComponent = 4,
             ratioPerPercentVar = FALSE,
             plotLabelRepel = TRUE) {
        barplotPercentVar(pca)
        for (i in 1:(visualizedComponent - 1)) {
            for (j in (i + 1):visualizedComponent) {
                for (annot in annot2Plot) {
                    pca2d(
                        pca,
                        colorBy = sampleAnnot[annot],
                        pointSize = 4,
                        comp = c(i, j),
                        main = "Principal Component Analysis",
                        ratioPerPercentVar  = ratioPerPercentVar,
                        colorScale = colorScales[[annot]],
                        plotLabelRepel = plotLabelRepel,
                    )
                }
                pca2d(
                    pca,
                    pointSize = 1.2,
                    comp = c(i, j),
                    plotVars = TRUE,
                    outlierLabel = TRUE,
                    ratioPerPercentVar = ratioPerPercentVar,
                    colorScale = colorScales[[annot]],
                    main = "Correlation circle"
                )
                #change outlierLabelThres (between 0 and 1 to display names of more or less genes)
            }
        }
    }

getPlotPathPerComp<-function(plotNames, comps){
    expand<-expand.grid("resPerComparison",comps,plotNames)
    expand<-expand[order(expand[,2]),]
    apply(expand,1,paste0,collapse="/")
}

