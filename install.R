if (!requireNamespace("BiocManager", quietly = TRUE)){
	install.packages("BiocManager")
	BiocManager::install()	
}

library("BiocManager")

install("BiocParallel", update=FALSE)
install("foreach", update=FALSE)
install("data.table", update=FALSE)
install("parallel", update=FALSE)
install("DESeq2", update=FALSE)
install("ggplot2", update=FALSE)
install("ggrepel", update=FALSE)
install("ggbeeswarm", update=FALSE)
install("rjson", update=FALSE)
install("RJSONIO", update=FALSE)
install("stringr", update=FALSE)
install("doParallel", update=FALSE)
install("sva", update=FALSE)
install("WGCNA", update=FALSE)
install("ComplexHeatmap", update=FALSE)
install("gage", update=FALSE)
install("reticulate", update=FALSE)
install("devtools", update=FALSE)

#oob install
install.packages("devtools")
devtools::install_github("https://github.com/DimitriMeistermann/oob")
reticulate::py_install(c("igraph","numpy","leidenalg","trimap"))