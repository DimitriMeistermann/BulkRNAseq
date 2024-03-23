if (!requireNamespace("BiocManager", quietly = TRUE)){
	install.packages("BiocManager")
	BiocManager::install()	
}

if (!requireNamespace("devtools", quietly = TRUE))
	install.packages("devtools")

library("BiocManager")

install("BiocParallel", update=FALSE)
install("doParallel", update=FALSE)
install("foreach", update=FALSE)
install("parallel", update=FALSE)
install("stringr", update=FALSE)
install("DESeq2", update=FALSE)
install("sva", update=FALSE)
install("ggrepel", update=FALSE)
install("ggbeeswarm", update=FALSE)
install("rjson", update=FALSE)
install("RJSONIO", update=FALSE)
install("gage", update=FALSE)
install("pathview", update=FALSE)
install("qualpalr",update = FALSE)
install("WGCNA",update = FALSE)

#homemade package install
devtools::install_url("https://github.com/DimitriMeistermann/oob/archive/refs/tags/march2024.zip")
devtools::install_url("https://github.com/DimitriMeistermann/GSDS/archive/refs/tags/march2024.zip")
