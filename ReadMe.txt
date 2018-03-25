#############################
# Bulk RNA-seq pipeline v2.3#
#############################

#Installation and launch:
	1- Put all pipeline files in a new directory
	2- Put data file (expression table and sample annotation table) in the same format as examples (.tsv with header and row names) in the data directory.
	3- Run installDGE.R
	4- Change all necessary parameters in section #Config param#
	5- Run the script

#Warning
	-Sample Names must not contain special character like '(/-', but can contain dot. Condition names must not contain '_' (underscore)
	-Do not open expression file with Excel: Excel converts some gene names into date.
	-Don't hesitate to remove abnormal sample (see figs/DistribCountPerSample.pdf and results/SamplesAbstract.tsv).
	-With fdrtool package, Q-value has a minumim of 1.09048e-14, that can explain the clumping effect on the 2nd page of Volcano plots.
	-For each condition you must have a minimum n of 2, otherwise you can run the script up to 'save.image("rsave/Step2.R.RData")'

#Session info
	R version 3.3.1 (2016-06-21)
	Platform: x86_64-w64-mingw32/x64 (64-bit)
	Running under: Windows >= 8 x64 (build 9200)

	locale:
	[1] LC_COLLATE=French_France.1252  LC_CTYPE=French_France.1252    LC_MONETARY=French_France.1252 LC_NUMERIC=C                  
	[5] LC_TIME=French_France.1252    

	attached base packages:
	 [1] grid      stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

	other attached packages:
	 [1] RColorBrewer_1.1-2         limma_3.30.8               circlize_0.3.9             plyr_1.8.4                 genefilter_1.56.0         
	 [6] fdrtool_1.2.15             pathview_1.14.0            org.Hs.eg.db_3.4.0         rgl_0.97.0                 ComplexHeatmap_1.12.0     
	[11] pvclust_2.0-0              topGO_2.26.0               SparseM_1.74               GO.db_3.4.0                AnnotationDbi_1.36.1      
	[16] graph_1.52.0               UsingR_2.0-5               Hmisc_4.0-2                Formula_1.2-1              survival_2.40-1           
	[21] lattice_0.20-34            HistData_0.8-1             MASS_7.3-45                ggplot2_2.2.1              RedeR_1.22.0              
	[26] igraph_1.0.1               gplots_3.0.1               glasso_1.8                 gage_2.24.0                DESeq2_1.14.1             
	[31] SummarizedExperiment_1.4.0 GenomicRanges_1.26.2       GenomeInfoDb_1.10.2        IRanges_2.8.1              S4Vectors_0.12.1          
	[36] flashClust_1.01-2          BiocInstaller_1.24.0       Biobase_2.34.0             BiocGenerics_0.20.0       

	loaded via a namespace (and not attached):
	 [1] colorspace_1.3-2     rjson_0.2.15         class_7.3-14         modeltools_0.2-21    mclust_5.2.2         htmlTable_1.8       
	 [7] XVector_0.14.0       GlobalOptions_0.0.10 base64enc_0.1-3      flexmix_2.3-13       mvtnorm_1.0-5        splines_3.3.1       
	[13] robustbase_0.92-7    geneplotter_1.52.0   knitr_1.15.1         jsonlite_1.2         annotate_1.52.1      cluster_2.0.5       
	[19] kernlab_0.9-25       png_0.1-7            shiny_1.0.0          httr_1.2.1           backports_1.0.5      assertthat_0.1      
	[25] Matrix_1.2-8         lazyeval_0.2.0       acepack_1.4.1        htmltools_0.3.5      tools_3.3.1          gtable_0.2.0        
	[31] Rcpp_0.12.9          trimcluster_0.1-2    Biostrings_2.42.1    gdata_2.17.0         fpc_2.1-10           stringr_1.1.0       
	[37] mime_0.5             gtools_3.5.0         XML_3.98-1.5         dendextend_1.4.0     DEoptimR_1.0-8       zlibbioc_1.20.0     
	[43] scales_0.4.1         KEGGgraph_1.32.0     memoise_1.0.0        gridExtra_2.2.1      rpart_4.1-10         latticeExtra_0.6-28 
	[49] stringi_1.1.2        RSQLite_1.1-2        checkmate_1.8.2      caTools_1.17.1       BiocParallel_1.8.1   shape_1.4.2         
	[55] prabclus_2.2-6       matrixStats_0.51.0   bitops_1.0-6         htmlwidgets_0.8      magrittr_1.5         R6_2.2.0            
	[61] DBI_0.5-1            whisker_0.3-2        foreign_0.8-67       KEGGREST_1.14.0      RCurl_1.95-4.8       nnet_7.3-12         
	[67] tibble_1.2           KernSmooth_2.23-15   viridis_0.3.4        GetoptLong_0.1.5     locfit_1.5-9.1       data.table_1.10.4   
	[73] Rgraphviz_2.18.0     digest_0.6.11        diptest_0.75-7       xtable_1.8-2         httpuv_1.3.3         munsell_0.4.3       

#Credit and thanks
	Pipeline writen by Dimitri Meistermann, University of Nantes, PHD student in computational biology
	at CRTI (UMR 1064) and LS2N (UMR 6241).
	mail: dimitri.meistermann@univ-nantes.fr
	PHD supervised by Jérémie Bourdon (UMR 6241) and Laurent David (UMR 1064).
	
	special thanks to Hayat Hage (University of Nantes)