
Bulk RNA-seq pipeline v3.2
---------------------------------
# Installation and launch:
1. Put all pipeline files in a new directory
2. Put data file (expression table and sample annotation table) in the same format as examples (.tsv with header and row names) in the data directory.
3. Run installDGE.R
4. Change all necessary parameters in section #Config param#
5. Run the script

# Warning
- Sample Names must not contain special character like '(/-', but can contain dot. Condition names must not contain '_' (underscore)
- Do not open expression file with Excel: Excel converts some gene names into date.
- Don't hesitate to remove abnormal sample (see figs/DistribCountPerSample.pdf and results/SamplesAbstract.tsv).
- With fdrtool package, Q-value has a minumim of 1.09048e-14, that can explain the clumping effect on the 2nd page of Volcano plots.
- For each condition you must have a minimum n of 2, otherwise you can run the script up to 'save.image("rsave/Step2.R.RData")'

# Session info
	R version 3.6.0 (2019-04-26)
	Platform: x86_64-w64-mingw32/x64 (64-bit)
	Running under: Windows >= 8 x64 (build 9200)

	Matrix products: default

	locale:
	[1] LC_COLLATE=French_France.1252  LC_CTYPE=French_France.1252    LC_MONETARY=French_France.1252 LC_NUMERIC=C                  
	[5] LC_TIME=French_France.1252    

	attached base packages:
	 [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

	other attached packages:
	 [1] ggrepel_0.8.1               circlize_0.4.6              RColorBrewer_1.1-2          limma_3.40.0               
	 [5] plyr_1.8.4                  genefilter_1.66.0           fdrtool_1.2.15              GO.db_3.8.2                
	 [9] pathview_1.24.0             org.Hs.eg.db_3.8.2          rgl_0.100.19                ComplexHeatmap_2.0.0       
	[13] pvclust_2.0-0               reactome.db_1.68.0          AnnotationDbi_1.46.0        fgsea_1.10.0               
	[17] Rcpp_1.0.1                  ggplot2_3.1.1               gage_2.34.0                 DESeq2_1.24.0              
	[21] SummarizedExperiment_1.14.0 DelayedArray_0.10.0         BiocParallel_1.18.0         matrixStats_0.54.0         
	[25] Biobase_2.44.0              GenomicRanges_1.36.0        GenomeInfoDb_1.20.0         IRanges_2.18.0             
	[29] S4Vectors_0.22.0            BiocGenerics_0.30.0         flashClust_1.01-2           shiny_1.3.2                
	[33] BiocManager_1.30.4         

	loaded via a namespace (and not attached):
	 [1] colorspace_1.4-1        rjson_0.2.20            htmlTable_1.13.1        XVector_0.24.0          GlobalOptions_0.1.0    
	 [6] base64enc_0.1-3         clue_0.3-57             rstudioapi_0.10         bit64_0.9-7             splines_3.6.0          
	[11] geneplotter_1.62.0      knitr_1.22              Formula_1.2-3           jsonlite_1.6            annotate_1.62.0        
	[16] cluster_2.0.9           png_0.1-7               graph_1.62.0            compiler_3.6.0          httr_1.4.0             
	[21] backports_1.1.4         assertthat_0.2.1        Matrix_1.2-17           lazyeval_0.2.2          later_0.8.0            
	[26] acepack_1.4.1           htmltools_0.3.6         tools_3.6.0             gtable_0.3.0            glue_1.3.1             
	[31] GenomeInfoDbData_1.2.1  dplyr_0.8.0.1           fastmatch_1.1-0         Biostrings_2.52.0       crosstalk_1.0.0        
	[36] xfun_0.6                stringr_1.4.0           mime_0.6                miniUI_0.1.1.1          XML_3.98-1.19          
	[41] zlibbioc_1.30.0         scales_1.0.0            promises_1.0.1          KEGGgraph_1.44.0        curl_3.3               
	[46] yaml_2.2.0              memoise_1.1.0           gridExtra_2.3           rpart_4.1-15            latticeExtra_0.6-28    
	[51] stringi_1.4.3           RSQLite_2.1.1           checkmate_1.9.3         manipulateWidget_0.10.0 shape_1.4.4            
	[56] rlang_0.3.4             pkgconfig_2.0.2         bitops_1.0-6            lattice_0.20-38         purrr_0.3.2            
	[61] labeling_0.3            htmlwidgets_1.3         bit_1.1-14              tidyselect_0.2.5        magrittr_1.5           
	[66] R6_2.4.0                Hmisc_4.2-0             DBI_1.0.0               pillar_1.3.1            foreign_0.8-71         
	[71] withr_2.1.2             survival_2.44-1.1       KEGGREST_1.24.0         RCurl_1.95-4.12         nnet_7.3-12            
	[76] tibble_2.1.1            crayon_1.3.4            GetoptLong_0.1.7        locfit_1.5-9.1          data.table_1.12.2      
	[81] blob_1.1.1              Rgraphviz_2.28.0        digest_0.6.18           webshot_0.5.1           xtable_1.8-4           
	[86] httpuv_1.5.1            munsell_0.5.0      

# Credit and thanks
	Pipeline writen by Dimitri Meistermann, University of Nantes, PHD student in computational biology
	at CRTI (UMR 1064) and LS2N (UMR 6241).
	mail: dimitri.meistermann@univ-nantes.fr
	PHD supervised by Jérémie Bourdon (UMR 6241) and Laurent David (UMR 1064).
	
	special thanks to Hayat Hage (University of Nantes)
