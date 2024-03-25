
Bulk RNA-Seq pipeline v5
--------------------------

# Installation and launch:
1. Download the workflow (for example: `git clone https://github.com/DimitriMeistermann/BulkRNAseq.git`)
2. Put data files (expression table and sample annotation table) in the same format as examples (.tsv with header and row names).
	
	**[optional]**: You can provide a "*comparison to do*" file and a "*color scales*" file. See the section *Inputs* in *BulkRNAseq.Rmd* for more details.

3. Run *install.R*
4. Change all necessary parameters in the *Configuration of parameters* section in the script *BulkRNAseq.Rmd*.
5. Run the notebook *BulkRNAseq.Rmd*.

If you want to resume the script after an interruption, run the *Setup* code chunk, and then load the appropriate checkpoint  from the *outRobj* folder (in RStudio, go to file/open file).

# Tips for analyzing bulk RNA-Seq
## A word on experimental design
The experimental design is related to what variables affect gene expression but are known before the start of the experiment. For example, with the variables **genotype [KO, WT]**,  **sex [M, F]** , and **presenceOfMarker [true, false]**, the experimental design is:
$geneExpression \sim genotype + sex$
In this case the condition columns must be *genotype* and *sex*.

Differentially expressed genes are taking account of experimental design: If you want differentially expressed genes between *WT* and *KO*, you will have differentially expressed genes, all things being equal (sex is taken into account).

We assume here that *presenceOfMarker* is the measurement of a protein presence. It is not a part of experimental design as its value for each sample is not known before the experiment is carried out.
*presenceOfMarker* can be used in the script by being put in the `otherInterestingColumn` argument.

For more complex experimental designs, the script has to be slightly modified. [See this website for more details](https://rstudio-pubs-static.s3.amazonaws.com/329027_593046fb6d7a427da6b2c538caf601e1.html).

DESeq2 is theoretically compatible with quantitative variables, however a good way to deal with it is to convert it into factors, with by example the R function *cut*.

## A word on batch correction
Batch effects are technical bias between two or more group of samples. The most common case is when the sequencing was done on two different plates. It is important to think carefully about the experimental design to avoid case where batch effect correction is not precise, or even impossible.

- Tip 1: try to split your group of samples between batches. If one group of sample is only in one batch the batch effect will be harder to retrieve.
- Tip 2: try to put technical replicates between batches to control the quality of batch correction. If you have enough technical replicates you can even retrieve batch effect only from technical replicates, [example here](https://gitlab.univ-nantes.fr/E137833T/Castel_et_al_2020/-/blob/master/Code/1_bulk_rna_seq_preprocessing.R)

For enabling batch correction, you just have to set the argument *batchColumn* to the column name that contains batches attribution in sample annotations. If you want to disable it, set it to `NULL`.

## A word on file formats
Most of text results are in TSV (tabulated separated values). The point is to not worry about the separator between French and English version of CSV. It is also enabling quick paste of data into a spreadsheet. You can open these files with Excel in windows (right click, open with, choose another program, retrieve Excel.exe, be sure to have checked the box to memorize the choice).

Most of figures are in PDF. Hence figure are easy to open on every platform and vectorized (no problem with figure resolution and editable with Adobe Illustrator or Inkscape).

## A word on distribution of p-values
For most large set of p-values, an histogram of p-values is plotted to show their distribution. 

![pvalhistogramms](ReadMe/pvalHist.png)

- ***Case a:***
P-values are following an uniform distribution: samples are presumably from the same population. **Every significant result could be a false positive**. 
- ***Case b:***
A lot of p-values are more closed to 0. Samples are not coming from the same sample. In this case the samples are not coming from the same population: differentially expressed genes are real !
- ***Case c:***
Every other shape of the distribution should raise extreme caution: the model can not handle the data. It could be the case for different reasons: poor quality of data, lack of statistical power...

# Brief summary of the workflow

## Step 1: loading data, building model and generation of counts tables
- Loading homemade function from the repository veneR, NB: you can clone the repository locally if you want to use the script without an internet connection. Loading and verifying data given by the user, setting of the seed and creation of output folders.
- Quality control of genes / samples.
- Generation of color scales and the comparison matrix from the experimental design.
- Species-specific data download (download the *org.sp.eg.db* package if not already installed).
- Compute DESeq2 model
- Normalization and transformation of counts
- Batch correction
- Computing of CPM (Count Per Million)

## Step 2: Unsupervised and multivariate analyses
- Computing of the overdispersion plot
- Principal Component Analysis and Regression (and batch correction control if applicable).
- Correlation heatmap
- 2D TRIMAP from over dispersed genes, sample clustering from an density based clustering on an alternative 10D UMAP.
- gene module attribution by a hierarchical clustering of overdispersed genes from a distance of biweight midcorrelation (bicor).
- Computing modules activation scores and module memberships.
- Computing AUROCS for each gene for each sample cluster.
- Plotting the super heatmap.

## Step 3: Differential gene expression analyses
- For each comparison, computing of results and plot of DE analyses. By default those operations are multithreaded.
- Plotting heatmaps of DE genes.
- Rich Upset plot

## Step 4: Functional enrichment
- Downloading gene set databases
- Over representation analysis (ORA) of DE genes
- Gene Set Differential Scoring (GSDS, experimental) of DE genes
- Gene Set Enrichment Analysis (GSEA) of Principal Components
- ORA enrichment of gene modules.

## Step 5: Building of the web app

# Credits and thanks
Workflow written by Dimitri Meistermann
current mail: dimitri.meistermann@helsinki.fi.

Special thanks to Hayat Hage to have written the very first lines of this script.
