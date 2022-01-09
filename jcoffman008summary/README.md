#jcoffman008summary

This directory was designed by Laura Drepanos in Summer 2021 to analyze jcoffman008 using PCA and DESEQ output files. jcoffman008 was an experiment conducted by Ian Gans in the Coffman lab that collected RNAseq data on 4 biological replicates of wildtype and klf9 -/- fish at 3 hours before, at, and 2.5 hours after the onset of light. This analysis is specific to this experiment as it relies upon results from differential expression analysis on data collected at 3 separate times. 

## Inserting data
The PCA and DESEQ output files required to run the scripts in this directory are not included in the github repository to protect the Coffman Lab's experimental data. If you have this data and would like to use it to run the scripts in this directory, create a "data" directory and populate it with...
1. DESeq_out.tsv from DESeq run of KO vs. WT at 3 hours before lights on
2. DESeq_out.tsv from DESeq run of KO vs. WT at lights on
3. DESeq_out.tsv from DESeq run of KO vs. WT at 2.5 hours after lights on
4. DESeq_out.tsv from DESeq run of KO vs. WT at all times with time as a covariant. 
5. rlog.tsv from DESeq run of KO vs. WT at all times with time as a covariant.
6. Z_normalized.txt output file from rpca (of all jcoffman008 samples)
7. pca_eigenvalues.txt output file from rpca 
8. meaningful_pc_loading_scores.txt file from rpca 
9. design_meaningful.txt file from rpca
10. (if running GO summary) downloaded GO terms file from http://aug2020.archive.ensembl.org/biomart/martview/9fc906878735da071619afa6abdf59f6
	- data must be of same ensembl release (101) as used for alignment of RNAseq data 
	- select to download accession term, term name, and Gene stable ID 


## Running Scripts 

### Gene Summary 
PURPOSE: Summarize the transcriptional behavior of a user-defined gene in jcoffman008, with particular focus upon the response to the klf9-knockout, time of day, and replicate number. 

USAGE: Rscript runjcoffman008GeneSummary.R "genename"
	- modify input_jcoffman008GeneSummary.json to modify file paths, normalization method for expression plot ("rlog" vs. "zscore"), number of principal components to report upon/ summarize gene behavior within.  

OUTPUT: This script generates an html report "jcoffman008GeneSummary.html". The first tab displays all plots featured in the report, whereas each 		individual tab can be selected for a larger view of these plots and more explanation of the expression behavior of the gene. 

	Additionally,this script generates facetplots that summarize the distribution of samples along the top 4(default) principal components to show experiment-level patterns of variation that contextualize the the individual gene loadings along these components. The first time the script is run these plots are exported as a png that is used for future runs to avoid repeating this step that is identical regardless of the gene being summarized. This figure is highlighted in the second tab of the report.

	The third tab of the report provides the position of the user-defined gene within the distribution of loadings for the top 4(default) principal components to suggest the effect of the replicate number, time of day, and presence of klf9 in the genome on expression of this gene. 

	The fourth tab of the report provides normalized transcript counts of the user-defined gene in each jcoffman008 sample, demonstrating the extent of variation in expression between samples and allowing the user to identify trends separate from summary-level statistics. 

	The fifth tab specifically focuses on the question of if the klf9-knockout had an effect on the expression of the user-defined gene, and if this effect was consistent across the three times relative to the onset of light. The first plot demonstrates the log-fold change in expression of the gene at the three separate times in response to the klf9-knockout in the context of all genes to underscore the significance of this differential expression. The second plot only shows the log fold change of the user-defined gene to elucidate nuanced changes in the effect of the knockout at different times. This plot reports the p-value and transcript count (base Mean) to display the confidence of the log fold change estimate. Lastly. the MA plot in this tab focuses on the time-independent effect of the klf9 knockout upon the user defined gene by plotting its log2foldchange in the analysis of all klf9-knockout vs. wildtype samples (with time as a covariant) against its base Mean transcript count. 


### Sorted Lists 
PURPOSE: getSortedLists.R automates the process of filtering DESEQ and PCA output otherwise completed in Excel, yielding sorted lists that can easily be entered into GOrilla or another downstream analysis tool to understand processes effected by the klf9-knockout in jcoffman008. Additionally, having this script tracks the cutoffs and modifications from the DESEQ and rpca output files in downstream analysis to make this process reproducible and give more background information to the results of this analysis. 

Note that genes from DESEQ output are not filtered out on the basis of p-value to avoid false negatives (genes that do respond to the presence of klf9 but show a high p-value in DESEQ output). Instead, genes are filtered on the basis of low count (since these genes can easily have false significance and have a minimal effect if only expressed minimally) and logFoldChange standard error (since high variability creates noise that could yield potential additional false significance). No filtering is performed on the lists of genes sorted by PC3, for rpca output already contains a filtered list that eliminates low-expression genes. 

USAGE: Rscript getSortedLists.R 
	- modify input_getSortedLists.json to edit logFoldChange SE and baseMean cutoffs and filepaths 

OUTPUT: This script creates and populates the sorted_lists subdirectory with csvs that are filtered and sorted in both descending and ascending order to identify upregulated and downregulated genes from each output. Genes sorted by PC3 loadings as well as logfoldchange in the analysis of all KO vs. WT samples with time as a covariant can be entered into downstream analysis tools to highlight the biological effects of the klf9-knockout that are consistent across time. Alternatively, the sorted lists from time-specific runs would only yield the biological processes effected by the klf9-knockout at that time. A subdirectory within sorted_lists titled histograms provides the distribution and cutoffs of log2FoldChange Standard Error and base Mean to defend the filtering applied on the basis of expression level and variation. 

### GO Summary 
PURPOSE: Summarize the transcriptional behavior of genes in a user-defined Gene Ontology group in jcoffman008, highlighting the response to the klf9 knockout at the three separate times

USAGE: (sample GO): Rscript jcoffman008GOSummary.R "GO:0007601" 
      - add tag "zscore" after GO term if you wish to see the Z score of transcript counts of genes in this GO rather than log fold change 
      		(technically any tag after the GO term that isn't "log2fc" will result in the script showing z score )
	- edit input_jcoffman008GOsummary.json to change filepaths
OUTPUT: creates the GO_summaries directory and populates it with an additional png each time it is run with a unique GO term. If "znorm", is in the file name, the png will demonstrate the z-normalized transcript count of all genes in all jcoffman008 samples, separated by time and further by klf9 knockouts vs. wildtype samples. Otherwise, the png will be a figure of 3 violin plots demonstrating the log fold change of all genes in the GO in the three time-separated DESEQ analyses. The boxplot juxtaposed over the violin plots provide the 5-number summary of the same data. 

### REVIGO TreeMap Visualizations 
PURPOSE: Uses Revigo clustering algorithm to summarize and visualize set of enriched Gene Ontology groups. File mostly written by Revigo webtool- I changed the settings in the treemap function call to enhance the readabiltiy of the visualization.

USAGE: Run Revigo webtool (http://revigo.irb.hr/) (provided list of GO groups and FDRs from results of GSEA or GOrilla), which yields a file nearly identical to RevigoTreeMap.R this with revigo.data object populated specific to the data-- replace revigo.data in RevigoTreeMap.R with those lines then run Rscript RevigoTreeMap.R

OUTPUT: Generates revigo_treemap.pdf with TreeMap visualization. 

