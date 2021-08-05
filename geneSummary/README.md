#geneSummary

This directory was designed by Laura Drepanos in Summer 2021 to make the procedure used for analyzing the PCA output of jcoffman008 applicable to any data that has been run through the MDIBL Bioinformatics Core's PCA script (developed by Bianca Massacci, https://github.com/mdibl/biocore_automated-pca). The goal of this script is to analyze any gene in the experiment in the context of the results of PCA to understand potential responses of the genes to design variables (after establishing association between Principal Components and design variables).

## Inserting data
The PCA output files required to run the scripts in this directory are not included in the github repository to protect the Coffman Lab's experimental data. If you have this data and would like to use it to run the scripts in this directory, create a "data" directory and populate it with the following rpca output files: 
1. ...Z_normalized.txt 
2. ...pca_eigenvalues.txt 
3. ...meaningful_pc_loading_scores.txt  
4. ...design_meaningful.txt



## Running Scripts 

PURPOSE: Use PCA output to summarize the transcriptional behavior of a user-defined gene in the provided experimental data

USAGE: Rscript runGeneSummary.R "genename"
	- modify input_geneSummary.json to enter file path names to access data and design variables to determine faceting in the visualizations 
		note that "treatment variable" requires the header of the column with the variable that represents the treatment (i.e. genetics, cortisol treatment...), and "control" and "treatment" require the different options for cells within that column (i.e. "wildtype" and "klf9 gene knockout"). "separator_variable" should be the secondary variable of interest. for example, someone analyzing a dataset with data from multiple experiments may be curious to see if which experiment a sample belongs to is driving one of the principal components and if the gene of interest has variable expression between these experiments. 

OUTPUT: This script generates an html report "geneSummary.html". The first tab displays all plots generated in the script. The second tab provides facet plots demonstrating the distribution of samples (separated by the design variables indicated in the .json file) to suggest which design variable each principal component may be representing. The third tab places the user-defined gene in the context of these principal components to suggest how the design variables that drive the principal components impact expression of this gene. The final tab more simply shows the normalized transcript counts of this gene in each sample, which may more clearly demonstrate the association with design variables if the principal components did not reveal a clear association with design variables in tab 2. 
