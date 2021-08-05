#runs geneSummary.Rmd using parameters set in input_geneSummary.json
#usage: Rscript runGeneSummary.R "genename"
#about: summarizes the behavior of a gene of interest using the pca output of an experiment
#output: geneSummary.html

options(stringsAsFactors = FALSE)
library(jsonlite)

args = base::commandArgs(trailingOnly = TRUE)
if(length(args)<1){
	genename="smc1a" #default gene to summarize
}else{
	genename=args[1]
}

json = read_json("input_geneSummary.json")

rmarkdown::render('geneSummary.Rmd',
                  params = list(
                      #the gene to summarize in this experiment is specified as an argument to this script to facilitate the process of examining different genes
                      genename=genename,

					  #experiment_name determines text in the report file as well as the name of output files, such as the png of PCA loadings
					  experiment_name=json$experiment_name,

					  #treatment_variable is the name of the column in the design file that specifies the treatment in the experiment
					  treatment_variable=json$design_variables$treatment_variable,
					  #control and treatment are values found within the treatment_variable column in the design file
					  control=json$design_variables$control,
					  treatment= json$design_variables$treatment,

					  #separator_variable is another factor aside from treatment that divides the samples, such as the time they were collected
					  separator_variable=json$design_variables$separator_variable,

					  #the number of PCs you wish to report on. must be <= the number of meaningful PCs as determined by rpca
					  num_meaningful_pcs=json$design_variables$num_meaningful_pcs,

					  #note: the script assumes that the design file has a "replicate" column specifying the replicate number of each sample
					  
					  #the files that result from rpca that are needed to run this script
					  filepath_Z_normalized= json$filepath$Z_normalized,
					  filepath_loading_scores=json$filepath$loading_scores,
					  filepath_pca_eigenvalues=json$filepath$pca_eigenvalues,
					  filepath_design_meaningful=json$filepath$design_meaningful))



