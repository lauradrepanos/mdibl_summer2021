#runs jcoffman008genesummary.Rmd using parameters set in input_jcoffman008GeneSummary.json
#usage: Rscript runjcoffman008GeneSummary.R "genename"
#output: jcoffman008genesummary.html

options(stringsAsFactors = FALSE)
library(jsonlite)

args = base::commandArgs(trailingOnly = TRUE)
if(length(args)<1){
	genename="smc1a" #default gene to summarize
}else{
	genename=args[1]
}

json = read_json("input_jcoffman008GeneSummary.json")

rmarkdown::render('jcoffman008GeneSummary.Rmd',
                  params = list(
                      genename=genename,
          
					  experiment_name= json$experiment_name,
					 
					  filepath_Z_normalized= json$filepath$Z_normalized,
					  filepath_loading_scores=json$filepath$loading_scores,
					  filepath_pca_eigenvalues=json$filepath$pca_eigenvalues,
					  filepath_design_meaningful=json$filepath$design_meaningful,

					  filepath_rlog=json$filepath$rlog,
					  filepath_deseq=json$filepath$deseq,
					  filepath_deseq3hb=json$filepath$deseq3hb,
					  filepath_deseqlightson=json$filepath$deseqlightson,
					  filepath_deseq2.5ha=json$filepath$deseq2.5ha,

					  treatment_variable=json$design_variables$treatment_variable,
					  control=json$design_variables$control,
					  treatment= json$design_variables$treatment,
					  separator_variable=json$design_variables$separator_variable,
					  num_meaningful_pcs=json$design_variables$num_meaningful_pcs,
					  normalization=json$design_variables$normalization))

					  



