# GI_mapping

* include background on file name formatting
* define the goals of the package, include a figure

## To Do

### Phoebe
* write fxn to convert Python/bash-formatted input to R file path input?
* make annotations in a separate script - do for all library files, and take user input re: having RNAseq data or getting from DepMap
* change sample/rep labels from days to "plasmid", "early", "late"? Or define what those are in the config file and then import that info as parameters to R scripts?
* address cases w/ multiple early TP reps
* convert rule 1 output from .html file to .txt file ... unless the output is just a report?
* change counter_efficient.R output to make variable names just be the sample name (not counts_sampleName)
* should I filter out count = 0 before calculating CPM outliers?? 
* within pgRNA_counts_QC.Rmd
  * check CPM calculation to determine if it's correct, share w/ Daniel
* figure out how to knit Rmd with parameters => do this within the Snakemake pipeline
* figure out if "normalized" count from MAGeCK = CPM - if so, can just do LFC calculations in my own script
* try out formatR formatting and options: https://bookdown.org/yihui/rmarkdown-cookbook/opts-tidy.html


#### Completed:
* write a function to print kable output ... include some kind of nrow() cutoff?
* calculate coverage in pgRNA_counts_QC Rmd
* write a function to save plots and tbls
