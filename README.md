# GI_mapping

* include background on file name formatting
* define the goals of the package, include a figure
* write out all requirements for formatting files, etc. (or make a readthedocs?)

## To Do

### Phoebe
* add user input re: filtering counts, reps w/ low corr, reps w/ high SSMD, running MAGeCK or calculating LFC, GI score method to be used
* add ifelse statement to confirm if RNAseq data is available - if it's not available, don't run those adjustments or include in the + control filter 
* add optional rule to use MAGeCK to calculate LFCs - run based on config contents
* re-write updated conda env yaml files
* at what stage can I take the mean across reps? Alice says keep the reps in, but add the mean as another "rep" - decide what is best to do for this
* stick with average lm method for now, add in binned median as well
* name function arguments
* make sure I get the same output I got from my original analysis
* think about Jesse's comment from my committee meeting re: calculating SL GI scores for low-scoring pairs
* consider changing target_type from gene_ctrl and ctrl_gene to gene1_only and gene2_only (or something similar) in case people have a different library design
* convert id => pgRNA_id in all relevant files (annotation for last step of pgRNA counting pipeline?)
* make Rproj? If I can figure out how to do this in a useful way with the command line - or should I use renv instead of conda?
  * https://rstudio.github.io/renv/articles/collaborating.html
* remove really big DFs from memory once they are no longer needed
* save RData somehow (but make sure I'm not saving big unnecessary datasets)
  * save image from Snakemake runs so I can make sure I'm editing using the same versions of everything?
* write fxn to convert Python/bash-formatted input to R file path input?
* update counts QC Rmd to get library size dynamically
* make it so that positive controls are also expressed?
* figure out when/how to make results files (for Snakemake and Rmd)
* compare my calculated LFC values to those from MAGeCK
* figure out when to filter out low read count pgRNAs - in the counts_QC or calculate_LFC scripts?
* within pgRNA_counts_QC.Rmd
  * check CPM calculation to determine if it's correct, share w/ Daniel
* figure out if "normalized" count from MAGeCK = CPM - if so, can just do LFC calculations in my own script
* make annotations in a separate script - do for all library files, and take user input re: having RNAseq data or getting from DepMap
* set the default parameter so that there is an appropriate error message if no parameter is supplied?
* consider using `rmarkdown::render("MyDocument.Rmd", params = "ask")` (see: https://bookdown.org/yihui/rmarkdown/params-knit.html) to get user input rather than relying solely on the config file
* change sample/rep labels from days to "plasmid", "early", "late"? Or define what those are in the config file and then import that info as parameters to R scripts?
* address cases w/ multiple early TP reps
* change counter_efficient.R output to make variable names just be the sample name (not counts_sampleName)
* should I filter out count = 0 before calculating CPM outliers??
* try out formatR formatting and options: https://bookdown.org/yihui/rmarkdown-cookbook/opts-tidy.html
* if I end up copying MAGeCK source code, include this info: https://github.com/davidliwei/mageck/blob/master/COPYING


#### Completed:
* write a function to print kable output ... include some kind of nrow() cutoff?
* calculate coverage in pgRNA_counts_QC Rmd
* write a function to save plots and tbls
* save functions as separate scripts & call from each Rmd using source()
* figure out how to knit Rmd with parameters => do this within the Snakemake pipeline
* update id => pgRNA_id
* pre-process annotations file to get a saveable/shareable one => update get_pgRNA_annotations.Rmd
* get rid of "d." in saved variable names...except for RDS files?
* convert rule 1 output from .html file to .txt file
