# GI_mapping

The goal of this package is to map genetic interactions (GIs) among duplicated gene pairs. Users can input a dual gRNA counts table generated via [pgMAP](https://github.com/fredhutch/pgmap_pipeline) or another approach. The GI mapping pipeline performs quality control analyses, quantifies the growth effects of genetic perturbations and drug treatments, and identifies synthetic lethal interactions. Synthetic lethal interactions are quantified as shown below: in the left-hand panel, the sum of growth effects of each single-targeting perturbation is compared to the growth effect of the dual-targeting perturbation for each gene pair. Single- and dual non-targeting controls are used to model a distribution of no interactions, which is compared to the distribution of dual-targeting gRNAs for each gene pair to obtain a multiple testing-adjusted p-value, shown in the right-hand panel. 

This pipeline is deployed using Snakemake, and all output can be found in the `results/` folder. These results include tab-delimited files containing the growth effects of gene knockouts at the pgRNA and target level, as well as HTML-formatted reports containing relevant figures and reports-based output. 


## To Do
* at what stage can I take the mean across reps? Alice says keep the reps in, but add the mean as another "rep" - decide what is best to do for this
* name function arguments
* remove really big DFs from memory once they are no longer needed
* compare my calculated LFC values to those from MAGeCK
* make annotations in a separate script - do for all library files, and take user input re: having RNAseq data or getting from DepMap
* change sample/rep labels from days to "plasmid", "early", "late"? Or define what those are in the config file and then import that info as parameters to R scripts?
* address cases w/ multiple early TP reps
* convert rule 1 output from .html file to .txt file
* change counter_efficient.R output to make variable names just be the sample name (not counts_sampleName)
* should I filter out count = 0 before calculating CPM outliers??
* write out requirements for formatting files, naming experiments etc. 
