# GI_mapping

The goal of this package is to map genetic interactions (GIs) among duplicated gene pairs. Users can input a dual gRNA counts table generated via [pgMAP](https://github.com/fredhutch/pgmap_pipeline) or another approach. The GI mapping pipeline performs quality control analyses, quantifies the growth effects of genetic perturbations and drug treatments, and identifies synthetic lethal interactions. Synthetic lethal interactions are quantified as shown below: in the left-hand panel, the sum of growth effects of each single-targeting perturbation is compared to the growth effect of the dual-targeting perturbation for each gene pair. Single- and dual non-targeting controls are used to model a distribution of no interactions, which is compared to the distribution of dual-targeting gRNAs for each gene pair to obtain a multiple testing-adjusted p-value, shown in the right-hand panel. 

![GI mapping approach](https://github.com/pcrparrish/GI_mapping/blob/main/resources/GI_mapping_stats.png?raw=true)

This pipeline is deployed using Snakemake, and all output can be found in the `results/` folder. These results include tab-delimited files containing the growth effects of gene knockouts at the pgRNA and target level, as well as HTML-formatted reports containing relevant figures and reports-based output. 
