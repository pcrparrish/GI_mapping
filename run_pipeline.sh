#!/bin/bash

source activate snakemake_env
echo "Snakemake environment activated"

snakemake --snakefile workflow/Snakefile \
  --use-conda --conda-frontend mamba \
  --reason --cores all -R pgRNA_counts_QC
