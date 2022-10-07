#!/bin/bash

source activate snakemake_env

snakemake --snakefile workflow/Snakefile \
  --use-conda --conda-frontend mamba \
  --reason --cores all
