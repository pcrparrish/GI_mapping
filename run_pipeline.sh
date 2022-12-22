#!/bin/bash

source activate snakemake_env
echo "Snakemake environment activated"

echo "Running Snakemake pipeline"
snakemake --snakefile workflow/Snakefile \
  --use-conda --conda-frontend mamba \
  --reason --cores all --latency-wait 30 \
  -R pgRNA_counts_QC

echo -e "Exporting pipeline DAG to pdf"
snakemake --snakefile workflow/Snakefile \
  --dag > "workflow/report/dag.dot"

dot -Tpdf "workflow/report/dag.dot" > "workflow/report/pipeline_dag.pdf"

rm "workflow/report/dag.dot"
