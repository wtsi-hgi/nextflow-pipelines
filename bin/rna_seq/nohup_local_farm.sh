#!/usr/bin/env bash

# farm5 config
export http_proxy=http://wwwcache.sanger.ac.uk:3128
export https_proxy=http://wwwcache.sanger.ac.uk:3128
export HTTP_PROXY=http://wwwcache.sanger.ac.uk:3128
export HTTPS_PROXY=http://wwwcache.sanger.ac.uk:3128
export PATH=/software/singularity-v3.5.1/bin/:$PATH
eval "$(conda shell.bash hook)"
conda activate nextflow

rm -f hs_err_pid* && rm -f timeline* && rm -f trace* && rm -rf report* && rm -f bsub.o && rm -f bsub.e && rm -f .nextflow.log && rm -f nohup.log && rm -f nohup.out && nohup ./nextflow-pipelines/bin/rna_seq/nextflow.sh > nohup.log 2>&1 &
