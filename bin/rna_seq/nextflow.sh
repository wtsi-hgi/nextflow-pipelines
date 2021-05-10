#!/usr/bin/env bash
echo starting nextflow

#export NXF_OPTS="-Xms3G -Xmx3G -Dnxf.pool.maxThreads=2000"
export NXF_OPTS="-Xms5G -Xmx5G"
#export NXF_VER=20.04.0-edge
nextflow run ./nextflow-pipelines/pipelines/rna_seq.nf --nf_ci_loc $PWD -c ./nextflow-pipelines/nextflow.config -profile farm4_singularity_gn5 -resume 
