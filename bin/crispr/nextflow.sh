#!/usr/bin/env bash
echo starting nextflow

#export NXF_OPTS="-Xms7G -Xmx7G -Dnxf.pool.type=sync -Dnxf.pool.maxThreads=4000"
export NXF_OPTS="-Xms3G -Xmx3G"
nextflow run ./nextflow-pipelines/pipelines/crispr.nf -c ./nextflow-pipelines/nextflow.config -profile farm5_lsf_singularity -resume 
