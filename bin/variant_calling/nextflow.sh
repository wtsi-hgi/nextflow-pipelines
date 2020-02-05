#!/usr/bin/env bash
echo starting nextflow

export NXF_VER=19.12.0-edge
export NXF_OPTS="-Xms6G -Xmx6G -Dnxf.pool.type=sync -Dnxf.pool.maxThreads=4000"
#export NXF_OPTS="-Xms2G -Xmx2G"
nextflow run ./nextflow-pipelines/pipelines/variant_calling.nf -c ./nextflow-pipelines/nextflow.config -profile farm4_singularity_gn5 -resume 
