#!/usr/bin/env bash
echo starting nextflow

export NXF_VER=19.12.0-edge
export NXF_OPTS='-Xms2G -Xmx2G -Dnxf.pool.maxThreads=2000'
nextflow run ./nextflow-pipelines/pipelines/scrna.nf -c ./nextflow-pipelines/nextflow.config -profile farm4_singularity_gn5 -resume 
