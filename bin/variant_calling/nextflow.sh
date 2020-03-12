#!/usr/bin/env bash
echo starting nextflow

#export NXF_OPTS="-Xms7G -Xmx7G -Dnxf.pool.type=sync -Dnxf.pool.maxThreads=4000"
export NXF_OPTS="-Xms5G -Xmx5G"
nextflow run ./nextflow-pipelines/pipelines/as_vqsr.nf -c ./nextflow-pipelines/nextflow.config -profile farm4_singularity_gn5 -resume 
