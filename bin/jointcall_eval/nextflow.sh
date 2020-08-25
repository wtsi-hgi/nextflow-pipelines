#!/usr/bin/env bash
echo starting nextflow

#export NXF_VER=19.12.0-edge
#export NXF_OPTS="-Xms3G -Xmx3G -Dnxf.pool.maxThreads=2000"
export NXF_OPTS="-Xms4G -Xmx4G"
nextflow run ./nextflow-pipelines/pipelines/jointcall_eval.nf \
	 -c ./nextflow-pipelines/nextflow.config -profile farm4_singularity_gn5 -resume 
