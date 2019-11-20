#!/usr/bin/env bash
echo starting nextflow

export _JAVA_OPTIONS='-XX:StartFlightRecording=duration=60s,filename=myrecording.jfr'
export NXF_OPTS='-Xms2G -Xmx2G -Dnxf.pool.maxThreads=2000'
nextflow run ./nextflow-pipelines/pipelines/crispr.nf -c ./nextflow-pipelines/nextflow.config -profile farm4_singularity_gn5 -resume 
