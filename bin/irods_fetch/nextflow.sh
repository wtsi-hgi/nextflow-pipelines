#!/usr/bin/env bash
echo starting nextflow

export NXF_OPTS="-Xms3G -Xmx3G -Dnxf.pool.maxThreads=2000"
nextflow run ./nextflow-pipelines/pipelines/irods_fetch.nf -c ./nextflow-pipelines/nextflow.config -profile farm4_singularity_gn5 -resume

