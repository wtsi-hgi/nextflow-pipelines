#!/usr/bin/env bash
echo starting nextflow

export NXF_OPTS='-Xms4G -Xmx4G -Dnxf.pool.maxThreads=2000'
nextflow run ./nextflow-pipelines/pipelines/sv-pipeline.nf -c ./nextflow_wr_openstack.config -resume 
