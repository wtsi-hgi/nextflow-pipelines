#!/usr/bin/env bash
echo starting nextflow

export NXF_OPTS='-Xms4G -Xmx4G -Dnxf.pool.maxThreads=2000'
#nextflow run ./nextflow-pipelines/pipelines/scrna.nf -c ./nextflow_wr/nextflow_wr_openstack.config -resume 
nextflow run ./nextflow-pipelines/pipelines/scrna_genotype.nf -c ./nextflow_wr/nextflow_wr_openstack.config -resume 
