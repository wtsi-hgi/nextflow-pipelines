#!/usr/bin/env bash
echo starting nextflow

export NXF_OPTS="-Xms8G -Xmx8G -Dnxf.pool.maxThreads=2000"
nextflow run ./pipelines/crispr.nf -profile farm4_singularity_gn5 -resume 
