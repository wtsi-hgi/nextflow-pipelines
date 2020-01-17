#!/usr/bin/env bash
echo starting nextflow

export NXF_OPTS="-Xms2G -Xmx2G -Dnxf.pool.maxThreads=2000"
/software/mistral_2.13.4_RC5_x86_64/mistral_launch.sh nextflow run ./nextflow-pipelines/pipelines/rna_seq.nf -c ./nextflow-pipelines/nextflow.config -profile farm4_singularity_gn5 -resume 
