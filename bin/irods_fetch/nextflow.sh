#!/usr/bin/env bash
echo starting nextflow

export NXF_OPTS="-Xms8G -Xmx8G -Dnxf.pool.maxThreads=2000"
nextflow run ./nextflow-pipelines/pipelines/irods_fetch.nf -profile farm4_singularity_gn5 --samplefile ./data/single_file_samples_to_read_1_nov_2019.csv --studyid 5679 --runtag 5679 -resume

echo nextflow started locally with nohup
