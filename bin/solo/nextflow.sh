#!/usr/bin/env bash
echo starting nextflow

export NXF_OPTS='-Xms4G -Xmx4G'
export DATE=$(date +"%Y-%m-%d %T")
export DATE=$(echo $DATE | tr -d '[:space:]' | sed 's/:/-/g') 
echo $DATE
/software/hgi/installs/anaconda3/envs/nextflow20/bin/nextflow run ./nextflow-pipelines/pipelines/solo.nf --current_date \'${DATE}\' -c ./nextflow-pipelines/nextflow.config -profile farm5_lsf_singularity -resume
