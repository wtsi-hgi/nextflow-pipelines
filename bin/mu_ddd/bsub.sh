#!/usr/bin/env bash

if [ -z "${LSF_GROUP}" ]; then
    echo "you must define env variable LSF_GROUP to submit bsub"
    echo "e.g.: export LSF_GROUP=your_group"
    exit 1
    #GROUP='hgi'
else 
    GROUP=${LSF_GROUP}
fi
## on farm5
export http_proxy=http://wwwcache.sanger.ac.uk:3128
export https_proxy=http://wwwcache.sanger.ac.uk:3128
export HTTP_PROXY=http://wwwcache.sanger.ac.uk:3128
export HTTPS_PROXY=http://wwwcache.sanger.ac.uk:3128
export PATH=/software/singularity-v3.5.1/bin/:$PATH
eval "$(conda shell.bash hook)"
export CONDA_ENVS_DIRS=/lustre/scratch118/humgen/resources/conda_envs
conda activate $CONDA_ENVS_DIRS/nextflow

echo starting bsub
rm -f hs_err_pid* && rm -f timeline* && rm -f trace* && rm -rf report* && rm -f bsub.o && rm -f bsub.e && rm -f .nextflow.log && bsub -G $GROUP -R'select[mem>6000] rusage[mem=6000] span[hosts=1]' -M 6000 -n 2 -o bsub.o -e bsub.e -q normal ./nextflow-pipelines/bin/mu_ddd/nextflow.sh > bjob.id
echo finished bsub
cat bjob.id
