#!/usr/bin/env bash

LSF_GROUP=${LSF_GROUP:-hgi}
export http_proxy=http://wwwcache.sanger.ac.uk:3128
export https_proxy=http://wwwcache.sanger.ac.uk:3128
export HTTP_PROXY=http://wwwcache.sanger.ac.uk:3128
export HTTPS_PROXY=http://wwwcache.sanger.ac.uk:3128
export PATH=/software/singularity-v3.5.1/bin/:$PATH
eval "$(conda shell.bash hook)"
conda activate nextflow20

echo starting bsub
rm -f hs_err_pid* && rm -f timeline* && rm -f trace* && rm -rf report* && rm -f bsub.o && rm -f bsub.e && rm -f .nextflow.log && bsub -G $LSF_GROUP -R'select[mem>5000] rusage[mem=5000] span[hosts=1]' -M 5000 -n 1 -o bsub.o -e bsub.e -q normal ./nextflow-pipelines/bin/solo/nextflow.sh > bjob.id
echo finished bsub
cat bjob.id
