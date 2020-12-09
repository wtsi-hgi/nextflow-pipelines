#!/usr/bin/env bash

export http_proxy=http://wwwcache.sanger.ac.uk:3128
export https_proxy=http://wwwcache.sanger.ac.uk:3128
export HTTP_PROXY=http://wwwcache.sanger.ac.uk:3128
export HTTPS_PROXY=http://wwwcache.sanger.ac.uk:3128
export PATH=/software/singularity-v3.5.1/bin/:$PATH
#export PATH=/software/hgi/installs/anaconda3/condabin:$PATH
#export PATH=/software/hgi/installs/anaconda3/envs/nextflow/bin:$PATH
#export CONDA_ENVS_DIRS=/software/hgi/installs/anaconda3/envs
#export CONDA_PKGS_DIRS=/software/hgi/installs/anaconda3/pkgs
#export CAPSULE_CACHE_DIR=/lustre/scratch114/projects/gains_team282/Nextflow-RNAseq-gains/capsules
#export CAPSULE_CACHE_DIR=/software/hgi/installs/anaconda3/capsule_cache
eval "$(conda shell.bash hook)"
#conda activate $CONDA_ENVS_DIRS/nextflow20
conda activate nextflow20

echo starting bsub
rm -f hs_err_pid* && rm -f timeline* && rm -f trace* && rm -rf report* && rm -f bsub.o && rm -f bsub.e && rm -f .nextflow.log && bsub -G hgi -R'select[mem>5000] rusage[mem=5000] span[hosts=1]' -M 5000 -n 1 -o bsub.o -e bsub.e -q normal ./nextflow-pipelines/bin/scrna/nextflow.sh > bjob.id
echo finished bsub
cat bjob.id
