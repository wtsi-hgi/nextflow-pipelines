#!/usr/bin/env bash
## on farm5
export PATH=/software/singularity-v3.5.1/bin/:$PATH
export PATH=/software/hgi/installs/anaconda3/condabin:$PATH
export PATH=/software/hgi/installs/anaconda3/envs/nextflow/bin:$PATH
export CONDA_ENVS_DIRS=/software/hgi/installs/anaconda3/envs
export CONDA_PKGS_DIRS=/software/hgi/installs/anaconda3/pkgs
#export CAPSULE_CACHE_DIR=/lustre/scratch114/projects/gains_team282/Nextflow-RNAseq-gains/capsules
export CAPSULE_CACHE_DIR=/software/hgi/installs/anaconda3/capsule_cache
eval "$(conda shell.bash hook)"
conda activate $CONDA_ENVS_DIRS/nextflow

echo starting bsub
rm -f hs_err_pid* && rm -f timeline* && rm -f trace* && rm -rf report* && rm -f bsub.o && rm -f bsub.e && rm -f .nextflow.log && bsub -G hgi -R'select[mem>8000] rusage[mem=8000] span[hosts=1]' -M 8000 -n 2 -o bsub.o -e bsub.e -q yesterday ./nextflow-pipelines/bin/crispr/nextflow.sh > bjob.id
echo finished bsub
cat bjob.id
