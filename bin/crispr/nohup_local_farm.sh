#!/usr/bin/env bash

#farm 4 old config:
#export PATH=/software/hgi/installs/conda/condabin:$PATH
#export LSB_DEFAULTGROUP=hgi
#export GOPATH=/lustre/scratch115/projects/interval_wgs/nextflow/wr_latest/wr/gopath
#export PATH=/lustre/scratch115/projects/interval_wgs/nextflow/wr0.19.9:$PATH
#export CAPSULE_CACHE_DIR=/software/hgi/installs/conda/capsule_cache
#export CONDA_ENVS_DIRS=/software/hgi/installs/conda/conda_envs
#export CONDA_PKGS_DIRS=/software/hgi/installs/conda/conda_dirs
#export PATH=/software/hgi/installs/conda/conda_envs/nextflow/bin:$PATH
#export PATH=/software/singularity-v3.2.0/bin:$PATH
#export PATH=/software/hgi/installs/basespace_cli/:$PATH
# conda activate $CONDA_ENVS_DIRS/nextflow

## on farm5
export PATH=/software/singularity-v3.4.0/bin/:$PATH
export PATH=/software/hgi/installs/anaconda3/condabin:$PATH
export PATH=/software/hgi/installs/anaconda3/envs/nextflow/bin:$PATH
export CONDA_ENVS_DIRS=/software/hgi/installs/anaconda3/envs
export CONDA_PKGS_DIRS=/software/hgi/installs/anaconda3/pkgs
#export CAPSULE_CACHE_DIR=/lustre/scratch114/projects/gains_team282/Nextflow-RNAseq-gains/capsules
export CAPSULE_CACHE_DIR=/software/hgi/installs/anaconda3/capsule_cache
##### conda activate $CONDA_ENVS_DIRS/nextflow

rm -f hs_err_pid* && rm -f timeline* && rm -f trace* && rm -rf report* && rm -f bsub.o && rm -f bsub.e && rm -f .nextflow.log && rm -f nohup.log && rm -f nohup.out && nohup ./nextflow-pipelines/bin/crispr/nextflow.sh > nohup.log 2>&1 &
