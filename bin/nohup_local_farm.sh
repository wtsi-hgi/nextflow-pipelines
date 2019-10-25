#!/usr/bin/env bash
export PATH=/software/hgi/installs/conda/condabin/:$PATH
export LSB_DEFAULTGROUP=hgi
export GOPATH=/lustre/scratch115/projects/interval_wgs/nextflow/wr_latest/wr/gopath
export PATH=/lustre/scratch115/projects/interval_wgs/nextflow/wr0.19.9:$PATH
export CAPSULE_CACHE_DIR=/software/hgi/installs/conda/capsule_cache
export CONDA_ENVS_DIRS=/software/hgi/installs/conda/conda_envs
export CONDA_PKGS_DIRS=/software/hgi/installs/conda/conda_dirs
export PATH=/software/hgi/installs/conda/conda_envs/nextflow/bin:$PATH
export PATH=/software/singularity-v3.2.0/bin:$PATH
export PATH=/software/hgi/installs/basespace_cli/:$PATH
# conda activate $CONDA_ENVS_DIRS/nextflow

rm -f hs_err_pid* && rm -f timeline* && rm -f trace* && rm -rf report* && rm -f bsub.o && rm -f bsub.e && rm -f .nextflow.log && rm -f nohup.log && rm -f nohup.out && nohup ./bin/nextflow.sh > nohup.log 2>&1 &
