#!/usr/bin/env bash
export GOPATH=/lustre/scratch115/projects/interval_wgs/nextflow/wr_latest/wr/gopath
export PATH=/lustre/scratch115/projects/interval_wgs/nextflow/wr0.19.9:$PATH
export LSB_DEFAULTGROUP=hgi
export CAPSULE_CACHE_DIR=/software/hgi/envs/
export CONDA_ENVS_DIRS=/software/hgi/envs/
export CONDA_PKGS_DIRS=/software/hgi/envs/
# conda activate /software/hgi/envs

echo starting bsub
rm -f hs_err_pid* && rm -f timeline* && rm -f trace* && rm -rf report* && rm -f bsub.o && rm -f bsub.e && rm -f .nextflow.log && bsub -G hgi -R'select[mem>8000] rusage[mem=8000] span[hosts=1]' -M 8000 -n 2 -o bsub.o -e bsub.e -q basement ./bin/nextflow.sh > bjob.id
echo finished bsub
cat bjob.id
