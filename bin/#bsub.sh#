#!/usr/bin/env bash
export LSB_DEFAULTGROUP=hgi
export CAPSULE_CACHE_DIR=/lustre/scratch115/projects/interval_wgs/nextflow/cache
export GOPATH=/lustre/scratch115/projects/interval_wgs/nextflow/wr_latest/wr/gopath
export PATH=/lustre/scratch115/projects/interval_wgs/nextflow/wr0.19.9:$PATH
export CONDA_ENVS_DIRS=/lustre/scratch115/realdata/mdt2/projects/bioaid/mercury_gn5/bioaid/.conda
export CONDA_PKGS_DIRS=/lustre/scratch115/realdata/mdt2/projects/bioaid/mercury_gn5/bioaid/.conda
export PATH=/software/singularity-v3.2.0/bin:$PATH

echo starting bsub
rm -f hs_err_pid* && rm -f timeline* && rm -f trace* && rm -rf report* && rm -f bsub.o && rm -f bsub.e && rm -f .nextflow.log && bsub -G hgi -R'select[mem>8000] rusage[mem=8000] span[hosts=1]' -M 8000 -n 2 -o bsub.o -e bsub.e -q yesterday ./bin/nextflow.sh > bjob.id
echo finished bsub
cat bjob.id
