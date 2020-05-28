#!/usr/bin/env bash

git checkout rna_seq_dddmouse
git pull --recurse-submodules 
git submodule sync
git submodule update --init --recursive --remote
git add .gitmodules && git add nextflow-pipelines && git commit -m "modified submodule URL" && git push -u origin

# farm5 config
export http_proxy=http://wwwcache.sanger.ac.uk:3128
export https_proxy=http://wwwcache.sanger.ac.uk:3128
export HTTP_PROXY=http://wwwcache.sanger.ac.uk:3128
export HTTPS_PROXY=http://wwwcache.sanger.ac.uk:3128
export PATH=/software/singularity-v3.5.1/bin/:$PATH
#export PATH=/software/hgi/installs/anaconda3/condabin:$PATH
#export PATH=/software/hgi/installs/anaconda3/envs/nextflow/bin:$PATH
#export CONDA_PKGS_DIRS=/software/hgi/installs/anaconda3/pkgs
#export CAPSULE_CACHE_DIR=/lustre/scratch114/projects/gains_team282/Nextflow-RNAseq-gains/capsules
#export CAPSULE_CACHE_DIR=/software/hgi/installs/anaconda3/capsule_cache
eval "$(conda shell.bash hook)"
export CONDA_ENVS_DIRS=/lustre/scratch118/humgen/resources/conda_envs
conda activate $CONDA_ENVS_DIRS/nextflow

rm -f hs_err_pid* && rm -f timeline* && rm -f trace* && rm -rf report* && rm -f bsub.o && rm -f bsub.e && rm -f .nextflow.log && rm -f nohup.log && rm -f nohup.out && nohup ./nextflow-pipelines/bin/mosdepth/nextflow.sh > nohup.log 2>&1 &
