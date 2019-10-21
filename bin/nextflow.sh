#!/usr/bin/env bash
echo starting nextflow

export NXF_VER=19.09.0-edge
export NXF_OPTS="-Xms6G -Xmx6G -Dnxf.pool.type=sync -Dnxf.pool.maxThreads=2000"
export NXF_JAVA_HOME="/lustre/scratch115/realdata/mdt2/projects/bioaid/mercury_gn5/bioaid/env"
export JAVA_HOME="/lustre/scratch115/realdata/mdt2/projects/bioaid/mercury_gn5/bioaid/env"
export JAVA_CMD="/lustre/scratch115/realdata/mdt2/projects/bioaid/mercury_gn5/bioaid/env/bin/java"
nextflow run ./pipelines/crispr.nf -profile farm4_singularity_gn5 -resume 
#nextflow run main.nf -profile farm4_singularity_gn5 -resume 
