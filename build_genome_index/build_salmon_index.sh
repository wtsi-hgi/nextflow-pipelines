#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate /lustre/scratch118/humgen/resources/conda/star

salmon index \
-t Homo_sapiens.GRCh38.cdna.all.fa \
-i salmon_index_Homo_sapiens.GRCh38.cdna.all \
-k 31

