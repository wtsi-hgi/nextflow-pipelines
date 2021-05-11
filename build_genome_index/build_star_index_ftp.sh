#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate /lustre/scratch118/humgen/resources/conda/star

mkdir star_index_Homo_sapiens.GRCh38.97_75bp_ftp
STAR \
        --runMode genomeGenerate \
        --runThreadN 4 \
        --sjdbGTFfile Homo_sapiens.GRCh38.97.gtf \
        --sjdbOverhang 75 \
        --genomeDir star_index_Homo_sapiens.GRCh38.97_75bp_ftp \
        --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa

#mkdir star_index_Homo_sapiens.GRCh38.97_100bp
#STAR \
#        --runMode genomeGenerate \
#        --runThreadN 4 \
#        --sjdbGTFfile Homo_sapiens.GRCh38.97.gtf \
#        --sjdbOverhang 100 \
#        --genomeDir star_index_Homo_sapiens.GRCh38.97_100bp \
#        --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa
