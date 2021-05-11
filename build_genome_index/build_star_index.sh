#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate /lustre/scratch118/humgen/resources/conda/star

mkdir star_index_Homo_sapiens.GRCh38.99_75bp
STAR \
        --runMode genomeGenerate \
        --runThreadN 4 \
        --sjdbGTFfile Homo_sapiens.GRCh38.99.gtf \
        --sjdbOverhang 75 \
        --genomeDir star_index_Homo_sapiens.GRCh38.99_75bp \
        --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa

mkdir star_index_Homo_sapiens.GRCh38.99_100bp
STAR \
        --runMode genomeGenerate \
        --runThreadN 4 \
        --sjdbGTFfile Homo_sapiens.GRCh38.99.gtf \
        --sjdbOverhang 100 \
        --genomeDir star_index_Homo_sapiens.GRCh38.99_100bp \
        --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa
