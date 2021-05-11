#!/bin/bash

# (april 7h 2020)
# get dna primary assembly
wget ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# get cdna 
wget ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
# get gtf annotation

# choose gtf
wget ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz
cp /lustre/scratch119/humgen/projects/interval_wgs/analysis/data/gencode31/gencode.v31.annotation.gtf.gz .

# unzip
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

gunzip Homo_sapiens.GRCh38.97.gtf.gz
gunzip gencode.v31.annotation.gtf.gz
