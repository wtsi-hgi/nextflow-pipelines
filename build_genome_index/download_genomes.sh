#!/bin/bash

# http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/

# (april 7h 2020)
# get dna primary assembly
wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# get cdna 
wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
# get gtf annotation
wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
# get ncrna
wget ftp://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz

gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.99.gtf.gz
gunzip Homo_sapiens.GRCh38.ncrna.fa.gz

# combine ncrna + cdna
cp Homo_sapiens.GRCh38.ncrna.fa Homo_sapiens.GRCh38.ncrna_plus_cdna.fa
cat Homo_sapiens.GRCh38.cdna.all.fa >> Homo_sapiens.GRCh38.ncrna_plus_cdna.fa
