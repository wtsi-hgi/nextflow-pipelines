#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {stop("arguments needed", call.=FALSE)}

library(tximport)
library(magrittr)
library(readr)

experiment_df = read_tsv(args[0])
txi_gene_counts = read.csv(args[1])
txi_transcript_counts = read.csv(args[2])
txi_tpm_counts = read.csv(args[3]) # tx.salmon.scale <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM"
write.csv(tx.salmon.scale$counts, "txi_lengthScaledTPM_gene_counts.csv")

#######
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = txi_gene_counts,
                              colData = experiment_df,
                              design = ~ experiment_df$model_formula[0])

design(ddsMF) <- formula(~ type + condition)
ddsMF <- DESeq(ddsMF)
Again, we access the results using the results function.



save.image("./deseq2.rdata")
