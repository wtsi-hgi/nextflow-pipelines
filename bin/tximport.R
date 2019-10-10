#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

library(ensembldb)
library(args[1]) # EnsDb.Hsapiens.v91
edb <- args[1] 
organism(edb)

Tx <- transcripts(edb, return.type="DataFrame")
tx2gene <- Tx[,c(1,7)]
head(tx2gene)

library(readr)
library(tximport)
files <- file.path(dir, "salmon", samples$run, "quant.sf.gz")
names(files) <- paste0("sample", 1:6)
all(file.exists(files))

# gene-level summarization 
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
head(txi$counts)

# We can avoid gene-level summarization by setting txOut=TRUE,
# giving the original transcript level estimates as a list of matrices.
txi.tx <- tximport(files, type = "salmon", txOut = TRUE)

saveRDS(sce, file="./tximport.rds")
