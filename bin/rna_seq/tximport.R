#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)}

local_libs = "/lustre/scratch118/humgen/resources/conda/star/lib/R/library"
.libPaths(local_libs)
.libPaths()

library(AnnotationHub)
library(tximport)
library(magrittr)
library(readr)
queries = query(AnnotationHub(), c(args[1], "Homo sapiens"))
edb = queries[[1]]
edb
Tx <- transcripts(edb, return.type="DataFrame")
head(Tx)
tx2gene <- Tx[,c(9,7)]
## tx_id         gene_id
##           <character>     <character>
##1      ENST00000387314 ENSG00000210049
head(tx2gene)

files_df = read.csv(args[2], header = FALSE)
files <- file.path(files_df[,1])

names(files) <- as.character(files)
all(file.exists(files))

# gene-level summarization 
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
head(txi$counts)
write.csv(txi$counts, "txi_gene_counts.csv")

# We can avoid gene-level summarization by setting txOut=TRUE,
# giving the original transcript level estimates as a list of matrices.
txi.tx <- tximport(files, type = "salmon", txOut = TRUE)
write.csv(txi.tx$counts, "txi_transcript_counts.csv")
                                        #
# count table from TPM 
# accounts for transcript length changes across samples and library size differences.
tx.salmon.scale <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                            countsFromAbundance = "lengthScaledTPM")
write.csv(tx.salmon.scale$counts, "txi_lengthScaledTPM_gene_counts.csv")

save.image("./tximport.rdata")
