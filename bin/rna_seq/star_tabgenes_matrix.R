#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(AnnotationHub)
library(tximport)
library(magrittr)
library(readr)
files_df = read.csv(args[2], header = FALSE)
files <- file.path(files_df[,1])
names(files) <- as.character(files) %>% gsub(".quant.sf","",.)
all(file.exists(files))
#################

library(gplots)
library(biomaRt)
library(reshape2)
library(ggplot2)
library(broom)
library(magrittr)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(stringr)
library(tibble)
library(readr)
library(openxlsx)
library(dendextend)
library(RColorBrewer)
library(genefilter)
library(ggbiplot)
library(DESeq2)
library(ggrepel)
select = dplyr :: select

dir.create('./outputs', showWarnings = FALSE)

rename = dplyr::rename
select = dplyr::select

${samplename}.ReadsPerGene.out.tab

write.csv(tx.salmon.scale$counts, "txi_lengthScaledTPM_gene_counts.csv")
