#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {stop("At least one argument must be supplied", call.=FALSE)}
variance_threshold = as.numeric(args[2])
pcent_samples_threshold = as.numeric(args[3])
print(paste0('variance_threshold: ', as.character(variance_threshold)))
print(paste0('pcent_samples_threshold: ', as.character(pcent_samples_threshold)))

library(genefilter)
library(reshape2)
library(ggplot2) 
library(magrittr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(readr)
library(openxlsx)
library(dendextend)
library(RColorBrewer)
library(purrr)
library(readr)
library(rsnps)
library(purrrlyr)
library(data.table)
library(tidyr)
library(stringr)

select = dplyr::select
rename = dplyr::rename
select = dplyr::select

bed = read_tsv(args[1])

as_df = bed %>% unite('ID',`#chr`:gene_id,sep='___') %>% as.data.frame
row.names(as_df) <- as_df$ID

print(paste0('number of genes before filtering for missingness or 0 counts: ', as.character(dim(as_df)[1])))
as_df %<>% select(-ID) %>% 
  purrr::discard(~sum(is.na(.x) | .x == 0)/length(.x) >= pcent_samples_threshold)
print(paste0('number of genes after filterting: ', as.character(dim(as_df)[1])))

print(paste0('number of genes before filtering with variance cutoff: ', as.character(dim(as_df)[1])))
as_df <- as.data.frame(varFilter(data.matrix(as_df), var.func=IQR,
                                 var.cutoff=variance_threshold, filterByQuantile=TRUE))
print(paste0('number of genes after filterting: ', as.character(dim(as_df)[1])))

IDs_cols = c('#chr','start','end','gene_id')
as_df$ID = row.names(as_df)
as_df %<>% 
separate(ID, into = IDs_cols, sep='___') %>%
select(`#chr`,start, end ,gene_id, everything())

write_tsv(as_df, 'expression_data.filtered.bed')
