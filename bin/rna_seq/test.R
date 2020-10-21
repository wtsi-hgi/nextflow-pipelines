#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)}

print(args[1])

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
select = dplyr :: select

dir.create('./outputs', showWarnings = FALSE)

rename = dplyr::rename
select = dplyr::select
#data = read_tsv(args[1]) 
data = read.csv(args[1]) %>% rename(ENS  = X)
print(colnames(data))
