#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(tximport)
library(magrittr)
library(readr)
library(gplots)
library(biomaRt)
library(reshape2)
library(ggplot2)
library(broom)
library(magrittr)
library(dplyr)
library(tidyr)
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
library(purrr)

select = dplyr :: select
rename = dplyr::rename
select = dplyr::select

dir.create('./outputs', showWarnings = FALSE)

tabs = list.files('.') %>% keep(~ str_detect(.x,'.ReadsPerGene.out.tab'))
tabs_names = map(tabs, ~ gsub('.ReadsPerGene.out.tab','',.x))
tabs_read = pmap(list(tabs, tabs_names), ~ read_tsv(.x, skip=5, col_names = c('Gene','x2','x3','Count')) %>%
               select(Gene,Count) %>% mutate(sample = .y))

reduced = tabs_read %>% purrr::reduce(bind_rows)
mat = reduced %>% spread(sample, Count)
write_tsv(mat, 'star_tabgenes_matrix.tsv')
