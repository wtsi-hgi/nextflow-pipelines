#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## added chron lib there (in /lustre, outside the singularity container)
Sys.setenv(LANG="en_US.UTF-8")
lustre_libs = "/lustre/scratch118/humgen/resources/rlibs3.6.0"
.libPaths(lustre_libs)

library(magrittr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(readr)
select = dplyr :: select
rename = dplyr :: rename
mutate = dplyr :: mutate

#####################
input_spreadsheet = './inputs/samples_google_spreadsheet.tsv'
last_sync = './sync_status/samples_metainfo.tsv'

last_sync = read_tsv(last_sync) %>% select(sanger_sample_id, sample_status) %>%
    filter(!is.na(sanger_sample_id))

print(input_spreadsheet)
df = read_tsv(input_spreadsheet) %>% select(sanger_sample_id) %>%
    filter(!is.na(sanger_sample_id)) %>% full_join(last_sync) %>%
    mutate(sample_status = ifelse(is.na(sample_status),'na',sample_status)) %>%
    filter(sample_status != 'Sequenced + Cell Ranger') %>%
    distinct(sanger_sample_id, .keep_all = "True") %>%
    arrange(sample_status, sanger_sample_id) 

write_tsv(df, './sync_status/to_unsync.tsv') 
#        tibble(csv_file = list.files('.')) %>% rowwise %>%
 #       filter(str_detect(csv_file, '.metrics_summary.csv')) %>%
