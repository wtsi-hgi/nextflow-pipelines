library(tidyverse)
setwd("~/Documents/workspace/hgi_projects/ibdx10/")
ibd_exomes <- readr::read_tsv("data/ss_wh_info_12_oct_2019.txt")
passing_files_per_sample <- ibd_exomes %>% filter(manual_qc == 1) %>% group_by(supplier_name, sanger_sample_id, id_study_lims) %>% summarise(n_files = n())
samples_to_read <- passing_files_per_sample %>% filter(n_files == 1) %>% select(supplier_name, sanger_sample_id, id_study_lims)
multifile_samples <- passing_files_per_sample %>% filter(n_files > 1) %>% select(supplier_name, sanger_sample_id, id_study_lims)

unpassed_files_per_sample <- ibd_exomes %>% filter(is.na(manual_qc) | (manual_qc != 1)) %>% group_by(supplier_name, sanger_sample_id, id_study_lims) %>% summarise(n_files = n())
unpassed_sample_count <- length(unique(unpassed_files_per_sample$supplier_name))
