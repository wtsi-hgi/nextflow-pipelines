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
setwd("/Users/gn5/interval_rnaseq/inputs/")

# imeta qu -z seq -d study = 'HG_Transcriptome sequencing in the INTERVAL cohort' and target = 1 and manual_qc = 1
# ----
#  collection: /seq/illumina/runs/29/29810/lane4/plex93
# dataObj: 29810_4#93.cram
# ----
#  collection: /seq/illumina/runs/29/29810/lane4/plex94
# dataObj: 29810_4#94.cram

# csv comes from seq Untitled.sql
pheno = read_delim('5591.csv', 
                   quote = "'", delim=',', col_names = FALSE) %>% 
  as_tibble %>% setNames(c('manual_qc',
'study.id_study_lims','study.name','sample.common_name','sample.description','samplename',
'sample.sanger_sample_id','sample.supplier_name','sample.donor_id','sample.gender',
'iseq_flowcell.last_updated','iseq_product_metrics.num_reads','iseq_product_metrics.id_run',
'iseq_product_metrics.position','iseq_product_metrics.tag_index')) %>% mutate(
  irods_path = paste0('/seq/illumina/runs',
                      '/',substr(iseq_product_metrics.id_run, start = 1, stop = 2),
                      '/',iseq_product_metrics.id_run,
                      '/lane',iseq_product_metrics.position,
                      '/plex',iseq_product_metrics.tag_index,
                      '/',iseq_product_metrics.id_run,'_',iseq_product_metrics.position,'#',
                      iseq_product_metrics.tag_index,'.cram'))

write.xlsx(pheno, 'pheno.xlsx')
pheno %>% head # %>% View
pheno$irods_path

batch7 = read_delim('INTERVAL_IDs_Batch7.txt', delim=',', col_names = FALSE) %>%
  as_tibble %>% setNames(c('samplename'))

batch7 %<>% left_join(pheno)
write.xlsx(batch7, 'batch7_pheno.xlsx')
write_csv(batch7, 'batch7_pheno.csv') 

batch7 %>% head # View
