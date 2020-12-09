#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## added chron lib there (in /lustre, outside the singularity container)
Sys.setenv(LANG="en_US.UTF-8")
lustre_libs = "/lustre/scratch118/humgen/resources/rlibs3.6.0"
.libPaths(lustre_libs)

library(lubridate)
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
library(hdf5r)
library(cowplot)
library(Seurat)
library(Hmisc)
library(chron)
select = dplyr :: select
rename = dplyr :: rename
mutate = dplyr :: mutate

# inputs:
#####################
minimal_spreadsheet = args[1]
imeta_info = args[2]
sync_status = args[3]

print(minimal_spreadsheet)
print(imeta_info)
print(sync_status)

####################
# custom function to convert date_of_sample to season
date_to_season <- function(input_date, input_format = "%d-%m-%Y") {
    as_date = as.Date(input_date, input_format)

    begin_winter = as.Date(paste0('21-12-',as.character(years(as_date))), "%d-%m-%Y")
    end_winter = as.Date(paste0('19-3-',as.character(years(as_date))), "%d-%m-%Y")  
    is_in_winter = ((as.numeric(as_date - begin_winter) >= 0) | (as.numeric(as_date - end_winter) < 0))
    
    begin_summer = as.Date(paste0('21-6-',as.character(years(as_date))), "%d-%m-%Y")
    end_summer = as.Date(paste0('22-9-',as.character(years(as_date))), "%d-%m-%Y")  
    is_in_summer = ((as.numeric(as_date - begin_summer) >= 0) && (as.numeric(as_date - end_summer) < 0))
    
    begin_spring = as.Date(paste0('20-3-',as.character(years(as_date))), "%d-%m-%Y")
    end_spring = as.Date(paste0('20-6-',as.character(years(as_date))), "%d-%m-%Y")  
    is_in_spring = ((as.numeric(as_date - begin_spring) >= 0) && (as.numeric(as_date - end_spring) < 0))
    
    begin_autumn = as.Date(paste0('23-9-',as.character(years(as_date))), "%d-%m-%Y")
    end_autumn = as.Date(paste0('20-12-',as.character(years(as_date))), "%d-%m-%Y")  
    is_in_autumn = ((as.numeric(as_date - begin_autumn) >= 0) && (as.numeric(as_date - end_autumn) < 0))
    
    if (is_in_winter) {season = 'winter'} 
    else if (is_in_summer) {season = 'summer'} 
    else if (is_in_spring) {season = 'spring'} 
    else if (is_in_autumn) {season = 'autumn'} 
    else {season = 'get_season_error'} 
    
  return(season)
}

#convert string epithelial_immune_ratio to float
to_float <- function(input) {
  if (is.na(input) | input == 'NA') {return(NA)}  
    
  split_ = str_split(input, ':')[[1]]

    if (length(split_) <2) {return(NA)}    
    else if (as.numeric(split_[1]) == 0) {
        
        return(round(as.numeric(split_[2])/as.numeric(0.0001),3))
    } 
    else {
        return(round(as.numeric(split_[2])/as.numeric(split_[1]),3))
    } 
}

to_ratio <- function(input) {
  if (is.na(input) | input == 'NA') {return(NA)}  
    
  split_ = str_split(input, ':')[[1]]

    if (length(split_) <2) {return(NA)}    
    else if (as.numeric(split_[1]) == 0) {
        
        return(paste0(as.character(as.numeric(split_[2])),":","0.0001"))
    } 
    else {
        return(paste0(as.character(as.numeric(split_[2])),":",as.character(as.numeric(split_[1]))))
    } 
}
#####################

# script:
#####################

## load in cellranger csv summaries
csv_files =
        tibble(csv_file = list.files('.')) %>% rowwise %>%
        filter(str_detect(csv_file, '.metrics_summary.csv')) %>%
        mutate(sanger_sample_id = gsub(".metrics_summary.csv","",csv_file))

csv_files$content = map(csv_files$csv_file, ~ read_csv(.x))
csv_files %<>% unnest(content)

minimal_spreadsheet = read_tsv(minimal_spreadsheet) %>% filter(!is.na(sanger_sample_id))

columns_to_remove =
c("frozen_banked","qc_complete","patient_number","supplier_sample_name","biopsy_type_comments","ethnicity","patient_comments","smoking_comments","medications","ibd_bioresource","viability","digest_time","time_to_processing","total_time","extra_biopsy_OCT_embedded","extra_biopsy_organoid","experimental_comments","other_comments","n_cells_input","n_cells_output","cell_efficiency","mean_reads_per_cell","median_genes_per_cell","fraction_reads_per_cell","pilot","cell_ranger_complete","extra_biopsy_followup","final_tissue_protocol","protocol_version")
columns_rm = colnames(minimal_spreadsheet)[colnames(minimal_spreadsheet) %in% columns_to_remove]
print(columns_rm)

imeta_info = read_tsv(imeta_info) %>% filter(!is.na(sanger_sample_id)) %>% select(-X10) %>%
    mutate(id_run = paste0('id_run_',as.character(id_run)),
           lane = paste0('lane_',as.character(lane)),
           library_id = paste0('library_id_',as.character(library_id)),
           study_id = paste0('study_id_',as.character(study_id))) %>%
    dplyr::rename(imeta_study = study)
sync_status = read_tsv(sync_status) %>% filter(!is.na(sanger_sample_id))

minimal_spreadsheet %<>%
  select(-one_of(columns_rm)) %>%
  left_join(imeta_info) %>%
  left_join(sync_status) %>%
  filter(!is.na(disease_status),!is.na(biopsy_type)) %>%  
  filter(frozen_processed != 'yes', !is.na(frozen_processed), frozen_processed != '') %>%  
  filter(disease_status != 'uc') %>%  
  rowwise %>%
  mutate(season_sample_collected = date_to_season(date_of_sample)) %>%
  mutate(immune_epithelial_ratio_float = to_float(epithelial_immune_ratio)) %>%
  mutate(immune_epithelial_ratio = to_ratio(epithelial_immune_ratio)) %>%
  mutate(inflammation_status_grouped =
             ifelse(inflammation_status %in% c('moderate','severe'), 'moderate_to_severe', inflammation_status)) %>%
  ungroup %>%
  mutate(
       smoker_at_time_of_biopsy = ifelse(smoking_status == 'ex-smoker', 'no', smoking_status),
        month_sample_collected = lubridate::month(date_of_sample), 
        biopsy_type_original = biopsy_type,
         biopsy_type = gsub("neoti","ti",biopsy_type_original) %>%
           gsub("r","rectum",.) %>%
           gsub("ti","TI",.) %>%
           Hmisc::capitalize(.),
         disease_status_pretty = gsub("cd","Crohn's disease", disease_status) %>% Hmisc::capitalize(.),
         disease_status_pretty = gsub("uc","Ulcerative colitis", disease_status_pretty) %>% Hmisc::capitalize(.),
         disease_status_pretty = gsub("Uc","Ulcerative colitis", disease_status_pretty) %>% Hmisc::capitalize(.),
         chromium_time = chron::times(chromium_time),
         collection_time = chron::times(collection_time),
         time_to_chromium_processing =
             chron::times(chromium_time) - chron::times(collection_time),
         hours_to_chromium_processing =
             hours(time_to_chromium_processing) +
             minutes(time_to_chromium_processing)/60 +
             seconds(time_to_chromium_processing)/3600,
         sequenced = ifelse(is.na(total_reads) ,"Not sequenced","Sequenced"),
         cellranger_synced = ifelse(is.na(sync_status) | sync_status != 'synced' ,"No","Yes"))  %>%
  mutate(sample_status = ifelse(cellranger_synced == "Yes","Sequenced + Cell Ranger",sequenced) %>%
           Hmisc::capitalize(.)) %>%
  select(-sequenced, -sync_status, -cellranger_synced) %>%
       left_join(csv_files) %>%
       select(-csv_file)

print(minimal_spreadsheet)

write_tsv(minimal_spreadsheet, "samples_metainfo.tsv")
write.xlsx(as.data.frame(minimal_spreadsheet), "samples_metainfo.xlsx")

final_only = minimal_spreadsheet %>%
    filter(!is.na(protocol), protocol %in% c("tissue_v2","blood_final")) %>%
    ungroup %>%
    group_by(biopsy_type,disease_status_pretty, sample_status) %>%
    mutate(n_group = dplyr::n()) %>%
    distinct(biopsy_type, disease_status_pretty, sample_status, .keep_all = TRUE)

plt = ggplot2::ggplot(final_only , aes(
  biopsy_type,fill = sample_status,label = n_group,y = n_group)) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single")) +
  facet_wrap(~ disease_status_pretty) +
  geom_text(position = position_dodge2(width = 0.9, preserve = "single"), vjust = -0.1) +
  ylim(0, max(final_only$n_group)+5) +
  theme_bw(base_size = 12) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Tissue",
  y = "Number of biopsies",
  title = "protocol == tissue_v2 | blood_final",
  fill = "") +
  theme(legend.position="bottom")

out_file <- "samples_status.pdf"
pdf(file = out_file, height = 4, width = 5.5)
print(plt)
dev.off()

plot2 = minimal_spreadsheet %>% 
    mutate(biopsy_disease = paste0(biopsy_type,'_',disease_status_pretty)) %>% 
    filter(biopsy_disease != 'Blood_Healthy', biopsy_disease != 'Rectum_Crohn\'s disease') %>%
    filter(!is.na(protocol), protocol %in% c("tissue_v2","blood_final")) %>%
    ungroup %>%
    group_by(biopsy_type,disease_status_pretty, sample_status) %>%
    mutate(n_group = dplyr::n()) %>%
    distinct(biopsy_type, disease_status_pretty, sample_status, .keep_all = TRUE)

plt2 = ggplot2::ggplot(plot2 , aes(
  biopsy_type,fill = sample_status,label = n_group,y = n_group)) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single")) +
  facet_wrap(~ disease_status_pretty, scales='free_x') +
  geom_text(position = position_dodge2(width = 0.9, preserve = "single"), vjust = -0.1) +
  ylim(0, max(final_only$n_group)+5) +
  theme_bw(base_size = 12) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Tissue",
  y = "Number of biopsies",
  title = "protocol == tissue_v2 | blood_final",
  fill = "") +
  theme(legend.position="bottom")

out_file <- "samples_status_dropped.pdf"
pdf(file = out_file, height = 4, width = 5.5)
print(plt2)
dev.off()



