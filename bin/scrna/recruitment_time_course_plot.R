#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## added chron lib there (in /lustre, outside the singularity container)
Sys.setenv(LANG="en_US.UTF-8")
lustre_libs = "/lustre/scratch118/humgen/resources/rlibs3.6.0"
.libPaths(lustre_libs)

library(dplyr)
library(reshape2)
library(ggplot2)
library(chron)
library(ggrepel)
library(tidyverse)

# inputs:
#####################
minimal_spreadsheet = args[1]
out_file <- "recruitment_time_course_plot.pdf"

# script:
#####################
minimal_spreadsheet = read.csv(minimal_spreadsheet, sep = "\t")
minimal_spreadsheet$date_of_sample=as.Date(minimal_spreadsheet$date_of_sample, "%d-%m-%Y")
minimal_spreadsheet=minimal_spreadsheet %>% filter(!is.na(protocol), protocol %in% c("tissue_v2","blood_final")) #%>%
minimal_spreadsheet=minimal_spreadsheet[!(minimal_spreadsheet$disease_status=="healthy" & minimal_spreadsheet$biopsy_type=="Blood"),]
for_timecourse = minimal_spreadsheet[,c('date_of_sample', 'disease_status', "biopsy_type")]

df <- subset(for_timecourse, !is.na(date_of_sample)) %>%
    dplyr::arrange(date_of_sample)

df_list <- list()
i <- 1
disease_status <- unique(df$disease_status)
biopsy_type <- unique(df$biopsy_type)
for (biopsy in biopsy_type) {
    df_tmp <- subset(df,
            biopsy_type == biopsy
        ) %>%
        dplyr::mutate(value = 1, disease_status = "total_samples") %>%
        dplyr::mutate(rolling_total = cumsum(value)) %>%
        dplyr::select(
            date_of_sample,
            rolling_total,
            biopsy_type,
            disease_status
        ) %>%
        as.data.frame
    df_list[[i]] <- df_tmp
    i <- i + 1
    for (disease in disease_status) {
        df_tmp <- subset(df,
                biopsy_type == biopsy &
                disease_status == disease
            ) %>%
            dplyr::mutate(value = 1) %>%
            dplyr::mutate(rolling_total = cumsum(value)) %>%
            dplyr::select(
                date_of_sample,
                rolling_total,
                biopsy_type,
                disease_status
            ) %>%
            as.data.frame
        df_list[[i]] <- df_tmp
        i <- i + 1
    }
}
df_plt <- data.frame(data.table::rbindlist(df_list)) %>%
    dplyr::mutate(
        tempvar = paste(biopsy_type, disease_status, sep = "-"),
        rowid = paste(biopsy_type, disease_status, date_of_sample, sep = "-"),
        rolling_total_text = ""
    ) %>%
    dplyr::filter(
        !(tempvar %in% c("Blood-total_samples", "Rectum-total_samples"))
    )

# Only plot every 5
n <- 5
df_agg <- df_plt
df_tmp <- df_agg %>%
    dplyr::group_by(biopsy_type, disease_status) %>%
    dplyr::filter(dplyr::row_number() %% 5 == 1) %>%
    dplyr::mutate(
        rolling_total_text = rolling_total
    ) %>%
    as.data.frame
df_plt <- subset(df_agg, !(rowid %in% df_tmp$rowid))
df_plt <- rbind(df_plt, df_tmp) %>%
    dplyr::arrange(date_of_sample) %>%
    as.data.frame

# Be sure the first and last are included
df_agg <- df_plt
df_tmp <- df_agg %>%
    dplyr::arrange(date_of_sample) %>%
    dplyr::group_by(biopsy_type, disease_status) %>%
    dplyr::filter(dplyr::row_number()==1 | dplyr::row_number()==dplyr::n()) %>%
    dplyr::mutate(
        rolling_total_text = rolling_total
    ) %>%
    as.data.frame
df_plt <- subset(df_agg, !(rowid %in% df_tmp$rowid))
df_plt <- rbind(df_plt, df_tmp) %>%
    dplyr::arrange(date_of_sample) %>%
    as.data.frame

# Clear out the TI instances where the total == cd
filt <- df_plt$biopsy_type == "TI" &
    df_plt$disease_status == "total_samples" &
    df_plt$rolling_total < 20
df_plt$rolling_total_text[filt] <- ""

# Clean up disease status name
df_plt$disease_status <- factor(
    df_plt$disease_status,
    levels = c("total_samples", "cd", "healthy"),
    labels = c("Total", "Crohn's disease", "Healthy")
)

# Plot
plt3 <- ggplot2::ggplot(df_plt, ggplot2::aes(
    date_of_sample,
    y = rolling_total,
    group = disease_status
))
plt3 <- plt3 + ggplot2::geom_line(
    ggplot2::aes(color = disease_status),
    alpha = 0.75
)
plt3 <- plt3 + ggplot2::geom_point(
    ggplot2::aes(color = disease_status),
    alpha = 0.25
)
plt3 <- plt3 + ggplot2::facet_wrap(
    ~ biopsy_type,
    nrow = 3,
    strip.position = "left",
    scales = "free_y"
 )
plt3 <- plt3 + ggrepel::geom_text_repel(
    ggplot2::aes(label = rolling_total_text),
    force = 3
)
plt3 <- plt3 + ggplot2::theme_bw(base_size = 12)
plt3 <- plt3 + ggplot2::scale_color_brewer(palette = "Dark2")
plt3 <- plt3 + ggplot2::labs(
    x = "Date sample collected",
    y = "Cumulative number of biopsies",
    color = "Sample type"
)
plt3 <- plt3 + ggplot2::theme(
    legend.position = "bottom",
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 0.95, vjust = 0.2))
plt3


# write output pdf figure
pdf(file = out_file, height = 5, width = 5)
print(plt3)
dev.off()
