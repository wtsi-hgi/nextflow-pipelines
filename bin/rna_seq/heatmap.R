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
data = read_tsv(args[1]) 
print(data)

meths = c("euclidean", "maximum", "manhattan", "canberra","binary","minkowski")
meths_num = 6
meths_hclust = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
meths_hclust_num = 2

as_matrix = data %>% select(-ENSEMBL_ID) %>% as.matrix
row.names(as_matrix) = data$ENSEMBL_ID
as_matrix %>% dim
# as_matrix = varFilter(as_matrix, var.func=IQR, var.cutoff=0.90, filterByQuantile=TRUE)
as_matrix = varFilter(as_matrix, var.func=IQR, var.cutoff=0.97, filterByQuantile=TRUE)
as_matrix %>% dim
as_matrix %>% class

# center and scale vars..
# directly from values
#as_mat %<>% mutate_all(funs(scale(.) %>% as.vector)) # vars= everything()) #c("y","z"))

#.. or from their rank
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
as_mat_raw = as_matrix
as_matrix %<>% t %>% as.data.frame
as_matrix %<>%
  mutate_all(funs(range01(rank(.,na.last = F), na.rm = T) %>% as.vector))
 #, .vars = vars(-ENSEMBL_ID)) 
row.names(as_matrix) = colnames(as_mat_raw)

as_mat_raw[1:5,1:5]
as_matrix[1:5,1:5]


# dist orders rows
dend_r <- t(as_matrix) %>% dist(method = meths[meths_num]) %>% 
  hclust(method = meths_hclust[meths_hclust_num] ) %>% 
  as.dendrogram %>% 
  ladderize %>%
  color_branches(k=3) 

dend_c <- as_matrix %>% 
  dist(method = meths[meths_num]) %>% 
  hclust(method = meths_hclust[meths_hclust_num] ) %>% 
  as.dendrogram %>% 
  ladderize %>%
  color_branches(k=3)

cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cols = c(rev(cool), rev(warm))
mypalette <- colorRampPalette(cols)(255)


dend_cols = tibble(ID_key = row.names(as_matrix)) %>% left_join(
  tibble(ID_key = labels(dend_c), 
                   col = get_leaves_branches_col(dend_c))) # %>% dplyr::left_join(pheno %>% filter(position == 2))
	
dend_cols %>% dim	



pdf('outputs/salmon_heatmap_toppc.pdf', width=16, height = 10)
gplots::heatmap.2(as_matrix %>% as.matrix, 
main = "Salmon: made with the 3% most variable ensembl gene IDs",
                  srtCol = 35,
                  Rowv = dend_c,
                  Colv = dend_r,
                  trace="none", hline = NA, tracecol = "darkgrey",         
                  margins =c(4,10),      
                  key.xlab = "",
                  denscol = "grey",
                  density.info = "density", cexRow = 1,
                  col = mypalette, cexCol = 0.15 #, 
#RowSideColors = dend_cols$pos_colors,
#labRow = paste0(dend_cols$ID_key, ' - ',dend_cols$sanger_sample_id)
		  )
#legend("left", 
#legend=levels(as.factor(dend_cols$groups)), 
#fill=unique(dend_cols$pos_colors)
#, bg = 'white', title="samples groups:", cex=1)

#legend("bottomleft", 
#  legend=levels(as.factor(dend_cols$col)), 
#  fill=levels(as.factor(dend_cols$col)), 
# bg = 'white', title="unbiased groups:", cex=1)
dev.off()

print(dend_cols)

# PCA
pdf(paste0('outputs/','salmon_PCA_unbiased_toppc','.pdf'), height = 10, width =10 ) 
ir.pca <- prcomp(as_matrix)
g <- ggbiplot(ir.pca, choices = 1:2,  #obs.scale = 1, var.scale = 1, 
              labels.size = 3, 
              size = 11, 
              labels = dend_cols$ID_key,
              #labels = paste0(dend_cols$ID_key, ' - ',dend_cols$sanger_sample_id),
              # groups = dend_cols$col, 
              ellipse = FALSE, circle = TRUE, var.axes = F) + 
  ggtitle('Salmon: made with the 3% most variable ensembl gene IDs') 
g <- g + scale_color_discrete(name = '') + 
  geom_point(aes(colour = dend_cols$col),size = 3,alpha = 0.5)
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')  + 
 scale_color_manual(name = 'unbiased groups:',
  labels=levels(as.factor(dend_cols$col)), 
  values=levels(as.factor(dend_cols$col))
  ) +
theme_bw()
print(g);dev.off()
