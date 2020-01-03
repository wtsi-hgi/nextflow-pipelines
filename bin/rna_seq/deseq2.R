#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# experiment_formula = "~ g:a"
# keep_cols = c("g;a")
# comparisons = c("gKO.a_0_vs_gWT.a_0;gKO.a_F_vs_gWT.a_F;gKO.a_T_vs_gWT.a_T;gKO.a_TF_vs_gWT.a_TF")
# keep_cols = strsplit(keep_cols, ";")[[1]]
# comparisons = strsplit(comparisons, ";")[[1]]

library(AnnotationHub)
library(tximport)
library(magrittr)
library(readr)
queries = query(AnnotationHub(), c(args[1], "Homo sapiens"))
edb = queries[[1]]
edb
Tx <- transcripts(edb, return.type="DataFrame")
head(Tx)
tx2gene <- Tx[,c(9,7)]
## tx_id         gene_id
##           <character>     <character>
##1      ENST00000387314 ENSG00000210049
head(tx2gene)

files_df = read.csv(args[2], header = FALSE)
files <- file.path(files_df[,1])

names(files) <- as.character(files) %>% gsub(".quant.sf","",.)
all(file.exists(files))

# gene-level summarization 
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
head(txi$counts)
write.csv(txi$counts, "txi_gene_counts.csv")

# We can avoid gene-level summarization by setting txOut=TRUE,
# giving the original transcript level estimates as a list of matrices.
txi.tx <- tximport(files, type = "salmon", txOut = TRUE)
write.csv(txi.tx$counts, "txi_transcript_counts.csv")
                                        #
# count table from TPM 
# accounts for transcript length changes across samples and library size differences.
tx.salmon.scale <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                            countsFromAbundance = "lengthScaledTPM")
write.csv(tx.salmon.scale$counts, "txi_lengthScaledTPM_gene_counts.csv")

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

BiocManager::install("caret", lib=".")
library(caret, lib=".")
BiocManager::install("pheatmap", lib=".")
library("pheatmap",lib=".")
BiocManager::install("apeglm", lib=".")
library(apeglm, lib=".")


dir.create('./outputs', showWarnings = FALSE)

rename = dplyr::rename
select = dplyr::select

experiment_df = read_tsv(args[3]) %>% as.data.frame
experiment_formula =  experiment_df$model_formula[1] #"~ g:a"
keep_cols = experiment_df$model_cols[1] #"~ g:a" c("g;a")
keep_cols = strsplit(keep_cols, ";")[[1]]
comparisons = experiment_df$contrast_list[1] # c("gKO.a_0_vs_gWT.a_0;gKO.a_F_vs_gWT.a_F;gKO.a_T_vs_gWT.a_T;gKO.a_TF_vs_gWT.a_TF")
comparisons = strsplit(comparisons, ";")[[1]]

# experiment_df %<>% select(g,r,a) %>% mutate_if(is.character, as.factor)
# experiment_df %<>% select(g,r,a) %>% mutate_if(is.numeric, as.factor)

# experiment_df %<>% dplyr::group_by(g,a) %>% dplyr::mutate(r=row_number()) %>% ungroup %>% as.data.frame
experiment_df %<>% mutate_if(is.character, as.factor)
experiment_df %<>% mutate_if(is.numeric, as.factor)
row.names(experiment_df) = experiment_df$sample
experiment_df = experiment_df[,keep_cols]


x = dummyVars(as.formula(experiment_formula), data = experiment_df, fullRank = TRUE)
experiment_df_matrix = data.frame(predict(x, newdata = experiment_df))
print(experiment_df_matrix)

colnames(txi$counts) == row.names(experiment_df)
    
dds <- DESeqDataSetFromTximport(txi,
                                   colData = experiment_df,
                                   design = as.matrix(experiment_df_matrix))
print(design(dds))
print(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
print(dds)
dds <- DESeq(dds)
print(resultsNames(dds))

for (comparison in comparisons) {
    comparison_list = list(
        strsplit(comparison, "_vs_")[[1]][1],
        strsplit(comparison, "_vs_")[[1]][2])

res <- results(dds,
               name=comparison,
               contrast = comparison_list)

# LFCs on the log2 scale
# apeglm provides empirical Bayes shrinkage estimators for effect sizes for a variety of GLM models;
#resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm") 
#resOrdered <- resLFC[order(resLFC$pvalue),]

resOrdered <- res[order(res$pvalue),]
print(summary(resOrdered))
#How many adjusted p-values were less than 0.1?
print(sum(resOrdered$padj < 0.1, na.rm=TRUE))
print(mcols(resOrdered)$description)

write.csv(as.data.frame(resOrdered), 
          file=paste0("./outputs/",comparison,".csv"))

padj_max = 0.1
resSig <- subset(resOrdered, padj < padj_max)
write.csv(as.data.frame(resSig), 
          file=paste0("./outputs/",comparison,".padj_inf_",as.character(padj_max),".csv"))

#It is more useful visualize the MA-plot for the shrunken log2 fold changes
#which remove the noise associated with log2 fold changes from low count genes
#without requiring arbitrary filtering thresholds.
pdf(paste0("./outputs/",comparison ,".MA_plot.pdf"))
print(plotMA(resOrdered, ylim=c(-2,2)))
dev.off()

# pdf(paste0("./outputs/",comparison, "MA_plotCounts_min.pdf"))
# print(plotCounts(dds, gene=which.min(resOrdered$padj), intgroup="g"))
# dev.off()

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
print(head(assay(vsd), 3))

df <- as.data.frame(experiment_df[,keep_cols]) %>% mutate_all(as.character)
row.names(df) = row.names(experiment_df)
colnames(df) = keep_cols
print(df)

colnames(assay(vsd)[select,]) == row.names(df)
pdf(paste0("./outputs/",comparison,".heatmap_vsd.pdf"))
print(
pheatmap(assay(vsd)[select,],
         main = "vsd heatmap",
         cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df))
dev.off()

pdf(paste0("./outputs/",comparison,".heatmap_rld.pdf"))
print(
pheatmap(assay(rld)[select,],
         main = "rld heatmap",
         cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df))
dev.off()


# pcaData <- plotPCA(vsd, intgroup=keep_cols, returnData=TRUE)
# percentVar <- round(100 * attr(pcaData, "percentVar"))

# nba_plot = ggplot(pcaData, aes(PC1, PC2, color=g, shape=a)) +
#  geom_point(size=3) +
#  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#    coord_fixed()

# pdf(paste0("./outputs/","comparison",".PCA_ggplot_vsd.pdf"))
# print(
# nba_plot + 
#  geom_label_repel(aes(label = name),
#                  box.padding   = 0.35, 
#                  point.padding = 0.5,
#                  segment.color = 'grey50') +
#  theme_classic())
# dev.off()

pdf(paste0("./outputs/","comparison",".PCA_vsd.pdf"))
print(
plotPCA(vsd, intgroup=keep_cols))
dev.off()
}

save.image("./outputs/deseq2.rdata")
