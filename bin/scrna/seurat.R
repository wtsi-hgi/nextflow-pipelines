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
library(monocle3)


select = dplyr :: select
rename = dplyr :: rename
setwd("/home/ubuntu/rsynced_nov21") 
dir.create('./Seurat_outputs', showWarnings = FALSE)

paths = read_csv('paths.csv', col_names = TRUE) %>%
	filter(raw=='filtered') # raw of filtered

#### samples to process
to_process = paths$matrixdir
sample_ids = paths$sampleid

to_process 
sample_ids

#### loop through samples
for (i in 1:length(to_process)) {
#### read 10X data

#sample_id = strsplit(found_path$filepath[i], '/')[[1]][length(strsplit(found_path$filepath[i], '/')[[1]])-1]
sample_id = sample_ids[i]
print(sample_id)
dir.create(paste0('Seurat_outputs/' ,sample_id), showWarnings = FALSE)

pbmc.data <-  Read10X(to_process[i])

#### raw stats
n_genes_raw = dim(pbmc.data)[1]
n_cells_raw = dim(pbmc.data)[2]

#### QC filters
# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells = 3,
   min.features  = 200, project = "10X_PBMC", assay = "RNA")
mito.genes <- grep(pattern = "^MT-", x = rownames(pbmc@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@assays[["RNA"]][mito.genes, ])/Matrix::colSums(pbmc@assays[["RNA"]])

mean.percent.mito = paste0(round(mean(percent.mito) * 100,2), '%')
mean.percent.mito

# Set mitogene threshold to 20%
#Seurat v2 function, but shows compatibility in Seurat v3
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
#in case the above function does not work simply do: # pbmc$percent.mito <- percent.mito
# We filter out cells that have unique gene counts (nFeature_RNA) over 2,500 or less than
# 200 Note that > and < are used to define a 'gate'.
#-Inf and Inf should be used if you don't want a lower or upper threshold.
pbmc <- subset(x = pbmc,
  subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito >  -Inf & percent.mito < 0.2 )

#### stats
n_genes = dim(pbmc[["RNA"]]@counts)[1]
n_cells = dim(pbmc[["RNA"]]@counts)[2]
# for each cell (column) count the number of expressed/non-0 genes, then median or mean across cells
median_n_genes_per_cell = apply(pbmc@assays[["RNA"]][,], 2, function(x) sum(x != 0 )) %>% median
mean_n_genes_per_cell = apply(pbmc@assays[["RNA"]][,], 2, function(x) sum(x != 0 )) %>% mean

##### outputs
n_genes_raw
n_cells_raw
n_genes
n_cells
mito.genes
percent.mito
median_n_genes_per_cell
mean_n_genes_per_cell

## normalize, scale and cluster
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nCounts_RNA", "percent.mito"))
pbmc <- FindVariableFeatures(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR,
                      x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
pbmc <- RunPCA(object = pbmc,  npcs = 30, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5, algorithm = 1)

n_cells_per_cluster = Idents(pbmc) %>% table
n_clusters = length(n_cells_per_cluster)

##### outputs
n_cells_per_cluster
n_clusters

pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters
pdf(paste0('Seurat_outputs/' ,sample_id,"/TSNEPlot.pdf"), width = 6, height = 6)
print(TSNEPlot(object = pbmc))
dev.off()
print(paste0('Seurat_outputs/' ,sample_id,"/TSNEPlot.pdf"))

out_table = tibble(
  sample_id = sample_id,
  n_genes_raw = n_genes_raw,
  n_cells_raw = n_cells_raw,
  n_genes = n_genes,
  n_cells = n_cells,
  mito.genes = paste0(mito.genes, collapse = ','),
  mean.percent.mito = mean.percent.mito,
  median_n_genes_per_cell = round(median_n_genes_per_cell,1),
  mean_n_genes_per_cell = round(mean_n_genes_per_cell,1),
  n_clusters = n_clusters,
  n_cells_per_cluster = paste0(paste(names(n_cells_per_cluster), n_cells_per_cluster, sep=": "), collapse=', ')
)

metrics_summary = as_tibble(read_csv(paste0(to_process[i],'/../metrics_summary.csv')))
out_table = bind_cols(out_table, metrics_summary )

write_tsv(out_table, paste0('Seurat_outputs/' ,sample_id,"/stats.tsv"))
write.xlsx(out_table, paste0('Seurat_outputs/' ,sample_id,"/stats.xlsx"))

## add diff genes each cluster
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25) %>% arrange(cluster, p_val)
write.xlsx(pbmc.markers,  paste0('Seurat_outputs/' ,sample_id,"/clusters_markers_FindAllMarkers.xlsx"))
