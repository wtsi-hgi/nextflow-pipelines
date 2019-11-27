
metrics_summary = as_tibble(read_csv(paste0(to_process[i],'/../metrics_summary.csv')))
out_table = bind_cols(out_table, metrics_summary )

write_tsv(out_table, paste0('Seurat_outputs/' ,sample_id,"/stats.tsv"))
write.xlsx(out_table, paste0('Seurat_outputs/' ,sample_id,"/stats.xlsx"))

## add diff genes each cluster
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25) %>% arrange(cluster, p_val)
write.xlsx(pbmc.markers,  paste0('Seurat_outputs/' ,sample_id,"/clusters_markers_FindAllMarkers.xlsx"))

if (i == 1) {binded = out_table} else {binded = binded %>% bind_rows(out_table)}
print(i)
}

write_tsv(binded, paste0('Seurat_outputs/' ,"stats.tsv"))
write.xlsx(binded, paste0('Seurat_outputs/',"stats.xlsx"))
