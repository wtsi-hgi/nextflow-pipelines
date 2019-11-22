nextflow.preview.dsl=2
params.runtag = 'UkB_scRNA_fase2_4pooled'
params.read2 = 'discard' // used by count_crispr_reads

include iget_cellranger from '../modules/crispr/irods_cellranger.nf' params(run: true, outdir: params.outdir)
