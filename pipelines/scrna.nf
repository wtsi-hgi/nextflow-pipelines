nextflow.preview.dsl=2
// branch unassign
params.runtag = 'UkB_scRNA_fase2_4pooled'
params.run_cellsnp = true
params.run_vireo = true
params.run_seurat = false
params.run_seurat_on_raw = false // run seurat on raw_feature_bc_matrix (in addition to filtered_feature_bc_matrix)


params.cellsnp_vcf_candidate_snps = "$baseDir/../assets/scrna/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
Channel.fromPath(params.cellsnp_vcf_candidate_snps)
    .ifEmpty { exit 1, "cellsnp_vcf_candidate_snps missing: ${params.cellsnp_vcf_candidate_snps}" }
    .set {ch_cellsnp_vcf_candidate_snps}


include iget_cellranger from '../modules/scrna/irods_cellranger.nf' params(run: true, outdir: params.outdir)
include cellsnp from '../modules/scrna/cellsnp.nf' params(run: true, outdir: params.outdir)
include vireo from '../modules/scrna/vireo.nf' params(run: true, outdir: params.outdir)
include split_vireo_barcodes from '../modules/scrna/split_vireo_barcodes.nf' params(run: true, outdir: params.outdir)
include seurat from '../modules/scrna/seurat.nf' params(run: true, outdir: params.outdir)
// include seurat_rbindmetrics from '../modules/scrna/seurat_rbindmetrics.nf' params(run: true, outdir: params.outdir)


workflow run_seurat {
    get: cellranger_data_raw
    get: cellranger_data_filtered
    get: cellranger_data_metrics_summary
    
    main:
    if (params.run_seurat_on_raw)
	input_seurat = cellranger_data_filtered.map{samplename,dir -> [samplename,dir,"filtered"]} // filtered_feature_bc_matrix
	.mix(cellranger_data_raw.map{samplename,dir -> [samplename,dir,"raw"]}) // raw_feature_bc_matrix
	.combine(cellranger_data_metrics_summary, by: 0) // add metrics_summary.csv file of each sample
    else
	input_seurat = cellranger_data_filtered.map{samplename,dir -> [samplename,dir,"filtered"]} // filtered_feature_bc_matrix
	    .combine(cellranger_data_metrics_summary, by: 0) // add metrics_summary.csv file of each sample

    seurat(input_seurat)
    
    emit: seurat.out.tsneplot_pdf
    emit: seurat.out.stats_xslx
    emit: seurat.out.diffexp_xlsx
    emit: seurat.out.seurat_rdata
}


workflow {

    // 1.A: from irods cellranger:
    Channel.fromPath("${baseDir}/../../inputs/study_5631_phase2pooled.csv")
	.splitCsv(header: true)
	.map { row -> tuple("${row.study_id}", "${row.run_id}", "${row.samplename}", "${row.well}", "${row.sanger_sample_id}",
			    "${row.supplier_sample_name}", "${row.pooled}","${row.n_pooled}", "${row.cellranger}") }
	.map { a,b,c,d,e,f,g,h, i -> [c,b,f,h] }
	.set{ch_samplename_runid_sangersampleid_npooled}

    ch_samplename_runid_sangersampleid = ch_samplename_runid_sangersampleid_npooled
	.map { a,b,c,d -> [a,b,c] }

    ch_samplename_npooled = ch_samplename_runid_sangersampleid_npooled
	.map { a,b,c,d -> [a,d] }

    iget_cellranger(ch_samplename_runid_sangersampleid)
    
    if (params.run_cellsnp)
	cellsnp(iget_cellranger.out.cellranger_sample_bam_barcodes, ch_cellsnp_vcf_candidate_snps.collect())
    
    if (params.run_vireo) {
	vireo(cellsnp.out.cellsnp_output_dir.combine(ch_samplename_npooled, by: 0))

        split_vireo_barcodes(vireo.out.vireo_output_dir.
			     combine(ch_samplename_npooled, by: 0).
			     combine(iget_cellranger.out.cellranger_filtered, by: 0))
	
        split_vireo_barcodes.out.cellranger_deconv_dirs
	    .transpose()
	    .map { samplename,deconv_dir -> tuple(deconv_dir.getName().replaceAll(~/cellranger_deconv_/, ""),deconv_dir) }
            .set{ch_cellranger_filtered_deconv}
	
	//ch_cellranger_filtered_deconv.view()
	//iget_cellranger.out.cellranger_raw.view()
	//iget_cellranger.out.cellranger_metrics_summary.view()
	
	
	run_seurat(iget_cellranger.out.cellranger_raw,
		   ch_cellranger_filtered_deconv,
		   iget_cellranger.out.cellranger_metrics_summary
	     	       .map{samplename, metrics_file -> tuple(["${samplename}_0","${samplename}_1","${samplename}_2","${samplename}_3","${samplename}_4"],
							  metrics_file)}
		       .transpose()
	)
    }
    
    //if (params.run_seurat)
//	run_seurat(iget_cellranger.out[2],ch_cellranger_filtered_deconv, iget_cellranger.out[4])
	//run_seurat(iget_cellranger.out[2],iget_cellranger.out[3],iget_cellranger.out[4])

}
