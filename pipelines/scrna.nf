nextflow.preview.dsl=2
params.runtag = 'cellranger_fetch'
params.current_date = 'current_date'
params.run_velocyto = false
params.run_cellsnp = false
params.run_vireo = false
params.run_seurat = false
params.run_seurat_on_raw = false // run seurat on raw_feature_bc_matrix (in addition to filtered_feature_bc_matrix)
params.run_solo= false
params.run_push_to_git = false
params.on_complete_uncache_irods_search = true
params.on_complete_remove_workdir_failed_tasks = true
params.min_reads = 500   // used by crams_to_fastq_gz

// params.cellsnp_vcf_candidate_snps = "$baseDir/../assets/scrna/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
// Channel.fromPath(params.cellsnp_vcf_candidate_snps)
//    .ifEmpty { exit 1, "cellsnp_vcf_candidate_snps missing: ${params.cellsnp_vcf_candidate_snps}" }
//    .set {ch_cellsnp_vcf_candidate_snps}
Channel.fromPath("${baseDir}/../../inputs/solo_params.json")
    .set{ch_solo_params_json}

// copied from /nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf 
Channel.fromPath("${baseDir}/../../inputs/refdata-cellranger-GRCh38-3.0.0_genes.gtf")
    .set{ch_cellranger_gtf}

include iget_cram_cellranger from '../modules/scrna/iget_cram_cellranger_fromsample.nf' params(run: true, outdir: params.outdir)
include baton_imeta from '../modules/scrna/baton_imeta.nf' params(run: true, outdir: params.outdir)
include push_to_git from '../modules/scrna/push_to_git.nf' params(run: true, outdir: params.outdir)
include info_stats from '../modules/scrna/info_stats.nf' params(run: true, outdir: params.outdir)

include velocyto from '../modules/scrna/velocyto.nf' params(run: true, outdir: params.outdir)
include cellsnp from '../modules/scrna/cellsnp.nf' params(run: true, outdir: params.outdir)
include vireo from '../modules/scrna/vireo.nf' params(run: true, outdir: params.outdir)
include seurat from '../modules/scrna/seurat.nf' params(run: true, outdir: params.outdir)
// include seurat_rbindmetrics from '../modules/scrna/seurat_rbindmetrics.nf' params(run: true, outdir: params.outdir)

include crams_to_fastq_gz from '../modules/scrna/crams_to_fastq_anyflag.nf' params(run:true, outdir: params.outdir,    min_reads: params.min_reads)
include fastqc from '../modules/scrna/fastqc.nf' params(run: true, outdir: params.outdir)
include multiqc from '../modules/scrna/multiqc.nf' params(run: true, outdir: params.outdir,   runtag : params.runtag)

include solo from '../modules/solo/solo.nf' params(run: true, outdir: params.outdir)

workflow run_seurat {
    take: cellranger_data_raw
    take: cellranger_data_filtered
    take: cellranger_data_metrics_summary
    
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

    Channel.fromPath("${baseDir}/../../inputs/samples_google_spreadsheet.tsv")
	.set{ch_google_spreadsheet_tsv}

    ch_google_spreadsheet_tsv
	.splitCsv(header: true, sep: '\t')
        .take(-1)
	.map { row -> tuple(row.sanger_sample_id,row.biopsy_type,row.disease_status)}
        .filter { it[0] != '' }
        .filter { it[1] != '' }
        .filter { it[2] != '' }
	// .filter { it[0] =~ /5892STDY9430443|5892STDY9430444|scrnacdb9430539/ } 
	//.filter { it[0] =~ /5892STDY8356878|5892STDY8357455|5892STDY8357550|5892STDY8357646|5892STDY8357647|5892STDY8644304|OTARscRNA8355918|Crohns_Disease_Collection_Study8727199|OTARscRNA8356110|Crohns_Disease_Collection_Study8727393|5892STDY8644400|5892STDY8644401|OTARscRNA8966185|OTARscRNA8966193|OTARscRNA9294490|OTARscRNA9294497|OTARscRNA9294498|OTARscRNA9294502|OTARscRNA9294504/ } 
	.set{ch_sanger_sample_id_biopsy}

    //// get cellranger data
    iget_cram_cellranger(ch_sanger_sample_id_biopsy)
    iget_cram_cellranger.out.sync_status
	.map{a,b -> b}
	.set{ch_sync_status}

    // solo doublets
    if (params.run_solo) {
	iget_cram_cellranger.out.cellranger_filtered
	    .map{samplename,b -> tuple(samplename, file("${b}.h5"))}
	    .set{to_solo}
	
	solo(to_solo, ch_solo_params_json.collect())
    }
    //// get imeta data
    baton_imeta(ch_sanger_sample_id_biopsy.map{a,b,c->a})

    ch_imeta_info_tsv = baton_imeta.out.imeta.map{a,b->b}
	.collectFile(name: 'samples_imeta_spreadsheet.tsv', newLine: true, keepHeader: true,
    		     storeDir: "${params.outdir}/", sort: true)

    ch_sync_status_tsv = ch_sync_status
	.collectFile(name: 'sync_status.tsv', newLine: true, keepHeader: true,
    		     storeDir: "${params.outdir}/", sort: true)

    file("${baseDir}/../../results/metrics_summary/").mkdirs()
    ch_metrics_summary = iget_cram_cellranger.out.cellranger_metrics_summary
	.map{a,b -> tuple(a,b.copyTo("${baseDir}/../../results/metrics_summary/${a}.metrics_summary.csv"))}
	.map{a,b -> file("${baseDir}/../../results/metrics_summary/${a}.metrics_summary.csv")}
    //ch_metrics_summary.view()

//    // bams to fastq to fastqc to multiQC
//    iget_cram_cellranger.out.cellranger_sample_bam_barcodes
//	.map{sanger_sample_id,bam,bai,raw_barcode -> tuple(sanger_sample_id,'batch',bam)}
//	.set{ch_bams}
//    crams_to_fastq_gz(ch_bams)
//    crams_to_fastq_gz.out[0].set{ch_samplename_batch_fastqs}
//    fastqc(ch_samplename_batch_fastqs
//	   .map{ samplename, batch, fastq -> tuple( samplename, fastq ) })
//    multiqc(fastqc.out.collect())
    
    //// join and plot infos 
    info_stats(ch_google_spreadsheet_tsv, ch_imeta_info_tsv, ch_sync_status_tsv, ch_metrics_summary.collect())

    iget_cram_cellranger.out.iget_command
	.map { samplename,command_file -> command_file }
	.splitCsv(header: false, sep: ' ')
	.map { get,kr,irods_loc,samplename -> "${samplename},${irods_loc}"}
	.collectFile(name: 'cellranger_irods_locations.csv', newLine: true , seed: "sanger_sample_id,cellranger_irods_path",
		     storeDir: "${params.outdir}/", sort: true)
	.set{ch_cellranger_irods_locations_csv}

    if (params.run_velocyto)
	velocyto(
	iget_cram_cellranger.out.cellranger_full_dir
	   // .filter { it =~ /5892STDY8356878|5892STDY8357455|5892STDY8357550|5892STDY8357646|5892STDY8357647|5892STDY8644304|OTARscRNA8355918|Crohns_Disease_Collection_Study8727199|OTARscRNA8356110|Crohns_Disease_Collection_Study8727393|5892STDY8644400|5892STDY8644401|OTARscRNA8966185|OTARscRNA8966193|OTARscRNA9294490|OTARscRNA9294497|OTARscRNA9294498|OTARscRNA9294502|OTARscRNA9294504/ } 
		 , ch_cellranger_gtf.collect())

    //// sync info tables and pdf to git
	push_to_git(info_stats.out[0]
		    .mix(info_stats.out[1])
		    .mix(info_stats.out[2])
		    .mix(info_stats.out[4])
		    .mix(info_stats.out[5])
		    .mix(info_stats.out[6])
		    .mix(info_stats.out[7])
		    .mix(info_stats.out[8])
		    .mix(ch_cellranger_irods_locations_csv)
		    .mix(info_stats.out[3]), params.current_date)
	push_to_git.out.stdout.view()
//

//
//    ch_samplename_runid_sangersampleid = ch_samplename_runid_sangersampleid_npooled
//	.map { a,b,c,d -> [a,b,c] }
//
//    ch_samplename_npooled = ch_samplename_runid_sangersampleid_npooled
//	.map { a,b,c,d -> [a,d] }
//
//    
//    if (params.run_cellsnp)
//	cellsnp(iget_cellranger.out.cellranger_sample_bam_barcodes, ch_cellsnp_vcf_candidate_snps.collect())
//    
//    if (params.run_vireo)
//	vireo(cellsnp.out.cellsnp_output_dir.combine(ch_samplename_npooled, by: 0))
//
//    if (params.run_seurat)
//	run_seurat(iget_cellranger.out[2],iget_cellranger.out[3],iget_cellranger.out[4])
    // list work dirs to remove (because they are Irods searches, so need to always rerun on each NF run):
    // these are removed on workflow.onComplete if (params.on_complete_uncache_irods_search), see below.
    
    //baton_imeta.out.work_dir_to_remove // set always remove
    push_to_git.out.work_dir_to_remove // set always remove
	.mix(info_stats.out.work_dir_to_remove) // set always remove
	.mix(iget_cram_cellranger.out.work_dir_to_remove) // set to remove only if cellranger data not found
	.filter { it != "dont_remove" }
	.collectFile(name: 'irods_work_dirs_to_remove.csv', newLine: true, sort: true,
		     storeDir:params.outdir)
    
}
workflow.onError {
    log.info "Pipeline execution stopped with the following message: ${workflow.errorMessage}" }

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Command line: $workflow.commandLine"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    
    if (params.on_complete_uncache_irods_search) {
	log.info "You have selected \"on_complete_uncache_irods_search = true\"; will therefore attempt to remove Irods work dirs to forcefully uncache them even if successful."
	if (! file("${params.outdir}/irods_work_dirs_to_remove.csv").isEmpty()) {
	    log.info "file ${params.outdir}/irods_work_dirs_to_remove.csv exists and not empty ..."
	    file("${params.outdir}/irods_work_dirs_to_remove.csv")
		.eachLine {  work_dir ->
		if (file(work_dir).isDirectory()) {
		    log.info "removing work dir $work_dir ..."
		    file(work_dir).deleteDir()   
		} } } }
    
    if (params.on_complete_remove_workdir_failed_tasks) {
	log.info "You have selected \"on_complete_remove_workdir_failed_tasks = true\"; will therefore remove work dirs of all tasks that failed (.exitcode file not 0)."
	// work dir and other paths are hardcoded here ... :
	def proc = "bash ./nextflow-pipelines/bin/scrna/del_work_dirs_failed.sh ${workDir}".execute()
	def b = new StringBuffer()
	proc.consumeProcessErrorStream(b)
	log.info proc.text
	log.info b.toString() }
}
