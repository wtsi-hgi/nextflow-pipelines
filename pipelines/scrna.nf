nextflow.preview.dsl=2
params.runtag = 'UkB_scRNA_fase2_4pooled'
// params.read2 = 'discard' 

include iget_cellranger from '../modules/scrna/irods_cellranger.nf' params(run: true, outdir: params.outdir)

workflow {

    // 1.A: from irods cellranger:
    Channel.fromPath("${baseDir}/../../inputs/study_5631_phase2pooled.csv")
	.splitCsv(header: true)
	.map { row -> tuple("${row.study_id}", "${row.run_id}", "${row.samplename}", "${well}", "${row.sanger_sample_id}",
			    "${row.supplier_sample_name}", "${row.pooled}", "${row.cellranger}") }
	.map { a,b,c,d,e,f,g,h -> [c,b,f] }
	.set{ch_samplename_runid_sangersampleid}

    ch_samplename_runid_sangersampleid.view()
    
    iget_cellranger(ch_samplename_runid_sangersampleid)

}

