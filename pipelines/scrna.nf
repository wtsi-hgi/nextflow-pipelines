nextflow.preview.dsl=2
params.runtag = 'UkB_scRNA_fase2_4pooled'
// params.read2 = 'discard' 

params.cellsnp_vcf_candidate_snps = "$baseDir/../assets/scrna/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
Channel.fromPath(params.cellsnp_vcf_candidate_snps)
    .ifEmpty { exit 1, "cellsnp_vcf_candidate_snps missing: ${params.cellsnp_vcf_candidate_snps}" }
    .set {ch_cellsnp_vcf_candidate_snps}


include iget_cellranger from '../modules/scrna/irods_cellranger.nf' params(run: true, outdir: params.outdir)
include cellsnp from '../modules/scrna/cellsnp.nf' params(run: true, outdir: params.outdir)
// include vireo from '../modules/scrna/vireo.nf' params(run: true, outdir: params.outdir)

workflow {

    // 1.A: from irods cellranger:
    Channel.fromPath("${baseDir}/../../inputs/study_5631_phase2pooled.csv")
	.splitCsv(header: true)
	.map { row -> tuple("${row.study_id}", "${row.run_id}", "${row.samplename}", "${row.well}", "${row.sanger_sample_id}",
			    "${row.supplier_sample_name}", "${row.pooled}", "${row.cellranger}") }
	.map { a,b,c,d,e,f,g,h -> [c,b,f] }
	.set{ch_samplename_runid_sangersampleid}

    iget_cellranger(ch_samplename_runid_sangersampleid)
    
    // cellsnp(iget_cellranger.out[1], ch_cellsnp_vcf_candidate_snps.collect())
    
    //vireo(cellsnp.out[0])

}
