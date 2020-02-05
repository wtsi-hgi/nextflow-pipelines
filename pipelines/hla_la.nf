nextflow.preview.dsl=2
params.run_hla_la = true

include copy_number_v2 from '../modules/hla_la/hla_la.nf' params(run: true, outdir: params.outdir)

workflow {

// accession_number,merge_location
// EGAN00001672988,/seq/illumina/library_merge/20636762.HXV2.paired308.4af907808c/20636762.HXV2.paired308.4af907808c.cram
// .map{row -> tuple(file(row[0]).getName().replaceAll(~/.noXY.recode.vcf/, "").replaceAll(~/.noXY.cn.vcf/, ""), file(row[0]))}
    Channel.fromPath("${baseDir}/../../inputs/EGAN_to_cram_19534_drop_dups.csv")
	.splitCsv(header: true)
	.take(1)
	.map { row -> tuple(row.accession_number, row.irods_cram)}
	.set{ch_eganid_irodscram}

    if (params.run_hla_la)
	hla_la(ch_eganid_irodscram)
}
