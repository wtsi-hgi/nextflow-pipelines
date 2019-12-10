nextflow.preview.dsl=2
params.runtag = 'copy_number_v2'
params.run_copy_number = true

include copy_number_v2 from '../modules/sv-pipeline/copy_number.nf' params(run: true, outdir: params.outdir)


workflow {

    Channel.fromPath("${baseDir}/../../inputs/samples_copy_number_v2.csv")
	.splitCsv(header: true)
	.map { row -> tuple(row.samplename, row.egan_id, row.hist_root_file, row.samplename_gt_vcf)}
	.set{ch_copy_number_v2}

    
    if (params.run_copy_number)
	copy_number_v2(ch_copy_number_v2)

}
