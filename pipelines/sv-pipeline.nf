nextflow.preview.dsl=2
params.runtag = 'remove_XY'
params.run_vcf_remove_chrXY = false
params.run_copy_number = false

include copy_number_v2 from '../modules/sv-pipeline/copy_number.nf' params(run: true, outdir: params.outdir)
include vcf_remove_chrXY from '../modules/sv-pipeline/vcf_remove_chrXY.nf' params(run: true, outdir: params.outdir)


workflow {

    //Channel.fromPath("${baseDir}/../../inputs/copy_number_input_v2.csv")
    Channel.fromPath("${baseDir}/../../inputs/copy_number_input_v3.csv")
	.splitCsv(header: true)
	.map { row -> tuple(row.samplename, row.EGAN_id, file(row.root_file), file(row.gt_vcf))}
	.set{ch_copy_number_v2}

    if (params.run_copy_number)
	copy_number_v2(ch_copy_number_v2)



    Channel.fromPath("/home/ubuntu/data2/results/copy_number/*.cn.vcf")
	.take(4)
	.map{gt_vcf -> tuple(gt_vcf.getName().replaceAll(~/.cn.vcf/, ""), gt_vcf)}
	.set{ch_vcf_remove_chrXY}
    
    if (params.run_vcf_remove_chrXY)
	vcf_remove_chrXY(ch_vcf_remove_chrXY)


}
