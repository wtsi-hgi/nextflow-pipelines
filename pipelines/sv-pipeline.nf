nextflow.preview.dsl=2
params.runtag = 'split_chr'
params.run_vcf_remove_chrXY = false
params.run_copy_number = false
params.run_split_all_chr = true

include copy_number_v2 from '../modules/sv-pipeline/copy_number.nf' params(run: true, outdir: params.outdir)
include vcf_remove_chrXY from '../modules/sv-pipeline/vcf_remove_chrXY.nf' params(run: true, outdir: params.outdir)
include vcf_split_all_chr from '../modules/sv-pipeline/vcf_split_all_chr.nf' params(run: true, outdir: params.outdir)


workflow {

    //Channel.fromPath("${baseDir}/../../inputs/copy_number_input_v2.csv")
    Channel.fromPath("${baseDir}/../../inputs/copy_number_input_v3.csv")
	.splitCsv(header: true)
	.map { row -> tuple(row.samplename, row.EGAN_id, file(row.root_file), file(row.gt_vcf))}
	.set{ch_copy_number_v2}

    if (params.run_copy_number)
	copy_number_v2(ch_copy_number_v2)

    if (params.run_split_all_chr) {
	Channel.fromPath("/home/ubuntu/data2/to_split_by_chr.txt") // lines
	.splitCsv(header: false)
	.take(2)
        .map{row -> tuple(file(row[0]).getName().replaceAll(~/.noXY.recode.vcf/, "").replaceAll(~/.noXY.cn.vcf/, ""), file(row[0]))}
	    .set{ch_to_split}
	
        vcf_split_all_chr(ch_to_split)
    }
    
    if (params.run_vcf_remove_chrXY) {
	Channel.fromPath("/home/ubuntu/data2/to_rm_XY.txt") // 11867 lines
	.splitCsv(header: false)
	.take(-1)
        .map{row -> tuple(file(row[0]).getName().replaceAll(~/.cn.vcf/, ""), file(row[0]))}
	.set{ch_vcf_remove_chrXY}
    
        vcf_remove_chrXY(ch_vcf_remove_chrXY)
    }
}
