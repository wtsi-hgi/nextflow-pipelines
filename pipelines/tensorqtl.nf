nextflow.preview.dsl=2

params.inputs_dir = "/lustre/scratch114/projects/ukbb_scrna/tensorqtl/test_data"

// choose input genotype data in either vcf or plink format:
params.convert_vcf_format = true // if true vcf must be specified:
params.vcf = '' // else use:
params.plink_prefix_path = "${params.inputs_dir}/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered"

params.prefix ="GEUVADIS.445_samples"
params.expression_bed = file("${params.inputs_dir}/GEUVADIS.445_samples.expression.bed.gz"),
params.covariated_file = file("${params.inputs_dir}/GEUVADIS.445_samples.covariates.txt")

params.filter_genotypes_data  = true
params.filter_expression_data = true
params.run_tensorqtl = true


include convert_vcf_format from '../modules/tensorqtl/convert_vcf_format.nf' params(run:true, outdir: params.outdir)
include filter_genotypes_data from '../modules/tensorqtl/filter_genotypes_data.nf' params(run:true, outdir: params.outdir)
include filter_expression_data from '../modules/tensorqtl/filter_expression_data.nf' params(run:true, outdir: params.outdir)
include tensorqtl from '../modules/tensorqtl/tensorqtl.nf' params(run:true, outdir: params.outdir)

workflow {
    
    if(params.convert_vcf_format) {
	convert_vcf_format(params.vcf, params.prefix)
	params.plink_prefix_path = convert_vcf_format.out.plink_prefix_path
	params.bed = convert_vcf_format.out.bed
	params.bim = convert_vcf_format.out.bim
	params.fam = convert_vcf_format.out.fam
    } else {
	params.bed = "${params.plink_prefix_path}.bed"
	params.bim = "${params.plink_prefix_path}.bim"
	params.fam = "${params.plink_prefix_path}.fam"
    }
    
    if (params.run_tensorqtl) {
	
	Channel.value(tuple(
	    params.prefix,
	    params.plink_prefix_path, 
	    file(params.bed),file(params.bim),file(params.fam),
	    file(params.expression_bed),
	    file(params.covariated_file)))
	    .set{to_tensorqtl}
	
	to_tensorqtl.view()
	tensorqtl(to_tensorqtl)
    }
}
