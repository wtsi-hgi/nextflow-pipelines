nextflow.preview.dsl=2

// required input parameters:
params.inputs_dir = "/lustre/scratch114/projects/ukbb_scrna/tensorqtl/test_data" // dir where input data is located
params.plink_prefix = "GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered" // .bed,.bim,.fam must be in dir params.inputs_dir
params.output_prefix ="GEUVADIS.445_samples" // prefix of tensorqtl outputs
params.expression_bed = "GEUVADIS.445_samples.expression.bed.gz"
params.covariates_file = "GEUVADIS.445_samples.covariates.txt"

// must also choose input genotype data in either vcf or plink format:
params.convert_vcf_format = false // if true, path to vcf file must be specified:
params.vcf = '' // vcf file name. vcf.gz compressed . File must be in dir params.inputs_dir

// choose which tasks to run:
params.filter_genotype_data  = true // if true, must set params.MAF_threshold:
params.MAF_threshold = 0.1 // 0.1%, minor allele frequency minimum threshold
params.filter_expression_data = true // if true, must set:
params.variance_threshold = 10 // percetange variance distribution minimum cutoff
params.pcent_samples = 40 // required mininum percentage of samples expression level > 0 (run on all genes separately)
params.run_tensorqtl = true
params.convert_plink_to_vcf = true // convert unfiltered .bed,.bim,.fam to vcf format

include convert_vcf_format from '../modules/tensorqtl/convert_vcf_format.nf' params(run:true, outdir: params.outdir)
include convert_plink_format from '../modules/tensorqtl/convert_plink_format.nf' params(run:true, outdir: params.outdir)
include filter_genotype_data from '../modules/tensorqtl/filter_genotype_data.nf' params(run:true, outdir: params.outdir)
include filter_expression_data from '../modules/tensorqtl/filter_expression_data.nf' params(run:true, outdir: params.outdir)
include tensorqtl from '../modules/tensorqtl/tensorqtl.nf' params(run:true, outdir: params.outdir)

workflow {
    
    if(params.convert_vcf_format) {
	convert_vcf_format(Channel.fromPath("${params.inputs_dir}/${params.vcf}"), params.plink_prefix)
	ch_bed_bim_fam = convert_vcf_format.out.bed_bim_fam
    } else {
	Channel.value(tuple(
	    file("${params.inputs_dir}/${params.plink_prefix}.bed"),
	    file("${params.inputs_dir}/${params.plink_prefix}.bim"),
	    file("${params.inputs_dir}/${params.plink_prefix}.fam")))
	    .set{ch_bed_bim_fam}
    }
    
    if(params.convert_vcf_format) {
	convert_plink_format(ch_bed_bim_fam, params.plink_prefix)
    }
    
    if (params.filter_genotype_data) {
	filter_genotype_data(ch_bed_bim_fam, params.plink_prefix, params.MAF_threshold)

	ch_plink_prefix = "${params.plink_prefix}.filtered"
	ch_bed_bim_fam = filter_genotype_data.out.bed_bim_fam
    } else {
	ch_plink_prefix = params.plink_prefix 
    }
    
    if (params.filter_expression_data) {
	filter_expression_data(Channel.fromPath("${params.inputs_dir}/${params.expression_bed}"),
			       params.variance_threshold,
			       params.pcent_samples
	)

	ch_expression_bed = filter_expression_data.out.expression_data
    } else {
	ch_expression_bed = Channel.fromPath("${params.inputs_dir}/${params.expression_bed}")
    }
    
    if (params.run_tensorqtl) {
	tensorqtl(
	    params.output_prefix, // string
	    ch_plink_prefix,  // string
	    ch_bed_bim_fam, // tuple of 3 files, .bed .bim and .fam, filenames prefix must match with params.plink_prefix
	    ch_expression_bed, // file
	    Channel.fromPath("${params.inputs_dir}/${params.covariates_file}"))
    }
}
