nextflow.preview.dsl = 2

params.restrict_vcf_to_chr = true // optional, speed up by e.g. restricting to chromosome 22
params.restrict_vcf_to_region = true // optional, subset input vcf with bed file

params.run_rtg_vcfeval = false
params.run_bcftools_stats = false
params.run_mendelian_errors = false

params.genome = "/lustre/scratch114/projects/interval_wes/giab_wes/gatk_haplotype/hs38DH.fa"

params.restrict_region_bed = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/intervals/Agilent_no_overlaps/S04380110_Padded+1_merged.bed"
params.restrict_chr = "chr22"

params.tag = "hail"
params.vcf = "/lustre/scratch114/projects/interval_wes/giab_wes/hail_dataproc/interval_wes.split_multi.vcf.gz" 
params.vcf_sample = "Sample_Diag-excap51-HG002-EEogPU"
params.vcf_trio = "Sample_Diag-excap51-HG004-EEogPU,Sample_Diag-excap51-HG003-EEogPU,Sample_Diag-excap51-HG002-EEogPU" // in order mother,father,child 

params.rtg_vcf_baseline = "/lustre/scratch114/projects/interval_wes/giab_wes/rtg_vcfeval/refs/HG002_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_CHROM1-22_v3.2.2_highconf_chr.hg38.split.chr22.vcf.gz"
params.rtg_sample = "INTEGRATION"
params.rtg_genome_template = "/lustre/scratch114/projects/interval_wes/giab_wes/gatk_haplotype/hs38DH.fa.sdf"
params.rtg_region = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/intervals/Agilent_no_overlaps/S04380110_Padded+1_merged.bed"

//params.vcf = "/lustre/scratch114/projects/interval_wes/giab_wes/gatk_vqsr/giab_gatk.vcf.gz" 
//params.vcf_sample = "HG002"
//params.vcf_trio = "HG004,HG003,HG002" // in order mother,father,child 

include subset_sample from '../modules/jointcall_eval/subset_sample.nf' params(run: true, outdir: params.outdir)
include subset_sample_and_region from '../modules/jointcall_eval/subset_sample_and_region.nf' params(run: true, outdir: params.outdir)
include subset_trio from '../modules/jointcall_eval/subset_trio.nf' params(run: true, outdir: params.outdir)
include subset_trio_and_region from '../modules/jointcall_eval/subset_trio_and_region.nf' params(run: true, outdir: params.outdir)
include subset_chr from '../modules/jointcall_eval/subset_chr.nf' params(run: true, outdir: params.outdir)

include bcftools_stats from '../modules/jointcall_eval/bcftools_stats.nf' params(run: true, outdir: params.outdir)
include rtg_vcfeval from '../modules/jointcall_eval/rtg_vcfeval.nf' params(run: true, outdir: params.outdir)
include mendelian_errors from '../modules/jointcall_eval/mendelian_errors.nf' params(run: true, outdir: params.outdir)

workflow {

    if (params.restrict_vcf_to_chr) {
	subset_chr(params.vcf, params.restrict_chr)
	ch_subset_vcf = subset_chr.out.vcf
    } else {
	ch_subset_vcf = params.vcf
    }
    
    if (params.restrict_vcf_to_region) {
	subset_sample_and_region(ch_subset_vcf, params.sample, params.restrict_region_bed)
	ch_subset_vcf_tbi = subset_sample_and_region.out.vcf_tbi
    } else {
	subset_sample(ch_subset_vcf, params.sample)
	ch_subset_vcf_tbi = subset_sample.out.vcf_tbi
    }
    
    if (params.run_bcftools_stats) {
	bcftools_stats(ch_subset_vcf_tbi) }

    if (params.run_rtg_vcfeval) {
	rtg_vcfeval(ch_subset_vcf_tbi,
		    params.sample,
		    params.rtg_sample,
		    params.rtg_vcf_baseline,
		    params.rtg_genome_template,
		    params.rtg_region) }
    
    if (params.run_mendelian_errors) {
	if (params.restrict_vcf_to_region) {
	    subset_trio_and_region(ch_subset_vcf, params.vcf_trio, params.restrict_region_bed)
	    ch_subset_trio_vcf_tbi = subset_trio_and_region.out.vcf_tbi
	} else {
	    subset_trio(ch_subset_vcf, params.vcf_trio)
	    ch_subset_trio_vcf_tbi = subset_trio.out.vcf_tbi
	}
	mendelian_errors(ch_subset_trio_vcf_tbi) }

}
