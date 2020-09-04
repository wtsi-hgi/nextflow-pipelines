nextflow.preview.dsl = 2

params.restrict_vcf_to_chr = true // optional, speed up by e.g. restricting to chromosome 22
params.restrict_vcf_to_region = true // optional, subset input vcf with bed file

params.run_bcftools_stats = true
params.run_vcfstats = true
params.run_rtg_vcfeval = true
params.run_mendelian_errors = true

params.genome = "/lustre/scratch114/projects/interval_wes/giab_wes/gatk_haplotype/hs38DH.fa"

params.restrict_chr = "chr22"
params.restrict_region_bed = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/intervals/Agilent_no_overlaps/S04380110_Padded+1_merged.bed"


/////////////////////////////////////////////
params.tag = "hail_HG002"
// need index .tbi of vcf as well:
params.vcf = "/lustre/scratch114/projects/interval_wes/giab_wes/hail_dataproc/interval_wes.split_multi.vcf.gz" 
params.vcf_sample = "Sample_Diag-excap51-HG002-EEogPU"
params.vcf_trio = "Sample_Diag-excap51-HG004-EEogPU,Sample_Diag-excap51-HG003-EEogPU,Sample_Diag-excap51-HG002-EEogPU" // in order mother,father,child 

// need index .tbi of vcf as well:
params.rtg_vcf_baseline = "/lustre/scratch114/projects/interval_wes/giab_wes/rtg_vcfeval/refs/HG002_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_CHROM1-22_v3.2.2_highconf_chr.hg38.split.chr22.vcf.gz"
params.rtg_sample = "INTEGRATION"
params.rtg_genome_template = "/lustre/scratch114/projects/interval_wes/giab_wes/gatk_haplotype/hs38DH.fa.sdf"
params.rtg_region = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/intervals/Agilent_no_overlaps/S04380110_Padded+1_merged.bed"
/////////////////////////////////////////////


// Illumina Dragen
/////////////////////////////////////////////
//params.tag = "illumina_HG002"
//params.vcf = "/lustre/scratch115/projects/interval_wes/illumina_dragen/openstack+GIAB/interval_wes_4070_QC_passed_diploid.vcf.gz"
//params.vcf_sample = "Sample_Diag-excap51-HG002-EEogPU"
//params.vcf_trio = "Sample_Diag-excap51-HG004-EEogPU,Sample_Diag-excap51-HG003-EEogPU,Sample_Diag-excap51-HG002-EEogPU" // in order mother,father,child 
/////////////////////////////////////////////

// GATK vcf
/////////////////////////////////////////////
//params.tag = "gatk_HG002"
//params.vcf = "/lustre/scratch114/projects/interval_wes/giab_wes/gatk_vqsr/giab_gatk.vcf.gz" 
//params.vcf_sample = "HG002"
//params.vcf_trio = "HG004,HG003,HG002" // in order mother,father,child 
/////////////////////////////////////////////

include subset_chr from '../modules/jointcall_eval/subset_chr.nf' params(run: true, outdir: params.outdir)

include subset_sample from '../modules/jointcall_eval/subset_sample.nf' params(run: true, outdir: params.outdir)
include subset_sample_and_region from '../modules/jointcall_eval/subset_sample_and_region.nf' params(run: true, outdir: params.outdir)
include subset_trio from '../modules/jointcall_eval/subset_trio.nf' params(run: true, outdir: params.outdir)
include subset_trio_and_region from '../modules/jointcall_eval/subset_trio_and_region.nf' params(run: true, outdir: params.outdir)

include bcftools_stats from '../modules/jointcall_eval/bcftools_stats.nf' params(run: true, outdir: params.outdir)
include vcfstats from '../modules/jointcall_eval/vcfstats.nf' params(run: true, outdir: params.outdir)
include rtg_vcfeval from '../modules/jointcall_eval/rtg_vcfeval.nf' params(run: true, outdir: params.outdir)
include mendelian_errors from '../modules/jointcall_eval/mendelian_errors.nf' params(run: true, outdir: params.outdir)

workflow {
    ch_genome = Channel.fromPath(params.genome)
    
    ch_input_vcf_tbi = Channel.from(params.vcf)
	.map{vcf -> tuple(file("${vcf}"), file("${vcf}.tbi"))}
    ch_input_vcf_tbi.view()
    
    if (params.restrict_vcf_to_chr) {
        ch_input_vcf_tbi.view()
	subset_chr(ch_input_vcf_tbi, params.restrict_chr)
	ch_subset_vcf_tbi = subset_chr.out.vcf_tbi
    } else {
	ch_subset_vcf_tbi = ch_input_vcf_tbi
    }
    ch_subset_vcf_tbi.view()



    
    if (params.restrict_vcf_to_region) {
	ch_restrict_region_bed = Channel.fromPath(params.restrict_region_bed)
	subset_sample_and_region(ch_subset_vcf_tbi, params.vcf_sample, ch_restrict_region_bed)
	ch_subset_sample_vcf_tbi = subset_sample_and_region.out.vcf_tbi
    } else {
	subset_sample(ch_subset_vcf, params.vcf_sample)
	ch_subset_sample_vcf_tbi = subset_sample.out.vcf_tbi
    }
    ch_subset_sample_vcf_tbi.view()   
    
    if (params.run_bcftools_stats) {
	bcftools_stats(ch_subset_sample_vcf_tbi, params.tag) }

    if (params.run_vcfstats) {
	vcfstats(ch_subset_sample_vcf_tbi, params.tag) }
    
    if (params.run_rtg_vcfeval) {

	ch_rtg_region = Channel.fromPath(params.rtg_region)
	ch_rtg_genome_template = Channel.fromPath(params.rtg_genome_template)
	
	ch_rtg_vcf_baseline_tbi = Channel.from(params.rtg_vcf_baseline)
	    .map{rtg_vcf -> tuple(file("${rtg_vcf}"), file("${rtg_vcf}.tbi"))}

	rtg_vcfeval(ch_subset_sample_vcf_tbi,
		    params.vcf_sample,
		    params.rtg_sample,
		    ch_rtg_vcf_baseline_tbi,
		    ch_rtg_genome_template,
    		    ch_rtg_region,
                    params.tag)
    }
    
    if (params.run_mendelian_errors) {
	if (params.restrict_vcf_to_region) {
	    subset_trio_and_region(ch_subset_vcf_tbi, params.vcf_trio, ch_restrict_region_bed)
	    ch_subset_trio_vcf = subset_trio_and_region.out.vcf_tbi
	} else {
	    subset_trio(ch_subset_vcf_tbi, params.vcf_trio)
	    ch_subset_trio_vcf = subset_trio.out.vcf_tbi
	}
	mendelian_errors(ch_subset_trio_vcf, params.vcf_trio, params.tag) }

}
