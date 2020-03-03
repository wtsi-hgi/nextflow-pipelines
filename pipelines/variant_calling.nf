nextflow.preview.dsl=2
params.run_strip = true
params.run_vep = true
params.run_concat = true
params.run_vqsr = false

params.vcfs_dir = "/lustre/scratch118/humgen/hgi/projects/ibdx10/variant_calling/joint_calling/vcfs_concatenated"
Channel.fromPath("${params.vcfs_dir}/*.vcf.gz")
	.set{ch_vcfs_gz}
Channel.fromPath("${params.vcfs_dir}/*.vcf.gz.csi")
	.set{ch_vcfs_gz_csi}

include strip_vcf from '../modules/variant_calling/strip_vcf.nf' params(run: true, outdir: params.outdir)
include vep_vcf from '../modules/variant_calling/vep_vcf.nf' params(run: true, outdir: params.outdir)
include concat_vcfs from '../modules/variant_calling/concat_vcfs.nf' params(run: true, outdir: params.outdir)
include vqsr_vcf from '../modules/variant_calling/vqsr_vcf.nf' params(run: true, outdir: params.outdir)

workflow {
    ch_vcfs_gz
	.map{vcf -> tuple(vcf.getSimpleName(),vcf)}
	.combine(
	ch_vcfs_gz_csi
	    .map{csi -> tuple(csi.getSimpleName(),csi)}, by: 0)
	.take(8) // replace with take(-1) to select all inputs
	.set{ch_name_vcf_csi}
//   ch_name_vcf_csi.view()

    if (params.run_strip) {
	strip_vcf(ch_name_vcf_csi)
//	strip_vcf.out.name_vcf_csi.view()
	
	if (params.run_vep) {
	    vep_vcf(strip_vcf.out.name_vcf_csi)
//	vep_vcf.out.name_vcf_csi.view()
	    
	    if (params.run_concat) {
		concat_vcfs(vep_vcf.out.name_vcf_csi.collect())
		
		if (params.run_vqsr) {
		    vqsr_vcf(concat_vcfs.out.concat_vcf)
		}
	    }
	}
    }
}
