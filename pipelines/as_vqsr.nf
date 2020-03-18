nextflow.preview.dsl=2
params.run_vqsr = true
params.apply_vqsr = true


params.vcfs_dir = "/lustre/scratch118/humgen/hgi/projects/interval_wes/joint_calls/output_vcf/stripped_vcf"
Channel.fromPath("${params.vcfs_dir}/interval_wes_stripped.vcf.gz")
	.set{ch_vcfs_gz}
Channel.fromPath("${params.vcfs_dir}/interval_wes_stripped.vcf.gz.tbi")
	.set{ch_vcfs_gz_tbi}

include vqsr_vcf from '../modules/variant_calling/as_vqsr_vcf_1st_step.nf' params(run: true, outdir: params.outdir)
include vqsr_vcf_apply from '../modules/variant_calling/as_vqsr_vcf_2nd_step.nf' params(run: true, outdir: params.outdir)

workflow {
    ch_vcfs_gz
	.map{vcf -> tuple(vcf.getSimpleName(),vcf)}
	.combine(
	ch_vcfs_gz_tbi
	    .map{tbi -> tuple(tbi.getSimpleName(),tbi)}, by: 0)
	.take(2) // replace with take(-1) to select all inputs
	.set{ch_name_vcf_tbi}

    



   // ch_name_vcf_tbi.view()
   if (params.run_vqsr) {
     vqsr_vcf(ch_name_vcf_tbi.map{name,vcf,tbi -> tuple(file(name),file(vcf),file(tbi))})
    //vqsr_vcf.out.name_vcf_csi.view()
    vqsr_vcf.out.recal.view()
    vqsr_vcf.out.tranches.view()
   }
    ch_name_vcf_tbi.combine(
    vqsr_vcf.out.recal).combine(
    vqsr_vcf.out.tranches).set{vqsr_step1_out}

  vqsr_step1_out.view()
    
   if (params.apply_vqsr) {
     vqsr_vcf_apply(vqsr_step1_out.map{name,vcf,tbi, snp_recal, indel_recal, snp_tranch, indel_tranch  -> tuple(file(name),file(vcf),file(tbi), file(snp_recal), file(indel_recal),file(snp_tranch), file(indel_tranch))})
   }

//vqsr_vcf(concat_vcfs.out.concat_vcf.map{vcf,csi,tbi -> tuple(vcf,tbi)})
  //  vqsr_vcf(ch_name_vcf_tbi)

      println "End of workflow."
}