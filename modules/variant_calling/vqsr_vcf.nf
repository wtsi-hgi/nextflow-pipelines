params.run = true

process vqsr_vcf {
    memory '6G'
    tag "$name"
    cpus 1
    //conda '/lustre/scratch118/humgen/hgi/projects/ibdx10/variant_calling/joint_calling/ibd_concat_nextflow/bcftools'
    //scratch '/tmp'
    //stageInMode 'copy'
    //stageOutMode 'copy'
    time '700m'
    queue 'normal'
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    publishDir "${params.outdir}/vqsr_vcf/$name/", mode: 'symlink', overwrite: true, pattern: "*.vqsr.vcf.gz"
    publishDir "${params.outdir}/vqsr_vcf/$name/", mode: 'symlink', overwrite: true, pattern: "*.vqsr.vcf.gz.csi"
    
    maxRetries 2

    when:
    params.run
     
    input:
    tuple file(vcf), file(csi)
    
    output:
    tuple file("*.vqsr.vcf.gz"), file("*.vqsr.vcf.gz.csi"), emit: vcf_csi 

    script:
    def simplename = vcf.getSimpleName()
""" 
sleep 30
export DIR=\$PWD 
/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/output_vcf/stripped_vcf/VQSR_indel.sh ${vcf} && \\
/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/output_vcf/stripped_vcf/VQSR_snp.sh  ${vcf} && \\
/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/output_vcf/stripped_vcf/VQSR_indel_apply.sh ${vcf} && \\ 
/lustre/scratch118/humgen/hgi/projects/wtsi_joint_exomes/output_vcf/stripped_vcf/VQSR_snp_apply.sh ${vcf} \\
"""
}

// bgzip ${simplename}.vqsr.vcf
// bcftools index ${simplename}.vqsr.vcf.gz
