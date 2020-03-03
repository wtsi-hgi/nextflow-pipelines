params.run = true

process vep_vcf {
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
    publishDir "${params.outdir}/vep_vcf/$name/", mode: 'symlink', overwrite: true, pattern: "*.sorted.vcf.gz"
    
    maxRetries 2

    when:
    params.run
     
    input:
    tuple val(name), file(vcf), file(csi)
    
    output:
    tuple val(name), file("*.stripped.vcf.gz"), file("*.stripped.vcf.gz.csi"), emit: name_vcf_csi 

    script:
    def simplename = vcf.getSimpleName()
""" 
bcftools view -G -o ${simplename}.stripped.vcf.gz -O z $vcf
bcftools index ${simplename}.stripped.vcf.gz
"""
}

