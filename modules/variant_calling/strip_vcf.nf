params.run = true

process strip_vcf {
    memory '6G'
    tag "$name"
    cpus 1
    conda '/lustre/scratch118/humgen/hgi/projects/ibdx10/variant_calling/joint_calling/ibd_concat_nextflow/bcftools'
    //scratch '/tmp'
    //stageInMode 'copy'
    //stageOutMode 'copy'
    time '700m'
    queue 'normal'
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    publishDir "${params.outdir}/strip_vcf/$name/", mode: 'symlink', overwrite: true, pattern: "*.stripped.vcf.gz"
    publishDir "${params.outdir}/strip_vcf/$name/", mode: 'symlink', overwrite: true, pattern: "*.stripped.vcf.gz.csi"
    publishDir "${params.outdir}/strip_vcf_nvariants/", mode: 'symlink', overwrite: true, pattern: "${simplename}.zero_variants.txt"
    
    maxRetries 2

    when:
    params.run
     
    input:
    tuple val(name), file(vcf), file(csi)
    
    output:
    tuple val(name), file("*.stripped.vcf.gz"), file("*.stripped.vcf.gz.csi"), emit: name_vcf_csi optional true
    tuple val(name), file("${simplename}.zero_variants.txt"), emit: name_nvariants_txt

    script:
    def simplename = vcf.getSimpleName()
""" 
sleep 10
N_VARIANTS=\$(zcat $vcf | grep -v '^#' |  wc -l)
if [ \"\${N_VARIANTS}\" -eq "0" ]
then
    echo \"N variants in $vcf is 0\" > ${simplename}.zero_variants.txt && sleep 5
else
    echo \"N variants in $vcf is not 0\" > ${simplename}.zero_variants.txt && sleep 5
    bcftools view -G -o ${simplename}.stripped.vcf.gz -O z $vcf
    bcftools index ${simplename}.stripped.vcf.gz
fi
"""
}

