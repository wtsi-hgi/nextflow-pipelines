params.run = true

process hla_la {
    memory '4G'
    tag "$samplename"
    cpus 1
    disk '20 GB'
    scratch '/tmp'
    stageInMode 'copy'
    stageOutMode 'rsync'
    time '1000m'
    container "hla-la-1.0.1"
    maxForks 2
    // containerOptions = "--bind /home/ubuntu"
    // errorStrategy 'terminate'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    publishDir "${params.outdir}/copy_number_splitchr/chr1/", mode: 'symlink', overwrite: true, pattern: "${samplename}.cn.chr1.recode.vcf"
    maxRetries 1

    when:
    params.run
     
    input: 
    tuple val(samplename), file(cn_vcf)
    
    output: 
    tuple val(samplename), file("${samplename}.cn.chr*.recode.vcf"), emit: samplename_cn_vcf

    script:
    """ 
    """
}
