params.run = true

process bcftools_stats {
    //tag ""
    memory = '4G'
    time '240m'
    cpus 1
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    maxForks 12
    publishDir "${params.outdir}/${tag}/bcftools/", mode: 'symlink', overwrite: true
    conda '/lustre/scratch118/humgen/resources/conda_envs/rtg_tools'

    when:
    params.run

    input: 
    set file(vcf), file(tbi)
    val(tag)

    output: 
    file("${vcf}.bcftools.stats") //, emit: vcf_tbi

    script:
    """
export PATH=/lustre/scratch118/humgen/resources/conda_envs/rtg_tools/bin:\$PATH

bcftools stats $vcf > ${vcf}.bcftools.stats
    """
}
