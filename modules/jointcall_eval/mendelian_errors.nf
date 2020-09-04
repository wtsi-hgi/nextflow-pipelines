params.run = true

process mendelian_errors {
    tag "${trio}"
    memory = '4G'
    time '240m'
    cpus 1
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    maxForks 12
    publishDir "${params.outdir}/${tag}/bcftools_trio_mendelian/", mode: 'symlink', overwrite: true
    conda '/lustre/scratch118/humgen/resources/conda_envs/rtg_tools'
    
    when:
    params.run

    input: 
    set file(vcf), file(tbi)
    val(trio)
    val(tag)

    output: 
    file("bcftools.trio_mendelian.txt")

    script:
    """
export PATH=/lustre/scratch118/humgen/resources/conda_envs/rtg_tools/bin:\$PATH

bcftools +mendelian \\
  $vcf -t $trio -c -o bcftools.trio_mendelian.txt
    """
}
