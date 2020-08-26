params.run = true

process bcftools_stats {
    tag "${sample}"
    memory = '4G'
    time '240m'
    cpus 1
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    maxForks 12
    publishDir "${params.outdir}/subset_sample/${sample}/", mode: 'symlink', overwrite: true
    conda '/lustre/scratch118/humgen/resources/conda_envs/rtg_tools'

    when:
    params.run

    input: 
    set file(vcf), file(tbi)
    val(sample)

    output: 
    tuple file("subset_sample.vcf.gz"), file('subset_sample.vcf.gz.tbi'), emit: vcf_tbi

    script:
    """
export PATH=/lustre/scratch118/humgen/resources/conda_envs/rtg_tools/bin:\$PATH

bcftools view -Oz -o subset_sample.vcf.gz \\
-s $sample $vcf
    """
}
