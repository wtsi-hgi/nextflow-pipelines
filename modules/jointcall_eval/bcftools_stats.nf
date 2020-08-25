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

    when:
    params.run

    input: 
    file vcf, val sample

    output: 
    tuple file("subset_sample.vcf.gz"), file('subset_sample.vcf.tbi'), emit: vcf_tbi

    script:
    """
bcftools view -Oz -o subset_sample.vcf.gz \\
-s $sample
    """
}
