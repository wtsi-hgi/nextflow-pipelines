params.run = true

process subset_trio {
    tag "${trio}"
    memory = '4G'
    time '240m'
    cpus 1
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    maxForks 12
    publishDir "${params.outdir}/subset_trio/${trio}/", mode: 'symlink', overwrite: true

    when:
    params.run

    input: 
    file vcf, val trio

    output: 
    tuple file("subset_trio.vcf.gz"), file('subset_trio.vcf.tbi'), emit: vcf_tbi

    script:
    """
bcftools view -Oz -o subset_trio.vcf.gz \\
-s $trio
    """
}
