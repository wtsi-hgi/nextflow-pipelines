params.run = true

process subset_chr {
    tag "${chr}"
    memory = '4G'
    time '240m'
    cpus 1
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    maxForks 12
    publishDir "${params.outdir}/subset_chr/${chr}/", mode: 'symlink', overwrite: true

    when:
    params.run

    input: 
    file vcf, val chr

    output: 
    tuple file("subset_chr.vcf.gz"), file('subset_chr.vcf.tbi'), emit: vcf_tbi

    script:
    """
bcftools view -Oz -o subset_chr.vcf.gz \\
-r $chr
    """
}
