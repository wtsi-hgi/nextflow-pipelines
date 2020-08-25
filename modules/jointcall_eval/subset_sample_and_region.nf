params.run = true

process subset_sample {
    tag "${sample}"
    memory = '4G'
    time '240m'
    cpus 1
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    maxForks 12
    publishDir "${params.outdir}/subset_sample_and_region/${sample}/", mode: 'symlink', overwrite: true

    when:
    params.run

    input: 
    file vcf, val sample, val region_bed

    output: 
    tuple file("subset_sample_region.vcf.gz"), file('subset_sample_region.vcf.tbi'), emit: vcf_tbi

    script:
    """
bcftools view -Oz -o subset_sample_region.vcf.gz \\
-s $sample -R $region_bed
    """
}
