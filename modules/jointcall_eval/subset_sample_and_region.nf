params.run = true

process subset_sample_and_region {
    tag "${sample}"
    memory = '4G'
    time '240m'
    cpus 1
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    maxForks 12
    publishDir "${params.outdir}/subset_sample_and_region/${sample}/", mode: 'symlink', overwrite: true
    conda '/lustre/scratch118/humgen/resources/conda_envs/rtg_tools'

    when:
    params.run

    input: 
    set file(vcf), file(tbi)
    val(sample)
    val(region_bed)

    output: 
    tuple file("subset_sample_region.vcf.gz"), file('subset_sample_region.vcf.gz.tbi'), emit: vcf_tbi

    script:
    """
export PATH=/lustre/scratch118/humgen/resources/conda_envs/rtg_tools/bin:\$PATH

bcftools view -Oz -o subset_sample_region.vcf.gz \\
-s $sample -R $region_bed $vcf

tabix -p subset_sample_region.vcf.gz
    """
}
