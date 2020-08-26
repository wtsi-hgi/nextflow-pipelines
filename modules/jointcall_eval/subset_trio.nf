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
    conda '/lustre/scratch118/humgen/resources/conda_envs/rtg_tools'

    when:
    params.run

    input: 
    set file(vcf), file(tbi)
    val(trio)

    output: 
    tuple file("subset_trio.vcf.gz"), file('subset_trio.vcf.gz.tbi'), emit: vcf_tbi

    script:
    """
export PATH=/lustre/scratch118/humgen/resources/conda_envs/rtg_tools/bin:\$PATH

bcftools view -Oz -o subset_trio.vcf.gz \\
-s $trio $vcf

tabix -p subset_trio.vcf.gz
    """
}
