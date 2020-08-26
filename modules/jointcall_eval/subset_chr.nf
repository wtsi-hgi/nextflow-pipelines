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
    conda '/lustre/scratch118/humgen/resources/conda_envs/rtg_tools'
    
    when:
    params.run

    input: 
    set file(vcf), file(tbi)
    val(chr)

    output: 
    tuple file("subset_chr.vcf.gz"), file('subset_chr.vcf.gz.tbi'), emit: vcf_tbi

    script:
    """
export PATH=/lustre/scratch118/humgen/resources/conda_envs/rtg_tools/bin:\$PATH

bcftools view -Oz -o subset_chr.vcf.gz \\
-r $chr $vcf

tabix -p vcf subset_chr.vcf.gz
    """
}
