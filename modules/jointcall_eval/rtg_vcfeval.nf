params.run = true

process rtg_vcfeval {
    tag "${tag} ${sample}"
    memory = '4G'
    time '240m'
    cpus 4
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    maxForks 12
    publishDir "${params.outdir}/${tag}/rtg_vcfeval/", mode: 'symlink', overwrite: true
    conda '/lustre/scratch118/humgen/resources/conda_envs/rtg_tools'

    when:
    params.run

    input: 
    set file(vcf), file(tbi)
    val(sample)
    val(rtg_sample)
    set file(rtg_vcf), file(rtg_tbi)
    file(rtg_genome_template)
    file(rtg_region)
    val(tag)
    
    output: 
    file("${tag}_${sample}")

    script:
    """
export PATH=/lustre/scratch118/humgen/resources/conda_envs/rtg_tools/bin:\$PATH

rtg vcfeval \\
--baseline $rtg_vcf \\
--calls $vcf \\
--output ${tag}_${sample} \\
--evaluation-regions $rtg_region \\
--bed-regions $rtg_region \\
--template $rtg_genome_template \\
--sample "${rtg_sample},${sample}" \\
--threads=4 \\
--output-mode \"split\"
    """
}
