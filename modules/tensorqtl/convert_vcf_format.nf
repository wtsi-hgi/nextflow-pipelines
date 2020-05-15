params.run = true

process convert_vcf_format {
    tag "${params.plink_prefix}"
    queue 'normal'
    maxForks 2
    conda '/lustre/scratch118/humgen/resources/conda_envs/tensorqtl'
    memory = '55G'
    time '700m'
    cpus 1
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/convert_vcf_format/", mode: 'symlink', pattern: "${plink_prefix}.*", overwrite: true

    when:
    params.run

    input: 
    set file(vcf_gz)
    set val(plink_prefix)

    output: 
    tuple file("${plink_prefix}.bed"), file("${plink_prefix}.bim"), file("${plink_prefix}.fam"), emit: bed_bim_fam

    script:
    """
export PATH=/lustre/scratch118/humgen/resources/conda_envs/tensorqtl/bin:\$PATH

plink2 --make-bed --vcf ${vcf_gz} --out ${plink_prefix}
    """
}
