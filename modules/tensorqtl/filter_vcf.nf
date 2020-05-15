params.run = true

process filter_vcf {
    tag "${prefix}"
    queue 'gpu-normal'
    clusterOptions '-gpu num=1'
    // clusterOptions '-gpu "num=1:mode=exclusive_process"'
    maxForks 2
    conda '/lustre/scratch118/humgen/resources/conda_envs/singlecell_eqtl'
    memory = '55G'
    time '700m'
    cpus 1
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/tensorqtl/", mode: 'symlink', pattern: "${prefix}.cis_qtl.txt.gz", overwrite: true

    when:
    params.run

    input: 
    set val(prefix), val(plink_prefix_path), file(bed), file(bim), file(fam)
    set val(MAF_threshold)

    output: 
    tuple val(prefix), file("${prefix}.cis_qtl.txt.gz"), emit: tensorqtl_parquet

    script:
    """
export PATH=/lustre/scratch118/humgen/resources/conda_envs/singlecell_eqtl/bin:\$PATH

tabix -p ${vcf_gz}
# add MAF
bcftools +fill-tags ${vcf_gz} -Oz -o sub2.vcf.gz
# filter MAF
bcftools view -i 'MAF[0]>${MAF_threshold}'
    """
}
