params.run = true

process filter_genotype_data {
    tag "${prefix}"
    queue 'normal'
    maxForks 2
    conda '/lustre/scratch118/humgen/resources/conda_envs/tensorqtl'
    memory = '55G'
    time '700m'
    cpus 1
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/filter_genotype_data/", mode: 'symlink', pattern: "${plink_prefix}.filtered.*", overwrite: true

    when:
    params.run

    input: 
    set file(bed), file(bim), file(fam)
    val(plink_prefix)
    val(MAF_threshold)

    output: 
    tuple file("${plink_prefix}.filtered.bed"), file("${plink_prefix}.filtered.bim"), file("${plink_prefix}.filtered.fam"), emit: bed_bim_fam

    script:
    """
export PATH=/lustre/scratch118/humgen/resources/conda_envs/tensorqtl/bin:\$PATH

# convert to vcf.gz
# plink2 --bfile ${plink_prefix} --keep-allele-order --recode vcf --out ${plink_prefix}.vcf
# bgzip -c ${plink_prefix}.vcf ${plink_prefix}.vcf.gz
# tabix -p vcf ${plink_prefix}.vcf.gz

# apply filters
# add MAF
# bcftools +fill-tags ${plink_prefix}.vcf.gz -Oz -o sub2.vcf.gz
# filter MAF
# bcftools view -i 'MAF[0]>${MAF_threshold}'

# clean up
# plink2 --make-bed --vcf ${plink_prefix}.vcf.gz --out ${plink_prefix}
# rm ${plink_prefix}.vcf*


# for testing pipeline 
cp ${plink_prefix}.bed ${plink_prefix}.filtered.bed
cp ${plink_prefix}.bim ${plink_prefix}.filtered.bim
cp ${plink_prefix}.fam ${plink_prefix}.filtered.fam
    """
}
