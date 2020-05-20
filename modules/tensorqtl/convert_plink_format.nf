params.run = true

process convert_plink_format {
    tag "${plink_prefix}"
    queue 'normal'
    maxForks 2
    conda '/lustre/scratch118/humgen/resources/conda_envs/tensorqtl'
    memory = '55G'
    time '700m'
    cpus 1
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/convert_plink_format/", mode: 'symlink', pattern: "${plink_prefix}.vcf.gz", overwrite: true

    when:
    params.run

    input: 
    tuple file(bed), file(bim), file(fam)
    val(plink_prefix)

    output: 
    tuple file("${plink_prefix}.vcf.gz"), emit: vcf_gz

    script:
    """
export PATH=/lustre/scratch118/humgen/resources/conda_envs/tensorqtl/bin:\$PATH

plink2 --bfile ${plink_prefix} --keep-allele-order --recode vcf --out ${plink_prefix}.pre
bgzip -c ${plink_prefix}.pre.vcf > ${plink_prefix}.pre.vcf.gz
tabix -p vcf ${plink_prefix}.pre.vcf.gz

# add chr to chr names
for i in {1..22} X Y MT
do
  echo \"\$i chr\$i\" >> chr_name_conv.txt
done
bcftools annotate --rename-chrs chr_name_conv.txt ${plink_prefix}.pre.vcf.gz \\
  -Oz -o ${plink_prefix}.vcf.gz

rm ${plink_prefix}.pre*
    """
}
