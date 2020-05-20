params.run = true

process convert_vcf_format {
    tag "${plink_prefix}"
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

tabix -p vcf ${vcf_gz}

# add chr to chr names
for i in {1..22} X Y MT
do
  echo \"\$i chr\$i\" >> chr_name_conv.txt
done
bcftools annotate --rename-chrs chr_name_conv.txt ${vcf_gz} \\
  -Oz -o ${plink_prefix}.pre.vcf.gz
plink2 --make-bed --vcf ${plink_prefix}.pre.vcf.gz --out ${plink_prefix}
sed -i s/^/chr/g ${plink_prefix}.bim
sed -i s/^chr23/chrX/g ${plink_prefix}.bim

rm ${vcf_gz}.tbi
rm ${plink_prefix}.pre.vcf.gz*
    """
}
