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
    val(GT_miss_threshold)

    output: 
    tuple file("${plink_prefix}.filtered.bed"), file("${plink_prefix}.filtered.bim"), file("${plink_prefix}.filtered.fam"), emit: bed_bim_fam

    script:
    """
export PATH=/lustre/scratch118/humgen/resources/conda_envs/tensorqtl/bin:\$PATH

# convert to vcf.gz
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

tabix -p vcf ${plink_prefix}.vcf.gz

# apply filters
# add MAF
bcftools +fill-tags ${plink_prefix}.vcf.gz -Oz -o to_filter.vcf.gz
# keep only biallelic SNPs with: -m2 -M2 -v snps
# keep only sites with GT missingness above threshold
# also filter MAF with: -i 'MAF[0]>${MAF_threshold}' 
bcftools view -m2 -M2 -v snps -i 'MAF[0]>${MAF_threshold} && AN>=2*N_SAMPLES*${GT_miss_threshold}' \\
    to_filter.vcf.gz -Oz -o ${plink_prefix}.vcf.gz 

plink2 --make-bed --vcf ${plink_prefix}.vcf.gz --out ${plink_prefix}.filtered
sed -i s/^/chr/g ${plink_prefix}.filtered.bim
sed -i s/^chr23/chrX/g ${plink_prefix}.filtered.bim

rm ${plink_prefix}.vcf*


# for testing pipeline 
# cp ${plink_prefix}.bed ${plink_prefix}.filtered.bed
# cp ${plink_prefix}.bim ${plink_prefix}.filtered.bim
# cp ${plink_prefix}.fam ${plink_prefix}.filtered.fam
    """
}
