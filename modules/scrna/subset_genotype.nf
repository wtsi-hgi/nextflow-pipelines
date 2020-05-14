params.run = true

process 'subset_genotype' {

    disk '20 GB'
    scratch '/tmp'
    stageInMode 'symlink'
    stageOutMode 'rsync'

    cpus = 2
    time '8000m'
    memory = 20.GB 
// {  100.GB + 50.GB * (task.attempt-1) }
    // queue 'basement'

    // maxForks 2
    tag "$samplename"
    container "lifebitai.bcftools"
    containerOptions = "--bind /"


    errorStrategy = { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    
    publishDir "${params.outdir}/subset_genotype/", mode: 'symlink', pattern: "${samplename}.subset.vcf.gz"

    when:
    params.run 

    input:
    set val(samplename), file(cellsnp_vcf), file(sample_subset_file), file(donor_vcf)
    
    output:
    tuple val(samplename), file("${samplename}.subset.vcf.gz"), emit: samplename_subsetvcf

  script:
   """
tabix -p vcf ${donor_vcf} 
# tabix -p vcf ${cellsnp_vcf}
bcftools view ${donor_vcf} -R ${cellsnp_vcf} -S ${sample_subset_file} -Oz -o ${samplename}.subset.vcf.gz
rm ${donor_vcf}.tbi
# rm ${cellsnp_vcf}.tbi
   """
}
