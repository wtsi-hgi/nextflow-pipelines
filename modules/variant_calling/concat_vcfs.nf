params.run = true

process concat_vcfs {
    memory '6G'
    tag "concat"
    cpus 1
    conda '/lustre/scratch118/humgen/hgi/projects/ibdx10/variant_calling/joint_calling/ibd_concat_nextflow/bcftools'
    scratch '/tmp'
    stageInMode 'copy'
    stageOutMode 'copy'
    time '700m'
    queue 'normal'
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    publishDir "${params.outdir}/concat/", mode: 'symlink', overwrite: true, pattern: "*.genome.sorted.vcf.gz"
    publishDir "${params.outdir}/concat/", mode: 'symlink', overwrite: true, pattern: "*.genome.sorted.vcf.gz.csi"
    publishDir "${params.outdir}/concat/", mode: 'symlink', overwrite: true, pattern: "vcf_files_sorted"
    maxRetries 2

    when:
    params.run
     
    input:
    tuple val(batch), file(vcf_files)
    
    output:
    tuple val(batch), file("*.genome.vcf.gz"), file("*.genome.vcf.gz.csi"), emit: batch_vcf 
    tuple val(batch), file("vcf_files_sorted"), emit: batch_vcflist 

    script:
""" 
ls *.vcf.gz | grep -v tbi | grep -v csi > vcf_files
cat vcf_files | sort > vcf_files_sorted

bcftools concat -n -O z -f vcf_files_sorted -o concat.vcf.gz
bcftools sort -O z concat.vcf.gz -o genome.sorted.vcf.gz
bcftools index genome.sorted.vcf.gz
"""
}

