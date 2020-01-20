params.run = true

process concat_vcfs {
    memory '3G'
    tag "$vcfs_location $name"
    cpus 2
    // disk '20 GB'
    time '700m'
    queue 'normal'
    container "graphtyper"
    containerOptions = "--bind /lustre"
    // errorStrategy 'terminate'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    publishDir "${params.outdir}/vcfs_concat/", mode: 'symlink', overwrite: true, pattern: "*.vcf.gz"
    maxRetries 3

    when:
    params.run
     
    input:
    val(vcfs_location)
    val(name)
    
    output: 
    tuple file("${name}.vcf.gz"), emit: vcf_gz

    script:
""" 
bcftools concat \$(find $vcfs_location -name '*.vcf.gz' | sort | paste -sd ' ') --allow-overlaps | bcftools sort -o ${name}.vcf.gz -O z
bcftools index ${name}.vcf.gz
"""
}

