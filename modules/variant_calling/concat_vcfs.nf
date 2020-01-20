params.run = true

process index_cram {
    memory '3G'
    tag "$cram_file"
    cpus 1
    // disk '20 GB'
    time '100m'
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
    tuple val(vcfs_location), val(name)
    
    output: 
    tuple file("${name}.vcf.gz"), emit: vcf_gz

    script:
""" 
bcftools concat \$(find $vcfs_location -name '*.vcf.gz' | sort | paste -sd ' ') --allow-overlaps | bcftools sort -o ${name}.vcf.gz -O z
bcftools index ${name}.vcf.gz
"""
}

