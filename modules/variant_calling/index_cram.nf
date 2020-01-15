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
    publishDir "${params.outdir}/cram_index/", mode: 'symlink', overwrite: true, pattern: "${cram_file}.crai"
    maxRetries 3

    when:
    params.run
     
    input:
    file(cram_file)

    output: 
    tuple file("${cram_file}.crai"), emit: indexes

    script:
""" 
samtools index $cram_file
"""
}

