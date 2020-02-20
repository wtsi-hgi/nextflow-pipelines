params.run = true

process sect_concat_vcfs {
    memory '60G'
    tag "$batch"
    cpus 2
    // disk '20 GB'
    time '700m'
    queue 'normal'
    container "graphtyper"
    containerOptions = "--bind /lustre"
    // errorStrategy 'terminate'
    errorStrategy { task.attempt <= 0 ? 'retry' : 'ignore' }
    //publishDir "${params.outdir}/setc_vcfs_concat/", mode: 'symlink', overwrite: true, pattern: "${name}.vcf.gz"
    maxRetries 0

    when:
    params.run
     
    input:
    tuple val(batch), file(vcf_files)
    file(bed)
    
    output:
    stdout
    //tuple file("${name}.vcf.gz"), file("${name}.vcf.gz.csi"), emit: vcf_gz
    //tuple file("to_concat.list"), file("to_concat_non_empty.list"), emit: concat_list

    script:
""" 
ls -ltra
"""
}

