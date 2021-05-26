params.run = true

params.ref_dir = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/"

process sort_cram {
    memory '10G'
    tag "$cram_file"
    //cpus 1
    disk '20 GB'
    //time '100m'
    //queue 'normal'
    container  = 'file:///software/hgi/containers/gatk-4.1.4.1.sif'
    containerOptions = "--bind /lustre --bind ${params.ref_dir}:/ref --bind /tmp:/tmp"
    // errorStrategy 'terminate'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    //publishDir "${params.outdir}/cram_index/", mode: 'symlink', overwrite: true, pattern: "${cram_file}.crai"
    maxRetries 3

    when:
    params.run
     
    input:
    path cram_file

    output:
    path "${cram_file}.sorted"
    //tuple file("${cram_file}.sorted"), emit: indexes

    script:
""" 
/gatk/gatk --java-options "-Xms4g -Xmx4g  -XX:+UseSerialGC" SortSam -I ${cram_file} -O ${cram_file}.sorted -SO queryname -R /ref/hs38DH.fa --TMP_DIR /tmp --MAX_RECORDS_IN_RAM 300000
"""
}

