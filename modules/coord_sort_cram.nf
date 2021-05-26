params.run = true

params.ref_dir = "/lustre/scratch118/humgen/resources/ref/Homo_sapiens/HS38DH/"

process coord_sort_cram {
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
    path cram_file_sorted_dups

    output:
    tuple path("${cram_file_sorted_dups}.coord"), path("${cram_file_sorted_dups}.coord.bai")
    //tuple file("${cram_file}.sorted"), emit: indexes

    script:
""" 
/gatk/gatk --java-options "-Xms4g -Xmx4g  -XX:+UseSerialGC" SortSam -I ${cram_file_sorted_dups} -O ${cram_file_sorted_dups}.coord -SO "coordinate" -R /ref/hs38DH.fa --CREATE_INDEX true --TMP_DIR /tmp
"""
}

