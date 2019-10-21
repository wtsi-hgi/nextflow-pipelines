params.run = true
params.runtag = 'runtag'
params.read2 = 'discard'

process count_crispr_reads {
    tag "read_counts $samplename"
    // container "nfcore-rnaseq"
    publishDir "${params.outdir}/read_counts", mode: 'symlink'
    memory = '10G'
    cpus 1
    time '300m'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    when:
    params.run

    input:
    set val(samplename), file(fastq_files), val(guide_library) 
    file(library_files)

    output:
    file("${samplename}.counts.txt")
    file("${samplename}.mapping.txt")
    set val(samplename), stdout

    shell:
    """
    if [ \"${params.read2}\"  == \"discard\" ]; then
    rm -f *_2.fastq.gz
    fi

    python3 ${workflow.projectDir}/../bin/crispr/count_reads.py \$(ls *.fastq.gz) \"${guide_library}\"
    """
}
