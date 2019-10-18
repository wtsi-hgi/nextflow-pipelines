params.run = true
params.runtag = 'runtag'
params.read2 = 'discard'

process counts_crispr_reads {
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
    set val(samplename), file(fastsq)
    file(library)

    output:
    file("${samplename}.counts.txt")
    file("${samplename}.mapping.txt")

    shell:
    """
    if [ \"${params.read2}\"  == \"discard\" ]; then
    rm *_2.fastq.gz
    fi

    python3 ${workflow.projectDir}/bin/crispr/read_counts.py 
    """
}
