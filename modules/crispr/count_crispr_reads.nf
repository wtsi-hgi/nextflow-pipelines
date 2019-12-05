params.run = true
params.runtag = 'runtag'
params.read2 = 'discard'

process count_crispr_reads {
    tag "read_counts $samplename"
    // container "nfcore-rnaseq"
    // containerOptions = "--bind /lustre"
    publishDir "${params.outdir}/read_counts", mode: 'symlink',
        saveAs: { filename ->
            if (filename ==~ /.*\.unmapped_sequences\.txt/) "unmapped_sequences/$filename"
            else filename
        }
    memory = '10G'
    cpus 1
    time '300m'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    when:
    params.run

    input:
    set val(samplename), file(fastq_files), val(guide_library), val(includeG), val(start_match_sample), val(end_match_sample), val(start_match_lib), val(end_match_lib) 
    file(library_files)

    output:
    set val(guide_library), file("*.counts.txt")
    file("${samplename}.mapping.txt")
    set val(samplename), val(guide_library), stdout
    set val(samplename), file("*.unmapped_sequences.txt")

    shell:
    """
    if [ \"${params.read2}\"  == \"discard\" ]; then
    rm -f *_2.fastq.gz
    fi

    python3 ${workflow.projectDir}/../bin/crispr/count_reads.py \$(ls *.fastq.gz) \"${guide_library}\" ${samplename}.counts.txt ${samplename}.mapping.txt ${includeG} ${start_match_sample} ${end_match_sample} ${start_match_lib} ${end_match_lib}
    """
}
