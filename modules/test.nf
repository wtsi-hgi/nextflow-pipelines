params.run = true
params.runtag = 'runtag'

process test {
    cache 'false'
    tag "$aligner"
    container "nfcore-rnaseq"
    publishDir "${params.outdir}/combined", mode: 'symlink'
    label 'merge_feature'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    when:
    params.run

    input:
    set val(aligner), file(files), val(files_list)  //_collectd, list_of_files //from ch_merge_fc_byaligner

    output: 
    set val(aligner), stdout
    // file(list_of_files)
    //file '*-fc-genecounts.txt'

    //shell:
    //suffix=['star':'.star.gene.fc.txt', 'hisat2':'.hisat2.gene.fc.txt']
    //aligner = metafile.baseName   // not strictly necessary
    //outputname = "${params.runtag}-${aligner}-fc-genecounts.txt"
    //thesuffix  = suffix[aligner] ?: '.txt'
    script:
    """
ls ${files}
    """
}
