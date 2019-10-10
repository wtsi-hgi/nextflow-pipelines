params.run = true
params.runtag = 'runtag'

process merge_featureCounts {
    cache 'deep'
    tag "$aligner"
    container "nfcore-rnaseq"
    publishDir "${params.outdir}/combined", mode: 'symlink'
    label 'merge_feature'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    when:
    params.run

    input:
    file metafile //from ch_merge_fc_byaligner

    output:
    file '*-fc-genecounts.txt'

    shell:
    suffix=['star':'.star.gene.fc.txt', 'hisat2':'.hisat2.gene.fc.txt']
    aligner = metafile.baseName   // not strictly necessary
    outputname = "${params.runtag}-${aligner}-fc-genecounts.txt"
    thesuffix  = suffix[aligner] ?: '.txt'
    '''
    export PATH=/opt/conda/envs/nf-core-rnaseq-1.3/bin:\$PATH

    python3 !{workflow.projectDir}/bin/merge_featurecounts.py        \\
      --rm-suffix !{thesuffix}                                       \\
      -c 1 --skip-comments --header                                  \\
      -o !{outputname} -I !{metafile}
    '''
}
