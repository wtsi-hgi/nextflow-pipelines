params.run = true
params.runtag = 'runtag'

process merge_featureCounts {
    tag "$aligner"
    container "nfcore-rnaseq"
    publishDir "${params.outdir}/combined", mode: 'symlink'
    containerOptions = "--bind /lustre"
    label 'merge_feature'
    memory = '100G'
    cpus 2
    time '600m'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    when:
    params.run

    input:
    file(collected_fc_gene_txt)

    output:
    file '*-fc-genecounts.txt'
    file("fofn_gene_featurecount.txt")

    shell:
    suffix=['star':'.star.gene.fc.txt', 'hisat2':'.hisat2.gene.fc.txt']
    aligner = "star"  // not strictly necessary
    outputname = "${params.runtag}-${aligner}-fc-genecounts.txt"
    thesuffix  = suffix[aligner] ?: '.txt'
    '''
    export PATH=/opt/conda/envs/nf-core-rnaseq-1.3/bin:$PATH

    ls . | grep gene.fc.txt\$ > fofn_gene_featurecount.txt

    python3 !{workflow.projectDir}/../bin/rna_seq/merge_featurecounts.py        \\
      --rm-suffix !{thesuffix}                                       \\
      -c 1 --skip-comments --header                                  \\
      -o !{outputname} -I fofn_gene_featurecount.txt
    '''
}
