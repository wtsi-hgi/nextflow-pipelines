params.run = true
params.runtag = 'runtag'

process merge_salmoncounts {
    tag ""
    scratch '/tmp'
    stageInMode 'copy'
    stageOutMode 'rsync'
    container "nfcore-rnaseq"
    publishDir "${params.outdir}/combined", mode: 'rellink'
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    containerOptions = "--bind /lustre"
    maxRetries 6
    // memory = '10G'
    memory = {  80.GB + 20.GB * (task.attempt-1) }
    time '400m'

    input:
    file (all_quant_sf)
    file (all_quant_genes_sf)

    when:
    params.run

    output:
    set file('*transcounts.txt'), file('*transtpm.txt'), file('*genecounts.txt'), file('*genetpm.txt')
    file("fofn_quant_sf_salmon.txt")
    file("fofn_quant_genes_sf_salmon.txt")
    
    script:
    def outtranscount = "${params.runtag}-salmon-transcounts.txt"
    def outgenescount = "${params.runtag}-salmon-genecounts.txt"
    def outtranstpm   = "${params.runtag}-salmon-transtpm.txt"
    def outgenestpm   = "${params.runtag}-salmon-genetpm.txt"
    """
    export PATH=/opt/conda/envs/nf-core-rnaseq-1.3/bin:\$PATH

    ls . | grep .quant.sf\$ > fofn_quant_sf_salmon.txt
    ls . | grep .quant.genes.sf\$ > fofn_quant_genes_sf_salmon.txt

    python3 $workflow.projectDir/../bin/rna_seq/merge_featurecounts.py           \\
      --rm-suffix .quant.genes.sf                                     \\
      -c -1 --skip-comments --header                                  \\
      -o $outgenescount -I fofn_quant_genes_sf_salmon.txt

    python3 $workflow.projectDir/../bin/rna_seq/merge_featurecounts.py           \\
      --rm-suffix .quant.sf                                           \\
      -c -1 --skip-comments --header                                  \\
      -o $outtranscount -I fofn_quant_sf_salmon.txt

    python3 $workflow.projectDir/../bin/rna_seq/merge_featurecounts.py           \\
      --rm-suffix .quant.genes.sf                                     \\
      -c -2 --skip-comments --header                                  \\
      -o $outgenestpm -I fofn_quant_genes_sf_salmon.txt

    python3 $workflow.projectDir/../bin/rna_seq/merge_featurecounts.py           \\
      --rm-suffix .quant.sf                                           \\
      -c -2 --skip-comments --header                                  \\
      -o $outtranstpm -I fofn_quant_sf_salmon.txt
    """
}
