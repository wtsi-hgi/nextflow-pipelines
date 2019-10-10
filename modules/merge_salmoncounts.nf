params.run = true
params.runtag = 'runtag'

process merge_salmoncounts {
    tag "${input_trans}/${input_genes}"
    publishDir "${params.outdir}/combined", mode: 'symlink'
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    time '120m'

    input:
    file input_trans from ch_salmon_trans.map { it.toString() }.collectFile(name: 'trans.meta', newLine: true)
    file input_genes from ch_salmon_genes.map { it.toString() }.collectFile(name: 'genes.meta', newLine: true)

    when:
    params.run

    output:
    set file('*counts.txt'), file('*tpm.txt')

    script:
    def outtranscount = "${params.runtag}-salmon-transcounts.txt"
    def outgenescount = "${params.runtag}-salmon-genecounts.txt"
    def outtranstpm   = "${params.runtag}-salmon-transtpm.txt"
    def outgenestpm   = "${params.runtag}-salmon-genetpm.txt"
    """
    python3 $workflow.projectDir/bin/merge_featurecounts.py           \\
      --rm-suffix .quant.genes.sf                                     \\
      -c -1 --skip-comments --header                                  \\
      -o $outgenescount -I $input_genes
    python3 $workflow.projectDir/bin/merge_featurecounts.py           \\
      --rm-suffix .quant.sf                                           \\
      -c -1 --skip-comments --header                                  \\
      -o $outtranscount -I $input_trans
    python3 $workflow.projectDir/bin/merge_featurecounts.py           \\
      --rm-suffix .quant.genes.sf                                     \\
      -c -2 --skip-comments --header                                  \\
      -o $outgenestpm -I $input_genes
    python3 $workflow.projectDir/bin/merge_featurecounts.py           \\
      --rm-suffix .quant.sf                                           \\
      -c -2 --skip-comments --header                                  \\
      -o $outtranstpm -I $input_trans
    """
}
