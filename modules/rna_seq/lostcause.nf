params.run = true
params.runtag = 'runtag'

process lostcause {

    publishDir "${params.outdir}/combined", mode: 'symlink'

    input:
    file (inputs) //from ch_lostcause.collectFile{ ['lostcause.txt', it.text] }

    output:
    file('*.lostcause_mqc.txt') //into ch_multiqc_lostcause

    script:
    def outputname = "${params.runtag}.lostcause_mqc.txt"
    """
    echo -e "# plot_type: 'table'\n# section_name: 'Lost samples'" > $outputname
    echo -e "Sample\tProcess\tMessage" >> $outputname
    cat $inputs | sort >> $outputname
    """
}
