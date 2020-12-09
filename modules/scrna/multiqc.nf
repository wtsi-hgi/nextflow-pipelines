params.run = true
params.runtag = 'runtag'

process multiqc {
    container "nfcore-rnaseq"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    memory = '10G'
    cpus 2
    time '300m'

    publishDir "${params.outdir}/multiqc", mode: 'copy',
      saveAs: {filename ->
          if (filename.indexOf("multiqc.html") > 0) "combined/$filename"
          else if (filename.indexOf("_data") > 0) "$filename"
          else null
      }

    when:
    params.run

    input:
    file (fastqc:'fastqc/*') //from ch_multiqc_fastqc.collect().ifEmpty([])
    //file ('lostcause/*') //from ch_multiqc_lostcause.collect().ifEmpty([])
    //file ('mapsummary/*') //from ch_multiqc_mapsum.collect().ifEmpty([])
    //file ('featureCounts/*') //from ch_multiqc_fc_aligner.collect().ifEmpty([])
    //file ('featureCounts_biotype/*') //from ch_multiqc_fcbiotype_aligner.collect().ifEmpty([])
    //file ('star/*') //from ch_alignment_logs_star.collect().ifEmpty([])
    //file ('salmon/*') //from ch_alignment_logs_salmon.collect().ifEmpty([])

    output:
    file "*_multiqc.html"
    file "*_data"

    script:
    def filename = "${params.runtag}_multiqc.html"
    def reporttitle = "${params.runtag}"
    """
    export PATH=/opt/conda/envs/nf-core-rnaseq-1.3/bin:\$PATH

    multiqc . -f --title "$reporttitle" --filename "$filename" -m fastqc
    """
}
// multiqc . -f --title "$reporttitle" --filename "$filename" -m custom_content -m featureCounts -m star -m fastqc -m salmon
