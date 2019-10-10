params.run = true 

params.fcextra = ""   
params.singleend = false 
params.forward_stranded = false
params.reverse_stranded = true
params.unstranded = false

process featureCounts {
    tag "${samplename}"
    container "nfcore-rnaseq"
    memory = '15G'
    time '300m'
    cpus 1
    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5
    publishDir "${params.outdir}/featureCounts/", mode: 'symlink',
        saveAs: {filename ->
            if (filename.indexOf(".biotype_counts_mqc.txt") > 0) "biotype_counts/$filename"
            else if (filename.indexOf(".gene.featureCounts.txt.summary") > 0) "gene_count_summaries/$filename"
            else if (filename.indexOf(".gene.featureCounts.txt") > 0) "gene_counts/$filename"
            else "$filename"
        }

    when:
    params.run
    
    input:
    set val(aligner), val(samplename), file(thebam) //from ch_featurecounts
    file gtf //from ch_gtf_featureCounts  //.collect()
    file biotypes_header

    output:
    set val(aligner), file("*.gene.fc.txt") //into ch_merge_fc
    set val(aligner), file("*.gene.fc.txt.summary") //into ch_multiqc_fc
    set val(aligner), file("*.biotype_counts*mqc.txt") //into ch_multiqc_fcbiotype

    script:
    def extraparams = params.fcextra.toString() - ~/^dummy/
    def fc_direction = 0
    def tag = "${samplename}.${aligner}"

    def pairedend = params.singleend ? "" : "-p"
    if (params.forward_stranded && !params.unstranded) {
        fc_direction = 1
    } else if (params.reverse_stranded && !params.unstranded){
        fc_direction = 2
    }
    outfile = "${tag}.gene.fc.txt"
    """
    export PATH=/opt/conda/envs/nf-core-rnaseq-1.3/bin:$PATH

    featureCounts -T ${task.cpus} -a $gtf -g gene_id          \\
      -o ${outfile} $pairedend                                \\
      -s $fc_direction ${extraparams} $thebam
    cut -f 1,7 ${outfile} > reduced.${outfile}   #  This
    mv reduced.${outfile} ${outfile}             #  reduces the file size from ~ 30M to ~1M
    featureCounts -T ${task.cpus} -a $gtf -g gene_id  \\
      -o ${tag}.biotype.fc.txt $pairedend                     \\
      -s $fc_direction ${extraparams} $thebam
    cut -f 1,7 ${tag}.biotype.fc.txt |                        \\
        tail -n +3 | cat $biotypes_header - >> ${tag}.biotype_counts_mqc.txt
    """
}
