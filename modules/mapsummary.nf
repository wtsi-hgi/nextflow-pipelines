params.mito_name = 'MT'

process mapsummary {
    tag "${samplename}"
    container "nfcore-rnaseq"
    publishDir "${params.outdir}/mapsummary/", mode: 'symlink'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    memory = '8G'
    cpus 1
    time '300m'

    input:
    set val(samplename), file(thestats) // from ch_mapsummary

    output:
    file "*_mqc.txt" //into ch_multiqc_mapsum

    script:
    """
    export PATH=/opt/conda/envs/nf-core-rnaseq-1.3/bin:$PATH

    python $baseDir/bin/mito.py -m ${params.mito_name} -t $thestats > ${samplename}_mqc.txt
    """
}
