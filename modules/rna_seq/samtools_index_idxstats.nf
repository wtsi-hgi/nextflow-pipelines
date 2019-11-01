params.run = true

process samtools_index_idxstats {
    tag "${samplename}"
    container "nfcore-rnaseq"
    memory = '8G'
    cpus 1
    time '300m'
    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5
    publishDir "${params.outdir}/idxstats/", mode: 'symlink', pattern: "*.idxstats"

    when:
    params.run

    input:
    set val(aligner), val(samplename), file(thebam) //from ch_indexbam

    output:
    set val(samplename), file("*.idxstats") //into ch_mapsummary

    script:
    """
    export PATH=/opt/conda/envs/nf-core-rnaseq-1.3/bin:\$PATH

    samtools index $thebam
    samtools idxstats $thebam > ${samplename}.idxstats
    rm ${thebam}.bai
    """
}
