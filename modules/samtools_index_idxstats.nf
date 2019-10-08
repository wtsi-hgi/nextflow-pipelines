params.run = true

process samtools_index_idxstats {
    tag "${samplename}"
    container "nfcore-rnaseq"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    when:
    params.run
    
    publishDir "${params.outdir}/idxstats/", mode: 'symlink', pattern: "*.idxstats"

    input:
    set val(aligner), val(samplename), file(thebam) //from ch_indexbam

    output:
    set val(samplename), file("*.idxstats") //into ch_mapsummary

    script:
    """
    export PATH=/opt/conda/envs/nf-core-rnaseq-1.3/bin:$PATH

    samtools index $thebam
    samtools idxstats $thebam > ${samplename}.idxstats
    rm ${thebam}.bai
    """
}
