params.run = true

process mosdepth_from_local {
    tag "${study_id}"
    memory = '4G'
    time '240m'
    cpus 1
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    maxForks 12
    publishDir "${params.outdir}/samples/study_id_${study_id}/", mode: 'copy', pattern: "samples.tsv", overwrite: true
    publishDir "${params.outdir}/samples/study_id_${study_id}/", mode: 'copy', pattern: "samples_noduplicates.tsv", overwrite: true

    when:
    params.run

    input: 
    val study_id

    output: 
    tuple val(study_id), file('samples.tsv'), emit: samples_tsv
    tuple val(study_id), file('samples_noduplicates.tsv'), emit: samples_noduplicates_tsv

    script:
    """
    """
}
// awk removes duplicates in the sense that one sanger sample can have several run_id
