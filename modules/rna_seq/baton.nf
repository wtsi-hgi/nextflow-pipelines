params.run = true
params.dropqc = ""

process baton_study_id {
    tag "${study_id}"
    memory = '4G'
    time '240m'
    cpus 1
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    maxForks 12
    publishDir "${params.outdir}/", mode: 'copy', pattern: "samples.tsv", overwrite: true

    when:
    params.run

    input: 
    val study_id

    output: 
    tuple val(study_id), file('samples.tsv'), emit: samples_tsv

    script:
    """
    bash $workflow.projectDir/../bin/rna_seq/baton.sh ${study_id}
    """
}
