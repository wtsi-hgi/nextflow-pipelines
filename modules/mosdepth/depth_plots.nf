params.run = true

process depth_plots {
    tag ""
    memory = '4G'
    time '240m'
    cpus 1
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    maxForks 12
    publishDir "${params.outdir}/depth_plots/", mode: 'copy', pattern: "samples.tsv", overwrite: true

    when:
    params.run

    input: 
    val study_id

    output: 
    tuple val(study_id), file('samples.tsv'), emit: samples_tsv

    script:
    """
    python3 $workflow.projectDir/../bin/mosdepth/depth_plots.py
    """
}
