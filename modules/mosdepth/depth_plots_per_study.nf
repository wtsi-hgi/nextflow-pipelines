params.run = true

process depth_plots_per_study {
    tag "$study_id"
    memory = '4G'
    time '240m'
    cpus 1
    conda '/lustre/scratch118/humgen/resources/conda_envs/mosdepth'
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    maxForks 12
    publishDir "${params.outdir}/per_study_depth_plots/${study_id}/", mode: 'copy', pattern: "*.png", overwrite: true

    when:
    params.run

    input: 
    set val(study_id), file(mosdepth_region_files)

    output: 
    tuple val(study_id), file('*.png'), emit: study_plots

    script:
    """
    export PATH=/lustre/scratch118/humgen/resources/conda_envs/mosdepth/bin:\$PATH
    python3 $workflow.projectDir/../bin/mosdepth/depth_plots_per_study.py
    """
}
