params.run = true 

process info_stats {
    tag "info_stats"
    container "singularity-rstudio-seurat-tximport"
    containerOptions = "--bind /tmp --bind /lustre"
    queue 'normal'
    time '400m'
    memory = '10G'
    cpus 2
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 1
    publishDir "${params.outdir}/info_stats/", mode: 'copy'
    scratch false 
    cache false

    when:
    params.run

    input:
    file(minimal_samples_spreadsheet_tsv)
    file(imeta_info_tsv)
    file(sync_status_tsv)
    file(metrics_summary_csvs)

    output:
    file("samples_status.pdf")
    file("samples_metainfo.tsv")
    file("samples_metainfo.xlsx")
    file("samples_status_dropped.pdf")
    file("samples_status_pre_pipelines.pdf")
    file("samples_status_dropped_pre_pipelines.pdf")
    file("samples_metainfo_supplier_name.tsv")
    file("samples_metainfo_supplier_name.xlsx")
    file("recruitment_time_course_plot.pdf")
    env(WORK_DIR), emit: work_dir_to_remove

    script:
    """
/usr/bin/Rscript $workflow.projectDir/../bin/scrna/info_stats.R $minimal_samples_spreadsheet_tsv $imeta_info_tsv $sync_status_tsv
/usr/bin/Rscript $workflow.projectDir/../bin/scrna/info_stats_supplier_name.R $minimal_samples_spreadsheet_tsv $imeta_info_tsv $sync_status_tsv
/usr/bin/Rscript $workflow.projectDir/../bin/scrna/recruitment_time_course_plot.R samples_metainfo_supplier_name.tsv

WORK_DIR=\$PWD
    """
}
