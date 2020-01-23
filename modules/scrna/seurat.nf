params.run = true 

process seurat {
    tag "seurat $samplename $raw_filtered"

    //// FCE
    disk '100 GB'
    scratch '/tmp'
    stageInMode 'symlink'
    stageOutMode 'rsync'
    cpus = 8
    time '8000m'
    tontainer "single_cell"
    containerOptions = "--bind /"
    memory = {  100.GB + 50.GB * (task.attempt-1) }
    ////// FCE 

    /// farm
    // containerOptions = "--bind /tmp --bind /lustre"
    // queue 'long'
    // time '1400m'
    // memory = '80G'
    // cpus 2
    // errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
    // maxRetries 2
    // scratch false
    /// farm

    container "singularity-rstudio-seurat-tximport"
    publishDir "${params.outdir}/seurat/$raw_filtered/", mode: 'symlink'

    when:
    params.run

    input:
    set val(samplename), file(cellranger_matrix_dir), val(raw_filtered), file(metrics_summary_csv)

    output:
    tuple val(samplename), val(raw_filtered), file("${samplename}_${raw_filtered}_TSNEPlot.pdf"), emit: tsneplot_pdf
    tuple val(samplename), val(raw_filtered), file("${samplename}_${raw_filtered}_stats.tsv"), file("${samplename}_${raw_filtered}_stats.xlsx"), emit: stats_xslx
    tuple val(samplename), val(raw_filtered), file("${samplename}_${raw_filtered}_clusters_markers_FindAllMarkers.xlsx"), emit: diffexp_xlsx
    tuple val(samplename), val(raw_filtered), file("${samplename}_${raw_filtered}_seuratimage.rdata"), emit: seurat_rdata

    script:
    """
    /usr/bin/Rscript $workflow.projectDir/../bin/scrna/seurat.R $samplename $cellranger_matrix_dir $raw_filtered $metrics_summary_csv
    """
}
