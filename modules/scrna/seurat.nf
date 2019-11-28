params.run = true 

process seurat {
    tag "seurat $samplename $raw_filtered"
    container "singularity-rstudio-seurat-tximport"
    containerOptions = "--bind /tmp --bind /lustre"
    queue 'long'
    time '1400m'
    memory = '80G'
    cpus 2
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/seurat/$raw_filtered/", mode: 'symlink'
    scratch false 

    when:
    params.run

    input:
    set val(samplename), file(cellranger_matrix_dir), val(raw_filtered), file(metrics_summary_csv)

    output:
    tuple( val(samplename), val(raw_filtered), file("${samplename}_${raw_filtered}_TSNEPlot.pdf")), emit: tsneplot_pdf
    set val(samplename), val(raw_filtered), file("${samplename}_${raw_filtered}_stats.tsv"), file("${samplename}_${raw_filtered}_stats.xlsx"), emit: stats_xslx
    set val(samplename), val(raw_filtered), file("${samplename}_${raw_filtered}_clusters_markers_FindAllMarkers.xlsx"), emit: diffexp_xlsx
    set val(samplename), val(raw_filtered), file("${samplename}_${raw_filtered}_seurat_image.rdata"), emit: seurat_rdata

    script:
    """
    /usr/bin/Rscript $workflow.projectDir/../bin/scrna/seurat.R $samplename $cellranger_matrix_dir $raw_filtered $metrics_summary_csv
    """
}
