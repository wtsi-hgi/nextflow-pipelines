params.run = true 

process seurat {
    tag "seurat $samplename $raw_filtered"
    container "singularity-rstudio-seurat-tximport"
    containerOptions = "--bind /tmp --bind /lustre"
    time '900m'
    memory = '80G'
    cpus 2
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    
    publishDir "${params.outdir}/seurat/$raw_filtered/", mode: 'symlink'

    when:
    params.run

    input:
    val(samplename), file(cellranger_matrix_dir), val(raw_filtered), file(metrics_summary_csv)

    output:
    val(samplename), val(raw_filtered), file("${samplename}_TSNEPlot.pdf")
    val(samplename), val(raw_filtered), file("${samplename}_stats.tsv"), file("${samplename}_stats.xlsx")
    val(samplename), val(raw_filtered), file("${samplename}_clusters_markers_FindAllMarkers.xlsx")
    val(samplename), val(raw_filtered), file("${samplename}_seurat_image.rdata")

    script:
    """
    /usr/bin/Rscript $workflow.projectDir/../bin/scrna/seurat.R $samplename $cellranger_matrix_dir $raw_filtered $metrics_summary_csv
    """
}
