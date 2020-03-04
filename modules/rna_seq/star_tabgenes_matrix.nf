params.run = true 

process star_tabgenes_matrix {
    tag "star_tabgenes_matrix"
    memory = '30G'
    container "singularity-rstudio-seurat-tximport"
    containerOptions = "--bind /tmp --bind /lustre"
    time '400m'
    cpus 1
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    
    publishDir "${params.outdir}/heatmap/", mode: 'symlink'

    when:
    params.run

    input:
    file (tagenes_files)

    output:
    tuple val("star"), file("star_tabgenes_matrix.tsv"), emit: star_matrix

    script:
    """
    /usr/bin/Rscript $workflow.projectDir/../bin/rna_seq/star_tabgenes_matrix.R
    """
}
