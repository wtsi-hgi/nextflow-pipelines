params.run = true 
params.ensembl_lib = "Ensembl 91 EnsDb"

process seurat {
    tag "seurat $params.ensembl_lib"
    container "singularity-rstudio-seurat-tximport"
    containerOptions = "--bind /tmp --bind /lustre"
    time '900m'
    memory = '80G'
    cpus 2
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    
    publishDir "${params.outdir}/seurat", mode: 'symlink'

    when:
    params.run

    input:
    val(samplename), file(cellranger_matrix_dir), val(raw_filtered), file(metrics_summary_csv)

    output:
    file("${samplename}_seurat_image.rdata")

    script:
    """
    /usr/bin/Rscript $workflow.projectDir/../bin/scrna/seurat.R $samplename $cellranger_matrix_dir $raw_filtered $metrics_summary_csv
    """
}
