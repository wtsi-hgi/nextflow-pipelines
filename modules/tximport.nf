params.run = true 

process tximport {
    tag "tximport $samplename"
    memory = '80G'
    container "singularity-rstudio-seurat-tximport"
    time '400m'
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    
    publishDir "${params.outdir}/tximport", mode: 'symlink'

    when:
    params.run

    input:
    file (quant_sf_all) 

    output:
    //file "${samplename}.quant.sf" // into ch_salmon_trans
    //file "${samplename}.quant.genes.sf" //into ch_salmon_genes
    // file "my_outs/${samplename}" optional true // into ch_alignment_logs_salmon

    script:
    """
    Rscript $workflow.projectDir/bin/tximport.R
    """

    // TODO: prepare columns for merging; extract correct column and transpose (paste) it.
    // Include the row names so merger can check identity.
    // The merge step will concatenate the rows and re-transpose to obtain final result.
}
