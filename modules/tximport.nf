params.run = true 
params.ensembl_lib = "Ensembl 91 EnsDb"

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
    file (quantf_list_txt)
    file (quant_sf_files)  // from collect()

    output:
    //file "${samplename}.quant.sf" // into ch_salmon_trans
    //file "${samplename}.quant.genes.sf" //into ch_salmon_genes
    // file "my_outs/${samplename}" optional true // into ch_alignment_logs_salmon

    script:
    """
    /usr/bin/Rscript --vanilla $workflow.projectDir/bin/tximport.R \"$params.ensembl_lib\" \"$quantf_list_txt\" 
    """

    // TODO: prepare columns for merging; extract correct column and transpose (paste) it.
    // Include the row names so merger can check identity.
    // The merge step will concatenate the rows and re-transpose to obtain final result.
}
