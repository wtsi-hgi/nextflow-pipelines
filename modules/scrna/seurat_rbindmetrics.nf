params.run = true 
params.ensembl_lib = "Ensembl 91 EnsDb"

process seurat {
    tag "seurat $params.ensembl_lib"
    memory = '80G'
    container "singularity-rstudio-seurat-seurat"
    containerOptions = "--bind /tmp --bind /lustre"
    time '400m'
    cpus 1
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    
    publishDir "${params.outdir}/seurat", mode: 'symlink'

    when:
    params.run

    input:
    file (quant_sf_files)  // from collect()

    output:
    file("fofn_quantfiles.txt")
    file("txi_gene_counts.csv")
    file("txi_transcript_counts.csv")
    file("txi_lengthScaledTPM_gene_counts.csv")
    file("seurat.rdata")
    //file "${samplename}.quant.sf" // into ch_salmon_trans
    //file "${samplename}.quant.genes.sf" //into ch_salmon_genes
    // file "my_outs/${samplename}" optional true // into ch_alignment_logs_salmon

    script:
    """
    ls . | grep .quant.sf\$ > fofn_quantfiles.txt

    /usr/bin/Rscript $workflow.projectDir/../bin/scrna/seurat.R \"$params.ensembl_lib\" fofn_quantfiles.txt 
    """

    // TODO: prepare columns for merging; extract correct column and transpose (paste) it.
    // Include the row names so merger can check identity.
    // The merge step will concatenate the rows and re-transpose to obtain final result.
}
