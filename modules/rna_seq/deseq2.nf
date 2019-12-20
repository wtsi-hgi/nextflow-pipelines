params.run = true 

process deseq2 {
    tag "$matrix_deseq2 $deseq2_tsv"
    memory = '80G'
    container "singularity-rstudio-seurat-tximport"
    containerOptions = "--bind /tmp --bind /lustre"
    time '400m'
    cpus 1
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    
    publishDir "${params.outdir}/DESeq2/", mode: 'symlink'

    when:
    params.run

    input:
    file(deseq2_tsv)
    file(txi_gene_counts_csv)
    file(txi_transcript_counts_csv)
    file(txi_lengthScaledTPM_gene_counts_csv)

    output:
    file("deseq2.rdata")

    script:
    """
    /usr/bin/Rscript $workflow.projectDir/../bin/rna_seq/deseq2.R $deseq2_tsv $txi_gene_counts_csv $txi_transcript_counts_csv $txi_lengthScaledTPM_gene_counts_csv
    """
}
