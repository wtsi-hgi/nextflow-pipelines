params.run = true 
params.ensembl_lib = "Ensembl 91 EnsDb"

process deseq2 {
    tag "$deseq2_tsv"
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
    file (quant_sf_files)  // from collect()
    file(deseq2_tsv)

    output:
    file("deseq2.rdata")

    script:
    """
    ls . | grep .quant.sf\$ > fofn_quantfiles.txt

    /usr/bin/Rscript $workflow.projectDir/../bin/rna_seq/deseq2.R  \"$params.ensembl_lib\" fofn_quantfiles.txt $deseq2_tsv
    """
}
