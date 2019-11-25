params.run = true

process 'vireo' {
    tag "vireo $samplename $run_id"
    conda "/software/hgi/installs/anaconda3/envs/hgi_singlecell"
    memory = '3G'
    time '120m'
    cpus 1
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/vireo/", mode: 'symlink'

    when:
    params.run 

    input:
    set val(samplename), val(bam_file), val(bai_file)
    
    output:
    set val(samplename), file("cellsnp_${samplename}") 

  script:
   """
cellSNP -s ${bam_file} -b $BARCODE -O cellsnp_${samplename} -R $REGION_VCF -p 20 --minMAF 0.1 --minCOUNT 20
   """
}
// https://github.com/single-cell-genetics/cellSNP
// Mode 1: pileup a list of SNPs for single cells in a big BAM/SAM file
// Require: a single BAM/SAM file, e.g., from cellranger, a VCF file for a list of common SNPs.
// This mode is recommended comparing to mode 2, if a list of common SNP is known, e.g., human (see Candidate SNPs below)
// As shown in the above command line, we recommend filtering SNPs with <20UMIs or <10% minor alleles for downstream donor deconvolution, by adding --minMAF 0.1 --minCOUNT 20
