params.run = true

process 'cellsnp' {
    tag "cellSNP $samplename"
    container "single_cell"

    errorStrategy = { task.attempt <= 4 ? 'retry' : 'ignore' }
    cpus = 2 
    memory = {  200.GB + 50.GB * (task.attempt-1) }
    maxRetries 4
    queue 'long'
    time '1400m'

    publishDir "${params.outdir}/cellsnp/", mode: 'symlink'

    when:
    params.run 

    input:
    set val(samplename), file(bam_file), file(bai_file), file(barcodes_tsv_gz)
    file(region_vcf)
    
    output:
    set val(samplename), file("cellsnp_${samplename}") 

  script:
   """
export PATH=/opt/conda/envs/conda_env/bin:/opt/conda/bin:\$PATH
zcat ${barcodes_tsv_gz} > barcodes.txt
cellSNP -s ${bam_file} -b barcodes.txt -O cellsnp_${samplename} -R ${region_vcf} -p 20 --minMAF 0.1 --minCOUNT 20
   """
}
// hg19: genome1K.phase3.SNP_AF5e2.chr1toX.hg19.vcf.gz (http://ufpr.dl.sourceforge.net/project/cellsnp/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz)


// https://github.com/single-cell-genetics/cellSNP
// Mode 1: pileup a list of SNPs for single cells in a big BAM/SAM file
// Require: a single BAM/SAM file, e.g., from cellranger, a VCF file for a list of common SNPs.
// This mode is recommended comparing to mode 2, if a list of common SNP is known, e.g., human (see Candidate SNPs below)
// As shown in the above command line, we recommend filtering SNPs with <20UMIs or <10% minor alleles for downstream donor deconvolution, by adding --minMAF 0.1 --minCOUNT 20
