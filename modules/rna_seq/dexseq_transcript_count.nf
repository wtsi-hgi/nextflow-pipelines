params.run = true

process 'dexseq_transcript_count' {
    tag "${samplename}"
    // container "nfcore-rnaseq"
    conda "/lustre/scratch118/humgen/resources/conda/star"
    time '700m'

    errorStrategy = { task.attempt <= 2 ? 'retry' : 'ignore' }
    cpus =  1 // {  2 * 2 * Math.min(2, task.attempt) }
    memory = 10.GB //{  80.GB + 40.GB * (task.attempt-1) }
    maxRetries 2
    
    publishDir "${params.outdir}/basic_dexseq_transcript_count/", mode: 'rellink', pattern: "${samplename}.dexseq.txt"

  input:
    set val(samplename), file(bam) //from ch_star // _reads_only
    file(dexseq_gff)

    when:
    params.run
    
  output:
    set val(samplename), file ("${samplename}.dexseq.txt")

  script:

  """
  export PATH=/lustre/scratch118/humgen/resources/conda/star/bin:\$PATH

python /lustre/scratch118/humgen/resources/rstudio_server_libs/DEXSeq/python_scripts/dexseq_count.py \\
    -p yes \\
    -s reverse \\
    -f bam \\
    -a 10 \\
    -r pos \\
    ${dexseq_gff} ${bam} ${samplename}.dexseq.txt
  """
}
// https://bioconductor.org/packages/3.11/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html

