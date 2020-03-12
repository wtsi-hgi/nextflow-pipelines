params.run = true

process 'dexseq_transcript_count' {
    tag "${samplename}"
    // container "nfcore-rnaseq"
    conda "/lustre/scratch118/humgen/resources/conda/star"
    time '700m'

    errorStrategy = { task.attempt <= 2 ? 'retry' : 'ignore' }
    cpus =   {  2 * 2 * Math.min(2, task.attempt) }
    memory = {  80.GB + 40.GB * (task.attempt-1) }
    maxRetries 3
    
    publishDir "${params.outdir}/star_transcriptomesam/$samplename/", mode: 'symlink', pattern: "*.bam"
    publishDir "${params.outdir}/star_transcriptomesam/$samplename/", mode: 'symlink', pattern: "*.out"
    publishDir "${params.outdir}/star_transcriptomesam/$samplename/", mode: 'symlink', pattern: "*.tab"

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

python /lustre/scratch118/humgen/resources/rstudio_server_libs/DEXSeq/python_scripts/dexseq_count.py ${dexseq_gff} ${bam} ${samplename}.dexseq.txt
  """
}

