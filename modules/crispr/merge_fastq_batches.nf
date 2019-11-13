params.run = true

process 'merge_fastq_batches' {
    tag "merge batch $samplename"
    memory = '15G'
    time '120m'
    cpus 1
    maxForks 100
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/merge_fastq_batches/", mode: 'symlink'

    when:
    params.run 

    input:
    set val(samplename), val(batches), file(fastqs)
    
  output:
    set val(samplename), file("${samplename}.fastq.gz")

  script:
   """
   cat *.fastq.gz > ${samplename}.fastq.gz
   """
}
