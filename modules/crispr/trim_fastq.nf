params.run = true

process 'fastx_trimmer' {
    tag "fastx_trimmer $samplename"
    memory = '10G'
    time '120m'
    container "fastxtoolkit"  //.img// FASTX Toolkit 0.0.13
    cpus 1
    maxForks 100
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/merge_fastq_batches/", mode: 'symlink'

    when:
    params.run 

    input:
    set val(samplename), val(star_trim), file(fastq_gz_input)
    
  output:
    set val(samplename), file("${samplename}.fastq.gz")

  script:
   """
   /usr/local/fastxToolkit-0.0.13/bin/fastx_trimmer -f ${start_trim} -z -i ${fastq_gz_input} -o tmp.fastq.gz

   rm -f ${samplename}.fastq.gz
   mv tmp.fastq.gz ${samplename}.fastq.gz
   """
}
