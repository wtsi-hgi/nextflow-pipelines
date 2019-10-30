params.run = true

process 'fastx_trimmer' {
    tag "fastx_trimmer $samplename"
    memory = '10G'
    time '120m'
    container "fastxtoolkit"  //.img// FASTX Toolkit 0.0.13
    containerOptions = "--bind /lustre"
    cpus 1
    maxForks 100
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/fastxtoolkit_trimmer/", mode: 'symlink'

    when:
    params.run 

    input:
    set val(samplename), val(batch), val(start_trim), file(fastq_gz_input)
    
  output:
    set val(samplename), val(batch), file("${samplename}.${batch}.fastq.gz")

  script:
   """
   gunzip -c ${fastq_gz_input} | /usr/local/fastxToolkit-0.0.13/bin/fastx_trimmer -v -f ${start_trim} -z -Q33 -o tmp.fastq.gz

   rm -f ${samplename}.${batch}.fastq.gz
   mv tmp.fastq.gz ${samplename}.${batch}.fastq.gz
   """
}
