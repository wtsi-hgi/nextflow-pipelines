params.run = true

process 'leafcutter_bam2junc' {
    tag "${samplename}"
    container "leafcutter"
    memory = '8G'
    cpus 1
    time '320m'
    errorStrategy { task.attempt <= 8 ? 'retry' : 'ignore' }
    maxRetries 8
    
    publishDir "${params.outdir}/leafcutter/bam2junc", mode: 'symlink', pattern: "*.junc"
    // publishDir "${params.outdir}/leafcutter/bam2junc", mode: 'copy', pattern: "*.bam.bed"

  input:
    set val(samplename), file (bamfile), file (baifile) //from star_aligned_with_bai

    when:
    params.run
    
  output:
    file ('*.junc')

  script:

  """
  export PATH=/home/leafcutter/scripts:/home/leafcutter/clustering:\$PATH

  echo Converting ${bamfile} to ${samplename}.junc
  sh /home/leafcutter/scripts/bam2junc.sh ${bamfile} ${samplename}.junc
  """
}
