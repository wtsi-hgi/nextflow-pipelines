params.run = true

process 'leafcutter_bam2junc' {
    tag "${samplename}"
    queue "long"
    time '2800m'
    container "leafcutter"


    errorStrategy = { task.attempt <= 4 ? 'retry' : 'ignore' }
    cpus =   {  2 * 2 * Math.min(2, task.attempt) }
    memory = {  10.GB + 20.GB * (task.attempt-1) }
    maxRetries 4
    
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

process 'leafcutter_bam2junc_regtools' {
    tag "${samplename}"
    queue "long"
    time '2800m'
    container "leafcutter"


    errorStrategy = { task.attempt <= 4 ? 'retry' : 'ignore' }
    cpus =   {  2 * 2 * Math.min(2, task.attempt) }
    memory = {  10.GB + 20.GB * (task.attempt-1) }
    maxRetries 4
    
    publishDir "${params.outdir}/leafcutter_regtools/bam2junc", mode: 'symlink', pattern: "*.junc"
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
  # pre regootls method: sh /home/leafcutter/scripts/bam2junc.sh ${bamfile} ${samplename}.junc
         
  regtools junctions extract -s 1 -a 8 -m 50 -M 500000 ${bamfile} -o ${samplename}.junc    
  """
}
