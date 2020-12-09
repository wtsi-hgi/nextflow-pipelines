params.run = true

process 'fastqc' {
    memory = '5G'
    cpus 2
    tag "fastqc ${samplename}"
    container "nfcore-rnaseq"
    errorStrategy 'retry'
    maxRetries 3
    time '120m'
    // Singularity nfcore-rnaseq.img:~> fastqc --version
    // FastQC v0.11.8
    
    // singularity pull --name nfcore-rnaseq.img docker://nfcore/rnaseq
    // have fastqc version FastQC v0.11.8, was pulled Thursday May 16th 2019
    // publishDir "${params.outdir}/STAR_2pass_bams/${samplename}/", mode: 'copy'
    publishDir "${params.outdir}/fastqc/", mode: 'symlink',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when:
    params.run 

    input:
    set val(samplename), file(reads)
    
    output:
    set file("*fastqc.zip"), file("*fastqc.html")

  script:
  """
  export PATH=/opt/conda/envs/nf-core-rnaseq-1.3/bin:\$PATH

  fastqc -t ${task.cpus} -q $reads
  """
}
