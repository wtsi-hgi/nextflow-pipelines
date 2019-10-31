params.run = true

process 'star_2pass_1st_pass' {
    tag "1st pass ${samplename}"
    container "nfcore-rnaseq"
    time '600m'

    errorStrategy = { task.attempt <= 2 ? 'retry' : 'ignore' }
    cpus =   {  2 * 2 * Math.min(2, task.attempt) }
    memory = {  80.GB + 40.GB * (task.attempt-1) }
    maxRetries 3
    
    publishDir "${params.outdir}/star_pass2_1stpass/$samplename", mode: 'symlink', pattern: "*.out"
    publishDir "${params.outdir}/star_pass2_1stpass/$samplename", mode: 'symlink', pattern: "*.tab"
    
    publishDir "${params.outdir}/star_pass2_1stpass_multiqc/", mode: 'copy',
        saveAs: { filename ->
            if (filename ==~ /.*\.out\.tab/) "STARcounts/$filename"
            else if (filename.indexOf(".bam") == -1) "STARlogs/$filename"
            else null
        }

  input:
    set val(samplename), file(reads) //from ch_star // _reads_only
    file genomeDir //from ch_star_index.collect()
    file gtf //from ch_gtf_star.collect()

    when:
    params.run
    
  output:
    set val(samplename), file("*Log.final.out")
    file "*.SJ.out.tab"
    set file("*.Log.final.out"), file("*.Log.out"), file("*.progress.out") //into ch_alignment_logs_star

  script:

  """
  export PATH=/opt/conda/envs/nf-core-rnaseq-1.3/bin:\$PATH

    # first pass 
    STAR --genomeDir ${genomeDir} \\
        --sjdbGTFfile $gtf \\
        --readFilesIn $reads --readFilesCommand zcat \\
        --runThreadN ${task.cpus} \\
        --outSAMtype BAM Unsorted \\
        --outFileNamePrefix ${samplename}.

  rm *.bam
  """
}
