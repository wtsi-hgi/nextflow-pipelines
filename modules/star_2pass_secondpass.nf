params.run = true

process 'star_2pass_2nd_pass' {
    tag "2nd pass ${samplename}"
    container "nfcore-rnaseq"
    time '600m'

    errorStrategy = { task.exitStatus == 130 && task.attempt <= 2 ? 'retry' : 'ignore' }
    cpus =   {  2 * 2 * Math.min(2, task.attempt) }
    memory = {  80.GB + 40.GB * (task.attempt-1) }
    maxRetries 3
    
    publishDir "${params.outdir}/star_pass2_2ndpass/$filename", mode: 'symlink', pattern: "*.bam"
    publishDir "${params.outdir}/star_pass2_2ndpass/$filename", mode: 'symlink', pattern: "*.bam.bai"
    publishDir "${params.outdir}/star_pass2_2ndpass/$filename", mode: 'symlink', pattern: "*.out"
    publishDir "${params.outdir}/star_pass2_2ndpass/$filename", mode: 'symlink', pattern: "*.tab"
    
    publishDir "${params.outdir}/star_pass2_2ndpass_multiqc/", mode: 'copy',
        saveAs: { filename ->
            if (filename ==~ /.*\.ReadsPerGene\.out\.tab/) "STARcounts/$filename"
            else if (filename.indexOf(".bam") == -1) "STARlogs/$filename"
            else null
        }

  input:
    set val(samplename), file(reads) //from ch_star // _reads_only
    file genomeDir //from ch_star_index.collect()
    file gtf //from ch_gtf_star.collect()
    file filtered_tab // merged and filtered junctions from 1st pass

    when:
    params.run
    
  output:
    set val(samplename), file ('*.bam'), file ('*.bai') //into star_aligned_with_bai
    set val(samplename), file("*Log.final.out"), file ('*.bam') //into star_aligned
    // file "*.ReadsPerGene.out.tab" into ch_merge_starcounts
    file "*.out" //into ch_alignment_logs_star
    file "*.SJ.out.tab"
    file "*.ReadsPerGene.out.tab" //into ch_merge_starcounts

  script:

  """
  export PATH=/opt/conda/envs/nf-core-rnaseq-1.3/bin:$PATH

  # 2nd pass
  STAR --genomeDir ${genomeDir} \
       --readFilesIn ${reads} \
       --runThreadN ${task.cpus} \
       --readFilesCommand zcat \
       --limitSjdbInsertNsj 10000000 \
       --sjdbFileChrStartEnd ${filtered_tab} \
       --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix ${samplename}.

  # Index the BAM file
  samtools index ${samplename}.Aligned.sortedByCoord.out.bam
  """
}
