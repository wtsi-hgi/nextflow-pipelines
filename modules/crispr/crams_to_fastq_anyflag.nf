params.run = true
params.min_reads = 500

process crams_to_fastq_gz {
    tag "crams to fastq_gz ${samplename} ${batch}"

    //container 'nfcore-rnaseq' 
    // has samtools Version: 1.9 (using htslib 1.9)
    
    container "samtools-1.6" // has: samtools 1.6 Using htslib 1.6
    containerOptions = "--bind /lustre/scratch117/core/sciops_repository/cram_cache --bind /lustre/scratch118/core/sciops_repository/cram_cache"
    // errorStrategy 'terminate'
    errorStrategy 'retry'
    maxRetries 3
    time '300m'
    cpus 1
    memory '2G'

//  if (params.scratch) {    // This is tricky; need to get job requirements correct to ensure space exists.
//     scratch true          // At the moment we don't use this. Perhaps with a retry regime ... but a lot of fuss to solve.
//  }                        // I've left it as a reminder it's an option (svd).
    publishDir "${params.outdir}/fastq12/", mode: 'symlink'
    
    when:
    params.run
    
    input: 
    set val(samplename), val(batch), file(crams) 
    output: 
    set val(samplename), val(batch), file("*.fastq.gz")
    file('*.lostcause.txt') optional true 
    file('numreads.txt') optional true 
    script:

        // 0.7 factor below: see https://github.com/samtools/samtools/issues/494
        // This is not confirmed entirely just yet.
        // def avail_mem = task.memory == null ? '' : "${ sprintf "%.0f", 0.7 * ( task.memory.toBytes() - 2000000000 ) / task.cpus}"
    def cramfile = "${batch}.${samplename}_merged.cram"
    """
    export REF_PATH=/lustre/scratch117/core/sciops_repository/cram_cache/%2s/%2s/%s:/lustre/scratch118/core/sciops_repository/cram_cache/%2s/%2s/%s:URL=http:://sf2-farm-srv1.internal.sanger.ac.uk::8000/%s

    export PATH=/opt/conda/envs/nf-core-rnaseq-1.3/bin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
    
    samtools merge -@ ${task.cpus} -f $cramfile ${crams}

    f1=${batch}.${samplename}_1.fastq.gz
    f2=${batch}.${samplename}_2.fastq.gz
    f0=${batch}.${samplename}.fastq.gz

    numreads=\$(samtools view -c -F 0x900 $cramfile)
    if (( numreads >= ${params.min_reads} )); then
                              # -O {stdout} -u {no compression}
                              # -N {always append /1 and /2 to the read name}
                              # -F 0x900 (bit 1, 8, filter secondary and supplementary reads)
      echo -n \$numreads > numreads.txt
      samtools collate    \\
          -O -u           \\
          -@ ${task.cpus} \\
          $cramfile pfx-${samplename} | \\
      samtools fastq      \\
          -N              \\
          -F 0x900        \\
          -@ ${task.cpus} \\
          -1 \$f1 -2 \$f2  -0 \$f0 \\
          -
      sleep 2
      find . -name \"*.fastq.gz\" -type 'f' -size -160k -delete
    else
      echo -e "${samplename}\\tcram\\tlowreads" > ${batch}.${samplename}.lostcause.txt
    fi
    """
}

// export REF_CACHE=/lustre/scratch117/core/sciops_repository/cram_cache/%2s/%2s/%s:/lustre/scratch118/core/sciops_repository/cram_cache/%2s/%2s/%s
