params.run = true
params.min_reads = 500

process crams_to_fastq_gz {
    tag "crams to fastq_gz ${samplename}"

    //container 'nfcore-rnaseq' 
    // has samtools Version: 1.9 (using htslib 1.9)
    
    container "samtools-1.6" // has: samtools 1.6 Using htslib 1.6
    containerOptions = "--bind /lustre/scratch117/core/sciops_repository/cram_cache --bind /lustre/scratch118/core/sciops_repository/cram_cache"
    // errorStrategy 'terminate'
    errorStrategy 'retry'
    maxRetries 6
    time '400m'
    cpus 1
    memory '4G'

//  if (params.scratch) {    // This is tricky; need to get job requirements correct to ensure space exists.
//     scratch true          // At the moment we don't use this. Perhaps with a retry regime ... but a lot of fuss to solve.
//  }                        // I've left it as a reminder it's an option (svd).
    publishDir "${params.outdir}/fastq12/", mode: 'symlink'
    
    when:
    params.run
    
    input: 
        set val(samplename), file(crams) 
    output: 
        set val(samplename), file("${samplename}_1.fastq.gz"), file("${samplename}_2.fastq.gz") optional true
        file('*.lostcause.txt') optional true 
        file('numreads.txt') optional true 
    script:

        // 0.7 factor below: see https://github.com/samtools/samtools/issues/494
        // This is not confirmed entirely just yet.
        // def avail_mem = task.memory == null ? '' : "${ sprintf "%.0f", 0.7 * ( task.memory.toBytes() - 2000000000 ) / task.cpus}"
    def cramfile = "${samplename}_merged.cram"
    """
    export REF_PATH=/lustre/scratch117/core/sciops_repository/cram_cache/%2s/%2s/%s:/lustre/scratch118/core/sciops_repository/cram_cache/%2s/%2s/%s:URL=http:://sf2-farm-srv1.internal.sanger.ac.uk::8000/%s

    export PATH=/opt/conda/envs/nf-core-rnaseq-1.3/bin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
    
    samtools merge -@ ${task.cpus} -f $cramfile ${crams}

    f1=${samplename}_1.fastq.gz
    f2=${samplename}_2.fastq.gz

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
          -1 \$f1 -2 \$f2 \\
          -
      sync \$f1 \$f2          # this line and next to tackle k8s weirdness (see k8s)
      sleep 1
    else
      echo -e "${samplename}\\tcram\\tlowreads" > ${samplename}.lostcause.txt
    fi
    """
}

// export REF_CACHE=/lustre/scratch117/core/sciops_repository/cram_cache/%2s/%2s/%s:/lustre/scratch118/core/sciops_repository/cram_cache/%2s/%2s/%s
