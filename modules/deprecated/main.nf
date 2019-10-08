Channel.fromPath('/lustre/scratch115/projects/interval_rna/inputs/*.cram').
map{ it -> [ it.toString().replaceAll(~/.*\/(.*).cram/, "\$1"), it ] }.
groupTuple(). //take(4).
take(4).
into{ch_cram_files;ch_cram_files1; ch_cram_files2}
// take(2). // 372 cram files in /warehouse/cellgeni/tic-109/batch7/data/study5591-tic329.tar
// [INT_RNA7710970, /lustre/scratch115/projects/interval_wgs/nextflow/interval_rnaseq/inputs/tar/INT_RNA7710970.cram]
// [INT_RNA7710909, /lustre/scratch115/projects/interval_wgs/nextflow/interval_rnaseq/inputs/tar/INT_RNA7710909.cram]
ch_cram_files1.view()


process crams_to_fastqgz {
    tag "crams to fastq_gz ${samplename}"

    //container 'nfcore-rnaseq' 
    // has samtools Version: 1.9 (using htslib 1.9)
    
    container "samtools-1.6" // has: samtools 1.6 Using htslib 1.6
    containerOptions = "--bind /lustre/scratch117/core/sciops_repository/cram_cache --bind /lustre/scratch118/core/sciops_repository/cram_cache"
    // errorStrategy 'terminate'
    errorStrategy 'retry'
    maxRetries 3
    time '60m'
    cpus 1
    memory '2G'

//  if (params.scratch) {    // This is tricky; need to get job requirements correct to ensure space exists.
//     scratch true          // At the moment we don't use this. Perhaps with a retry regime ... but a lot of fuss to solve.
//  }                        // I've left it as a reminder it's an option (svd).
    publishDir "${params.outdir}/fastq12/", mode: 'copy'

    input: 
        set val(samplename), file(crams) from ch_cram_files
    output: 
        set val(samplename), file("${samplename}_?.fastq.gz") optional true into ch_fastqs_irods
        file('*.lostcause.txt') optional true into ch_lostcause_cram
        file('numreads.txt') optional true into ch_numreads_crams
    script:

        // 0.7 factor below: see https://github.com/samtools/samtools/issues/494
        // This is not confirmed entirely just yet.
        // def avail_mem = task.memory == null ? '' : "${ sprintf "%.0f", 0.7 * ( task.memory.toBytes() - 2000000000 ) / task.cpus}"
    def cramfile = "${samplename}_merged.cram"
    """
    export REF_PATH=/lustre/scratch117/core/sciops_repository/cram_cache/%2s/%2s/%s:/lustre/scratch118/core/sciops_repository/cram_cache/%2s/%2s/%s:URL=http:://sf2-farm-srv1.internal.sanger.ac.uk::8000/%s
    export REF_CACHE=/lustre/scratch117/core/sciops_repository/cram_cache/%2s/%2s/%s:/lustre/scratch118/core/sciops_repository/cram_cache/%2s/%2s/%s

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

ch_fastqs_irods
    .into{ ch_rnaseq; ch_rnaseq2; ch_fastqc; ch_mixcr; ch_bracer; ch_tracer }




// /////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////

params.min_reads    = 500
params.min_pct_aln  = 5
params.pe_suffix_pattern  = '_{1,2}.fastq.gz'
params.se_suffix    = '.fastq.gz'
params.runtag = 'interval_rna'

// done for salmon14 from docker,
// in /lustre/scratch115/projects/interval_wgs/nextflow/salmon14_index
params.make_salmon_index = false

params.run_qc = true
params.make_star_first_pass = false
params.run_star = false

params.run_leafcutter_bam2junc = false
params.run_leafcutter_clustering = false
// if false set ch_first_pass_genome = Channel.fromPath('/lustre/scratch115/projects/interval_wgs/nextflow/hg_ibob_gmpr_rnaseq_5796/results/star2pass/genome2pass')

params.run_multiqc = false
params.run_salmon = false
params.run_salmon_mergecounts = false
params.save_bam = false

params.genome = 'GRCh38'
params.cdna = params.genome ? params.genomes[ params.genome ].cdna ?: false : false
Channel.fromPath(params.cdna).ifEmpty { exit 1, "cdna annot file not found" } .into { ch_cdna_salmon; ch_cdna_transgene }

process makeSalmonIndex {
    tag "$fasta"
    publishDir "/lustre/scratch115/projects/interval_wgs/nextflow/salmon14_index", mode: 'copy'
    containerOptions = "--bind /lustre/scratch115/projects/interval_wgs/nextflow/salmon14_index"
    container "salmon"
    cpus 2
    memory '10G'

    when:
    params.make_salmon_index

    input:
    file fasta from ch_cdna_salmon

    output:
    file "salmon"

    script:
    """
    mkdir salmon
    salmon index        \\
        -t $fasta       \\
        -p ${task.cpus} \\
        -i salmon
    """
}
process makeTransGeneMatrix {
    tag "$fasta"
    publishDir "/lustre/scratch115/projects/interval_wgs/nextflow/salmon14_index", mode: 'copy'
    containerOptions = "--bind /lustre/scratch115/projects/interval_wgs/nextflow/salmon14_index"
    container "salmon"
    cpus 2
    memory '10G'

    when:
    params.make_salmon_index

    input:
    file fasta from ch_cdna_transgene

    output:
    file "trans_gene*.txt"

    shell:
    '''
    perl -ne 'if (/^>(\\w+)(?:\\.\\d+)\\s+.*?gene:(\\w+)/){print "$1\\t$2\\n"}elsif(/^>(ERCC\\S+)/){print"$1\\t$1-gene\\n"}' \\
      !{fasta} > trans_gene.txt
    if (( $(grep -c ENST trans_gene.txt) < 1000 )); then
       echo 'Not enough Ensembl transcripts. This test is present because this script section is ugly.'
       echo 'Currently it makes an effort to recognise ERCC information.'
       echo 'If you want to run a gencode genome, update this section, make this file aware of gencode/Ensembl distinction'
       false
    fi
    '''
}


// /////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////



process 'fastqc' {
    memory = '5G'
    cpus 2
    tag "fastqc ${samplename}"
    container "nfcore-rnaseq"
    errorStrategy 'retry'
    maxRetries 3
    time '60m'
    // Singularity nfcore-rnaseq.img:~> fastqc --version
    // FastQC v0.11.8
    
    // singularity pull --name nfcore-rnaseq.img docker://nfcore/rnaseq
    // have fastqc version FastQC v0.11.8, was pulled Thursday May 16th 2019
    // publishDir "${params.outdir}/STAR_2pass_bams/${samplename}/", mode: 'copy'
    publishDir "${params.outdir}/fastqc/", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when:
    params.run_qc 

    input:
    set val(samplename), file(reads) from ch_fastqc
    
  output:
  file "*_fastqc.{zip,html}" into ch_multiqc_fastqc

  script:
  """
  export PATH=/opt/conda/envs/nf-core-rnaseq-1.3/bin:$PATH

  fastqc -t ${task.cpus} -q $reads
  """
}

ch_rnaseq
  .into { ch_hisat2; ch_star; ch_star2; ch_salmon }


params.salmon_index = params.genome ? params.genomes[ params.genome ].salmon ?: false : false

params.salmon_trans_gene = params.genome ? params.genomes[ params.genome ].salmon_trans_gene ?: false : false

ch_salmon_index = params.run_salmon
    ? Channel.fromPath(params.salmon_index)
       .ifEmpty { exit 1, "Salmon index not found: ${params.salmon_index}" }
    : Channel.empty()

ch_salmon_trans_gene = params.run_salmon
    ?  Channel.fromPath(params.salmon_trans_gene)
       .ifEmpty { exit 1, "Salmon index not found: ${params.salmon_trans_gene}" }
    : Channel.empty()

process salmon {
    tag "salmon $samplename"
    memory = '10G'
    container "salmon"
    time '400m'
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
// Singularity salmon.img:~> whereis salmon
// salmon: /usr/local/bin/salmon
// Singularity salmon.img:~> echo $PATH
// /home/salmon-0.14.1/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
    containerOptions = "--bind /lustre/scratch115/projects/interval_wgs/nextflow/salmon14_index"
    // container "salmon0113"
// [vagrant@localhost ~]$   singularity build salmon0113.sif docker://combinelab/salmon:0.11.3
// Docker image path: index.docker.io/combinelab/salmon:0.11.3

    
    publishDir "${params.outdir}/Salmon", mode: 'copy'

    when:
    params.run_salmon

    input:
    set val(samplename), file(reads) from ch_salmon
    file index from ch_salmon_index.collect()
    file trans_gene from ch_salmon_trans_gene.collect()

    output:
    file "${samplename}.quant.sf" into ch_salmon_trans
    file "${samplename}.quant.genes.sf" into ch_salmon_genes
    file "my_outs/${samplename}" optional true into ch_alignment_logs_salmon

    script:
    """
    salmon quant \\
        -i /lustre/scratch115/projects/interval_wgs/nextflow/salmon14_index/salmon \\
        -l ISR \\
        -p ${task.cpus} \\
        --seqBias \\
        --gcBias \\
        --posBias \\
        --no-version-check \\
        -q \\
        -o . \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        -g /lustre/scratch115/projects/interval_wgs/nextflow/salmon14_index/trans_gene.txt \\
        --useVBOpt \\
        --numBootstraps 100
    mv quant.sf ${samplename}.quant.sf
    mv quant.genes.sf ${samplename}.quant.genes.sf
    mkdir -p my_outs/${samplename}/libParams
    mkdir -p my_outs/${samplename}/aux_info
    ln -f aux_info/meta_info.json my_outs/${samplename}/aux_info/meta_info.json
    ln -f libParams/flenDist.txt  my_outs/${samplename}/libParams/flenDist.txt
    """

    // TODO: prepare columns for merging; extract correct column and transpose (paste) it.
    // Include the row names so merger can check identity.
    // The merge step will concatenate the rows and re-transpose to obtain final result.
}

process merge_salmoncounts {
    tag "${input_trans}/${input_genes}"
    publishDir "${params.outdir}/combined", mode: 'copy'
    label 'merge_feature'
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    time '120m'

    input:
    file input_trans from ch_salmon_trans.map { it.toString() }.collectFile(name: 'trans.meta', newLine: true)
    file input_genes from ch_salmon_genes.map { it.toString() }.collectFile(name: 'genes.meta', newLine: true)

    when:
    params.run_salmon_mergecounts

    output:
    set file('*counts.txt'), file('*tpm.txt')

    script:
    def outtranscount = "${params.runtag}-salmon-transcounts.txt"
    def outgenescount = "${params.runtag}-salmon-genecounts.txt"
    def outtranstpm   = "${params.runtag}-salmon-transtpm.txt"
    def outgenestpm   = "${params.runtag}-salmon-genetpm.txt"
    """
    python3 $workflow.projectDir/bin/merge_featurecounts.py           \\
      --rm-suffix .quant.genes.sf                                     \\
      -c -1 --skip-comments --header                                  \\
      -o $outgenescount -I $input_genes
    python3 $workflow.projectDir/bin/merge_featurecounts.py           \\
      --rm-suffix .quant.sf                                           \\
      -c -1 --skip-comments --header                                  \\
      -o $outtranscount -I $input_trans
    python3 $workflow.projectDir/bin/merge_featurecounts.py           \\
      --rm-suffix .quant.genes.sf                                     \\
      -c -2 --skip-comments --header                                  \\
      -o $outgenestpm -I $input_genes
    python3 $workflow.projectDir/bin/merge_featurecounts.py           \\
      --rm-suffix .quant.sf                                           \\
      -c -2 --skip-comments --header                                  \\
      -o $outtranstpm -I $input_trans
    """
}


// ch_star.
//    map { a,b -> b  }.
//    set { ch_star_reads_only }

params.dna = params.genome ? params.genomes[ params.genome ].dna ?: false : false // added
ch_dna_star = params.dna 
    ? Channel.fromPath(params.dna)
      .ifEmpty { exit 1, "file genome fasta not found: ${params.dna}" }
    : Channel.empty() // added

params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
ch_star_index = params.run_star
    ? Channel.fromPath(params.star_index)
      .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
    : Channel.empty()

ch_dna_star2 = params.dna 
    ? Channel.fromPath(params.dna)
      .ifEmpty { exit 1, "file genome fasta not found: ${params.dna}" }
    : Channel.empty() // added


// cell geni accept/reject channels;
 // Filter removes all 'aligned' channels that fail the check
params.fcextra = ""                          // feature counts extra parameters; currently for testing
params.singleend = false
params.forward_stranded = false
params.reverse_stranded = false
params.unstranded = false
forward_stranded = params.forward_stranded
reverse_stranded = params.reverse_stranded
unstranded = params.unstranded

params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
Channel.fromPath(params.gtf)
  .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
  .into { ch_gtf_star; ch_gtf_featureCounts; }
params.biotypes_header= "$baseDir/assets/biotypes_header.txt"
biotypes_header = file(params.biotypes_header)
output_docs = file("$baseDir/docs/output.md")
gene_biotype = params.gtf.matches(".*gencode.*") ? "gene_type" : "gene_biotype"

process 'STAR2_map2' {
    tag "${samplename}"
    // container = 'callings-nf'
    container "nfcore-rnaseq"
    memory = '40G'
    cpus 2
    time '600m'
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    // errorStrategy = { task.exitStatus == 130 ? 'retry' : 'ignore' }
    // cpus = {  8 * Math.min(2, task.attempt) }
    // memory = {  40.GB * task.attempt * 1.6 ** (task.attempt - 1) }
    // Note above grows to about 100GB on 2 retries.
    
    publishDir "${params.outdir}/star/", mode: 'copy'
    publishDir "${params.outdir}/star2pass2/", mode: 'copy', pattern: "*.bam"
    publishDir "${params.outdir}/star2pass2/", mode: 'copy', pattern: "*.bam.bai"

  input:
    set val(samplename), file(reads) from ch_star // _reads_only
    file genomeDir from ch_star_index.collect()
    // file genome_fasta3 from ch_dna_star.collect()
    file gtf from ch_gtf_star.collect()

    when:
    params.run_star
    
  output:
    set val(samplename), file("*Log.final.out"), file ('*.bam') into star_aligned
    set val(samplename), file ('*.bam'), file ('*.bai') into star_aligned_with_bai
    // file "*.ReadsPerGene.out.tab" into ch_merge_starcounts
    file "*.out" into ch_alignment_logs_star
    file "*.SJ.out.tab"
    file "*.ReadsPerGene.out.tab" into ch_merge_starcounts

  script:

  """
  export PATH=/opt/conda/envs/nf-core-rnaseq-1.3/bin:$PATH

    STAR --genomeDir ${genomeDir} \\
        --sjdbGTFfile $gtf \\
        --readFilesIn $reads --readFilesCommand zcat \\
        --runThreadN ${task.cpus} \\
        --twopassMode Basic \\
        --outWigType bedGraph \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within \\
        --runDirPerm All_RWX \\
        --quantMode GeneCounts \\
        --outFileNamePrefix ${samplename}.

  # Index the BAM file
  samtools index ${samplename}.Aligned.sortedByCoord.out.bam
  rm -f Log.out 
  rm -f log.out 
  """
}

process 'bam2junc' {
    tag "${samplename}"
    // container = 'callings-nf'
    container "leafcutter"
    memory = '5G'
    cpus 1
    time '120m'
    // errorStrategy 'terminate'
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    
    publishDir "${params.outdir}/leafcutter/bam2junc", mode: 'copy', pattern: "*.junc"
    // publishDir "${params.outdir}/leafcutter/bam2junc", mode: 'copy', pattern: "*.bam.bed"

  input:
    set val(samplename), file (bamfile), file (baifile) from star_aligned_with_bai

    when:
    params.run_leafcutter_bam2junc
    
  output:
    file ('*.junc') into star_bam2junc
    file ('*.junc') into star_bam2junc2
    // file "*.ReadsPerGene.out.tab" into ch_merge_starcounts

  script:

  """
  export PATH=/home/leafcutter/scripts:/home/leafcutter/clustering:$PATH

  echo Converting ${bamfile} to ${samplename}.junc
  sh /home/leafcutter/scripts/bam2junc.sh ${bamfile} ${samplename}.junc
  """
}

star_bam2junc2
  .map{it -> it.toString()}
  .collectFile(name: 'juncfiles.txt', newLine: true)
  .set { ch_juncfiles }

process 'leafcutter_clustering' {
    tag "${samplename}"
    // container = 'callings-nf'
    container "leafcutter"
    memory = '5G'
    cpus 1
    time '300m'
    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }
    maxRetries 5
    
    publishDir "${params.outdir}/leafcutter/clustering", mode: 'copy', pattern: "*.junc.clust.sorted.gz"
    publishDir "${params.outdir}/leafcutter/clustering", mode: 'copy', pattern: "clust_*"
    // publishDir "${params.outdir}/leafcutter/bam2junc", mode: 'copy', pattern: "*.bam.bed"

    input:
    file ('juncfiles.txt') from ch_juncfiles
    file (junc_files) from star_bam2junc.collect()

    when:
    params.run_leafcutter_clustering
    
  output:
    set file('clust_perind.counts.gz'), file('clust_perind_numers.counts.gz'), file('clust_pooled'),file('clust_refined'),file('clust_sortedlibs'), file('*.junc.clust.sorted.gz') into star_clustering
    // file "*.ReadsPerGene.out.tab" into ch_merge_starcounts

  script:
  """
  export PATH=/home/leafcutter/scripts:/home/leafcutter/clustering:$PATH

  python /home/leafcutter/clustering/leafcutter_cluster.py -j juncfiles.txt -m 50 -o clust -l 500000
  """
}


ch_star_accept = Channel.create()
ch_star_reject = Channel.create()


def star_filter(logs) {
    def percent_aligned = 100
    logs.eachLine { line ->
        matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/
        if (matcher.matches()) {
            percent_aligned = matcher[0][1]
        }
    }
    if (percent_aligned.toFloat() < params.min_pct_aln.toFloat()) {
        n_star_lowmapping++
        return false
    } else {
        return true
    }
}

  star_aligned
      .choice(ch_star_accept, ch_star_reject)
          { namelogbam -> star_filter(namelogbam[1]) ? 0 : 1 }

  ch_star_accept
  .map    { name, log, bam -> ["star", name, bam] }
  .into   { ch_fc_star; ch_bam_star }

              // { it -> [text: "${it[0]}\tSTAR\tlowmapping\n"] }
              // ^ This channel output will be merged with the it.text from a file
              // in other channels. This is slightly hacky. Dangersign.
              // This in pursuit of keeping track of where we lose samples.
              // If this is to be rejigged, then it is probably easiest to
              // implement the alignment check in shell code, and use file-based
              // logic similarly as in process irods and process crams_to_fastq.
  ch_star_reject
  .map    { it -> [text: "${it[0]}\tSTAR\tlowmapping\n"] }
  .set    { ch_lostcause_star }

// skip hisat 2 
ch_fc_star
  .set{ ch_featurecounts }

ch_bam_star
  .into{ ch_indexbam; ch_publishbam }

process featureCounts {
    tag "${samplename}"
    container "nfcore-rnaseq"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    publishDir "${params.outdir}/featureCounts/", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".biotype_counts_mqc.txt") > 0) "biotype_counts/$filename"
            else if (filename.indexOf(".gene.featureCounts.txt.summary") > 0) "gene_count_summaries/$filename"
            else if (filename.indexOf(".gene.featureCounts.txt") > 0) "gene_counts/$filename"
            else "$filename"
        }

    input:
    set val(aligner), val(samplename), file(thebam) from ch_featurecounts
    file gtf from ch_gtf_featureCounts.collect()
    file biotypes_header

    output:
    set val(aligner), file("*.gene.fc.txt") into ch_merge_fc
    set val(aligner), file("*.gene.fc.txt.summary") into ch_multiqc_fc
    set val(aligner), file("*.biotype_counts*mqc.txt") into ch_multiqc_fcbiotype

    script:
    def extraparams = params.fcextra.toString() - ~/^dummy/
    def fc_direction = 0
    def tag = "${samplename}.${aligner}"

    def pairedend = params.singleend ? "" : "-p"
    if (forward_stranded && !unstranded) {
        fc_direction = 1
    } else if (reverse_stranded && !unstranded){
        fc_direction = 2
    }
    outfile = "${tag}.gene.fc.txt"
    """
    export PATH=/opt/conda/envs/nf-core-rnaseq-1.3/bin:$PATH

    featureCounts -T ${task.cpus} -a $gtf -g gene_id          \\
      -o ${outfile} $pairedend                                \\
      -s $fc_direction ${extraparams} $thebam
    cut -f 1,7 ${outfile} > reduced.${outfile}   #  This
    mv reduced.${outfile} ${outfile}             #  reduces the file size from ~ 30M to ~1M
    featureCounts -T ${task.cpus} -a $gtf -g gene_id  \\
      -o ${tag}.biotype.fc.txt $pairedend                     \\
      -s $fc_direction ${extraparams} $thebam
    cut -f 1,7 ${tag}.biotype.fc.txt |                        \\
        tail -n +3 | cat $biotypes_header - >> ${tag}.biotype_counts_mqc.txt
    """
}

// featureCounts -T ${task.cpus} -a $gtf -g gene_id          \\
// Currently this prefers star if both star and hisat2 are run, otherwise takes hisat2
// This is for processes that we only want to run for one aligner, not both;
// For example when publishing bams, or pushing bams to multiqc.
def pick_aligner(String aligner) {
    return  aligner == 'star' || (!params.run_star && aligner == 'hisat2')
      ? true
      : false
}
// Note: we lose information about the used aligner currently.
process indexbam {
    tag "${samplename}"
    container "nfcore-rnaseq"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6

    when:
    pick_aligner(aligner)

    input:
    set val(aligner), val(samplename), file(thebam) from ch_indexbam

    output:
    set val(samplename), file("*.idxstats") into ch_mapsummary

    script:
    """
    export PATH=/opt/conda/envs/nf-core-rnaseq-1.3/bin:$PATH

    samtools index $thebam
    samtools idxstats $thebam > ${samplename}.idxstats
    """
}

ch_publishbam
  .until{ !params.save_bam }
  .subscribe {
      aligner     = it[0]
      samplename  = it[1]
      thebam      = it[2]
      bamname     = thebam.toString()
      if (pick_aligner(aligner)) {
        dir = "${params.outdir}/${aligner}-bams"
        thebam.copyTo("$dir/${samplename.md5()[0..1]}/${samplename}-${aligner}.bam")
      }
  }

params.mito_name = 'MT'
process mapsummary {
    tag "${samplename}"
    container "nfcore-rnaseq"
    publishDir "${params.outdir}/mapsummary/", mode: 'copy'
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6

    input:
    set val(samplename), file(thestats) from ch_mapsummary

    output:
    file "*_mqc.txt" into ch_multiqc_mapsum

    script:
    def mito_name = params.mito_name
    """
    export PATH=/opt/conda/envs/nf-core-rnaseq-1.3/bin:$PATH

    python $baseDir/bin/mito.py -m ${mito_name} -t $thestats > ${samplename}_mqc.txt
    """
}
ch_merge_fc
  .transpose()
  .groupTuple()
  .collectFile { id, files -> [ id, files.collect{ it.toString() }.join('\n') + '\n' ] }
  .set{ ch_merge_fc_byaligner }


process merge_featureCounts {
    tag "$aligner"
    container "nfcore-rnaseq"
    publishDir "${params.outdir}/combined", mode: 'link'
    label 'merge_feature'
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6

    input:
    file metafile from ch_merge_fc_byaligner

    output:
    file '*-fc-genecounts.txt'

    shell:
    suffix=['star':'.star.gene.fc.txt', 'hisat2':'.hisat2.gene.fc.txt']
    aligner = metafile.baseName   // not strictly necessary
    outputname = "${params.runtag}-${aligner}-fc-genecounts.txt"
    thesuffix  = suffix[aligner] ?: '.txt'
    '''
    export PATH=/opt/conda/envs/nf-core-rnaseq-1.3/bin:$PATH

    python3 !{workflow.projectDir}/bin/merge_featurecounts.py        \\
      --rm-suffix !{thesuffix}                                       \\
      -c 1 --skip-comments --header                                  \\
      -o !{outputname} -I !{metafile}
    '''
}

ch_lostcause_cram.
  mix(ch_lostcause_star).
  set  {ch_lostcause}

process lostcause {

    publishDir "${params.outdir}/combined", mode: 'link'

    input:
    file (inputs) from ch_lostcause.collectFile{ ['lostcause.txt', it.text] }

    output:
    file('*.lostcause_mqc.txt') into ch_multiqc_lostcause

    script:
    def outputname = "${params.runtag}.${workflow.runName}.lostcause_mqc.txt"
    """
    echo -e "# plot_type: 'table'\n# section_name: 'Lost samples'" > $outputname
    echo -e "Sample\tProcess\tMessage" >> $outputname
    cat $inputs | sort >> $outputname
    """
}

ch_multiqc_fc
  .filter{ pick_aligner(it[0]) }
  .map { it[1] }
  .set{ ch_multiqc_fc_aligner }

ch_multiqc_fcbiotype
  .filter{ pick_aligner(it[0]) }
  .map{ it[1] }
  .set{ ch_multiqc_fcbiotype_aligner }

process multiqc {
    container "nfcore-rnaseq"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    publishDir "${params.outdir}", mode: 'link',
      saveAs: {filename ->
          if (filename.indexOf("multiqc.html") > 0) "combined/$filename"
          else if (filename.indexOf("_data") > 0) "$filename"
          else null
      }

    when:
    params.run_qc && params.run_multiqc

    input:
    file ('lostcause/*') from ch_multiqc_lostcause.collect().ifEmpty([])
    file (fastqc:'fastqc/*') from ch_multiqc_fastqc.collect().ifEmpty([])
    file ('mapsummary/*') from ch_multiqc_mapsum.collect().ifEmpty([])
    file ('featureCounts/*') from ch_multiqc_fc_aligner.collect().ifEmpty([])
    file ('featureCounts_biotype/*') from ch_multiqc_fcbiotype_aligner.collect().ifEmpty([])
    file ('star/*') from ch_alignment_logs_star.collect().ifEmpty([])

    output:
    file "*_multiqc.html"
    file "*_data"

    script:
    def filename = "${params.runtag}_multiqc.html"
    def reporttitle = "${params.runtag}"
    """
    export PATH=/opt/conda/envs/nf-core-rnaseq-1.3/bin:$PATH

    multiqc . -f --title "$reporttitle" --filename "$filename" -m featureCounts -m star -m fastqc
    """
}
// multiqc . -f --title "$reporttitle" --filename "$filename" -m custom_content -m featureCounts -m star -m fastqc -m salmon

