Hosts several Nextflow pipelines used by Sanger hgi team, the largest one is 

### Nextflow-RNAseq

A bioinformatics analysis pipeline used for RNA sequencing data, written in the new nextflow DSL2 language syntax, leveraging nextflow modules.  
DSL2 documentation: https://www.nextflow.io/docs/edge/dsl2.html     

### Credits
- The original pipeline was forked from the amazing nf-core/rnaseq pipeline, originally developed at the National Genomics Infrastructure at SciLifeLab in Stockholm, Sweden by Phil Ewels (@ewels) and Rickard Hammarén (@Hammarn).
- this DSL2 rewrite is based on the [cellgeni/rnaseq fork](https://github.com/cellgeni/rnaseq) (from the Cellular Genetics program at the Wellcome Sanger Institute).






    
    //.filter { it[1] ==~ /^[cC].*/} //.filter { it[1] ==~ /^[cC].*/}
    
    //// from irods studyid and list of samplenames
    //iget_cram(
    //	Channel.fromPath("${baseDir}/../../inputs/samples.txt")
    //	    .flatMap{ it.readLines()}, "5933")
    crams_to_fastq_gz(iget_cram.out[0])
    ////

    //// from cram files:
    ////Channel.fromPath('/lustre/scratch115/projects/interval_rna/inputs/*.cram').
    ////map{ it -> [ it.toString().replaceAll(~/.*\/(.*).cram/, "\$1"), it ] }.
    ////groupTuple(sort: true). //take(4).
    ////set{ch_cram_files}
    ////crams_to_fastq_gz(ch_cram_files)
    ////
    crams_to_fastq_gz.out[0]
	.map{ samplename, fastq1, fastq2 -> tuple( samplename, tuple(fastq1, fastq2) ) }
	.set{ch_samplename_crams}
    
    fastqc(ch_samplename_crams)

    salmon(ch_samplename_crams, ch_salmon_index.collect(), ch_salmon_trans_gene.collect())

    merge_salmoncounts(salmon.out[0].collect(), salmon.out[1].collect())

    tximport(salmon.out[0].collect())
    heatmap(merge_salmoncounts.out[0].map{transcounts,transtpm,genecouts,genetpm-> genecouts})
    
    star_2pass_basic(ch_samplename_crams, ch_star_index.collect(), ch_gtf_star.collect())

    //star_2pass_1st_pass(ch_samplename_crams, ch_star_index.collect(), ch_gtf_star.collect())
    //star_2pass_merge_junctions(star_2pass_1st_pass.out[1].collect())
    //star_2pass_2nd_pass(ch_samplename_crams, ch_star_index.collect(), ch_gtf_star.collect(), star_2pass_merge_junctions.out)

    // star_out = star_2pass_2nd_pass.out // choose star_2pass_basic.out or star_2pass_2ndpass.out 
    star_out = star_2pass_basic.out


    //// 
    star_tabgenes_matrix(star_2pass_basic.out.samplename_readspergene_tab.collect())
    
    leafcutter_bam2junc(star_out[0])
    leafcutter_clustering(leafcutter_bam2junc.out.collect())

    filter_star_aln_rate(star_out[1].map{samplename,logfile,bamfile -> [samplename,logfile]}) // discard bam file, only STAR log required to filter
    
    filter_star_aln_rate.out.branch {
        filtered: it[1] == 'above_threshold'
        discarded: it[1] == 'below_threshold'}.set { star_filter }
    
    star_filter.filtered.combine(star_out[1], by:0) //reattach bam file
	.map{samplename,filter,logfile,bamfile -> ["star", samplename, bamfile]} // discard log file and attach aligner name
	.set{star_filtered} 
    
    samtools_index_idxstats(star_filtered)
    
    mapsummary(samtools_index_idxstats.out)
    
    featureCounts(star_filtered, ch_gtf_star.collect(), ch_biotypes_header.collect())

    merge_featureCounts(featureCounts.out[0].map{samplename, gene_fc_txt -> gene_fc_txt}.collect())

    crams_to_fastq_gz.out[1]
	.mix(star_filter.discarded.map{samplename, filter -> [text: "${samplename}\tSTAR\tlowmapping\n"]})
	.set{ch_lostcause }
    lostcause(ch_lostcause.collectFile({ ['lostcause.txt', it.text]},sort:true))

    featureCounts.out[1]
	.filter{ pick_aligner(it[0]) }
	.map { it[1] }
	.set{ ch_multiqc_fc_aligner }

    featureCounts.out[2]
	.filter{ pick_aligner(it[0]) }
	.map{ it[1] }
	.set{ ch_multiqc_fcbiotype_aligner }
    
    multiqc(lostcause.out.collect().ifEmpty([]),
	    fastqc.out.collect().ifEmpty([]),
	    mapsummary.out.collect().ifEmpty([]),
	    ch_multiqc_fc_aligner.collect().ifEmpty([]),
	    ch_multiqc_fcbiotype_aligner.collect().ifEmpty([]),
	    star_out[2].collect().ifEmpty([]),
	    salmon.out[2].collect().ifEmpty([]))
 
    if(params.run_deseq2)
	deseq2(salmon.out[0].collect(), Channel.fromPath(params.deseq2_tsv))
       